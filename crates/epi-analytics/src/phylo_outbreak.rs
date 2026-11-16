// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

use crate::math::{gamma_pdf, ln_gamma};
use epi_core::EpiError;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Joint genomic + epidemiological outbreak transmission tree reconstruction (`Outbreaker2` / `seqTrack` class)
/// combining aligned FASTA SNP / K2P distances with symptom-onset generation-interval likelihoods,
/// unsampled intermediate host inference (`kappa >= 2`), multiple index introductions, sequence QC gates,
/// and Newick / Graphviz DOT export.
///
/// References:
/// - Combined genomic + epidemiological transmission tree reconstruction: https://pmc.ncbi.nlm.nih.gov/articles/PMC8875843/
/// - Multi-introduction coalescent "history host" transmission tree inference: https://pmc.ncbi.nlm.nih.gov/articles/PMC10703282/

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutbreakIsolateInput {
    pub case_id: String,
    pub onset_day: f64,
    pub onset_label: String,
    pub location: String,
    pub quarantined: bool,
    pub aligned_sequence: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct OutbreakTreeConfig {
    /// Minimum required fraction of unambiguous A/C/G/T sites across alignment (default `0.90` = 90%).
    pub min_unambiguous_coverage: f64,
    /// Exclude sequences flagged `quarantined = true` in metadata while tracking them for re-audit.
    pub exclude_quarantined: bool,
    /// Mean generation / serial interval in days.
    pub generation_mean_days: f64,
    /// Standard deviation of generation / serial interval in days.
    pub generation_sd_days: f64,
    /// Expected genomic substitution rate in SNPs per genome per day (`mu * L`).
    pub mutation_rate_snps_per_day: f64,
    /// Within-host / bottleneck clock jitter in SNPs per transmission event.
    pub bottleneck_jitter_snps: f64,
    /// Maximum transmission chain hops `kappa` (`1` = direct, `2..=max_kappa` = `kappa - 1` unsampled hosts).
    pub max_kappa_hops: usize,
    /// Log-prior threshold for an independent index importation vs. local transmission link.
    pub import_log_likelihood_floor: f64,
}

impl Default for OutbreakTreeConfig {
    fn default() -> Self {
        Self {
            min_unambiguous_coverage: 0.90,
            exclude_quarantined: true,
            generation_mean_days: 11.7,
            generation_sd_days: 2.5,
            mutation_rate_snps_per_day: 0.045,
            bottleneck_jitter_snps: 0.35,
            max_kappa_hops: 3,
            import_log_likelihood_floor: -11.5,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwiseGenomicDistance {
    pub case_a: String,
    pub case_b: String,
    pub comparable_sites: usize,
    pub snp_count: usize,
    pub transitions: usize,
    pub transversions: usize,
    /// Kimura 2-Parameter (`K2P`) evolutionary distance per site.
    pub k2p_distance: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransmissionEdge {
    pub infector_id: Option<String>,
    pub infectee_id: String,
    pub infectee_location: String,
    pub infectee_onset_label: String,
    pub infectee_onset_day: f64,
    pub is_importation: bool,
    pub onset_lag_days: Option<f64>,
    pub snp_distance: Option<usize>,
    pub k2p_distance: Option<f64>,
    /// Inferred unsampled intermediate hosts (`kappa - 1`, `0` = direct link).
    pub unsampled_intermediates: usize,
    /// Posterior ancestor probability `P(infector -> infectee | D, t)` in `[0, 1]`.
    pub posterior_support: f64,
    pub edge_log_likelihood: f64,
    pub cluster_id: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CaseOffspringMetric {
    pub case_id: String,
    pub location: String,
    pub onset_label: String,
    pub map_secondary_cases: usize,
    /// Expected secondary cases summing posterior probabilities `sum_j P(i -> j)`.
    pub expected_ri: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutbreakTreeResult {
    pub alignment_length: usize,
    pub retained_cases: usize,
    pub quarantined_excluded_ids: Vec<String>,
    pub low_coverage_excluded_ids: Vec<(String, f64)>,
    pub total_importations: usize,
    pub total_unsampled_intermediates: usize,
    pub mean_case_ri: f64,
    pub max_superspreader_cases: usize,
    pub total_tree_log_likelihood: f64,
    pub edges: Vec<TransmissionEdge>,
    pub case_offspring: Vec<CaseOffspringMetric>,
    pub pairwise_distances: Vec<PairwiseGenomicDistance>,
    pub newick_tree: String,
}

fn is_unambiguous_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

fn is_transition(a: u8, b: u8) -> bool {
    matches!(
        (a, b),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

fn compute_pairwise_distance(
    id_a: &str,
    seq_a: &[u8],
    id_b: &str,
    seq_b: &[u8],
) -> PairwiseGenomicDistance {
    let len = seq_a.len().min(seq_b.len());
    let mut comp = 0usize;
    let mut snps = 0usize;
    let mut ti = 0usize;
    let mut tv = 0usize;

    for idx in 0..len {
        let ba = seq_a[idx].to_ascii_uppercase();
        let bb = seq_b[idx].to_ascii_uppercase();
        if is_unambiguous_base(ba) && is_unambiguous_base(bb) {
            comp += 1;
            if ba != bb {
                snps += 1;
                if is_transition(ba, bb) {
                    ti += 1;
                } else {
                    tv += 1;
                }
            }
        }
    }

    let k2p = if comp > 0 {
        let p = (ti as f64) / (comp as f64);
        let q = (tv as f64) / (comp as f64);
        let term1 = (1.0 - 2.0 * p - q).max(1e-6);
        let term2 = (1.0 - 2.0 * q).max(1e-6);
        (-0.5 * term1.ln() - 0.25 * term2.ln()).max(0.0)
    } else {
        0.0
    };

    PairwiseGenomicDistance {
        case_a: id_a.to_string(),
        case_b: id_b.to_string(),
        comparable_sites: comp,
        snp_count: snps,
        transitions: ti,
        transversions: tv,
        k2p_distance: k2p,
    }
}

fn log_poisson_pmf(k: usize, lambda: f64) -> f64 {
    let lam = lambda.max(1e-8);
    let kf = k as f64;
    kf * lam.ln() - lam - ln_gamma(kf + 1.0)
}

fn log_sum_exp(vals: &[f64]) -> f64 {
    let m = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !m.is_finite() {
        return f64::NEG_INFINITY;
    }
    let s: f64 = vals.iter().map(|&v| (v - m).exp()).sum();
    m + s.ln()
}

/// Reconstruct a joint genomic + symptom-onset outbreak transmission forest, enforcing QC thresholds,
/// inferring unsampled intermediate hosts (`kappa - 1`), and separating independent index introductions.
pub fn reconstruct_outbreak_tree(
    isolates: &[OutbreakIsolateInput],
    cfg: &OutbreakTreeConfig,
) -> Result<OutbreakTreeResult, EpiError> {
    if isolates.is_empty() {
        return Err(EpiError::Schema("No outbreak isolates supplied".to_string()));
    }

    let mut quarantined_excluded_ids = Vec::new();
    let mut low_coverage_excluded_ids = Vec::new();
    let mut retained: Vec<OutbreakIsolateInput> = Vec::new();

    let align_len = isolates
        .iter()
        .map(|s| s.aligned_sequence.len())
        .max()
        .unwrap_or(0);
    if align_len == 0 {
        return Err(EpiError::Schema(
            "Aligned FASTA sequences have zero length".to_string(),
        ));
    }

    for iso in isolates {
        if cfg.exclude_quarantined && iso.quarantined {
            quarantined_excluded_ids.push(iso.case_id.clone());
            continue;
        }
        let bytes = iso.aligned_sequence.as_bytes();
        if bytes.len() != align_len {
            return Err(EpiError::Schema(format!(
                "Isolate '{}' sequence length {} does not match alignment length {}",
                iso.case_id,
                bytes.len(),
                align_len
            )));
        }
        let unambig = bytes
            .iter()
            .filter(|&&b| is_unambiguous_base(b.to_ascii_uppercase()))
            .count();
        let cov = (unambig as f64) / (align_len as f64);
        if cov + 1e-9 < cfg.min_unambiguous_coverage {
            low_coverage_excluded_ids.push((iso.case_id.clone(), cov));
            continue;
        }
        retained.push(iso.clone());
    }

    if retained.len() < 2 {
        return Err(EpiError::Schema(
            "Fewer than 2 isolates passed quarantine and unambiguous site coverage QC".to_string(),
        ));
    }

    retained.sort_by(|a, b| {
        a.onset_day
            .partial_cmp(&b.onset_day)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.case_id.cmp(&b.case_id))
    });

    let n = retained.len();
    let mut dist_lookup: HashMap<(usize, usize), PairwiseGenomicDistance> = HashMap::new();
    let mut pairwise_distances = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            let d = compute_pairwise_distance(
                &retained[i].case_id,
                retained[i].aligned_sequence.as_bytes(),
                &retained[j].case_id,
                retained[j].aligned_sequence.as_bytes(),
            );
            dist_lookup.insert((i, j), d.clone());
            dist_lookup.insert(
                (j, i),
                PairwiseGenomicDistance {
                    case_a: d.case_b.clone(),
                    case_b: d.case_a.clone(),
                    ..d.clone()
                },
            );
            pairwise_distances.push(d);
        }
    }

    let mu_si = cfg.generation_mean_days.max(1.0);
    let sd_si = cfg.generation_sd_days.max(0.5);
    let base_shape = (mu_si * mu_si) / (sd_si * sd_si);
    let base_rate = mu_si / (sd_si * sd_si);
    let max_k = cfg.max_kappa_hops.clamp(1, 6);

    let mut expected_ri = vec![0.0_f64; n];
    let mut map_children_count = vec![0usize; n];
    let mut parent_idx: Vec<Option<usize>> = vec![None; n];
    let mut raw_edges = Vec::with_capacity(n);
    let mut total_tree_ll = 0.0;
    let mut total_unsampled = 0usize;

    for j in 0..n {
        let mut cand_parents: Vec<(Option<usize>, usize, f64)> =
            vec![(None, 0, cfg.import_log_likelihood_floor)];

        for i in 0..j {
            let dt = retained[j].onset_day - retained[i].onset_day;
            if dt < 0.5 {
                continue;
            }
            let d_ij = &dist_lookup[&(i, j)];

            let mut best_kappa = 1usize;
            let mut best_hop_ll = f64::NEG_INFINITY;
            for kappa in 1..=max_k {
                let kf = kappa as f64;
                let t_pdf = gamma_pdf(dt, kf * base_shape, base_rate).max(1e-30);
                let expected_snps = cfg.mutation_rate_snps_per_day.max(1e-5) * dt
                    + kf * cfg.bottleneck_jitter_snps.max(0.05);
                let g_ll = log_poisson_pmf(d_ij.snp_count, expected_snps);
                let hop_ll = t_pdf.ln() + g_ll - 0.65 * ((kappa - 1) as f64);
                if hop_ll > best_hop_ll {
                    best_hop_ll = hop_ll;
                    best_kappa = kappa;
                }
            }
            cand_parents.push((Some(i), best_kappa - 1, best_hop_ll));
        }

        let log_vals: Vec<f64> = cand_parents.iter().map(|c| c.2).collect();
        let lse = log_sum_exp(&log_vals);
        let probs: Vec<f64> = log_vals.iter().map(|&lv| (lv - lse).exp()).collect();

        for (c_idx, &(maybe_i, _, _)) in cand_parents.iter().enumerate() {
            if let Some(i) = maybe_i {
                expected_ri[i] += probs[c_idx];
            }
        }

        let mut best_c = 0usize;
        for c_idx in 1..cand_parents.len() {
            if cand_parents[c_idx].2 > cand_parents[best_c].2 {
                best_c = c_idx;
            }
        }

        let (chosen_parent, unsampled_hops, best_ll) = cand_parents[best_c];
        let post_p = probs[best_c];
        total_tree_ll += best_ll;

        if let Some(i) = chosen_parent {
            parent_idx[j] = Some(i);
            map_children_count[i] += 1;
            total_unsampled += unsampled_hops;
            let d_ij = &dist_lookup[&(i, j)];
            raw_edges.push(TransmissionEdge {
                infector_id: Some(retained[i].case_id.clone()),
                infectee_id: retained[j].case_id.clone(),
                infectee_location: retained[j].location.clone(),
                infectee_onset_label: retained[j].onset_label.clone(),
                infectee_onset_day: retained[j].onset_day,
                is_importation: false,
                onset_lag_days: Some(retained[j].onset_day - retained[i].onset_day),
                snp_distance: Some(d_ij.snp_count),
                k2p_distance: Some(d_ij.k2p_distance),
                unsampled_intermediates: unsampled_hops,
                posterior_support: post_p,
                edge_log_likelihood: best_ll,
                cluster_id: 0,
            });
        } else {
            raw_edges.push(TransmissionEdge {
                infector_id: None,
                infectee_id: retained[j].case_id.clone(),
                infectee_location: retained[j].location.clone(),
                infectee_onset_label: retained[j].onset_label.clone(),
                infectee_onset_day: retained[j].onset_day,
                is_importation: true,
                onset_lag_days: None,
                snp_distance: None,
                k2p_distance: None,
                unsampled_intermediates: 0,
                posterior_support: post_p,
                edge_log_likelihood: best_ll,
                cluster_id: 0,
            });
        }
    }

    let mut cluster_ids = vec![0usize; n];
    let mut next_cluster = 0usize;
    for j in 0..n {
        if let Some(i) = parent_idx[j] {
            cluster_ids[j] = cluster_ids[i];
        } else {
            next_cluster += 1;
            cluster_ids[j] = next_cluster;
        }
        raw_edges[j].cluster_id = cluster_ids[j];
    }

    let case_offspring: Vec<CaseOffspringMetric> = (0..n)
        .map(|i| CaseOffspringMetric {
            case_id: retained[i].case_id.clone(),
            location: retained[i].location.clone(),
            onset_label: retained[i].onset_label.clone(),
            map_secondary_cases: map_children_count[i],
            expected_ri: expected_ri[i],
        })
        .collect();

    let mean_ri = expected_ri.iter().sum::<f64>() / (n as f64);
    let max_super = map_children_count.iter().copied().max().unwrap_or(0);
    let newick_tree = build_newick_forest(&retained, &parent_idx, &raw_edges);

    Ok(OutbreakTreeResult {
        alignment_length: align_len,
        retained_cases: n,
        quarantined_excluded_ids,
        low_coverage_excluded_ids,
        total_importations: next_cluster,
        total_unsampled_intermediates: total_unsampled,
        mean_case_ri: mean_ri,
        max_superspreader_cases: max_super,
        total_tree_log_likelihood: total_tree_ll,
        edges: raw_edges,
        case_offspring,
        pairwise_distances,
        newick_tree,
    })
}

fn build_newick_subtree(
    node: usize,
    children_map: &[Vec<usize>],
    isolates: &[OutbreakIsolateInput],
    edges: &[TransmissionEdge],
) -> String {
    let kids = &children_map[node];
    let snp_len = edges[node].snp_distance.unwrap_or(0);
    let safe_id = isolates[node]
        .case_id
        .replace(['(', ')', ':', ';', ',', ' '], "_");
    if kids.is_empty() {
        format!("{safe_id}:{snp_len}")
    } else {
        let child_strs: Vec<String> = kids
            .iter()
            .map(|&ch| build_newick_subtree(ch, children_map, isolates, edges))
            .collect();
        format!("({safe_id}:0,{0}):{snp_len}", child_strs.join(","))
    }
}

fn build_newick_forest(
    isolates: &[OutbreakIsolateInput],
    parent_idx: &[Option<usize>],
    edges: &[TransmissionEdge],
) -> String {
    let n = isolates.len();
    let mut children_map = vec![Vec::new(); n];
    let mut roots = Vec::new();
    for j in 0..n {
        if let Some(i) = parent_idx[j] {
            children_map[i].push(j);
        } else {
            roots.push(j);
        }
    }
    if roots.len() == 1 {
        format!(
            "{};",
            build_newick_subtree(roots[0], &children_map, isolates, edges)
        )
    } else {
        let subtrees: Vec<String> = roots
            .into_iter()
            .map(|r| build_newick_subtree(r, &children_map, isolates, edges))
            .collect();
        format!("({})HistoryHost;", subtrees.join(","))
    }
}

/// Export reconstructed outbreak transmission tree as a Graphviz DOT directed network.
pub fn export_outbreak_graphviz_dot(res: &OutbreakTreeResult) -> String {
    let mut out = String::from("digraph OutbreakTransmissionTree {\n  rankdir=LR;\n  node [shape=box, style=\"rounded,filled\", fillcolor=\"#1e293b\", fontcolor=\"#f8fafc\"];\n");
    for edge in &res.edges {
        if let Some(ref inf) = edge.infector_id {
            let style = if edge.unsampled_intermediates > 0 {
                "dashed"
            } else {
                "solid"
            };
            out.push_str(&format!(
                "  \"{}\" -> \"{}\" [style={}, label=\"{} SNP | {:.0}d | p={:.2} | unsamp={}\"];\n",
                inf,
                edge.infectee_id,
                style,
                edge.snp_distance.unwrap_or(0),
                edge.onset_lag_days.unwrap_or(0.0),
                edge.posterior_support,
                edge.unsampled_intermediates
            ));
        } else {
            out.push_str(&format!(
                "  \"Import_Cluster_{}\" [shape=ellipse, fillcolor=\"#0f766e\"];\n  \"Import_Cluster_{}\" -> \"{}\" [style=dotted, label=\"Index p={:.2}\"];\n",
                edge.cluster_id, edge.cluster_id, edge.infectee_id, edge.posterior_support
            ));
        }
    }
    out.push_str("}\n");
    out
}

fn parse_iso_date_to_days(date_str: &str) -> Option<f64> {
    if let Ok(v) = date_str.trim().parse::<f64>() {
        return Some(v);
    }
    let parts: Vec<&str> = date_str.trim().split('-').collect();
    if parts.len() == 3 {
        let y: i64 = parts[0].parse().ok()?;
        let m: i64 = parts[1].parse().ok()?;
        let d: i64 = parts[2].parse().ok()?;
        let y_adj = y - if m <= 2 { 1 } else { 0 };
        let era = if y_adj >= 0 { y_adj } else { y_adj - 399 } / 400;
        let yoe = y_adj - era * 400;
        let doy = (153 * (m + if m > 2 { -3 } else { 9 }) + 2) / 5 + d - 1;
        let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
        return Some((era * 146097 + doe - 719468) as f64);
    }
    None
}

/// Parse multi-FASTA text + CSV/TSV symptom-onset metadata (`case_id,onset_date,location[,quarantined]`)
/// into `Vec<OutbreakIsolateInput>`.
pub fn parse_fasta_and_metadata(
    fasta_text: &str,
    metadata_csv: &str,
) -> Result<Vec<OutbreakIsolateInput>, EpiError> {
    let mut seq_map: HashMap<String, String> = HashMap::new();
    let mut cur_id: Option<String> = None;
    let mut cur_seq = String::new();

    for line in fasta_text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with(';') {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(prev_id) = cur_id.take() {
                seq_map.insert(prev_id, cur_seq.clone());
                cur_seq.clear();
            }
            let id_token = rest
                .split(|c: char| c.is_whitespace() || c == '|')
                .next()
                .unwrap_or("")
                .trim()
                .to_string();
            cur_id = Some(id_token);
        } else {
            cur_seq.extend(trimmed.chars().filter(|c| !c.is_whitespace()));
        }
    }
    if let Some(prev_id) = cur_id {
        seq_map.insert(prev_id, cur_seq);
    }
    if seq_map.is_empty() {
        return Err(EpiError::Schema(
            "No FASTA records found in sequence input".to_string(),
        ));
    }

    let mut lines = metadata_csv
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let header = lines
        .next()
        .ok_or_else(|| EpiError::Schema("Empty outbreak metadata CSV".to_string()))?;
    let delim = if header.contains('\t') { '\t' } else { ',' };

    let mut raw_rows = Vec::new();
    for line in lines {
        let cols: Vec<&str> = line.split(delim).map(|s| s.trim().trim_matches('"')).collect();
        if cols.len() < 2 {
            continue;
        }
        let cid = cols[0].to_string();
        let Some(seq) = seq_map.get(&cid) else {
            continue;
        };
        let onset_raw = cols[1].to_string();
        let abs_day = parse_iso_date_to_days(&onset_raw).ok_or_else(|| {
            EpiError::Schema(format!(
                "Cannot parse onset date '{}' for isolate '{}'",
                onset_raw, cid
            ))
        })?;
        let loc = cols.get(2).copied().unwrap_or("Unspecified").to_string();
        let quar = cols
            .get(3)
            .map(|s| {
                matches!(
                    s.to_ascii_lowercase().as_str(),
                    "true" | "1" | "yes" | "quarantined"
                )
            })
            .unwrap_or(false);
        raw_rows.push((cid, abs_day, onset_raw, loc, quar, seq.clone()));
    }

    if raw_rows.is_empty() {
        return Err(EpiError::Schema(
            "No FASTA isolate IDs matched metadata CSV rows".to_string(),
        ));
    }
    let min_day = raw_rows
        .iter()
        .map(|r| r.1)
        .fold(f64::INFINITY, f64::min);

    Ok(raw_rows
        .into_iter()
        .map(
            |(case_id, d, onset_label, location, quarantined, aligned_sequence)| {
                OutbreakIsolateInput {
                    case_id,
                    onset_day: d - min_day,
                    onset_label,
                    location,
                    quarantined,
                    aligned_sequence,
                }
            },
        )
        .collect())
}

/// Built-in 14-isolate multi-chain measles genotype D8 outbreak dataset (200-bp core variable window)
/// featuring 2 independent index introductions (`Cluster 1` & `Cluster 2`), a Clinic-A superspreading event,
/// an unsampled intermediate link (`kappa = 2`), 1 quarantined specimen, and 1 low-coverage (`< 90%`) specimen.
pub fn builtin_outbreak_fasta_and_metadata() -> (String, String) {
    let base_c1: Vec<u8> = b"ACGTACGTAGCTAGCTAGGCTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGGCTAGCT".to_vec();

    let mutate = |src: &[u8], edits: &[(usize, u8)]| -> String {
        let mut v = src.to_vec();
        for &(pos, ch) in edits {
            if pos < v.len() {
                v[pos] = ch;
            }
        }
        String::from_utf8(v).unwrap()
    };

    let s01 = mutate(&base_c1, &[]);
    let s02 = mutate(&base_c1, &[]);
    let s03 = mutate(&base_c1, &[(14, b'G')]);
    let s04 = mutate(&base_c1, &[(52, b'C')]);
    let s05 = mutate(&base_c1, &[(14, b'G'), (88, b'A')]);
    let s06 = mutate(&base_c1, &[(52, b'C')]);
    let s07 = mutate(&base_c1, &[(14, b'G'), (88, b'A'), (121, b'T')]);
    let s08 = mutate(&base_c1, &[(52, b'C'), (140, b'G'), (175, b'A')]);

    let c2_edits: &[(usize, u8)] = &[
        (3, b'C'),
        (11, b'A'),
        (19, b'T'),
        (27, b'C'),
        (35, b'A'),
        (43, b'G'),
        (61, b'T'),
        (69, b'C'),
        (77, b'A'),
        (95, b'G'),
        (103, b'T'),
        (115, b'C'),
        (131, b'A'),
        (149, b'G'),
        (163, b'T'),
        (187, b'C'),
    ];
    let base_c2 = mutate(&base_c1, c2_edits);
    let s09 = base_c2.clone();
    let s10 = mutate(base_c2.as_bytes(), &[(180, b'A')]);
    let s11 = mutate(base_c2.as_bytes(), &[(180, b'A'), (192, b'G')]);
    let s12 = mutate(base_c2.as_bytes(), &[(195, b'T')]);

    let s13 = mutate(&base_c1, &[(14, b'G')]);
    let mut low_cov = base_c1.clone();
    for slot in low_cov.iter_mut().take(50) {
        *slot = b'N';
    }
    let s14 = String::from_utf8(low_cov).unwrap();

    let fasta = format!(
        ">MSL-01\n{s01}\n>MSL-02\n{s02}\n>MSL-03\n{s03}\n>MSL-04\n{s04}\n>MSL-05\n{s05}\n>MSL-06\n{s06}\n>MSL-07\n{s07}\n>MSL-08\n{s08}\n>MSL-09\n{s09}\n>MSL-10\n{s10}\n>MSL-11\n{s11}\n>MSL-12\n{s12}\n>MSL-Q13\n{s13}\n>MSL-LC14\n{s14}\n"
    );

    let meta = "\
case_id,onset_date,location,quarantined
MSL-01,2026-08-01,Clinic-A,false
MSL-02,2026-08-12,Clinic-A,false
MSL-03,2026-08-13,School-North,false
MSL-04,2026-08-14,Clinic-A,false
MSL-05,2026-08-24,School-North,false
MSL-06,2026-08-26,Household-West,false
MSL-07,2026-09-04,School-North,false
MSL-08,2026-09-18,Community-East,false
MSL-09,2026-08-05,Airport-Import,false
MSL-10,2026-08-17,Transit-South,false
MSL-11,2026-08-29,Transit-South,false
MSL-12,2026-08-18,Hotel-Central,false
MSL-Q13,2026-08-15,Lab-Contam-Hold,true
MSL-LC14,2026-08-20,Degraded-Swab,false
"
    .to_string();

    (fasta, meta)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_outbreak_tree_qc_multi_intro_and_unsampled_hosts() {
        let (fa, meta) = builtin_outbreak_fasta_and_metadata();
        let isolates = parse_fasta_and_metadata(&fa, &meta).expect("parse builtin outbreak");
        assert_eq!(isolates.len(), 14);

        let cfg = OutbreakTreeConfig::default();
        let tree = reconstruct_outbreak_tree(&isolates, &cfg).expect("reconstruct tree");

        assert_eq!(tree.retained_cases, 12);
        assert_eq!(tree.quarantined_excluded_ids, vec!["MSL-Q13".to_string()]);
        assert_eq!(tree.low_coverage_excluded_ids.len(), 1);
        assert_eq!(tree.low_coverage_excluded_ids[0].0, "MSL-LC14");
        assert_eq!(tree.total_importations, 2);
        assert!(tree.total_unsampled_intermediates >= 1);
        assert!(tree.newick_tree.ends_with("HistoryHost;"));
    }
}
