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

use crate::lineage_logistics::{LineageCountRow, LineageCountTable};
use epi_core::EpiError;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// Depth weighting scheme for wastewater SNV deconvolution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DepthWeighting {
    /// Freyja standard weighting w_m = log2(depth_m + 1).
    Log2DepthPlusOne,
    /// Square-root read-depth weighting w_m = sqrt(depth_m).
    SqrtDepth,
    /// Uniform site weighting w_m = 1.
    Uniform,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WastewaterDemixConfig {
    pub depth_weighting: DepthWeighting,
    /// Minimum sequencing depth required to include an SNV site in optimization.
    pub min_site_depth: u32,
    /// Number of deterministic site-resampling bootstrap replicates for 95% CIs (0 disables).
    pub bootstrap_replicates: usize,
    /// Positive residual threshold (b_m - (A p)_m) flagging potential cryptic / uncatalogued mutations.
    pub cryptic_residual_cutoff: f64,
    /// Pseudo-sequencing total count per week when exporting to LineageCountTable.
    pub pseudo_total_reads_per_sample: f64,
}

impl Default for WastewaterDemixConfig {
    fn default() -> Self {
        Self {
            depth_weighting: DepthWeighting::Log2DepthPlusOne,
            min_site_depth: 10,
            bootstrap_replicates: 80,
            cryptic_residual_cutoff: 0.08,
            pseudo_total_reads_per_sample: 500.0,
        }
    }
}

/// Reference mutation barcode library mapping M signature SNVs across K target lineages.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BarcodeLibrary {
    pub lineage_names: Vec<String>,
    pub mutation_sites: Vec<String>,
    /// Row-major M x K matrix where barcode_matrix[m][k] in [0, 1] is expected mutation frequency.
    pub barcode_matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnvObservation {
    pub mutation_site: String,
    pub alt_allele_freq: f64,
    pub read_depth: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WastewaterSampleInput {
    pub sample_id: String,
    pub date_label: String,
    pub time_weeks: f64,
    pub sewershed: String,
    pub viral_copies_per_liter: f64,
    pub observations: Vec<SnvObservation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageAbundanceEstimate {
    pub lineage: String,
    pub abundance: f64,
    pub ci_low_95: f64,
    pub ci_high_95: f64,
    pub viral_load_weighted_copies: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrypticResidualSite {
    pub mutation_site: String,
    pub observed_af: f64,
    pub reconstructed_af: f64,
    pub residual_excess: f64,
    pub read_depth: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleDemixResult {
    pub sample_id: String,
    pub date_label: String,
    pub time_weeks: f64,
    pub sewershed: String,
    pub viral_copies_per_liter: f64,
    pub sites_used: usize,
    pub sites_total_library: usize,
    pub coverage_breadth: f64,
    pub weighted_l1_loss: f64,
    pub weighted_rmse: f64,
    pub lineages: Vec<LineageAbundanceEstimate>,
    pub cryptic_sites: Vec<CrypticResidualSite>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LongitudinalWastewaterResult {
    pub lineage_names: Vec<String>,
    pub samples: Vec<SampleDemixResult>,
}

impl LongitudinalWastewaterResult {
    /// Convert deconvolved longitudinal lineage proportions and viral loads into a `LineageCountTable`
    /// directly consumable by `fit_multinomial_lineage_logistics` for variant fitness and wastewater Rt.
    pub fn to_lineage_count_table(&self, pseudo_total_per_sample: f64) -> LineageCountTable {
        let base_scale = pseudo_total_per_sample.max(25.0);
        let mean_copies = if self.samples.is_empty() {
            1.0
        } else {
            (self
                .samples
                .iter()
                .map(|s| s.viral_copies_per_liter.max(1.0))
                .sum::<f64>()
                / (self.samples.len() as f64))
                .max(1.0)
        };

        let mut rows = Vec::with_capacity(self.samples.len());
        for s in &self.samples {
            let total_equiv = (s.viral_copies_per_liter.max(1.0) / mean_copies) * base_scale;
            let mut clade_counts = Vec::with_capacity(self.lineage_names.len());
            for name in &self.lineage_names {
                let p = s
                    .lineages
                    .iter()
                    .find(|l| l.lineage == *name)
                    .map(|l| l.abundance)
                    .unwrap_or(0.0);
                clade_counts.push((p * base_scale).max(0.0));
            }
            rows.push(LineageCountRow {
                period_label: s.date_label.clone(),
                total_cases: total_equiv.max(1.0),
                clade_counts,
            });
        }

        let shed = self
            .samples
            .first()
            .map(|s| s.sewershed.as_str())
            .unwrap_or("Sewershed");
        LineageCountTable {
            title: format!("Freyja L1 Deconvolved Wastewater Lineages ({shed})"),
            clade_names: self.lineage_names.clone(),
            rows,
        }
    }
}

/// Exact O(K log K) Euclidean projection of vector `v` onto the unit probability simplex {p >= 0, sum(p) = 1} (Duchi et al. 2008).
pub fn project_onto_probability_simplex(v: &[f64]) -> Vec<f64> {
    let k = v.len();
    if k == 0 {
        return Vec::new();
    }
    if k == 1 {
        return vec![1.0];
    }
    let mut u = v.to_vec();
    u.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let mut cssv = 0.0;
    let mut rho = 0_usize;
    for (j, &uj) in u.iter().enumerate() {
        cssv += uj;
        let cond = uj - (cssv - 1.0) / ((j + 1) as f64);
        if cond > 0.0 {
            rho = j;
        }
    }
    let sum_rho: f64 = u[..=rho].iter().sum();
    let theta = (sum_rho - 1.0) / ((rho + 1) as f64);
    v.iter().map(|&vi| (vi - theta).max(0.0)).collect()
}

fn site_weight(depth: u32, mode: DepthWeighting) -> f64 {
    let d = f64::from(depth);
    match mode {
        DepthWeighting::Log2DepthPlusOne => (d + 1.0).log2().max(0.1),
        DepthWeighting::SqrtDepth => d.sqrt().max(0.1),
        DepthWeighting::Uniform => 1.0,
    }
}

fn eval_weighted_l1(a_mat: &[Vec<f64>], b: &[f64], w: &[f64], p: &[f64]) -> f64 {
    let mut loss = 0.0;
    for m in 0..a_mat.len() {
        let pred: f64 = a_mat[m].iter().zip(p.iter()).map(|(&amk, &pk)| amk * pk).sum();
        loss += w[m] * (pred - b[m]).abs();
    }
    loss
}

/// Solve min_{p in Delta^K} sum_m w_m |(A p)_m - b_m| via continuation IRLS + Nesterov projected gradient + L1 subgradient polish.
pub fn solve_weighted_l1_simplex(a_mat: &[Vec<f64>], b: &[f64], w: &[f64], k_dim: usize) -> Vec<f64> {
    if k_dim == 0 {
        return Vec::new();
    }
    if k_dim == 1 || a_mat.is_empty() {
        let mut p = vec![0.0; k_dim];
        p[0] = 1.0;
        return p;
    }

    let m_rows = a_mat.len();
    let mut p = vec![1.0 / (k_dim as f64); k_dim];
    let mut irls_v = vec![1.0; m_rows];
    let mut grad = vec![0.0; k_dim];

    // Continuation Huber-IRLS stages tightening delta from 2e-2 down to 5e-5
    for &delta in &[2e-2_f64, 3e-3, 4e-4, 5e-5] {
        for _irls_outer in 0..6 {
            let mut lipschitz = 1e-6_f64;
            for m in 0..m_rows {
                let pred: f64 = a_mat[m].iter().zip(p.iter()).map(|(&amk, &pk)| amk * pk).sum();
                let abs_r = (pred - b[m]).abs().max(delta);
                irls_v[m] = w[m] / abs_r;
                let row_norm2: f64 = a_mat[m].iter().map(|&x| x * x).sum();
                lipschitz += irls_v[m] * row_norm2;
            }
            let step_alpha = 0.95 / lipschitz;

            // Nesterov FISTA inner projected gradient on smooth weighted quadratic
            let mut y = p.clone();
            let mut t_nesterov = 1.0_f64;
            for _inner in 0..35 {
                grad.fill(0.0);
                for m in 0..m_rows {
                    let pred: f64 =
                        a_mat[m].iter().zip(y.iter()).map(|(&amk, &yk)| amk * yk).sum();
                    let diff = irls_v[m] * (pred - b[m]);
                    for k in 0..k_dim {
                        grad[k] += diff * a_mat[m][k];
                    }
                }
                let candidate: Vec<f64> = (0..k_dim).map(|k| y[k] - step_alpha * grad[k]).collect();
                let p_next = project_onto_probability_simplex(&candidate);
                let t_next = 0.5 * (1.0 + (1.0 + 4.0 * t_nesterov * t_nesterov).sqrt());
                let beta = (t_nesterov - 1.0) / t_next;
                for k in 0..k_dim {
                    y[k] = p_next[k] + beta * (p_next[k] - p[k]);
                }
                p = p_next;
                t_nesterov = t_next;
            }
        }
    }

    // Exact L1 subgradient + simplex projection polish keeping strictly best L1 point
    let mut best_p = p.clone();
    let mut best_l1 = eval_weighted_l1(a_mat, b, w, &p);
    let w_sum = w.iter().sum::<f64>().max(1.0);
    for it in 1..=80 {
        grad.fill(0.0);
        for m in 0..m_rows {
            let pred: f64 = a_mat[m].iter().zip(p.iter()).map(|(&amk, &pk)| amk * pk).sum();
            let r = pred - b[m];
            let sgn = if r > 1e-12 {
                1.0
            } else if r < -1e-12 {
                -1.0
            } else {
                0.0
            };
            for k in 0..k_dim {
                grad[k] += w[m] * sgn * a_mat[m][k];
            }
        }
        let step = 0.08 / ((it as f64).sqrt() * w_sum);
        let cand: Vec<f64> = (0..k_dim).map(|k| p[k] - step * grad[k]).collect();
        p = project_onto_probability_simplex(&cand);
        let cur_l1 = eval_weighted_l1(a_mat, b, w, &p);
        if cur_l1 < best_l1 {
            best_l1 = cur_l1;
            best_p = p.clone();
        }
    }

    let mut sum_clean = 0.0;
    for val in &mut best_p {
        if *val < 1e-5 {
            *val = 0.0;
        }
        sum_clean += *val;
    }
    if sum_clean > 0.0 {
        for val in &mut best_p {
            *val /= sum_clean;
        }
    }
    best_p
}

fn splitmix64_next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Deconvolve a single wastewater sequencing sample against a reference SNV barcode library.
pub fn deconvolve_wastewater_sample(
    library: &BarcodeLibrary,
    sample: &WastewaterSampleInput,
    cfg: &WastewaterDemixConfig,
) -> Result<SampleDemixResult, EpiError> {
    let k_dim = library.lineage_names.len();
    let m_lib = library.mutation_sites.len();
    if k_dim == 0 || m_lib == 0 || library.barcode_matrix.len() != m_lib {
        return Err(EpiError::Schema(
            "Invalid wastewater barcode library dimensions".to_string(),
        ));
    }

    let site_lookup: BTreeMap<&str, usize> = library
        .mutation_sites
        .iter()
        .enumerate()
        .map(|(idx, s)| (s.as_str(), idx))
        .collect();

    let mut a_sub: Vec<Vec<f64>> = Vec::new();
    let mut b_sub: Vec<f64> = Vec::new();
    let mut w_sub: Vec<f64> = Vec::new();
    let mut d_sub: Vec<u32> = Vec::new();
    let mut name_sub: Vec<String> = Vec::new();

    for obs in &sample.observations {
        if obs.read_depth < cfg.min_site_depth {
            continue;
        }
        if let Some(&m_idx) = site_lookup.get(obs.mutation_site.trim()) {
            a_sub.push(library.barcode_matrix[m_idx].clone());
            b_sub.push(obs.alt_allele_freq.clamp(0.0, 1.0));
            w_sub.push(site_weight(obs.read_depth, cfg.depth_weighting));
            d_sub.push(obs.read_depth);
            name_sub.push(obs.mutation_site.trim().to_string());
        }
    }

    if a_sub.is_empty() {
        return Err(EpiError::Schema(format!(
            "Sample '{}' has zero SNV sites passing min_site_depth={}",
            sample.sample_id, cfg.min_site_depth
        )));
    }

    let p_hat = solve_weighted_l1_simplex(&a_sub, &b_sub, &w_sub, k_dim);
    let w_total: f64 = w_sub.iter().sum::<f64>().max(1e-12);
    let mut l1_sum = 0.0;
    let mut l2_sum = 0.0;
    let mut cryptic_sites = Vec::new();

    for m in 0..a_sub.len() {
        let recon: f64 = a_sub[m]
            .iter()
            .zip(p_hat.iter())
            .map(|(&amk, &pk)| amk * pk)
            .sum();
        let diff = b_sub[m] - recon;
        l1_sum += w_sub[m] * diff.abs();
        l2_sum += w_sub[m] * diff * diff;
        if diff >= cfg.cryptic_residual_cutoff && d_sub[m] >= cfg.min_site_depth.max(30) {
            cryptic_sites.push(CrypticResidualSite {
                mutation_site: name_sub[m].clone(),
                observed_af: b_sub[m],
                reconstructed_af: recon,
                residual_excess: diff,
                read_depth: d_sub[m],
            });
        }
    }
    cryptic_sites.sort_by(|x, y| {
        y.residual_excess
            .partial_cmp(&x.residual_excess)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let reps = cfg.bootstrap_replicates.min(400);
    let mut boot_draws: Vec<Vec<f64>> = vec![Vec::with_capacity(reps); k_dim];
    if reps >= 10 && a_sub.len() >= 4 {
        let mut rng_state = 0xD1B5_4A32_D192_ED03_u64
            ^ (sample.sample_id.bytes().fold(0_u64, |acc, b| {
                acc.wrapping_mul(131).wrapping_add(u64::from(b))
            }));
        let m_used = a_sub.len();
        let mut w_boot = vec![0.0_f64; m_used];
        for _ in 0..reps {
            w_boot.fill(0.0);
            for _ in 0..m_used {
                let pick = (splitmix64_next(&mut rng_state) as usize) % m_used;
                w_boot[pick] += w_sub[pick];
            }
            let p_b = solve_weighted_l1_simplex(&a_sub, &b_sub, &w_boot, k_dim);
            for k in 0..k_dim {
                boot_draws[k].push(p_b[k]);
            }
        }
    }

    let mut lineages = Vec::with_capacity(k_dim);
    for k in 0..k_dim {
        let est = p_hat[k];
        let (lo, hi) = if boot_draws[k].len() >= 10 {
            boot_draws[k].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let n_b = boot_draws[k].len();
            let lo_idx = ((0.025 * ((n_b - 1) as f64)).floor() as usize).min(n_b - 1);
            let hi_idx = ((0.975 * ((n_b - 1) as f64)).ceil() as usize).min(n_b - 1);
            (boot_draws[k][lo_idx].min(est), boot_draws[k][hi_idx].max(est))
        } else {
            (est, est)
        };
        lineages.push(LineageAbundanceEstimate {
            lineage: library.lineage_names[k].clone(),
            abundance: est,
            ci_low_95: lo,
            ci_high_95: hi,
            viral_load_weighted_copies: est * sample.viral_copies_per_liter.max(0.0),
        });
    }

    Ok(SampleDemixResult {
        sample_id: sample.sample_id.clone(),
        date_label: sample.date_label.clone(),
        time_weeks: sample.time_weeks,
        sewershed: sample.sewershed.clone(),
        viral_copies_per_liter: sample.viral_copies_per_liter,
        sites_used: a_sub.len(),
        sites_total_library: m_lib,
        coverage_breadth: (a_sub.len() as f64) / (m_lib as f64),
        weighted_l1_loss: l1_sum / w_total,
        weighted_rmse: (l2_sum / w_total).sqrt(),
        lineages,
        cryptic_sites,
    })
}

/// Deconvolve a longitudinal series of wastewater samples in chronological order.
pub fn deconvolve_wastewater_timeseries(
    library: &BarcodeLibrary,
    samples: &[WastewaterSampleInput],
    cfg: &WastewaterDemixConfig,
) -> Result<LongitudinalWastewaterResult, EpiError> {
    if samples.is_empty() {
        return Err(EpiError::Schema(
            "Wastewater sample series cannot be empty".to_string(),
        ));
    }
    let mut ordered = samples.to_vec();
    ordered.sort_by(|a, b| {
        a.time_weeks
            .partial_cmp(&b.time_weeks)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut results = Vec::with_capacity(ordered.len());
    for s in &ordered {
        results.push(deconvolve_wastewater_sample(library, s, cfg)?);
    }
    Ok(LongitudinalWastewaterResult {
        lineage_names: library.lineage_names.clone(),
        samples: results,
    })
}

/// Parse combined wastewater SNV CSV (`sample_id,date_label,time_weeks,sewershed,viral_copies_l,mutation_site,alt_af,read_depth,barcodes...`).
pub fn parse_wastewater_snv_csv(
    csv_text: &str,
) -> Result<(BarcodeLibrary, Vec<WastewaterSampleInput>), EpiError> {
    let mut lines = csv_text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let header = lines
        .next()
        .ok_or_else(|| EpiError::Schema("Wastewater SNV CSV is empty".to_string()))?;
    let h_cols: Vec<&str> = header.split(',').map(str::trim).collect();
    if h_cols.len() < 9 {
        return Err(EpiError::Schema(
            "Wastewater SNV CSV requires 8 sample/site columns followed by >=1 lineage barcode columns (e.g. bc_JN.1,bc_KP.2,bc_KP.3,bc_XEC)"
                .to_string(),
        ));
    }

    let lineage_names: Vec<String> = h_cols[8..]
        .iter()
        .map(|s| s.strip_prefix("bc_").unwrap_or(s).to_string())
        .collect();
    let k_dim = lineage_names.len();

    let mut site_barcodes: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut sample_map: BTreeMap<String, WastewaterSampleInput> = BTreeMap::new();
    let mut site_order: Vec<String> = Vec::new();
    let mut seen_sites: BTreeSet<String> = BTreeSet::new();

    for (idx, line) in lines.enumerate() {
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        if cols.len() < 8 + k_dim {
            return Err(EpiError::Schema(format!(
                "Wastewater CSV data row {} has {} columns, expected {}",
                idx + 2,
                cols.len(),
                8 + k_dim
            )));
        }
        let sid = cols[0].to_string();
        let date_label = cols[1].to_string();
        let time_weeks: f64 = cols[2]
            .parse()
            .map_err(|_| EpiError::Schema(format!("Invalid time_weeks '{}'", cols[2])))?;
        let sewershed = cols[3].to_string();
        let copies: f64 = cols[4].parse().unwrap_or(1.0e5);
        let site = cols[5].to_string();
        let alt_af: f64 = cols[6]
            .parse()
            .map_err(|_| EpiError::Schema(format!("Invalid alt_af '{}'", cols[6])))?;
        let depth: u32 = cols[7].parse().unwrap_or(100);

        if seen_sites.insert(site.clone()) {
            let mut bc_row = Vec::with_capacity(k_dim);
            for k in 0..k_dim {
                let v: f64 = cols[8 + k].parse().unwrap_or(0.0);
                bc_row.push(v.clamp(0.0, 1.0));
            }
            site_barcodes.insert(site.clone(), bc_row);
            site_order.push(site.clone());
        }

        let entry = sample_map
            .entry(sid.clone())
            .or_insert_with(|| WastewaterSampleInput {
                sample_id: sid,
                date_label,
                time_weeks,
                sewershed,
                viral_copies_per_liter: copies,
                observations: Vec::new(),
            });
        entry.observations.push(SnvObservation {
            mutation_site: site,
            alt_allele_freq: alt_af,
            read_depth: depth,
        });
    }

    let barcode_matrix: Vec<Vec<f64>> = site_order
        .iter()
        .map(|s| site_barcodes.get(s).cloned().unwrap_or_else(|| vec![0.0; k_dim]))
        .collect();

    let library = BarcodeLibrary {
        lineage_names,
        mutation_sites: site_order,
        barcode_matrix,
    };
    let samples: Vec<WastewaterSampleInput> = sample_map.into_values().collect();
    Ok((library, samples))
}

/// Built-in 10-week Willamette Valley sewershed SARS-CoV-2 tiled-amplicon SNV barcode & allele-frequency series.
pub fn builtin_wastewater_barcode_and_samples() -> (BarcodeLibrary, Vec<WastewaterSampleInput>) {
    let lineage_names = vec![
        "JN.1".to_string(),
        "KP.2".to_string(),
        "KP.3.1.1".to_string(),
        "XEC".to_string(),
    ];

    let site_defs: &[(&str, [f64; 4])] = &[
        ("S:L455S", [1.0, 1.0, 1.0, 1.0]),
        ("ORF1a:A2529V", [1.0, 0.0, 0.0, 0.0]),
        ("N:Q229K", [1.0, 0.05, 0.0, 0.0]),
        ("S:R346T", [0.0, 1.0, 0.0, 0.0]),
        ("S:F456L", [0.0, 1.0, 1.0, 1.0]),
        ("S:V1104L", [0.0, 1.0, 1.0, 0.0]),
        ("ORF1b:P314L_sub", [0.0, 0.95, 0.0, 0.0]),
        ("S:Q493E", [0.0, 0.0, 1.0, 1.0]),
        ("S:S31del", [0.0, 0.0, 1.0, 0.0]),
        ("ORF1a:T170I", [0.0, 0.0, 1.0, 0.0]),
        ("S:T22N", [0.0, 0.0, 0.0, 1.0]),
        ("S:F59S", [0.0, 0.0, 0.0, 1.0]),
        ("N:G204P_rec", [0.0, 0.0, 0.0, 1.0]),
        ("ORF3a:S253P", [0.85, 0.15, 0.0, 0.0]),
    ];

    let mutation_sites: Vec<String> = site_defs.iter().map(|(s, _)| (*s).to_string()).collect();
    let barcode_matrix: Vec<Vec<f64>> = site_defs.iter().map(|(_, r)| r.to_vec()).collect();
    let library = BarcodeLibrary {
        lineage_names,
        mutation_sites: mutation_sites.clone(),
        barcode_matrix: barcode_matrix.clone(),
    };

    let true_trajectories: &[(&str, f64, [f64; 4])] = &[
        ("2026-W35", 1.75e5, [0.56, 0.30, 0.11, 0.03]),
        ("2026-W36", 2.10e5, [0.46, 0.33, 0.16, 0.05]),
        ("2026-W37", 2.50e5, [0.36, 0.34, 0.22, 0.08]),
        ("2026-W38", 3.10e5, [0.27, 0.32, 0.29, 0.12]),
        ("2026-W39", 3.80e5, [0.19, 0.28, 0.36, 0.17]),
        ("2026-W40", 4.60e5, [0.13, 0.23, 0.41, 0.23]),
        ("2026-W41", 5.20e5, [0.09, 0.18, 0.43, 0.30]),
        ("2026-W42", 5.90e5, [0.06, 0.14, 0.42, 0.38]),
        ("2026-W43", 6.40e5, [0.04, 0.10, 0.40, 0.46]),
        ("2026-W44", 7.10e5, [0.03, 0.07, 0.36, 0.54]),
    ];

    let mut samples = Vec::with_capacity(true_trajectories.len());
    for (t_idx, &(wk_label, copies, shares)) in true_trajectories.iter().enumerate() {
        let mut obs = Vec::with_capacity(mutation_sites.len());
        for (m_idx, site_name) in mutation_sites.iter().enumerate() {
            let exact_af: f64 = (0..4).map(|k| barcode_matrix[m_idx][k] * shares[k]).sum();
            let jitter = (((t_idx * 17 + m_idx * 13) % 9) as f64 - 4.0) * 0.0035;
            let cryptic_bump = if t_idx == 9 && site_name == "S:S31del" {
                0.115
            } else {
                0.0
            };
            let depth = 180 + (((t_idx * 97 + m_idx * 53) % 920) as u32);
            obs.push(SnvObservation {
                mutation_site: site_name.clone(),
                alt_allele_freq: (exact_af + jitter + cryptic_bump).clamp(0.0, 1.0),
                read_depth: depth,
            });
        }
        samples.push(WastewaterSampleInput {
            sample_id: format!("WW-OR-CORV-{:02}", t_idx + 1),
            date_label: wk_label.to_string(),
            time_weeks: t_idx as f64,
            sewershed: "Corvallis-Willamette-WRF".to_string(),
            viral_copies_per_liter: copies,
            observations: obs,
        });
    }

    (library, samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simplex_projection_and_l1_wastewater_recovery() {
        let proj = project_onto_probability_simplex(&[0.7, -0.4, 0.5, 0.2]);
        let sum: f64 = proj.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
        assert!(proj.iter().all(|&p| p >= 0.0));

        let (lib, samples) = builtin_wastewater_barcode_and_samples();
        let series =
            deconvolve_wastewater_timeseries(&lib, &samples, &WastewaterDemixConfig::default())
                .expect("wastewater deconvolution succeeds");
        assert_eq!(series.samples.len(), 10);

        let first = &series.samples[0];
        let jn1_w0 = first
            .lineages
            .iter()
            .find(|l| l.lineage == "JN.1")
            .unwrap()
            .abundance;
        let xec_w0 = first
            .lineages
            .iter()
            .find(|l| l.lineage == "XEC")
            .unwrap()
            .abundance;
        assert!((jn1_w0 - 0.56).abs() < 0.035, "got JN.1={jn1_w0}");
        assert!((xec_w0 - 0.03).abs() < 0.035, "got XEC={xec_w0}");

        let last = series.samples.last().unwrap();
        let xec_end = last
            .lineages
            .iter()
            .find(|l| l.lineage == "XEC")
            .unwrap()
            .abundance;
        assert!(xec_end > 0.48, "expected XEC takeover > 0.48, got {xec_end}");
        assert!(!last.cryptic_sites.is_empty());

        let tbl = series.to_lineage_count_table(500.0);
        assert_eq!(tbl.clade_names.len(), 4);
        assert_eq!(tbl.rows.len(), 10);
    }
}
