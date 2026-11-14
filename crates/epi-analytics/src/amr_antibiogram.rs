use crate::math::normal_quantile;
use epi_core::EpiError;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// Clinical Laboratory Standards Institute (CLSI) categorical susceptibility interpretation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SirInterpretation {
    Susceptible,
    SusceptibleDoseDependent,
    Intermediate,
    Resistant,
}

impl SirInterpretation {
    pub fn from_code(s: &str) -> Option<Self> {
        match s.trim().to_ascii_uppercase().as_str() {
            "S" | "SUSCEPTIBLE" => Some(Self::Susceptible),
            "SDD" => Some(Self::SusceptibleDoseDependent),
            "I" | "INTERMEDIATE" => Some(Self::Intermediate),
            "R" | "RESISTANT" => Some(Self::Resistant),
            _ => None,
        }
    }

    pub fn is_non_susceptible(self) -> bool {
        matches!(self, Self::Intermediate | Self::Resistant)
    }
}

/// Isolate inclusion policy per CLSI M39 cumulative antibiogram guidelines.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IsolateDeduplicationPolicy {
    /// CLSI M39 standard: retain strictly the chronologically first diagnostic isolate per patient per species.
    ClsiM39FirstIsolatePerPatient,
    /// Unfiltered raw isolates (shown for diagnostic comparison to quantify repeat-culture bias).
    AllIsolatesUnfiltered,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AntimicrobialTestCall {
    pub drug_code: String,
    pub drug_name: String,
    pub drug_class: String,
    pub mic_ug_ml: Option<f64>,
    pub sir: SirInterpretation,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AmrIsolateRecord {
    pub isolate_id: String,
    pub patient_id: String,
    pub collection_day: u32,
    pub facility_unit: String,
    pub specimen_source: String,
    pub organism: String,
    pub susceptibility_calls: Vec<AntimicrobialTestCall>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AntibiogramConfig {
    pub deduplication: IsolateDeduplicationPolicy,
    /// Optional facility unit filter (e.g., "ICU" or empty for facility-wide).
    pub facility_unit_filter: String,
    /// Minimum isolates required by CLSI M39 for statistically reliable %S reporting (standard = 30).
    pub min_isolates_m39: usize,
    /// Empiric stewardship %S threshold below which an alert flag is raised (default = 80.0%).
    pub empiric_stewardship_threshold_pct: f64,
}

impl Default for AntibiogramConfig {
    fn default() -> Self {
        Self {
            deduplication: IsolateDeduplicationPolicy::ClsiM39FirstIsolatePerPatient,
            facility_unit_filter: String::new(),
            min_isolates_m39: 30,
            empiric_stewardship_threshold_pct: 80.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AntibiogramCell {
    pub organism: String,
    pub drug_code: String,
    pub drug_name: String,
    pub drug_class: String,
    pub isolates_tested: usize,
    pub susceptible_count: usize,
    pub sdd_count: usize,
    pub intermediate_count: usize,
    pub resistant_count: usize,
    /// Percent susceptible (%S) strictly = S / N_tested * 100 per CLSI M39.
    pub pct_susceptible: f64,
    pub wilson_ci_low_95: f64,
    pub wilson_ci_high_95: f64,
    pub mic50_ug_ml: Option<f64>,
    pub mic90_ug_ml: Option<f64>,
    pub meets_clsi_m39_n30: bool,
    pub empiric_caution_below_threshold: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganismResistanceSummary {
    pub organism: String,
    pub raw_isolates: usize,
    pub deduplicated_isolates: usize,
    pub duplicate_isolates_excluded: usize,
    pub meets_clsi_m39_n30: bool,
    /// Magiorakos consensus MDR: non-susceptible (I or R) to >=1 agent in >=3 antimicrobial classes.
    pub mdr_isolates: usize,
    pub pct_mdr: f64,
    /// Magiorakos consensus XDR: non-susceptible in all but <=2 tested antimicrobial classes (when >=4 classes tested).
    pub xdr_isolates: usize,
    pub pct_xdr: f64,
    /// Key sentinel phenotype count (MRSA for S. aureus, Carbapenem-R / CRE or 3GC-R for Gram-negatives).
    pub sentinel_phenotype_label: String,
    pub sentinel_resistant_count: usize,
    pub pct_sentinel_resistant: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CumulativeAntibiogramResult {
    pub policy: IsolateDeduplicationPolicy,
    pub facility_unit_filter: String,
    pub total_raw_isolates_matched: usize,
    pub total_included_isolates: usize,
    pub total_duplicates_excluded: usize,
    pub organisms: Vec<String>,
    pub drugs: Vec<(String, String, String)>,
    pub organism_summaries: Vec<OrganismResistanceSummary>,
    pub cells: Vec<AntibiogramCell>,
}

/// Wilson score 95% binomial confidence interval returned in percentage points [0, 100].
pub fn wilson_95_ci_pct(successes: usize, total: usize) -> (f64, f64) {
    if total == 0 {
        return (0.0, 0.0);
    }
    let n = total as f64;
    let p = (successes as f64) / n;
    let z = normal_quantile(0.975);
    let z2 = z * z;
    let denom = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denom;
    let half = (z * ((p * (1.0 - p) / n) + (z2 / (4.0 * n * n))).sqrt()) / denom;
    (
        ((center - half).clamp(0.0, 1.0)) * 100.0,
        ((center + half).clamp(0.0, 1.0)) * 100.0,
    )
}

fn mic_percentile(mics: &mut [f64], q: f64) -> Option<f64> {
    if mics.is_empty() {
        return None;
    }
    mics.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx = ((q.clamp(0.0, 1.0) * (mics.len() as f64)).ceil() as usize)
        .saturating_sub(1)
        .min(mics.len() - 1);
    Some(mics[idx])
}

/// Compute a CLSI M39 cumulative antibiogram, MIC50/MIC90 quantiles, Wilson 95% CIs, and Magiorakos MDR/XDR phenotypes.
pub fn compute_clsi_m39_antibiogram(
    records: &[AmrIsolateRecord],
    cfg: &AntibiogramConfig,
) -> Result<CumulativeAntibiogramResult, EpiError> {
    let unit_flt = cfg.facility_unit_filter.trim();
    let mut filtered: Vec<AmrIsolateRecord> = records
        .iter()
        .filter(|r| {
            unit_flt.is_empty()
                || unit_flt.eq_ignore_ascii_case("ALL")
                || r.facility_unit.eq_ignore_ascii_case(unit_flt)
        })
        .cloned()
        .collect();

    if filtered.is_empty() {
        return Err(EpiError::Schema(format!(
            "No clinical microbiology isolates matched facility unit filter '{}'",
            cfg.facility_unit_filter
        )));
    }

    // Sort chronologically so CLSI M39 retains strictly the earliest diagnostic isolate per (patient_id, organism)
    filtered.sort_by(|a, b| {
        a.collection_day
            .cmp(&b.collection_day)
            .then_with(|| a.isolate_id.cmp(&b.isolate_id))
    });

    let total_raw = filtered.len();
    let mut raw_per_org: BTreeMap<String, usize> = BTreeMap::new();
    for r in &filtered {
        *raw_per_org.entry(r.organism.clone()).or_insert(0) += 1;
    }

    let included: Vec<AmrIsolateRecord> = match cfg.deduplication {
        IsolateDeduplicationPolicy::ClsiM39FirstIsolatePerPatient => {
            let mut seen_pt_org: BTreeSet<(String, String)> = BTreeSet::new();
            let mut out = Vec::with_capacity(filtered.len());
            for rec in filtered {
                let key = (
                    rec.patient_id.trim().to_ascii_uppercase(),
                    rec.organism.trim().to_string(),
                );
                if seen_pt_org.insert(key) {
                    out.push(rec);
                }
            }
            out
        }
        IsolateDeduplicationPolicy::AllIsolatesUnfiltered => filtered,
    };

    let total_included = included.len();
    let total_dups = total_raw.saturating_sub(total_included);

    let mut org_set: BTreeSet<String> = BTreeSet::new();
    let mut drug_meta: BTreeMap<String, (String, String)> = BTreeMap::new();
    for rec in &included {
        org_set.insert(rec.organism.clone());
        for call in &rec.susceptibility_calls {
            drug_meta
                .entry(call.drug_code.clone())
                .or_insert_with(|| (call.drug_name.clone(), call.drug_class.clone()));
        }
    }

    let organisms: Vec<String> = org_set.into_iter().collect();
    let drugs: Vec<(String, String, String)> = drug_meta
        .into_iter()
        .map(|(code, (name, cls))| (code, name, cls))
        .collect();

    // Per-organism MDR / XDR / Sentinel summaries
    let mut organism_summaries = Vec::with_capacity(organisms.len());
    for org in &organisms {
        let org_isos: Vec<&AmrIsolateRecord> =
            included.iter().filter(|r| &r.organism == org).collect();
        let n_inc = org_isos.len();
        let n_raw = raw_per_org.get(org).copied().unwrap_or(n_inc);

        let mut mdr_c = 0_usize;
        let mut xdr_c = 0_usize;
        let mut sent_c = 0_usize;
        let is_staph = org.to_ascii_lowercase().contains("staphylococcus");
        let sentinel_label = if is_staph {
            "MRSA (Oxacillin-R)".to_string()
        } else if org.to_ascii_lowercase().contains("pseudomonas") {
            "Carbapenem-R P. aeruginosa".to_string()
        } else {
            "3GC / ESBL or CRE".to_string()
        };

        for iso in &org_isos {
            let mut tested_classes: BTreeSet<&str> = BTreeSet::new();
            let mut nonsusc_classes: BTreeSet<&str> = BTreeSet::new();
            let mut hit_sentinel = false;

            for c in &iso.susceptibility_calls {
                tested_classes.insert(c.drug_class.as_str());
                if c.sir.is_non_susceptible() {
                    nonsusc_classes.insert(c.drug_class.as_str());
                }
                if is_staph && c.drug_code == "OXA" && c.sir == SirInterpretation::Resistant {
                    hit_sentinel = true;
                } else if !is_staph
                    && (c.drug_code == "MEM" || c.drug_code == "CRO")
                    && c.sir == SirInterpretation::Resistant
                {
                    hit_sentinel = true;
                }
            }

            if nonsusc_classes.len() >= 3 {
                mdr_c += 1;
            }
            if tested_classes.len() >= 4
                && tested_classes.len().saturating_sub(nonsusc_classes.len()) <= 2
            {
                xdr_c += 1;
            }
            if hit_sentinel {
                sent_c += 1;
            }
        }

        let denom = (n_inc.max(1)) as f64;
        organism_summaries.push(OrganismResistanceSummary {
            organism: org.clone(),
            raw_isolates: n_raw,
            deduplicated_isolates: n_inc,
            duplicate_isolates_excluded: n_raw.saturating_sub(n_inc),
            meets_clsi_m39_n30: n_inc >= cfg.min_isolates_m39,
            mdr_isolates: mdr_c,
            pct_mdr: (mdr_c as f64) * 100.0 / denom,
            xdr_isolates: xdr_c,
            pct_xdr: (xdr_c as f64) * 100.0 / denom,
            sentinel_phenotype_label: sentinel_label,
            sentinel_resistant_count: sent_c,
            pct_sentinel_resistant: (sent_c as f64) * 100.0 / denom,
        });
    }

    // Build Organism x Antimicrobial matrix cells
    let mut cells = Vec::new();
    for org in &organisms {
        for (d_code, d_name, d_class) in &drugs {
            let mut s_cnt = 0_usize;
            let mut sdd_cnt = 0_usize;
            let mut i_cnt = 0_usize;
            let mut r_cnt = 0_usize;
            let mut mics: Vec<f64> = Vec::new();

            for iso in included.iter().filter(|r| &r.organism == org) {
                if let Some(call) = iso
                    .susceptibility_calls
                    .iter()
                    .find(|c| &c.drug_code == d_code)
                {
                    match call.sir {
                        SirInterpretation::Susceptible => s_cnt += 1,
                        SirInterpretation::SusceptibleDoseDependent => sdd_cnt += 1,
                        SirInterpretation::Intermediate => i_cnt += 1,
                        SirInterpretation::Resistant => r_cnt += 1,
                    }
                    if let Some(m) = call.mic_ug_ml {
                        if m.is_finite() && m > 0.0 {
                            mics.push(m);
                        }
                    }
                }
            }

            let n_tested = s_cnt + sdd_cnt + i_cnt + r_cnt;
            if n_tested == 0 {
                continue;
            }
            let pct_s = (s_cnt as f64) * 100.0 / (n_tested as f64);
            let (ci_lo, ci_hi) = wilson_95_ci_pct(s_cnt, n_tested);
            let mic50 = mic_percentile(&mut mics, 0.50);
            let mic90 = mic_percentile(&mut mics, 0.90);

            cells.push(AntibiogramCell {
                organism: org.clone(),
                drug_code: d_code.clone(),
                drug_name: d_name.clone(),
                drug_class: d_class.clone(),
                isolates_tested: n_tested,
                susceptible_count: s_cnt,
                sdd_count: sdd_cnt,
                intermediate_count: i_cnt,
                resistant_count: r_cnt,
                pct_susceptible: pct_s,
                wilson_ci_low_95: ci_lo,
                wilson_ci_high_95: ci_hi,
                mic50_ug_ml: mic50,
                mic90_ug_ml: mic90,
                meets_clsi_m39_n30: n_tested >= cfg.min_isolates_m39,
                empiric_caution_below_threshold: pct_s < cfg.empiric_stewardship_threshold_pct,
            });
        }
    }

    Ok(CumulativeAntibiogramResult {
        policy: cfg.deduplication,
        facility_unit_filter: if unit_flt.is_empty() {
            "ALL_UNITS".to_string()
        } else {
            unit_flt.to_string()
        },
        total_raw_isolates_matched: total_raw,
        total_included_isolates: total_included,
        total_duplicates_excluded: total_dups,
        organisms,
        drugs,
        organism_summaries,
        cells,
    })
}

/// Parse long-format clinical susceptibility CSV (`isolate_id,patient_id,collection_day,facility_unit,specimen_source,organism,drug_code,drug_name,drug_class,mic_ug_ml,sir`).
pub fn parse_amr_isolates_csv(csv_text: &str) -> Result<Vec<AmrIsolateRecord>, EpiError> {
    let mut iso_map: BTreeMap<String, AmrIsolateRecord> = BTreeMap::new();
    for (idx, raw_line) in csv_text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        if idx == 0 && cols[0].eq_ignore_ascii_case("isolate_id") {
            continue;
        }
        if cols.len() < 11 {
            return Err(EpiError::Schema(format!(
                "AMR CSV line {} requires 11 columns (isolate_id..sir)",
                idx + 1
            )));
        }
        let sir = SirInterpretation::from_code(cols[10]).ok_or_else(|| {
            EpiError::Schema(format!("Unrecognized SIR code '{}' on line {}", cols[10], idx + 1))
        })?;
        let day: u32 = cols[2].parse().unwrap_or(1);
        let mic = cols[9].parse::<f64>().ok();
        let iso_id = cols[0].to_string();

        let entry = iso_map
            .entry(iso_id.clone())
            .or_insert_with(|| AmrIsolateRecord {
                isolate_id: iso_id,
                patient_id: cols[1].to_string(),
                collection_day: day,
                facility_unit: cols[3].to_string(),
                specimen_source: cols[4].to_string(),
                organism: cols[5].to_string(),
                susceptibility_calls: Vec::new(),
            });
        entry.susceptibility_calls.push(AntimicrobialTestCall {
            drug_code: cols[6].to_string(),
            drug_name: cols[7].to_string(),
            drug_class: cols[8].to_string(),
            mic_ug_ml: mic,
            sir,
        });
    }
    if iso_map.is_empty() {
        return Err(EpiError::Schema(
            "No valid AMR isolate rows found in CSV".to_string(),
        ));
    }
    Ok(iso_map.into_values().collect())
}

/// Built-in 196-isolate regional hospital microbiology panel across 4 priority pathogens with realistic ICU repeat-culture bias.
pub fn builtin_hospital_amr_isolates() -> Vec<AmrIsolateRecord> {
    let mut records = Vec::new();
    let units = ["Inpatient Med-Surg", "ICU", "Emergency", "Outpatient Clinic"];

    // 1. Escherichia coli: 52 unique patients (>30 CLSI M39 reliable) + 12 repeat resistant ICU surveillance cultures
    for p in 1..=52_u32 {
        let unit = units[(p as usize) % 4];
        let is_esbl = p % 7 == 0;
        let is_fq_r = p % 4 == 0 || is_esbl;
        let is_sxt_r = p % 3 == 0;
        let calls = vec![
            AntimicrobialTestCall {
                drug_code: "CRO".to_string(),
                drug_name: "Ceftriaxone".to_string(),
                drug_class: "Cephalosporins".to_string(),
                mic_ug_ml: Some(if is_esbl { 32.0 } else { 0.25 }),
                sir: if is_esbl {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "TZP".to_string(),
                drug_name: "Piperacillin-Tazobactam".to_string(),
                drug_class: "Beta-Lactam-Inhibitors".to_string(),
                mic_ug_ml: Some(if p % 17 == 0 { 64.0 } else { 4.0 }),
                sir: if p % 17 == 0 {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "MEM".to_string(),
                drug_name: "Meropenem".to_string(),
                drug_class: "Carbapenems".to_string(),
                mic_ug_ml: Some(0.125),
                sir: SirInterpretation::Susceptible,
            },
            AntimicrobialTestCall {
                drug_code: "CIP".to_string(),
                drug_name: "Ciprofloxacin".to_string(),
                drug_class: "Fluoroquinolones".to_string(),
                mic_ug_ml: Some(if is_fq_r { 4.0 } else { 0.06 }),
                sir: if is_fq_r {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "GEN".to_string(),
                drug_name: "Gentamicin".to_string(),
                drug_class: "Aminoglycosides".to_string(),
                mic_ug_ml: Some(if p % 11 == 0 { 16.0 } else { 1.0 }),
                sir: if p % 11 == 0 {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "SXT".to_string(),
                drug_name: "TMP-SMX".to_string(),
                drug_class: "Folate-Pathway".to_string(),
                mic_ug_ml: Some(if is_sxt_r { 8.0 } else { 0.5 }),
                sir: if is_sxt_r {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
        ];
        records.push(AmrIsolateRecord {
            isolate_id: format!("ISO-EC-{p:03}"),
            patient_id: format!("PT-EC-{p:03}"),
            collection_day: p * 3,
            facility_unit: unit.to_string(),
            specimen_source: if p % 5 == 0 { "Blood" } else { "Urine" }.to_string(),
            organism: "Escherichia coli".to_string(),
            susceptibility_calls: calls.clone(),
        });
        // Serial repeat cultures on ESBL ICU patients (demonstrates why CLSI M39 deduplication prevents artificial %S drop)
        if is_esbl {
            for rep in 1..=2_u32 {
                records.push(AmrIsolateRecord {
                    isolate_id: format!("ISO-EC-{p:03}-DUP{rep}"),
                    patient_id: format!("PT-EC-{p:03}"),
                    collection_day: p * 3 + rep * 6,
                    facility_unit: "ICU".to_string(),
                    specimen_source: "Blood".to_string(),
                    organism: "Escherichia coli".to_string(),
                    susceptibility_calls: calls.clone(),
                });
            }
        }
    }

    // 2. Klebsiella pneumoniae: 38 unique patients (>30 CLSI M39 reliable) + 8 CRE repeat cultures
    for p in 1..=38_u32 {
        let is_cre = p % 12 == 0;
        let is_esbl = is_cre || p % 5 == 0;
        let calls = vec![
            AntimicrobialTestCall {
                drug_code: "CRO".to_string(),
                drug_name: "Ceftriaxone".to_string(),
                drug_class: "Cephalosporins".to_string(),
                mic_ug_ml: Some(if is_esbl { 64.0 } else { 0.5 }),
                sir: if is_esbl {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "TZP".to_string(),
                drug_name: "Piperacillin-Tazobactam".to_string(),
                drug_class: "Beta-Lactam-Inhibitors".to_string(),
                mic_ug_ml: Some(if is_cre || p % 6 == 0 { 128.0 } else { 8.0 }),
                sir: if is_cre || p % 6 == 0 {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "MEM".to_string(),
                drug_name: "Meropenem".to_string(),
                drug_class: "Carbapenems".to_string(),
                mic_ug_ml: Some(if is_cre { 16.0 } else { 0.25 }),
                sir: if is_cre {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "CIP".to_string(),
                drug_name: "Ciprofloxacin".to_string(),
                drug_class: "Fluoroquinolones".to_string(),
                mic_ug_ml: Some(if is_esbl || p % 4 == 0 { 4.0 } else { 0.125 }),
                sir: if is_esbl || p % 4 == 0 {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "GEN".to_string(),
                drug_name: "Gentamicin".to_string(),
                drug_class: "Aminoglycosides".to_string(),
                mic_ug_ml: Some(if is_cre || p % 9 == 0 { 16.0 } else { 1.0 }),
                sir: if is_cre || p % 9 == 0 {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
        ];
        records.push(AmrIsolateRecord {
            isolate_id: format!("ISO-KP-{p:03}"),
            patient_id: format!("PT-KP-{p:03}"),
            collection_day: p * 4,
            facility_unit: units[(p as usize + 1) % 4].to_string(),
            specimen_source: "Respiratory".to_string(),
            organism: "Klebsiella pneumoniae".to_string(),
            susceptibility_calls: calls.clone(),
        });
        if is_cre {
            for rep in 1..=2_u32 {
                records.push(AmrIsolateRecord {
                    isolate_id: format!("ISO-KP-{p:03}-DUP{rep}"),
                    patient_id: format!("PT-KP-{p:03}"),
                    collection_day: p * 4 + rep * 5,
                    facility_unit: "ICU".to_string(),
                    specimen_source: "Blood".to_string(),
                    organism: "Klebsiella pneumoniae".to_string(),
                    susceptibility_calls: calls.clone(),
                });
            }
        }
    }

    // 3. Pseudomonas aeruginosa: 24 unique patients (<30 CLSI M39 low-N flag triggered intentionally!)
    for p in 1..=24_u32 {
        let is_mdr = p % 5 == 0;
        let calls = vec![
            AntimicrobialTestCall {
                drug_code: "TZP".to_string(),
                drug_name: "Piperacillin-Tazobactam".to_string(),
                drug_class: "Beta-Lactam-Inhibitors".to_string(),
                mic_ug_ml: Some(if is_mdr { 128.0 } else { 8.0 }),
                sir: if is_mdr {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "MEM".to_string(),
                drug_name: "Meropenem".to_string(),
                drug_class: "Carbapenems".to_string(),
                mic_ug_ml: Some(if is_mdr || p % 7 == 0 { 16.0 } else { 1.0 }),
                sir: if is_mdr || p % 7 == 0 {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "CIP".to_string(),
                drug_name: "Ciprofloxacin".to_string(),
                drug_class: "Fluoroquinolones".to_string(),
                mic_ug_ml: Some(if is_mdr || p % 4 == 0 { 4.0 } else { 0.25 }),
                sir: if is_mdr || p % 4 == 0 {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "GEN".to_string(),
                drug_name: "Gentamicin".to_string(),
                drug_class: "Aminoglycosides".to_string(),
                mic_ug_ml: Some(if is_mdr { 16.0 } else { 2.0 }),
                sir: if is_mdr {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
        ];
        records.push(AmrIsolateRecord {
            isolate_id: format!("ISO-PA-{p:03}"),
            patient_id: format!("PT-PA-{p:03}"),
            collection_day: p * 5,
            facility_unit: "ICU".to_string(),
            specimen_source: "Respiratory".to_string(),
            organism: "Pseudomonas aeruginosa".to_string(),
            susceptibility_calls: calls,
        });
    }

    // 4. Staphylococcus aureus: 45 unique patients (>30 reliable; ~31% MRSA, 100% Vancomycin-S) + 7 MRSA repeat cultures
    for p in 1..=45_u32 {
        let is_mrsa = p % 3 == 0 && p <= 42;
        let calls = vec![
            AntimicrobialTestCall {
                drug_code: "OXA".to_string(),
                drug_name: "Oxacillin".to_string(),
                drug_class: "Penicillins".to_string(),
                mic_ug_ml: Some(if is_mrsa { 8.0 } else { 0.5 }),
                sir: if is_mrsa {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "VAN".to_string(),
                drug_name: "Vancomycin".to_string(),
                drug_class: "Glycopeptides".to_string(),
                mic_ug_ml: Some(1.0),
                sir: SirInterpretation::Susceptible,
            },
            AntimicrobialTestCall {
                drug_code: "CLI".to_string(),
                drug_name: "Clindamycin".to_string(),
                drug_class: "Lincosamides".to_string(),
                mic_ug_ml: Some(if is_mrsa || p % 5 == 0 { 4.0 } else { 0.25 }),
                sir: if is_mrsa || p % 5 == 0 {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
            AntimicrobialTestCall {
                drug_code: "SXT".to_string(),
                drug_name: "TMP-SMX".to_string(),
                drug_class: "Folate-Pathway".to_string(),
                mic_ug_ml: Some(if p % 15 == 0 { 8.0 } else { 0.5 }),
                sir: if p % 15 == 0 {
                    SirInterpretation::Resistant
                } else {
                    SirInterpretation::Susceptible
                },
            },
        ];
        records.push(AmrIsolateRecord {
            isolate_id: format!("ISO-SA-{p:03}"),
            patient_id: format!("PT-SA-{p:03}"),
            collection_day: p * 2,
            facility_unit: units[(p as usize + 2) % 4].to_string(),
            specimen_source: "Blood".to_string(),
            organism: "Staphylococcus aureus".to_string(),
            susceptibility_calls: calls.clone(),
        });
        if is_mrsa && p <= 21 {
            records.push(AmrIsolateRecord {
                isolate_id: format!("ISO-SA-{p:03}-DUP1"),
                patient_id: format!("PT-SA-{p:03}"),
                collection_day: p * 2 + 8,
                facility_unit: "ICU".to_string(),
                specimen_source: "Blood".to_string(),
                organism: "Staphylococcus aureus".to_string(),
                susceptibility_calls: calls,
            });
        }
    }

    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clsi_m39_deduplication_corrects_repeat_culture_bias_and_flags_low_n() {
        let recs = builtin_hospital_amr_isolates();
        let m39 = compute_clsi_m39_antibiogram(&recs, &AntibiogramConfig::default())
            .expect("CLSI M39 antibiogram succeeds");
        let raw_cfg = AntibiogramConfig {
            deduplication: IsolateDeduplicationPolicy::AllIsolatesUnfiltered,
            ..Default::default()
        };
        let raw_ab = compute_clsi_m39_antibiogram(&recs, &raw_cfg).expect("unfiltered succeeds");

        assert!(m39.total_duplicates_excluded >= 20);
        let ec_cro_m39 = m39
            .cells
            .iter()
            .find(|c| c.organism == "Escherichia coli" && c.drug_code == "CRO")
            .unwrap();
        let ec_cro_raw = raw_ab
            .cells
            .iter()
            .find(|c| c.organism == "Escherichia coli" && c.drug_code == "CRO")
            .unwrap();
        assert_eq!(ec_cro_m39.isolates_tested, 52);
        assert!(ec_cro_m39.meets_clsi_m39_n30);
        // Deduplicating repeat ICU ESBL cultures produces higher, unbiased %S than unfiltered cultures
        assert!(ec_cro_m39.pct_susceptible > ec_cro_raw.pct_susceptible);

        // P. aeruginosa has 24 unique patients (< 30), so CLSI M39 n>=30 flag must be false
        let pa_summary = m39
            .organism_summaries
            .iter()
            .find(|s| s.organism == "Pseudomonas aeruginosa")
            .unwrap();
        assert_eq!(pa_summary.deduplicated_isolates, 24);
        assert!(!pa_summary.meets_clsi_m39_n30);
    }
}
