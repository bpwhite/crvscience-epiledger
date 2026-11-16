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

use epi_core::{
    parse_flexible_date_to_mmwr, EpiError, SourceKind, TimeSeriesDataset, TimeSeriesPoint,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Deterministic CDC NSSP ESSENCE-style Concatenated Chief-Complaint & Discharge-Diagnosis (`CCDD`)
/// syndrome classification with clinical context-window negation (`NegEx`) and boolean `OR` / `ANDNOT` queries.
///
/// References:
/// - Evaluation and Revision of NSSP ESSENCE Syndromic Definitions: https://pmc.ncbi.nlm.nih.gov/articles/PMC11556616/
/// - Chief-Complaint and ICD-10-CM Exclusion Validation in NSSP: https://pmc.ncbi.nlm.nih.gov/articles/PMC11426074/
/// - All-Hazards Emergency Department Syndromic Surveillance: https://pmc.ncbi.nlm.nih.gov/articles/PMC5676514/

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum EssenceSyndromeId {
    IliRespiratory,
    CovidLikeIllness,
    AcuteGiEnteric,
    SuspectedOverdose,
    HeatRelatedIllness,
    FeverRashExanthem,
}

impl EssenceSyndromeId {
    pub const ALL: [Self; 6] = [
        Self::IliRespiratory,
        Self::CovidLikeIllness,
        Self::AcuteGiEnteric,
        Self::SuspectedOverdose,
        Self::HeatRelatedIllness,
        Self::FeverRashExanthem,
    ];

    pub fn code(self) -> &'static str {
        match self {
            Self::IliRespiratory => "ILI_NSSP_V2",
            Self::CovidLikeIllness => "CLI_NSSP_V2",
            Self::AcuteGiEnteric => "GI_ENTERIC_V1",
            Self::SuspectedOverdose => "DOSE_OVERDOSE_V2",
            Self::HeatRelatedIllness => "HRI_HEAT_V2",
            Self::FeverRashExanthem => "FEVER_RASH_V1",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::IliRespiratory => "Influenza-Like Illness (ILI)",
            Self::CovidLikeIllness => "COVID-Like Illness (CLI)",
            Self::AcuteGiEnteric => "Acute GI / Enteric",
            Self::SuspectedOverdose => "Suspected Unintentional Overdose (DOSE)",
            Self::HeatRelatedIllness => "Environmental Heat-Related Illness (HRI)",
            Self::FeverRashExanthem => "Fever + Rash / Exanthem",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdTriageEncounter {
    pub visit_id: String,
    pub visit_iso_date: String,
    pub facility_id: String,
    pub patient_age: u8,
    pub chief_complaint_text: String,
    /// Semicolon- or space-separated ICD-10-CM discharge diagnosis codes (may be empty at initial triage).
    pub icd10_cm_codes: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClassifiedEdEncounter {
    pub visit_id: String,
    pub visit_iso_date: String,
    pub facility_id: String,
    pub normalized_cc_text: String,
    pub matched_syndromes: Vec<EssenceSyndromeId>,
    pub matched_inclusion_terms: Vec<String>,
    pub suppressed_negated_terms: Vec<String>,
    pub custom_query_hit: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyndromeSummaryRow {
    pub syndrome: EssenceSyndromeId,
    pub code: String,
    pub label: String,
    pub matched_visits: usize,
    pub pct_of_ed_visits: f64,
    pub negated_false_positives_blocked: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EssenceBatchClassification {
    pub total_encounters: usize,
    pub encounters_matching_any_syndrome: usize,
    pub total_negations_blocked: usize,
    pub custom_query_expr: String,
    pub custom_query_matches: usize,
    pub syndrome_summaries: Vec<SyndromeSummaryRow>,
    pub classified_encounters: Vec<ClassifiedEdEncounter>,
}

const PRE_NEGATION_TRIGGERS: &[&str] = &[
    "denies",
    "denied",
    "deny",
    "no",
    "not",
    "without",
    "negative",
    "neg",
    "rule out",
    "r o",
    "free of",
    "absent",
    "resolved",
    "hx of",
    "history of",
];

const CLAUSE_BREAKERS: &[&str] = &[
    ";",
    ".",
    "but",
    "however",
    "except",
    "now",
    "presents",
    "presenting",
    "complains",
    "c o",
];

/// Normalize ED triage shorthand and punctuation while preserving clause boundary markers (`;` and `.`).
pub fn normalize_chief_complaint(raw: &str) -> String {
    let mut lowered = String::with_capacity(raw.len() + 16);
    for ch in raw.chars() {
        if ch.is_ascii_alphanumeric() || ch.is_whitespace() {
            lowered.push(ch.to_ascii_lowercase());
        } else if ch == ';' || ch == '.' {
            lowered.push(' ');
            lowered.push(ch);
            lowered.push(' ');
        } else if ch == '/' || ch == '-' || ch == ',' {
            lowered.push(' ');
        }
    }

    // Expand standard ED triage abbreviations token-by-token
    let mut expanded_tokens: Vec<&str> = Vec::new();
    let raw_toks: Vec<&str> = lowered.split_whitespace().collect();
    let mut i = 0;
    while i < raw_toks.len() {
        if i + 2 < raw_toks.len()
            && raw_toks[i] == "n"
            && raw_toks[i + 1] == "v"
            && raw_toks[i + 2] == "d"
        {
            expanded_tokens.extend(["nausea", "vomiting", "diarrhea"]);
            i += 3;
            continue;
        }
        if i + 1 < raw_toks.len() && raw_toks[i] == "n" && raw_toks[i + 1] == "v" {
            expanded_tokens.extend(["nausea", "vomiting"]);
            i += 2;
            continue;
        }
        match raw_toks[i] {
            "sob" => expanded_tokens.extend(["shortness", "of", "breath"]),
            "uri" => expanded_tokens.extend(["upper", "respiratory", "infection"]),
            "ili" => expanded_tokens.extend(["influenza", "like", "illness"]),
            "od" => expanded_tokens.push("overdose"),
            "ha" => expanded_tokens.push("headache"),
            "fev" | "tmax" => expanded_tokens.push("fever"),
            "vomit" | "emesis" => expanded_tokens.push("vomiting"),
            other => expanded_tokens.push(other),
        }
        i += 1;
    }
    expanded_tokens.join(" ")
}

/// Check if a candidate symptom phrase appears affirmatively (un-negated) or only under negation.
/// Returns `(affirmative_hit, negated_hit)`.
pub fn match_affirmative_or_negated(normalized_cc: &str, phrase: &str, neg_window_tokens: usize) -> (bool, bool) {
    let target_toks: Vec<&str> = phrase.split_whitespace().collect();
    if target_toks.is_empty() {
        return (false, false);
    }
    let cc_toks: Vec<&str> = normalized_cc.split_whitespace().collect();
    if cc_toks.len() < target_toks.len() {
        return (false, false);
    }

    let mut found_affirmative = false;
    let mut found_negated = false;

    for pos in 0..=(cc_toks.len() - target_toks.len()) {
        if cc_toks[pos..pos + target_toks.len()] == target_toks[..] {
            // Inspect up to `neg_window_tokens` preceding tokens up to the nearest clause breaker
            let win_start = pos.saturating_sub(neg_window_tokens.max(1));
            let mut active_slice_start = win_start;
            for back in (win_start..pos).rev() {
                if CLAUSE_BREAKERS.contains(&cc_toks[back]) {
                    active_slice_start = back + 1;
                    break;
                }
            }
            let prefix = cc_toks[active_slice_start..pos].join(" ");
            let prefix_padded = format!(" {prefix} ");
            let is_neg = PRE_NEGATION_TRIGGERS.iter().any(|trig| {
                let pat = format!(" {trig} ");
                prefix_padded.contains(&pat)
            });

            if is_neg {
                found_negated = true;
            } else {
                found_affirmative = true;
            }
        }
    }

    (found_affirmative, found_negated)
}

fn normalize_icd_list(icd_raw: &str) -> Vec<String> {
    icd_raw
        .split(|c: char| c == ';' || c == ',' || c.is_whitespace() || c == '|')
        .map(|s| s.trim().replace('.', "").to_ascii_uppercase())
        .filter(|s| !s.is_empty())
        .collect()
}

fn has_icd_prefix(codes: &[String], prefixes: &[&str]) -> Option<String> {
    for c in codes {
        for &p in prefixes {
            if c.starts_with(p) {
                return Some(format!("ICD:{c}"));
            }
        }
    }
    None
}

/// Evaluate a custom NSSP ESSENCE boolean query of the form:
/// `(^term1^ OR ^term2^) ANDNOT (^excl1^ OR ^excl2^)` against concatenated normalized CC + ICD string.
pub fn evaluate_custom_essence_query(normalized_ccdd: &str, query_expr: &str) -> bool {
    let q = query_expr.trim().to_ascii_lowercase();
    if q.is_empty() {
        return false;
    }

    let (inc_part, exc_part) = if let Some(idx) = q.find("andnot") {
        (&q[..idx], Some(&q[idx + 6..]))
    } else {
        (q.as_str(), None)
    };

    let extract_terms = |segment: &str| -> Vec<String> {
        let mut out = Vec::new();
        if segment.contains('^') {
            let parts: Vec<&str> = segment.split('^').collect();
            for k in (1..parts.len()).step_by(2) {
                let t = parts[k].trim().replace('.', "");
                if !t.is_empty() {
                    out.push(t);
                }
            }
        } else {
            for piece in segment.split("or") {
                let cleaned = piece
                    .trim()
                    .trim_matches(|c: char| c == '(' || c == ')')
                    .trim();
                if !cleaned.is_empty() {
                    out.push(cleaned.to_string());
                }
            }
        }
        out
    };

    let inc_terms = extract_terms(inc_part);
    if inc_terms.is_empty() {
        return false;
    }

    let padded = format!(" {normalized_ccdd} ");
    let inc_hit = inc_terms.iter().any(|t| padded.contains(t.as_str()));
    if !inc_hit {
        return false;
    }

    if let Some(exc_seg) = exc_part {
        let exc_terms = extract_terms(exc_seg);
        if exc_terms.iter().any(|t| padded.contains(t.as_str())) {
            return false;
        }
    }
    true
}

/// Classify a batch of ED triage encounters across all 6 CDC NSSP ESSENCE syndromes + optional custom query.
pub fn classify_ed_encounters_essence(
    encounters: &[EdTriageEncounter],
    custom_query_expr: &str,
    neg_window_tokens: usize,
) -> EssenceBatchClassification {
    let mut classified = Vec::with_capacity(encounters.len());
    let mut syn_counts: BTreeMap<EssenceSyndromeId, usize> = BTreeMap::new();
    let mut syn_neg_blocked: BTreeMap<EssenceSyndromeId, usize> = BTreeMap::new();
    for s in EssenceSyndromeId::ALL {
        syn_counts.insert(s, 0);
        syn_neg_blocked.insert(s, 0);
    }

    let mut any_matches = 0_usize;
    let mut total_neg_blocked = 0_usize;
    let mut custom_matches = 0_usize;

    for enc in encounters {
        let norm_cc = normalize_chief_complaint(&enc.chief_complaint_text);
        let icd_codes = normalize_icd_list(&enc.icd10_cm_codes);
        let ccdd_concat = format!("{} {}", norm_cc, icd_codes.join(" ").to_ascii_lowercase());

        let mut matched_syndromes = Vec::new();
        let mut matched_terms = Vec::new();
        let mut suppressed_terms = Vec::new();

        // Helper closure checking affirmative symptom list + global encounter exclusion phrases
        let mut check_terms = |syn: EssenceSyndromeId,
                               phrases: &[&str],
                               exclusions: &[&str]|
         -> Option<String> {
            if exclusions.iter().any(|ex| norm_cc.contains(ex)) {
                suppressed_terms.push(format!("{}:excluded_context", syn.code()));
                *syn_neg_blocked.entry(syn).or_insert(0) += 1;
                total_neg_blocked += 1;
                return None;
            }
            let mut hit_term: Option<String> = None;
            let mut only_neg = false;
            for &ph in phrases {
                let (aff, neg) = match_affirmative_or_negated(&norm_cc, ph, neg_window_tokens);
                if aff && hit_term.is_none() {
                    hit_term = Some(ph.to_string());
                } else if neg {
                    only_neg = true;
                    suppressed_terms.push(format!("negated({ph})"));
                }
            }
            if hit_term.is_none() && only_neg {
                *syn_neg_blocked.entry(syn).or_insert(0) += 1;
                total_neg_blocked += 1;
            }
            hit_term
        };

        // 1. CDC ILI: (fever AND (cough OR sore throat OR pharyngitis)) OR explicit influenza/ILI OR ICD J09-J11
        let ili_exclusions = &["immunization", "vaccine only", "flu shot"];
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["J09", "J10", "J11"]) {
            matched_syndromes.push(EssenceSyndromeId::IliRespiratory);
            matched_terms.push(icd_m);
        } else if let Some(direct_ili) = check_terms(
            EssenceSyndromeId::IliRespiratory,
            &["influenza", "flu like"],
            ili_exclusions,
        ) {
            matched_syndromes.push(EssenceSyndromeId::IliRespiratory);
            matched_terms.push(direct_ili);
        } else {
            let fev = check_terms(
                EssenceSyndromeId::IliRespiratory,
                &["fever", "febrile", "chills"],
                ili_exclusions,
            );
            let resp = check_terms(
                EssenceSyndromeId::IliRespiratory,
                &["cough", "sore throat", "pharyngitis"],
                ili_exclusions,
            );
            if let (Some(f_t), Some(r_t)) = (fev, resp) {
                matched_syndromes.push(EssenceSyndromeId::IliRespiratory);
                matched_terms.push(format!("{f_t}+{r_t}"));
            }
        }

        // 2. COVID-Like Illness (CLI): U071/J1282 OR (cough/shortness of breath/anosmia + fever/covid)
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["U071", "J1282"]) {
            matched_syndromes.push(EssenceSyndromeId::CovidLikeIllness);
            matched_terms.push(icd_m);
        } else if let Some(cli_hit) = check_terms(
            EssenceSyndromeId::CovidLikeIllness,
            &["covid", "sars cov", "loss of taste", "loss of smell", "anosmia"],
            &["covid test for travel", "vaccine"],
        ) {
            matched_syndromes.push(EssenceSyndromeId::CovidLikeIllness);
            matched_terms.push(cli_hit);
        } else {
            let dys = match_affirmative_or_negated(&norm_cc, "shortness of breath", neg_window_tokens).0
                || match_affirmative_or_negated(&norm_cc, "hypoxia", neg_window_tokens).0;
            let fev = match_affirmative_or_negated(&norm_cc, "fever", neg_window_tokens).0;
            if dys && fev {
                matched_syndromes.push(EssenceSyndromeId::CovidLikeIllness);
                matched_terms.push("fever+shortness_of_breath".to_string());
            }
        }

        // 3. Acute GI / Enteric: A08/A09/R112 OR vomiting/diarrhea/gastroenteritis (excluding pregnancy/chemo)
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["A08", "A09", "R197"]) {
            matched_syndromes.push(EssenceSyndromeId::AcuteGiEnteric);
            matched_terms.push(icd_m);
        } else if let Some(gi_hit) = check_terms(
            EssenceSyndromeId::AcuteGiEnteric,
            &["diarrhea", "vomiting", "gastroenteritis", "norovirus", "watery stool"],
            &["chemotherapy", "morning sickness", "hyperemesis gravidarum"],
        ) {
            matched_syndromes.push(EssenceSyndromeId::AcuteGiEnteric);
            matched_terms.push(gi_hit);
        }

        // 4. Suspected Unintentional Overdose (CDC DOSE): T40/T42/T43 (non-intentional) OR narcan/naloxone/overdose
        let od_exclusions = &["self harm", "suicide attempt", "intentional overdose", "detox clearance"];
        let mut od_intentional_excluded = false;
        if !od_exclusions.iter().any(|ex| norm_cc.contains(ex)) {
            if let Some(icd_m) = has_icd_prefix(&icd_codes, &["T400", "T401", "T402", "T404", "T405", "T436"]) {
                matched_syndromes.push(EssenceSyndromeId::SuspectedOverdose);
                matched_terms.push(icd_m);
            } else if let Some(od_hit) = check_terms(
                EssenceSyndromeId::SuspectedOverdose,
                &["overdose", "naloxone", "narcan", "fentanyl", "heroin", "pinpoint pupils"],
                od_exclusions,
            ) {
                matched_syndromes.push(EssenceSyndromeId::SuspectedOverdose);
                matched_terms.push(od_hit);
            }
        } else {
            od_intentional_excluded = true;
        }

        // 5. Heat-Related Illness (HRI): T67/X30 OR heat exhaustion / heat stroke / sunstroke (excluding "heat pack"!)
        let hri_exclusions = &["heat pack", "heating pad", "chili heat", "prickly heat"];
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["T67", "X30"]) {
            matched_syndromes.push(EssenceSyndromeId::HeatRelatedIllness);
            matched_terms.push(icd_m);
        } else if let Some(hri_hit) = check_terms(
            EssenceSyndromeId::HeatRelatedIllness,
            &[
                "heat exhaustion",
                "heat stroke",
                "sunstroke",
                "heat cramps",
                "working outside in heat",
                "overheated outdoors",
            ],
            hri_exclusions,
        ) {
            matched_syndromes.push(EssenceSyndromeId::HeatRelatedIllness);
            matched_terms.push(hri_hit);
        }

        // 6. Fever + Rash / Exanthem: B05/B06 OR (fever + rash / maculopapular / koplik)
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["B05", "B06"]) {
            matched_syndromes.push(EssenceSyndromeId::FeverRashExanthem);
            matched_terms.push(icd_m);
        } else {
            let rash = check_terms(
                EssenceSyndromeId::FeverRashExanthem,
                &["rash", "maculopapular", "exanthem", "koplik"],
                &["diaper rash", "poison oak", "ivy rash"],
            );
            let fev = match_affirmative_or_negated(&norm_cc, "fever", neg_window_tokens).0;
            if let Some(r_t) = rash {
                if fev || r_t == "koplik" {
                    matched_syndromes.push(EssenceSyndromeId::FeverRashExanthem);
                    matched_terms.push(format!("fever+{r_t}"));
                }
            }
        }

        if od_intentional_excluded {
            suppressed_terms.push("DOSE_OVERDOSE_V2:intentional_excluded".to_string());
            *syn_neg_blocked.entry(EssenceSyndromeId::SuspectedOverdose).or_insert(0) += 1;
            total_neg_blocked += 1;
        }

        for &syn in &matched_syndromes {
            *syn_counts.entry(syn).or_insert(0) += 1;
        }
        if !matched_syndromes.is_empty() {
            any_matches += 1;
        }

        let custom_hit = evaluate_custom_essence_query(&ccdd_concat, custom_query_expr);
        if custom_hit {
            custom_matches += 1;
        }

        classified.push(ClassifiedEdEncounter {
            visit_id: enc.visit_id.clone(),
            visit_iso_date: enc.visit_iso_date.clone(),
            facility_id: enc.facility_id.clone(),
            normalized_cc_text: norm_cc,
            matched_syndromes,
            matched_inclusion_terms: matched_terms,
            suppressed_negated_terms: suppressed_terms,
            custom_query_hit: custom_hit,
        });
    }

    let denom = (encounters.len().max(1)) as f64;
    let syndrome_summaries = EssenceSyndromeId::ALL
        .iter()
        .map(|&s| {
            let c = syn_counts.get(&s).copied().unwrap_or(0);
            SyndromeSummaryRow {
                syndrome: s,
                code: s.code().to_string(),
                label: s.label().to_string(),
                matched_visits: c,
                pct_of_ed_visits: (c as f64) * 100.0 / denom,
                negated_false_positives_blocked: syn_neg_blocked.get(&s).copied().unwrap_or(0),
            }
        })
        .collect();

    EssenceBatchClassification {
        total_encounters: encounters.len(),
        encounters_matching_any_syndrome: any_matches,
        total_negations_blocked: total_neg_blocked,
        custom_query_expr: custom_query_expr.to_string(),
        custom_query_matches: custom_matches,
        syndrome_summaries,
        classified_encounters: classified,
    }
}

/// Aggregate classified ED visits for a chosen syndrome into an MMWR-week `TimeSeriesDataset` for SQLite & EARS/Farrington.
pub fn build_syndrome_weekly_dataset(
    batch: &EssenceBatchClassification,
    syndrome: EssenceSyndromeId,
) -> TimeSeriesDataset {
    let mut week_buckets: BTreeMap<(i32, u8), (Option<epi_core::MmwrWeek>, usize)> = BTreeMap::new();
    for enc in &batch.classified_encounters {
        let Some(mmwr) = parse_flexible_date_to_mmwr(&enc.visit_iso_date) else {
            continue;
        };
        let entry = week_buckets
            .entry((mmwr.year, mmwr.week))
            .or_insert((Some(mmwr), 0));
        if enc.matched_syndromes.contains(&syndrome) {
            entry.1 += 1;
        }
    }

    let points: Vec<TimeSeriesPoint> = week_buckets
        .into_iter()
        .map(|((yr, wk), (mmwr_opt, cnt))| TimeSeriesPoint {
            date_label: format!("{yr}-W{wk:02}"),
            epiweek: mmwr_opt,
            count: cnt as f64,
            suppressed: false,
        })
        .collect();

    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for p in &points {
        for b in p.date_label.bytes() {
            hash ^= u64::from(b);
            hash = hash.wrapping_mul(0x0100_0000_01b3);
        }
        hash ^= p.count.to_bits();
        hash = hash.wrapping_mul(0x0100_0000_01b3);
    }

    TimeSeriesDataset {
        vintage_id: format!("nssp_essence_{}_v1", syndrome.code().to_ascii_lowercase()),
        title: format!("NSSP ESSENCE Weekly ED Visits: {}", syndrome.label()),
        source_kind: SourceKind::LocalCsvTsv,
        source_url: "local://nssp-essence-ccdd-classifier".to_string(),
        ingested_at_utc: "2026-11-13T00:00:00Z".to_string(),
        default_si_preset: "influenza".to_string(),
        checksum_hex: format!("{:08x}", (hash & 0xffff_ffff) as u32),
        points,
    }
}

/// Parse ED encounters CSV (`visit_id,visit_iso_date,facility_id,patient_age,chief_complaint_text,icd10_cm_codes`).
pub fn parse_ed_encounters_csv(csv_text: &str) -> Result<Vec<EdTriageEncounter>, EpiError> {
    let mut rows = Vec::new();
    for (idx, raw_line) in csv_text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        if idx == 0 && cols[0].eq_ignore_ascii_case("visit_id") {
            continue;
        }
        if cols.len() < 5 {
            return Err(EpiError::Schema(format!(
                "ED encounters CSV line {} requires >=5 columns",
                idx + 1
            )));
        }
        rows.push(EdTriageEncounter {
            visit_id: cols[0].to_string(),
            visit_iso_date: cols[1].to_string(),
            facility_id: cols[2].to_string(),
            patient_age: cols[3].parse().unwrap_or(35),
            chief_complaint_text: cols[4].to_string(),
            icd10_cm_codes: cols.get(5).copied().unwrap_or("").to_string(),
        });
    }
    if rows.is_empty() {
        return Err(EpiError::Schema(
            "No valid ED encounters found in CSV".to_string(),
        ));
    }
    Ok(rows)
}

/// Built-in 8-week 68-encounter ED syndromic triage dataset with clinical negations, exclusion idioms, and an ILI/GI wave.
pub fn builtin_ed_chief_complaint_encounters() -> Vec<EdTriageEncounter> {
    let mut out = Vec::new();
    let week_dates = [
        "2026-09-21",
        "2026-09-28",
        "2026-10-05",
        "2026-10-12",
        "2026-10-19",
        "2026-10-26",
        "2026-11-02",
        "2026-11-09",
    ];

    // Tricky negation & exclusion benchmark encounters (must NOT trigger false-positive syndrome counts)
    let tricky_negations: &[(&str, &str, &str)] = &[
        (
            "2026-09-22",
            "Lac right thumb; pt denies fever cough or sob",
            "S61.011A",
        ),
        (
            "2026-09-29",
            "Low back strain after lifting; applied heat pack at home without relief",
            "M54.50",
        ),
        (
            "2026-10-06",
            "Ankle inversion injury; no fever no chills negative for vomiting or diarrhea",
            "S93.401A",
        ),
        (
            "2026-10-13",
            "Here for routine flu shot immunization; denies cough or sore throat",
            "Z23",
        ),
        (
            "2026-10-20",
            "Intentional overdose self harm evaluation by psychiatry",
            "T43.202A",
        ),
        (
            "2026-10-27",
            "8 wks pregnant morning sickness hyperemesis gravidarum; denies fever",
            "O21.0",
        ),
        (
            "2026-11-03",
            "Contact poison oak rash on forearms after trail work; no fever",
            "L23.7",
        ),
        (
            "2026-11-10",
            "Rule out DVT left calf; denies sob denies chest pain no fever",
            "I82.402",
        ),
    ];

    for (idx, &(dt, cc, icd)) in tricky_negations.iter().enumerate() {
        out.push(EdTriageEncounter {
            visit_id: format!("ED-NEG-{:03}", idx + 1),
            visit_iso_date: dt.to_string(),
            facility_id: "OR-GOODSAM-ED".to_string(),
            patient_age: 29 + (idx as u8) * 4,
            chief_complaint_text: cc.to_string(),
            icd10_cm_codes: icd.to_string(),
        });
    }

    // Weekly background + escalating Weeks 6-8 ILI and Acute GI outbreak signal
    let ili_per_week = [1_usize, 1, 2, 2, 3, 5, 7, 10];
    let gi_per_week = [1_usize, 1, 1, 2, 2, 4, 5, 6];
    let mut seq = 1_usize;

    for (w_idx, &dt) in week_dates.iter().enumerate() {
        for k in 0..ili_per_week[w_idx] {
            let (cc, icd) = if k % 3 == 0 {
                ("High fever chills myalgias and productive cough x2 days", "J10.1")
            } else if k % 3 == 1 {
                ("Tmax 102.4 sore throat cough; denies vomiting", "R50.9;J02.9")
            } else {
                ("Fever sob and loss of smell after household covid exposure", "U07.1")
            };
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-GOODSAM-ED".to_string(),
                patient_age: 19 + ((seq * 7) % 62) as u8,
                chief_complaint_text: cc.to_string(),
                icd10_cm_codes: icd.to_string(),
            });
            seq += 1;
        }

        for k in 0..gi_per_week[w_idx] {
            let (cc, icd) = if k % 2 == 0 {
                ("Acute n/v/d watery stool after banquet; denies cough", "A08.4")
            } else {
                ("Intractable emesis and diarrhea x18h", "A09")
            };
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-RIVERBEND-ED".to_string(),
                patient_age: 12 + ((seq * 5) % 58) as u8,
                chief_complaint_text: cc.to_string(),
                icd10_cm_codes: icd.to_string(),
            });
            seq += 1;
        }

        // Add sentinel overdose, heat, and fever-rash alerts in selected weeks
        if w_idx % 2 == 1 {
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-RIVERBEND-ED".to_string(),
                patient_age: 34,
                chief_complaint_text:
                    "Found unresponsive with pinpoint pupils; EMS gave narcan 4mg intranasal"
                        .to_string(),
                icd10_cm_codes: "T40.411A".to_string(),
            });
            seq += 1;
        }
        if w_idx == 0 {
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-GOODSAM-ED".to_string(),
                patient_age: 46,
                chief_complaint_text:
                    "Dizziness and heat exhaustion after working outside in heat harvesting"
                        .to_string(),
                icd10_cm_codes: "T67.5XXA".to_string(),
            });
            seq += 1;
        }
        if w_idx >= 6 {
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-GOODSAM-ED".to_string(),
                patient_age: 7,
                chief_complaint_text:
                    "Unvaccinated child with high fever coryza and maculopapular rash spreading from face"
                        .to_string(),
                icd10_cm_codes: "B05.9".to_string(),
            });
            seq += 1;
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_negex_blocks_denied_symptoms_and_heat_pack_idiom() {
        let encounters = builtin_ed_chief_complaint_encounters();
        let batch = classify_ed_encounters_essence(
            &encounters,
            "(^fever^ OR ^cough^) ANDNOT (^immunization^)",
            5,
        );

        // First 8 encounters are all negated or excluded idioms -> zero syndrome matches!
        for neg_enc in batch.classified_encounters.iter().take(8) {
            assert!(
                neg_enc.matched_syndromes.is_empty(),
                "Negated/excluded encounter {} falsely matched {:?}",
                neg_enc.visit_id,
                neg_enc.matched_syndromes
            );
        }
        assert!(batch.total_negations_blocked >= 6);

        let ds = build_syndrome_weekly_dataset(&batch, EssenceSyndromeId::IliRespiratory);
        assert_eq!(ds.points.len(), 8);
        assert!(ds.points.last().unwrap().count > ds.points.first().unwrap().count);
    }
}
