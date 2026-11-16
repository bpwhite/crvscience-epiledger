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

use epi_analytics::renewal_rt::{compute_cori_rt, CoriRtConfig};
use epi_analytics::wastewater_demix::{
    builtin_wastewater_barcode_and_samples, deconvolve_wastewater_timeseries, WastewaterDemixConfig,
};
use epi_core::{
    builtin_measles_weekly_dataset, CellValue, SuppressionPolicy,
};
use epi_ingest::syndromic_essence::{
    build_syndrome_weekly_dataset, builtin_ed_chief_complaint_encounters,
    classify_ed_encounters_essence, EssenceSyndromeId,
};
use epi_store::{
    encode_pprl_record, link_and_deduplicate_pprl, release_differentially_private_series,
    verify_audit_chain, DpReleaseConfig, EpiStore, PprlConfig, PprlDemographicInput,
};
use rusqlite::Connection;
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};

const MAX_VALIDATION_FILE_BYTES: usize = 384 * 1024;
const PREPRINT_TELEMETRY_REL: &str =
    "../../../app-preprints/releases/epiledger-differential-privacy-epidemiology-validation-telemetry.json";

#[derive(Debug, Clone, Default, Serialize)]
struct T1DifferentialPrivacyAndAuditReport {
    vintage_id: String,
    series_points: usize,
    nchs_suppressed_small_cells_le9: usize,
    nchs_exact_cells: usize,
    dp_release_1_epsilon: f64,
    dp_release_1_alpha: f64,
    dp_release_1_theoretical_rmse: f64,
    dp_release_1_empirical_mae: f64,
    dp_release_1_empirical_rmse: f64,
    dp_cumulative_epsilon_after_two_releases: f64,
    dp_budget_cap: f64,
    budget_overrun_rejected: bool,
    audit_chain_entries_before_tamper: usize,
    audit_chain_valid_before_tamper: bool,
    tamper_detected: bool,
    tamper_broken_at_seq: Option<i64>,
    tamper_failure_reason: Option<String>,
}

#[derive(Debug, Clone, Default, Serialize)]
struct T2PprlAndNegExReport {
    pprl_perturbed_records: usize,
    pprl_ground_truth_entities: usize,
    pprl_linked_entities: usize,
    pprl_evaluated_pairs: usize,
    pprl_cartesian_pairs: usize,
    pprl_blocking_reduction_pct: f64,
    pprl_true_positives: usize,
    pprl_false_positives: usize,
    pprl_false_negatives: usize,
    pprl_pair_precision: f64,
    pprl_pair_recall: f64,
    pprl_pair_f1: f64,
    pprl_hospital_lab_unique_entities: usize,
    essence_total_encounters: usize,
    essence_negation_controls_tested: usize,
    essence_negation_false_positives: usize,
    essence_negation_accuracy: f64,
    essence_total_negations_blocked: usize,
    essence_ili_matched_visits: usize,
    essence_cli_matched_visits: usize,
}

#[derive(Debug, Clone, Default, Serialize)]
struct T3CoriRtAndWastewaterReport {
    measles_weeks_total: usize,
    measles_rt_valid_weeks: usize,
    measles_rt_peak_mean: f64,
    measles_rt_peak_cri_low_95: f64,
    measles_rt_peak_cri_high_95: f64,
    measles_rt_terminal_mean: f64,
    measles_rt_all_intervals_ordered: bool,
    respiratory_ili_weeks: usize,
    respiratory_used_daily_em_disaggregation: bool,
    respiratory_rt_terminal_mean: f64,
    respiratory_rt_terminal_cri_low_95: f64,
    respiratory_rt_terminal_cri_high_95: f64,
    wastewater_samples: usize,
    wastewater_lineages: Vec<String>,
    wastewater_mean_weighted_l1_loss: f64,
    wastewater_max_simplex_deviation: f64,
    wastewater_wk1_dominant_lineage: String,
    wastewater_wk1_dominant_share: f64,
    wastewater_wk10_dominant_lineage: String,
    wastewater_wk10_dominant_share: f64,
    wastewater_wk10_cryptic_sites_detected: Vec<String>,
    csv_output_path: String,
    csv_bytes: usize,
}

#[derive(Debug, Clone, Default, Serialize)]
struct ValidationSummary {
    schema_version: String,
    repository: String,
    version: String,
    preprint_id: String,
    concept_doi: String,
    validation_mode: String,
    t1_dp_and_audit: Option<T1DifferentialPrivacyAndAuditReport>,
    t2_pprl_and_negex: Option<T2PprlAndNegExReport>,
    t3_cori_rt_and_wastewater: Option<T3CoriRtAndWastewaterReport>,
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn validation_dir() -> PathBuf {
    let dir = repo_root().join("validation");
    std::fs::create_dir_all(&dir).expect("create validation directory");
    dir
}

fn summary_store() -> &'static Mutex<ValidationSummary> {
    static STORE: OnceLock<Mutex<ValidationSummary>> = OnceLock::new();
    STORE.get_or_init(|| {
        Mutex::new(ValidationSummary {
            schema_version: "1.0.0".to_string(),
            repository: "crvscience-epiledger".to_string(),
            version: "2.0.0".to_string(),
            preprint_id: "APP-PREPRINT-06".to_string(),
            concept_doi: "10.5281/zenodo.24222382".to_string(),
            validation_mode: "headless_real_rust_epiledger_workspace".to_string(),
            ..ValidationSummary::default()
        })
    })
}

fn write_and_mirror_summary(update: impl FnOnce(&mut ValidationSummary)) {
    let mut guard = summary_store().lock().expect("lock summary");
    update(&mut guard);
    let json = serde_json::to_string_pretty(&*guard).expect("serialize JSON");
    assert!(
        json.len() < MAX_VALIDATION_FILE_BYTES,
        "summary JSON {} bytes exceeds 384 KB",
        json.len()
    );
    let out_path = validation_dir().join("epiledger_validation_summary.json");
    std::fs::write(&out_path, json.as_bytes()).expect("write epiledger_validation_summary.json");

    let preprint_path = repo_root().join(PREPRINT_TELEMETRY_REL);
    if let Some(parent) = preprint_path.parent() {
        std::fs::create_dir_all(parent).expect("create app-preprints/releases dir");
    }
    std::fs::write(&preprint_path, json.as_bytes())
        .expect("mirror epiledger telemetry JSON to app-preprints/releases");
}

#[test]
fn t1_discrete_laplace_dp_cell_suppression_and_hash_chain_tamper_verification() {
    let ds = builtin_measles_weekly_dataset();
    let pid = std::process::id();
    let db_path = std::env::temp_dir().join(format!("epiledger_t1_audit_{pid}.sqlite3"));
    let _ = std::fs::remove_file(&db_path);

    let mut store = EpiStore::open_file(&db_path).expect("open temp EpiStore");
    store.commit_vintage(&ds).expect("commit builtin measles vintage");

    let mut suppressed_cnt = 0usize;
    let mut exact_cnt = 0usize;
    for pt in &ds.points {
        let cell = CellValue::from_count(
            Some(pt.count),
            pt.suppressed,
            SuppressionPolicy::NchsWonderLe9,
        );
        if cell == CellValue::Suppressed {
            suppressed_cnt += 1;
            assert!(pt.count > 0.0 && pt.count <= 9.0);
        } else {
            exact_cnt += 1;
        }
    }
    assert!(suppressed_cnt > 0 && exact_cnt > 0);

    let cfg1 = DpReleaseConfig {
        epsilon: 0.75,
        sensitivity_l1: 1.0,
        max_budget_cap: 2.0,
        deterministic_seed: Some(202604),
    };
    let r1 = release_differentially_private_series(&store, &ds, &cfg1)
        .expect("first DP release within cap");
    assert_eq!(r1.noisy_dataset.points.len(), ds.points.len());
    assert!((r1.cumulative_epsilon_for_vintage - 0.75).abs() < 1e-9);

    let cfg2 = DpReleaseConfig {
        epsilon: 1.00,
        sensitivity_l1: 1.0,
        max_budget_cap: 2.0,
        deterministic_seed: Some(202605),
    };
    let r2 = release_differentially_private_series(&store, &ds, &cfg2)
        .expect("second DP release (cum=1.75<=2.0) succeeds");
    assert!((r2.cumulative_epsilon_for_vintage - 1.75).abs() < 1e-9);

    // Budget cap enforcement: 1.75 + 0.50 = 2.25 > 2.00 must fail cleanly
    let cfg3 = DpReleaseConfig {
        epsilon: 0.50,
        sensitivity_l1: 1.0,
        max_budget_cap: 2.0,
        deterministic_seed: Some(202606),
    };
    let over_cap = release_differentially_private_series(&store, &ds, &cfg3);
    assert!(over_cap.is_err(), "DP budget overrun must be refused");

    let ver_clean = verify_audit_chain(&store).expect("verify clean chain");
    assert!(ver_clean.is_valid, "Clean hash chain must verify");
    assert!(ver_clean.total_entries >= 3);

    // Tamper payload on seq_no=1 outside EpiStore API to verify cryptographic SHA-256 detection
    {
        let raw_conn = Connection::open(&db_path).expect("open raw SQLite connection");
        raw_conn
            .execute(
                "UPDATE audit_hash_chain SET detail_summary = 'TAMPERED_RECORD' WHERE seq_no = 1",
                [],
            )
            .expect("tamper row 1");
    }

    let ver_tampered = verify_audit_chain(&store).expect("verify tampered chain");
    assert!(
        !ver_tampered.is_valid,
        "Tampered SHA-256 audit chain MUST fail verification"
    );
    assert_eq!(ver_tampered.broken_at_seq, Some(1));

    drop(store);
    let _ = std::fs::remove_file(&db_path);
    let _ = std::fs::remove_file(db_path.with_extension("sqlite3-wal"));
    let _ = std::fs::remove_file(db_path.with_extension("sqlite3-shm"));

    let report = T1DifferentialPrivacyAndAuditReport {
        vintage_id: ds.vintage_id,
        series_points: ds.points.len(),
        nchs_suppressed_small_cells_le9: suppressed_cnt,
        nchs_exact_cells: exact_cnt,
        dp_release_1_epsilon: r1.epsilon_spent,
        dp_release_1_alpha: r1.geometric_alpha,
        dp_release_1_theoretical_rmse: r1.theoretical_rmse,
        dp_release_1_empirical_mae: r1.empirical_mae,
        dp_release_1_empirical_rmse: r1.empirical_rmse,
        dp_cumulative_epsilon_after_two_releases: r2.cumulative_epsilon_for_vintage,
        dp_budget_cap: 2.0,
        budget_overrun_rejected: over_cap.is_err(),
        audit_chain_entries_before_tamper: ver_clean.total_entries,
        audit_chain_valid_before_tamper: ver_clean.is_valid,
        tamper_detected: !ver_tampered.is_valid,
        tamper_broken_at_seq: ver_tampered.broken_at_seq,
        tamper_failure_reason: ver_tampered.failure_reason,
    };

    write_and_mirror_summary(|s| {
        s.t1_dp_and_audit = Some(report);
    });
}

fn perturb_last_char_or_digit(s: &str, idx: usize) -> String {
    let mut chars: Vec<char> = s.chars().collect();
    if let Some(last) = chars.last_mut() {
        *last = if last.is_ascii_digit() {
            if *last == '0' { '1' } else { '0' }
        } else if idx % 2 == 0 {
            'Z'
        } else {
            'X'
        };
    }
    chars.into_iter().collect()
}

#[test]
fn t2_clk_bloom_pprl_precision_recall_f1_and_essence_negex() {
    let syl_a = ["BEL", "COR", "DEV", "FER", "GAL", "HAR", "KAS", "LOR", "MON", "NOR"];
    let syl_b = ["MAN", "SON", "TON", "BER", "DEN", "FEL", "HAM", "KIN", "LIN", "TER"];
    let g_pre = ["ALE", "BRI", "CAR", "DAN", "ELE", "FRE", "GAB", "HAN", "ISA", "JUL"];
    let g_suf = ["XANDER", "ANNA", "OLINE", "IEL", "NOR", "RICK", "RIEL", "NAH", "BEL", "IAN"];

    let cfg = PprlConfig::default();
    let n_entities = 40usize;
    let mut records = Vec::with_capacity(n_entities * 2);
    let mut truth_entity_by_rec: HashMap<String, usize> = HashMap::new();

    for e in 0..n_entities {
        let fam = format!("{}{}{}", syl_a[e % 10], syl_b[(e / 10) % 10], (b'A' + (e % 26) as u8) as char);
        let giv = format!("{}{}{}{}", g_pre[(e * 3) % 10], g_suf[(e * 7) % 10], (b'K' + (e % 11) as u8) as char, (b'A' + ((e * 5) % 17) as u8) as char);
        let yr = 1954 + (e * 17) % 58;
        let mo = 1 + (e * 7) % 12;
        let dy = 1 + (e * 11) % 28;
        let dob = format!("{yr:04}{mo:02}{dy:02}");
        let sex = if e % 2 == 0 { "F" } else { "M" };
        let zip = format!("{:05}", 97001 + (e * 37) % 899);

        let id_a = format!("ELR_{e:03}");
        let id_b = format!("ADT_{e:03}");
        truth_entity_by_rec.insert(id_a.clone(), e);
        truth_entity_by_rec.insert(id_b.clone(), e);

        records.push(PprlDemographicInput {
            record_id: id_a,
            source_system: "State-ELR".to_string(),
            encounter_date: "2026-10-05".to_string(),
            encounter_day: 278.0,
            family_name: fam.clone(),
            given_name: giv.clone(),
            dob_yyyymmdd: dob.clone(),
            sex: sex.to_string(),
            postal_zip: zip.clone(),
        });

        // Perturb secondary clinical record slightly (1-char family edit / ZIP typo / swapped minor field every 3rd patient)
        let fam_b = if e % 3 == 0 { perturb_last_char_or_digit(&fam, e) } else { fam };
        let zip_b = if e % 4 == 1 { perturb_last_char_or_digit(&zip, e) } else { zip };
        records.push(PprlDemographicInput {
            record_id: id_b,
            source_system: "Hospital-ADT".to_string(),
            encounter_date: "2026-10-08".to_string(),
            encounter_day: 281.0,
            family_name: fam_b,
            given_name: giv,
            dob_yyyymmdd: dob,
            sex: sex.to_string(),
            postal_zip: zip_b,
        });
    }

    let encoded_records: Vec<_> = records.iter().map(|r| encode_pprl_record(r, &cfg)).collect();
    let summary = link_and_deduplicate_pprl(&encoded_records, &cfg);
    let mut tp = 0usize;
    let mut fp = 0usize;
    let mut linked_truth_entities = std::collections::HashSet::new();
    for pair in &summary.candidate_pairs {
        if pair.is_linked_match {
            let ea = truth_entity_by_rec.get(&pair.record_a).copied();
            let eb = truth_entity_by_rec.get(&pair.record_b).copied();
            if ea.is_some() && ea == eb {
                tp += 1;
                linked_truth_entities.insert(ea.unwrap());
            } else {
                fp += 1;
            }
        }
    }
    let fn_cnt = n_entities.saturating_sub(linked_truth_entities.len());
    let precision = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
    let recall = if tp + fn_cnt > 0 { tp as f64 / (tp + fn_cnt) as f64 } else { 0.0 };
    let f1 = if precision + recall > 0.0 {
        2.0 * precision * recall / (precision + recall)
    } else {
        0.0
    };
    assert!(f1 >= 0.90, "Expected PPRL F1 >= 0.90, got {f1:.4} (P={precision:.4}, R={recall:.4}, TP={tp}, FP={fp}, FN={fn_cnt})");

    // Built-in hospital/lab reference linkage check
    let builtin_recs = epi_store::builtin_pprl_hospital_lab_records();
    let encoded_builtin: Vec<_> = builtin_recs.iter().map(|r| encode_pprl_record(r, &cfg)).collect();
    let builtin_sum = link_and_deduplicate_pprl(&encoded_builtin, &cfg);
    assert!(builtin_sum.duplicate_encounters_collapsed >= 1);

    // NSSP ESSENCE NegEx clinical negation filtering verification
    let ed_encounters = builtin_ed_chief_complaint_encounters();
    let batch = classify_ed_encounters_essence(
        &ed_encounters,
        "(^fever^ OR ^cough^) ANDNOT (^immunization^)",
        5,
    );
    let mut neg_fps = 0usize;
    for enc in batch.classified_encounters.iter().take(8) {
        if !enc.matched_syndromes.is_empty() {
            neg_fps += 1;
        }
    }
    assert_eq!(neg_fps, 0, "All 8 NegEx / exclusion controls must remain suppressed");
    let neg_acc = (8 - neg_fps) as f64 / 8.0;

    let ili_visits = batch
        .syndrome_summaries
        .iter()
        .find(|r| r.syndrome == EssenceSyndromeId::IliRespiratory)
        .map(|r| r.matched_visits)
        .unwrap_or(0);
    let cli_visits = batch
        .syndrome_summaries
        .iter()
        .find(|r| r.syndrome == EssenceSyndromeId::CovidLikeIllness)
        .map(|r| r.matched_visits)
        .unwrap_or(0);

    let report = T2PprlAndNegExReport {
        pprl_perturbed_records: records.len(),
        pprl_ground_truth_entities: n_entities,
        pprl_linked_entities: summary.unique_linked_entities,
        pprl_evaluated_pairs: summary.evaluated_candidate_pairs,
        pprl_cartesian_pairs: summary.total_cartesian_pairs,
        pprl_blocking_reduction_pct: summary.blocking_pair_reduction_pct,
        pprl_true_positives: tp,
        pprl_false_positives: fp,
        pprl_false_negatives: fn_cnt,
        pprl_pair_precision: precision,
        pprl_pair_recall: recall,
        pprl_pair_f1: f1,
        pprl_hospital_lab_unique_entities: builtin_sum.unique_linked_entities,
        essence_total_encounters: batch.total_encounters,
        essence_negation_controls_tested: 8,
        essence_negation_false_positives: neg_fps,
        essence_negation_accuracy: neg_acc,
        essence_total_negations_blocked: batch.total_negations_blocked,
        essence_ili_matched_visits: ili_visits,
        essence_cli_matched_visits: cli_visits,
    };

    write_and_mirror_summary(|s| {
        s.t2_pprl_and_negex = Some(report);
    });
}

#[test]
fn t3_cori_rt_posterior_cri_and_freyja_wastewater_deconvolution() {
    // 1. Built-in historical measles MMWR series + DP noisy count release
    let measles_ds = builtin_measles_weekly_dataset();
    let mut store = EpiStore::open_in_memory().expect("in-memory EpiStore");
    store.commit_vintage(&measles_ds).expect("commit measles");
    let dp_outcome = release_differentially_private_series(
        &store,
        &measles_ds,
        &DpReleaseConfig {
            epsilon: 0.75,
            sensitivity_l1: 1.0,
            max_budget_cap: 2.0,
            deterministic_seed: Some(202604),
        },
    )
    .expect("DP release for audit CSV");

    let raw_counts: Vec<f64> = measles_ds.points.iter().map(|p| p.count).collect();
    let measles_cfg = CoriRtConfig::default(); // SI mean 11.7d, sd 2.0d, tau=3 weeks
    let measles_rt = compute_cori_rt(&raw_counts, &measles_cfg);
    assert_eq!(measles_rt.estimates.len(), measles_ds.points.len());

    let mut valid_weeks = 0usize;
    let mut all_ordered = true;
    let mut peak_est = None::<epi_analytics::renewal_rt::CoriRtEstimate>;
    for est in measles_rt.estimates.iter().flatten() {
        valid_weeks += 1;
        if !(est.cri_low_95 <= est.median
            && est.median <= est.cri_high_95
            && est.cri_low_95 <= est.mean
            && est.mean <= est.cri_high_95)
        {
            all_ordered = false;
        }
        if peak_est.map(|p| est.mean > p.mean).unwrap_or(true) {
            peak_est = Some(*est);
        }
    }
    assert!(valid_weeks >= 20 && all_ordered);
    let peak = peak_est.expect("peak measles Rt");
    let term_measles = measles_rt
        .estimates
        .iter()
        .rev()
        .flatten()
        .next()
        .copied()
        .expect("terminal measles Rt");

    // Write validation/cori_rt_and_dp_audit.csv (< 384 KB)
    let csv_path = validation_dir().join("cori_rt_and_dp_audit.csv");
    {
        let mut wtr = csv::Writer::from_path(&csv_path).expect("create cori_rt_and_dp_audit.csv");
        wtr.write_record([
            "week_idx",
            "date_label",
            "mmwr_week",
            "raw_count",
            "nchs_wonder_cell",
            "dp_eps075_count",
            "dp_abs_error",
            "cori_rt_mean",
            "cori_rt_median",
            "cori_rt_cri_low_95",
            "cori_rt_cri_high_95",
        ])
        .expect("write CSV header");

        for (idx, (pt, dp_pt)) in measles_ds
            .points
            .iter()
            .zip(dp_outcome.noisy_dataset.points.iter())
            .enumerate()
        {
            let cell = CellValue::from_count(
                Some(pt.count),
                pt.suppressed,
                SuppressionPolicy::NchsWonderLe9,
            );
            let mmwr_lbl = pt
                .epiweek
                .map(|w| w.compact_label())
                .unwrap_or_else(|| "NA".to_string());
            let (m, med, lo, hi) = match measles_rt.estimates[idx] {
                Some(e) => (
                    format!("{:.4}", e.mean),
                    format!("{:.4}", e.median),
                    format!("{:.4}", e.cri_low_95),
                    format!("{:.4}", e.cri_high_95),
                ),
                None => ("NA".into(), "NA".into(), "NA".into(), "NA".into()),
            };
            wtr.write_record(&[
                (idx + 1).to_string(),
                pt.date_label.clone(),
                mmwr_lbl,
                format!("{:.0}", pt.count),
                cell.display_string(0),
                format!("{:.0}", dp_pt.count),
                format!("{:.1}", (dp_pt.count - pt.count).abs()),
                m,
                med,
                lo,
                hi,
            ])
            .expect("write audit row");
        }
        wtr.flush().expect("flush CSV");
    }
    let csv_bytes = std::fs::metadata(&csv_path).expect("CSV metadata").len() as usize;
    assert!(
        csv_bytes > 0 && csv_bytes < MAX_VALIDATION_FILE_BYTES,
        "cori_rt_and_dp_audit.csv size {csv_bytes} out of bounds"
    );

    // 2. Respiratory ED ILI weekly incidence via daily EM disaggregation (SI < 7.0d)
    let ed_encounters = builtin_ed_chief_complaint_encounters();
    let essence_batch = classify_ed_encounters_essence(
        &ed_encounters,
        "(^fever^ OR ^cough^) ANDNOT (^immunization^)",
        5,
    );
    let ili_ds = build_syndrome_weekly_dataset(&essence_batch, EssenceSyndromeId::IliRespiratory);
    let ili_counts: Vec<f64> = ili_ds.points.iter().map(|p| p.count).collect();
    let resp_cfg = CoriRtConfig {
        si_mean_days: 4.2,
        si_sd_days: 1.9,
        tau_weeks: 2,
        prior_a: 1.0,
        prior_b: 5.0,
    };
    let resp_rt = compute_cori_rt(&ili_counts, &resp_cfg);
    assert!(resp_rt.used_daily_em_disaggregation);
    let term_resp = resp_rt
        .estimates
        .iter()
        .rev()
        .flatten()
        .next()
        .copied()
        .expect("terminal respiratory Rt");
    assert!(
        term_resp.mean > 1.0,
        "Escalating ILI outbreak terminal Rt must exceed 1.0, got {}",
        term_resp.mean
    );

    // 3. Freyja L1 constrained-least-squares wastewater lineage deconvolution
    let (library, samples) = builtin_wastewater_barcode_and_samples();
    let demix = deconvolve_wastewater_timeseries(
        &library,
        &samples,
        &WastewaterDemixConfig::default(),
    )
    .expect("deconvolve builtin wastewater cohort");
    assert_eq!(demix.samples.len(), 10);

    let mut max_simplex_dev = 0.0_f64;
    let mut sum_l1 = 0.0_f64;
    for s in &demix.samples {
        let prob_sum: f64 = s.lineages.iter().map(|l| l.abundance).sum();
        max_simplex_dev = max_simplex_dev.max((prob_sum - 1.0).abs());
        sum_l1 += s.weighted_l1_loss;
    }
    assert!(max_simplex_dev < 1e-6);

    let mut s_lineages: BTreeMap<String, f64> = BTreeMap::new();
    for l in &demix.samples[0].lineages {
        s_lineages.insert(l.lineage.clone(), l.abundance);
    }
    let wk1_dom = demix.samples[0]
        .lineages
        .iter()
        .max_by(|a, b| a.abundance.total_cmp(&b.abundance))
        .expect("wk1 dominant");
    let wk10_dom = demix.samples[9]
        .lineages
        .iter()
        .max_by(|a, b| a.abundance.total_cmp(&b.abundance))
        .expect("wk10 dominant");
    let wk10_cryptic: Vec<String> = demix.samples[9]
        .cryptic_sites
        .iter()
        .map(|c| c.mutation_site.clone())
        .collect();
    assert!(!s_lineages.is_empty());
    assert!(wk10_cryptic.iter().any(|s| s == "S:S31del"));

    let report = T3CoriRtAndWastewaterReport {
        measles_weeks_total: measles_ds.points.len(),
        measles_rt_valid_weeks: valid_weeks,
        measles_rt_peak_mean: peak.mean,
        measles_rt_peak_cri_low_95: peak.cri_low_95,
        measles_rt_peak_cri_high_95: peak.cri_high_95,
        measles_rt_terminal_mean: term_measles.mean,
        measles_rt_all_intervals_ordered: all_ordered,
        respiratory_ili_weeks: ili_ds.points.len(),
        respiratory_used_daily_em_disaggregation: resp_rt.used_daily_em_disaggregation,
        respiratory_rt_terminal_mean: term_resp.mean,
        respiratory_rt_terminal_cri_low_95: term_resp.cri_low_95,
        respiratory_rt_terminal_cri_high_95: term_resp.cri_high_95,
        wastewater_samples: demix.samples.len(),
        wastewater_lineages: demix.lineage_names,
        wastewater_mean_weighted_l1_loss: sum_l1 / (demix.samples.len() as f64),
        wastewater_max_simplex_deviation: max_simplex_dev,
        wastewater_wk1_dominant_lineage: wk1_dom.lineage.clone(),
        wastewater_wk1_dominant_share: wk1_dom.abundance,
        wastewater_wk10_dominant_lineage: wk10_dom.lineage.clone(),
        wastewater_wk10_dominant_share: wk10_dom.abundance,
        wastewater_wk10_cryptic_sites_detected: wk10_cryptic,
        csv_output_path: "validation/cori_rt_and_dp_audit.csv".to_string(),
        csv_bytes,
    };

    write_and_mirror_summary(|s| {
        s.t3_cori_rt_and_wastewater = Some(report);
    });
}
