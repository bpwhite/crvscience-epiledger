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

use crate::sitrep_export::{render_standalone_sitrep_html, SitRepContext};
use epi_analytics::{
    analyze_contact_tracing_network, build_us_state_scan_zones, builtin_age_serosurvey_dataset,
    builtin_all_cause_weekly_mortality_dataset, builtin_emerging_wave_cases_and_deaths,
    builtin_hospital_amr_isolates, builtin_outbreak_fasta_and_metadata,
    builtin_respiratory_lineage_table, builtin_superspreading_outbreak_records,
    builtin_three_source_tb_pertussis_crc, builtin_tnd_vaccine_study_dataset,
    builtin_wastewater_barcode_and_samples, builtin_wildfire_pm25_dlnm_dataset,
    calibrate_and_project_seir, compute_aberrations, compute_clsi_m39_antibiogram,
    compute_cori_rt, compute_delay_adjusted_cfr, compute_disproportionality_2x2,
    compute_euromomo_excess_mortality, compute_kulldorff_spatial_clusters,
    compute_taylor_linearized_mean, compute_two_source_chapman, deconvolve_wastewater_timeseries,
    estimate_tnd_vaccine_effectiveness, export_bipartite_network_dot, export_outbreak_graphviz_dot,
    fit_catalytic_seroprevalence, fit_dlnm_quasi_poisson, fit_multinomial_lineage_logistics,
    fit_negative_binomial_offspring, fit_three_source_loglinear_crc, nowcast_parametric_delay,
    parse_amr_isolates_csv, parse_cfr_incidence_deaths_csv, parse_contact_tracing_csv,
    parse_crc_case_lists_csv, parse_dlnm_exposure_csv, parse_fasta_and_metadata,
    parse_lineage_counts_csv, project_competing_risk_hospital_surge,
    project_hospital_surge_from_seir,
    parse_serosurvey_csv, parse_wastewater_snv_csv, parse_weekly_mortality_csv,
    reconstruct_outbreak_tree, run_kulldorff_circular_scan, run_rolling_wis_backtest_and_ensemble,
    AberrationConfig, AberrationMethod, AntibiogramConfig, BacktestConfig, CatalyticModelKind,
    CoriRtConfig, DelayCfrConfig, DepthWeighting, DlnmConfig, ExcessMortalityConfig,
    ExposureBasisKind, HospitalSurgeConfig,
    IsolateDeduplicationPolicy, LineageBasisMode, LineageLogisticsConfig, LonelyPsuStrategy,
    MortalityBaselineMode, OffspringFitConfig, OutbreakTreeConfig, ParametricNowcastConfig,
    SeirConfig, SpatialScanInputZone, SurveyObservation, ThreeSourceOverlapCounts,
    WastewaterDemixConfig, ZeroTruncationMode,
};
use epi_core::{
    builtin_measles_weekly_dataset, builtin_pertussis_spatial_records, EpiError, SuppressionPolicy,
};
use epi_ingest::{
    build_syndrome_weekly_dataset, builtin_ed_chief_complaint_encounters,
    builtin_hl7_v2_elr_and_nssp_batch, classify_ed_encounters_essence, export_dataset_csv,
    export_timeseries_as_fhir_r4_bundle, import_delimited_timeseries, import_fhir_r4_bundle_file,
    import_hl7_v2_file, parse_ed_encounters_csv, parse_geojson_centroids_file,
    parse_hl7_v2_batch_str, parse_sas_xport_bytes, CdcSodaClient, CmuDelphiClient, ColumnMapping,
    EssenceSyndromeId, OpenFdaClient, SasXportDataset, SasXportValue, SodaDatasetKind,
    WhoGhoClient,
};
use epi_store::{
    builtin_pprl_hospital_lab_records, encode_pprl_record, generate_synthetic_pprl_scale_benchmark,
    link_and_deduplicate_pprl, list_audit_events, release_differentially_private_series,
    verify_audit_chain, DpReleaseConfig, EncodedPprlRecord, EpiStore, PprlConfig,
};
use std::path::Path;
use std::time::Instant;

pub fn extract_survey_obs_from_xpt(
    xpt: &SasXportDataset,
    stratum_col: &str,
    psu_col: &str,
    weight_col: &str,
    outcome_col: &str,
) -> Result<Vec<SurveyObservation>, EpiError> {
    let find_idx = |name: &str| -> Result<usize, EpiError> {
        xpt.variables
            .iter()
            .position(|v| v.name.eq_ignore_ascii_case(name.trim()))
            .ok_or_else(|| EpiError::Schema(format!("Variable '{name}' not found in SAS .XPT")))
    };

    let s_i = find_idx(stratum_col)?;
    let p_i = find_idx(psu_col)?;
    let w_i = find_idx(weight_col)?;
    let y_i = find_idx(outcome_col)?;

    let mut obs = Vec::with_capacity(xpt.rows.len());
    for row in &xpt.rows {
        let Some(SasXportValue::Numeric(s_v)) = row.get(s_i) else { continue };
        let Some(SasXportValue::Numeric(p_v)) = row.get(p_i) else { continue };
        let Some(SasXportValue::Numeric(w_v)) = row.get(w_i) else { continue };
        if *w_v <= 0.0 || !w_v.is_finite() {
            continue;
        }
        match row.get(y_i) {
            Some(SasXportValue::Numeric(y_v)) if y_v.is_finite() => {
                obs.push(SurveyObservation {
                    stratum: s_v.round() as i64,
                    psu: p_v.round() as i64,
                    weight: *w_v,
                    value: *y_v,
                    in_subpop: true,
                });
            }
            _ => {
                obs.push(SurveyObservation {
                    stratum: s_v.round() as i64,
                    psu: p_v.round() as i64,
                    weight: *w_v,
                    value: 0.0,
                    in_subpop: false,
                });
            }
        }
    }
    Ok(obs)
}

pub fn run_cli(args: &[String], db_path: &Path) -> Result<(), EpiError> {
    let mut store = EpiStore::open_file(db_path)?;

    if store.list_vintages()?.is_empty() {
        store.commit_vintage(&builtin_measles_weekly_dataset())?;
    }

    match args[0].as_str() {
        "vintages" => {
            let list = store.list_vintages()?;
            println!(
                "{:<30} {:<18} {:>7}  {:<10}  {}",
                "VINTAGE_ID", "SOURCE", "POINTS", "CHECKSUM", "TITLE"
            );
            for v in list {
                println!(
                    "{:<30} {:<18} {:>7}  {:<10}  {}",
                    v.vintage_id, v.source_kind, v.row_count, v.checksum_hex, v.title
                );
            }
        }
        "analyze" => {
            let vid = args
                .get(1)
                .cloned()
                .unwrap_or_else(|| "cdc_measles_2022_2026_v1".to_string());
            let ds = store
                .load_vintage(&vid)?
                .ok_or_else(|| EpiError::Schema(format!("Vintage '{vid}' not found")))?;
            let counts: Vec<f64> = ds.points.iter().map(|p| p.count).collect();

            let rt = compute_cori_rt(&counts, &CoriRtConfig::default());
            let ab = compute_aberrations(
                &counts,
                &AberrationConfig {
                    method: AberrationMethod::FarringtonQuasiPoisson,
                    ..Default::default()
                },
            );

            println!(
                "Dataset: {} ({} steps, {} Farrington alerts)",
                ds.title,
                ds.points.len(),
                ab.total_alerts
            );
            println!(
                "{:<14} {:>10} {:>12} {:>12} {:>22}",
                "PERIOD", "OBSERVED", "EXPECTED", "UPPER_LIM", "CORI_RT [95% CrI]"
            );
            let start = ds.points.len().saturating_sub(16);
            for i in start..ds.points.len() {
                let p = &ds.points[i];
                let a = &ab.points[i];
                let rt_str = match rt.estimates[i] {
                    Some(est) => format!(
                        "{:.2} [{:.2}, {:.2}]",
                        est.mean, est.cri_low_95, est.cri_high_95
                    ),
                    None => "--".to_string(),
                };
                let flag = if a.is_alert { " *ALERT*" } else { "" };
                println!(
                    "{:<14} {:>10.0} {:>12} {:>12} {:>22}{}",
                    p.date_label,
                    p.count,
                    a.expected_baseline
                        .map(|v| format!("{v:.1}"))
                        .unwrap_or_else(|| "--".into()),
                    a.upper_threshold
                        .map(|v| format!("{v:.1}"))
                        .unwrap_or_else(|| "--".into()),
                    rt_str,
                    flag
                );
            }
        }
        "nowcast" => {
            let vid = args
                .get(1)
                .cloned()
                .unwrap_or_else(|| "cdc_measles_2022_2026_v1".to_string());
            let mean_delay: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1.6);
            let ds = store
                .load_vintage(&vid)?
                .ok_or_else(|| EpiError::Schema(format!("Vintage '{vid}' not found")))?;
            let counts: Vec<f64> = ds.points.iter().map(|p| p.count).collect();

            let cfg = ParametricNowcastConfig {
                delay_mean_weeks: mean_delay,
                ..Default::default()
            };
            let nc = nowcast_parametric_delay(&counts, &cfg);

            println!(
                "Reporting-Delay Nowcast: {} (mean delay {:.1} wk, total imputed unreported = {:.1})",
                ds.title, cfg.delay_mean_weeks, nc.total_imputed_unreported
            );
            println!(
                "{:<14} {:>10} {:>12} {:>12} {:>14} {:>22}",
                "PERIOD", "RAW_OBS", "COMPLETE_%", "UNREPORTED", "NOWCAST_MEAN", "95%_PRED_INTERVAL"
            );
            let start = ds.points.len().saturating_sub(10);
            for i in start..ds.points.len() {
                let p = &ds.points[i];
                let npt = &nc.points[i];
                println!(
                    "{:<14} {:>10.0} {:>11.1}% {:>12.1} {:>14.1} {:>22}",
                    p.date_label,
                    npt.observed_count,
                    npt.completeness_fraction * 100.0,
                    npt.expected_unreported,
                    npt.nowcast_mean,
                    format!("[{:.1}, {:.1}]", npt.pi_low_95, npt.pi_high_95)
                );
            }
        }
        "seir-forecast" => {
            let vid = args
                .get(1)
                .cloned()
                .unwrap_or_else(|| "cdc_measles_2022_2026_v1".to_string());
            let horizon: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(16);
            let runs: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(240);
            let ds = store
                .load_vintage(&vid)?
                .ok_or_else(|| EpiError::Schema(format!("Vintage '{vid}' not found")))?;
            let counts: Vec<f64> = ds.points.iter().map(|p| p.count).collect();

            let cfg = SeirConfig {
                horizon_weeks: horizon,
                ensemble_runs: runs,
                ..Default::default()
            };
            let fc = calibrate_and_project_seir(&counts, true, None, &cfg);

            println!(
                "SEIR Scoped-Parallel Forecast ({} trajectories): {} (calibrated R0 = {:.2}, initial R_eff = {:.2}, peak Wk +{} = {:.1})",
                runs, ds.title, fc.calibrated_r0, fc.initial_r_eff, fc.peak_week_offset, fc.peak_reported_weekly
            );
            println!(
                "{:<10} {:>14} {:>22} {:>12} {:>12} {:>10}",
                "WEEK_OFF", "REPORTED_MEAN", "95%_STOCHASTIC_PI", "EXPOSED", "INFECTIOUS", "R_EFF"
            );
            for s in &fc.steps {
                println!(
                    "Wk +{:<6} {:>14.1} {:>22} {:>12.1} {:>12.1} {:>10.2}",
                    s.week_offset,
                    s.reported_weekly_mean,
                    format!("[{:.0}, {:.0}]", s.reported_pi_low_95, s.reported_pi_high_95),
                    s.exposed,
                    s.infectious,
                    s.effective_rt
                );
            }
        }
        "backtest-wis" => {
            let vid = args
                .get(1)
                .cloned()
                .unwrap_or_else(|| "cdc_measles_2022_2026_v1".to_string());
            let horizon: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
            let origins: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(14);
            let ds = store
                .load_vintage(&vid)?
                .ok_or_else(|| EpiError::Schema(format!("Vintage '{vid}' not found")))?;
            let counts: Vec<f64> = ds.points.iter().map(|p| p.count).collect();
            let bt_cfg = BacktestConfig {
                horizon_weeks: horizon,
                num_origins: origins,
                ..Default::default()
            };
            let bt =
                run_rolling_wis_backtest_and_ensemble(&counts, &bt_cfg, &SeirConfig::default())?;
            println!(
                "Rolling K=11 WIS Backtest & Inverse-WIS Ensemble: {} ({} origins x {}w horizon) | Best={}",
                ds.title, bt.origins_evaluated, bt.horizon_weeks, bt.best_model_label
            );
            println!(
                "{:<30} {:>9} {:>9} {:>9} {:>9} {:>9} {:>8} {:>8} {:>8}",
                "MODEL", "MEAN_WIS", "SHARP", "OVERPRED", "UNDERPRD", "LOG_WIS", "REL_WIS", "COV_95", "WEIGHT"
            );
            for sc in &bt.scorecards {
                let wt = if sc.final_ensemble_weight > 0.0 {
                    format!("{:.1}%", sc.final_ensemble_weight * 100.0)
                } else {
                    "--".to_string()
                };
                println!(
                    "{:<30} {:>9.2} {:>9.2} {:>9.2} {:>9.2} {:>9.3} {:>7.2}x {:>7.1}% {:>8}",
                    sc.model_label,
                    sc.mean_wis,
                    sc.mean_sharpness,
                    sc.mean_overprediction,
                    sc.mean_underprediction,
                    sc.mean_log_wis,
                    sc.relative_wis_vs_naive,
                    sc.coverage_95_pct,
                    wt
                );
            }
        }
        "export-sitrep" => {
            let vid = args
                .get(1)
                .cloned()
                .unwrap_or_else(|| "cdc_measles_2022_2026_v1".to_string());
            let out_html = args
                .get(2)
                .map(String::as_str)
                .unwrap_or("epiledger_sitrep.html");
            let mask = if args.get(3).map(String::as_str) == Some("mask") {
                SuppressionPolicy::NchsWonderLe9
            } else {
                SuppressionPolicy::Unmasked
            };
            let ds = store
                .load_vintage(&vid)?
                .ok_or_else(|| EpiError::Schema(format!("Vintage '{vid}' not found")))?;
            let counts: Vec<f64> = ds.points.iter().map(|p| p.count).collect();
            let rt_series = compute_cori_rt(&counts, &CoriRtConfig::default());
            let ab_series = compute_aberrations(&counts, &AberrationConfig::default());
            let nowcast = nowcast_parametric_delay(&counts, &ParametricNowcastConfig::default());
            let backtest = run_rolling_wis_backtest_and_ensemble(
                &counts,
                &BacktestConfig::default(),
                &SeirConfig::default(),
            )?;
            let zones = build_us_state_scan_zones(&builtin_pertussis_spatial_records(), true);
            let spatial_clusters = run_kulldorff_circular_scan(&zones, 0.25, 499, 5);

            let ctx = SitRepContext {
                dataset: &ds,
                suppression_policy: mask,
                rt_series: &rt_series,
                ab_series: &ab_series,
                nowcast: &nowcast,
                backtest: &backtest,
                spatial_clusters: &spatial_clusters,
            };
            let html = render_standalone_sitrep_html(&ctx);
            std::fs::write(out_html, html)?;
            println!(
                "Exported standalone HTML5 + inline vector SVG Epidemiological Situation Report to {out_html}"
            );
        }
        "dp-release" => {
            let vid = args
                .get(1)
                .cloned()
                .unwrap_or_else(|| "cdc_measles_2022_2026_v1".to_string());
            let eps: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.50);
            let cap: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(4.0);
            let ds = store
                .load_vintage(&vid)?
                .ok_or_else(|| EpiError::Schema(format!("Vintage '{vid}' not found")))?;
            let cfg = DpReleaseConfig {
                epsilon: eps,
                sensitivity_l1: 1.0,
                max_budget_cap: cap,
                deterministic_seed: None,
            };
            let out = release_differentially_private_series(&store, &ds, &cfg)?;
            store.commit_vintage(&out.noisy_dataset)?;
            if let Some(csv_out) = args.get(4) {
                export_dataset_csv(&out.noisy_dataset, csv_out, SuppressionPolicy::Unmasked)?;
            }
            println!(
                "Discrete Laplace DP Release: '{}' -> '{}' | eps={:.3} | cum_eps={:.3}/{:.2} (rem={:.3}) | Theory_RMSE={:.2} | Emp_RMSE={:.2} (MAE={:.2}) | AuditSHA256={}",
                vid,
                out.noisy_dataset.vintage_id,
                out.epsilon_spent,
                out.cumulative_epsilon_for_vintage,
                cap,
                out.remaining_budget,
                out.theoretical_rmse,
                out.empirical_rmse,
                out.empirical_mae,
                out.audit_entry_hash
            );
        }
        "audit-verify" => {
            let limit: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(20);
            let ver = verify_audit_chain(&store)?;
            println!(
                "Tamper-Evident SHA-256 Audit Chain Verification: VALID={} | TotalEntries={} | TipHash={}",
                ver.is_valid, ver.total_entries, ver.tip_hash_hex
            );
            if let Some(reason) = ver.failure_reason {
                println!("  FAILURE: {reason}");
            }
            for e in list_audit_events(&store, limit)? {
                println!(
                    "  #{:<4} {:<16} {:<26} prev={}.. sha256={}.. | {}",
                    e.seq_no,
                    e.action_kind,
                    e.target_ref,
                    &e.prev_hash_hex[0..10],
                    &e.entry_hash_hex[0..10],
                    e.detail_summary
                );
            }
        }
        "pprl-bench" => {
            let num_u: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(800);
            let use_block = args.get(2).map(|s| s != "cartesian").unwrap_or(true);
            let t0 = Instant::now();
            let raw = generate_synthetic_pprl_scale_benchmark(num_u);
            let mut cfg = PprlConfig::default();
            cfg.use_hash_blocking = use_block;
            let enc: Vec<EncodedPprlRecord> =
                raw.iter().map(|r| encode_pprl_record(r, &cfg)).collect();
            let sum = link_and_deduplicate_pprl(&enc, &cfg);
            let ms = t0.elapsed().as_secs_f64() * 1000.0;
            println!(
                "Hash-Blocked 512-Bit Bloom PPRL Scale Benchmark: {} encounters ({} true patients) | Blocking={} | CartesianPairs={} -> EvaluatedPairs={} ({:.1}% pruned) | LinkedPatients={} | Elapsed={:.2} ms",
                sum.total_input_records,
                num_u,
                use_block,
                sum.total_cartesian_pairs,
                sum.evaluated_candidate_pairs,
                sum.blocking_pair_reduction_pct,
                sum.unique_linked_entities,
                ms
            );
        }
        "lineage-fit" => {
            let table = if let Some(csv_path) = args.get(1).filter(|s| *s != "builtin") {
                let txt = std::fs::read_to_string(csv_path)?;
                parse_lineage_counts_csv(&txt, csv_path)?
            } else {
                builtin_respiratory_lineage_table()
            };
            let mode = match args.get(2).map(|s| s.to_ascii_lowercase()).as_deref() {
                Some("spline") => LineageBasisMode::NaturalCubicSpline,
                _ => LineageBasisMode::LogLinearFitness,
            };
            let forecast_wks: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(6);
            let cfg = LineageLogisticsConfig {
                basis_mode: mode,
                forecast_weeks: forecast_wks,
                ..Default::default()
            };
            let res = fit_multinomial_lineage_logistics(&table, &cfg)?;

            println!(
                "Multinomial Lineage Dynamics: {} | Mode={} | Ref={} | LogLik={:.1} | AIC={:.1} | Converged={} ({} iters)",
                table.title,
                res.basis_mode.label(),
                res.reference_clade,
                res.log_likelihood,
                res.aic,
                res.converged,
                res.iterations_used
            );
            println!(
                "{:<12} {:>8} {:>10} {:>10} {:>24} {:>12} {:>18}",
                "CLADE", "SEQ_N", "OBS_NOW", "FCST_END", "GROWTH s/wk [95% CI]", "GEN_ADV_%", "VARIANT_RT"
            );
            for c in &res.clade_summaries {
                let g_str = if c.is_reference {
                    "0.000 (ref)".to_string()
                } else {
                    format!(
                        "{:+.3} [{:+.3}, {:+.3}]",
                        c.growth_advantage_per_week, c.growth_ci_low_95, c.growth_ci_high_95
                    )
                };
                let rt_str = c
                    .latest_variant_rt
                    .map(|r| format!("{:.2} [{:.2},{:.2}]", r.mean, r.cri_low_95, r.cri_high_95))
                    .unwrap_or_else(|| "--".to_string());
                println!(
                    "{:<12} {:>8.0} {:>9.1}% {:>9.1}% {:>24} {:>+11.1}% {:>18}",
                    c.clade_name,
                    c.total_sequences,
                    c.latest_observed_share * 100.0,
                    c.forecast_end_share * 100.0,
                    g_str,
                    c.multiplicative_advantage_pct,
                    rt_str
                );
            }
        }
        "phylo-tree" => {
            let isolates = if let (Some(fa_path), Some(meta_path)) =
                (args.get(1).filter(|s| *s != "builtin"), args.get(2))
            {
                let fa_txt = std::fs::read_to_string(fa_path)?;
                let meta_txt = std::fs::read_to_string(meta_path)?;
                parse_fasta_and_metadata(&fa_txt, &meta_txt)?
            } else {
                let (fa_txt, meta_txt) = builtin_outbreak_fasta_and_metadata();
                parse_fasta_and_metadata(&fa_txt, &meta_txt)?
            };

            let cfg = OutbreakTreeConfig::default();
            let res = reconstruct_outbreak_tree(&isolates, &cfg)?;

            println!(
                "Joint Genomic + Symptom-Onset Outbreak Tree: {} retained ({} bp) | Introductions={} | UnsampledIntermediates={} | Mean_Ri={:.2} | LogLik={:.1}",
                res.retained_cases,
                res.alignment_length,
                res.total_importations,
                res.total_unsampled_intermediates,
                res.mean_case_ri,
                res.total_tree_log_likelihood
            );
            println!("\nNewick Outbreak Forest:\n{}", res.newick_tree);

            if let Some(nwk_out) = args.get(3) {
                std::fs::write(nwk_out, &res.newick_tree)?;
                println!("Wrote Newick tree to {nwk_out}");
            }
            if let Some(dot_out) = args.get(4) {
                let dot = export_outbreak_graphviz_dot(&res);
                std::fs::write(dot_out, dot)?;
                println!("Wrote Graphviz DOT transmission network to {dot_out}");
            }
        }
        "dlnm-fit" => {
            let ds = if let Some(csv_path) = args.get(1).filter(|s| *s != "builtin") {
                let txt = std::fs::read_to_string(csv_path)?;
                parse_dlnm_exposure_csv(&txt, csv_path, "exposure_units")?
            } else {
                builtin_wildfire_pm25_dlnm_dataset()
            };
            let max_lag: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(14);
            let exp_basis = match args.get(3).map(|s| s.to_ascii_lowercase()).as_deref() {
                Some("threshold") => ExposureBasisKind::QuadraticThreshold,
                _ => ExposureBasisKind::NaturalCubicSpline3Df,
            };
            let cfg = DlnmConfig {
                exposure_basis: exp_basis,
                max_lag,
                ..Default::default()
            };
            let res = fit_dlnm_quasi_poisson(&ds, &cfg)?;

            println!(
                "Gasparrini Cross-Basis DLNM: {} ({} steps, max_lag={}) | Ref={:.1} | MME={:.1} {} | P95_CumRR={:.3} [{:.3}, {:.3}] | AN={:.1} ({:.2}%)",
                res.title,
                res.effective_rows_fitted,
                res.max_lag,
                res.reference_exposure,
                res.minimum_risk_exposure,
                res.exposure_unit,
                res.p95_cumulative_rr,
                res.p95_rr_ci_low_95,
                res.p95_rr_ci_high_95,
                res.attributable_events,
                res.attributable_fraction_pct
            );
        }
        "excess-mortality" => {
            let ds = if let Some(csv_path) = args.get(1).filter(|s| *s != "builtin") {
                let txt = std::fs::read_to_string(csv_path)?;
                parse_weekly_mortality_csv(&txt, csv_path)?
            } else {
                builtin_all_cause_weekly_mortality_dataset()
            };
            let mode = match args.get(2).map(|s| s.to_ascii_lowercase()).as_deref() {
                Some("wmd") | Some("full") => MortalityBaselineMode::WorldMortalityFullYear,
                _ => MortalityBaselineMode::EuroMomoSpringAutumn,
            };
            let eval_wks: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(52);
            let cfg = ExcessMortalityConfig {
                baseline_mode: mode,
                evaluation_window_weeks: eval_wks,
                ..Default::default()
            };
            let res = compute_euromomo_excess_mortality(&ds, &cfg)?;

            println!(
                "Excess Mortality ({}) — {} | EvalExcess={:+.1} | MeanPScore={:+.1}% | YPLL75={:.0} | Alerts(>2Z)={} | Extreme(>4Z)={}",
                res.baseline_mode.label(),
                res.title,
                res.evaluation_cumulative_excess,
                res.evaluation_mean_p_score_pct,
                res.evaluation_ypll_75,
                res.substantial_alert_weeks,
                res.extreme_alert_weeks
            );
        }
        "import-hl7v2" => {
            let pprl_cfg = PprlConfig::default();
            let dedup = args
                .get(2)
                .map(|s| s != "raw" && s != "false")
                .unwrap_or(true);
            let out = if let Some(hl7_path) = args.get(1).filter(|s| *s != "builtin") {
                import_hl7_v2_file(hl7_path, &pprl_cfg, dedup)?
            } else {
                let raw = builtin_hl7_v2_elr_and_nssp_batch();
                parse_hl7_v2_batch_str(&raw, &pprl_cfg, dedup, "BuiltinBatch")?
            };
            if !out.dataset.points.is_empty() {
                store.commit_vintage(&out.dataset)?;
            }
            println!(
                "HL7 v2.5.1 ELR/NSSP Ingest: {} messages -> {} positive/syndromic hits -> {} PPRL-deduplicated incident cases across {} MMWR weeks (vintage '{}')",
                out.total_messages_parsed,
                out.positive_or_syndromic_hits,
                out.deduplicated_incident_cases,
                out.dataset.points.len(),
                out.dataset.vintage_id
            );
        }
        "pprl-link" => {
            let mut cfg = PprlConfig::default();
            if let Some(d_thr) = args.get(1).and_then(|s| s.parse().ok()) {
                cfg.dice_match_threshold = d_thr;
            }
            if let Some(win) = args.get(2).and_then(|s| s.parse().ok()) {
                cfg.episode_window_days = win;
            }
            let raw = builtin_pprl_hospital_lab_records();
            let enc: Vec<EncodedPprlRecord> =
                raw.iter().map(|r| encode_pprl_record(r, &cfg)).collect();
            let sum = link_and_deduplicate_pprl(&enc, &cfg);
            println!(
                "512-Bit HMAC-SHA256 Bigram Bloom PPRL: {} records -> {} unique linked patients -> {} incident episodes ({} duplicates collapsed, EM prior={:.2}%)",
                sum.total_input_records,
                sum.unique_linked_entities,
                sum.incident_episodes_retained,
                sum.duplicate_encounters_collapsed,
                sum.em_prior_match_rate * 100.0
            );
        }
        "sero-foi" => {
            let ds = if let Some(csv_path) = args.get(1).filter(|s| *s != "builtin") {
                let txt = std::fs::read_to_string(csv_path)?;
                parse_serosurvey_csv(&txt, csv_path)?
            } else {
                builtin_age_serosurvey_dataset()
            };
            let kind = match args.get(2).map(|s| s.to_ascii_lowercase()).as_deref() {
                Some("irreversible") | Some("irrev") => CatalyticModelKind::IrreversibleCatalytic,
                _ => CatalyticModelKind::ReversibleWaningCatalytic,
            };
            let fit = fit_catalytic_seroprevalence(&ds, kind)?;
            println!(
                "Sero-Catalytic Model ({}): FOI lambda={:.4}/yr [{:.4}, {:.4}] | Waning omega={:.4}/yr | Implied R0={:.2} | Weighted Immune={:.1}% | Empirical S_eff={:.1}% | AIC={:.1}",
                fit.model_kind.label(),
                fit.foi_lambda_per_year,
                fit.foi_ci_low_95,
                fit.foi_ci_high_95,
                fit.waning_omega_per_year,
                fit.implied_endemic_r0,
                fit.weighted_seropositive_fraction * 100.0,
                fit.effective_susceptible_fraction * 100.0,
                fit.aic
            );
        }
        "tnd-ve" => {
            let ds = builtin_tnd_vaccine_study_dataset();
            let res = estimate_tnd_vaccine_effectiveness(&ds)?;
            println!(
                "Test-Negative Design (TND) VE ({} cases, {} controls): Crude VE={:.1}% [{:.1}%, {:.1}%] -> Adjusted VE={:.1}% [{:.1}%, {:.1}%] (Bias removed={:+.1} pp)",
                res.total_cases,
                res.total_controls,
                res.crude_unadjusted_ve_pct,
                res.crude_ve_ci_low_95,
                res.crude_ve_ci_high_95,
                res.adjusted_overall_ve_pct,
                res.adjusted_ve_ci_low_95,
                res.adjusted_ve_ci_high_95,
                res.confounding_bias_removed_pp
            );
        }
        "spatial-scan" => {
            let records = builtin_pertussis_spatial_records();
            let clusters = compute_kulldorff_spatial_clusters(&records, true);
            println!("Single-State Kulldorff Poisson Excess LLR (Top 10, 2025 Pertussis):");
            for c in clusters.into_iter().take(10) {
                println!(
                    "  {:<18} (FIPS {})  Obs={:>5.0}  Exp={:>6.1}  RR={:>5.2}x  LLR={:>8.1}",
                    c.name, c.fips, c.observed_cases, c.expected_cases, c.relative_risk, c.kulldorff_llr
                );
            }
        }
        "zone-scan" => {
            let reps: usize = args
                .get(2)
                .and_then(|s| s.parse().ok())
                .unwrap_or(999);
            let zones = if let Some(geo_path) = args.get(1).filter(|s| *s != "builtin") {
                let feats = parse_geojson_centroids_file(geo_path)?;
                feats
                    .into_iter()
                    .map(|f| SpatialScanInputZone {
                        id: f.geoid,
                        name: f.name,
                        lat_deg: f.lat_deg,
                        lon_deg: f.lon_deg,
                        observed_cases: f.cases.unwrap_or(0.0),
                        population: f.population.unwrap_or(100_000.0),
                    })
                    .collect()
            } else {
                build_us_state_scan_zones(&builtin_pertussis_spatial_records(), true)
            };

            let t0 = Instant::now();
            let clusters = run_kulldorff_circular_scan(&zones, 0.25, reps, 6);
            let ms = t0.elapsed().as_secs_f64() * 1000.0;
            println!(
                "Kulldorff Scoped-Parallel Multi-Zone Scan (M={} Monte Carlo Replications in {:.2} ms):",
                reps, ms
            );
            for (rank, c) in clusters.iter().enumerate() {
                println!(
                    "#{:<2} Center: {:<14} (R={:>6.0} km, {} zones) | Obs={:.0}, Exp={:.1}, RR={:.2}x, LLR={:.1}, p={:.4} | Members: {}",
                    rank + 1,
                    c.center_name,
                    c.radius_km,
                    c.member_names.len(),
                    c.observed_cases,
                    c.expected_cases,
                    c.relative_risk,
                    c.kulldorff_llr,
                    c.monte_carlo_p_value,
                    c.member_names.join(", ")
                );
            }
        }
        "fetch-delphi" => {
            let region = args.get(1).map(String::as_str).unwrap_or("nat");
            let wks = args.get(2).map(String::as_str).unwrap_or("202340-202645");
            let ds = CmuDelphiClient::default().fetch_fluview_ilinet(
                &mut store, region, wks, true, 21_600,
            )?;
            println!(
                "Committed CMU Delphi FluView ILINet vintage '{}' ({} weeks)",
                ds.vintage_id,
                ds.points.len()
            );
        }
        "fetch-nndss" => {
            let disease = args.get(1).map(String::as_str).unwrap_or("Pertussis");
            let state_filter = args.get(2).map(String::as_str).unwrap_or("");
            let out = CdcSodaClient::default().fetch_with_lazy_cache(
                &mut store,
                SodaDatasetKind::NndssWeekly,
                disease,
                state_filter,
                800,
                3600 * 6,
            )?;
            println!(
                "Committed vintage '{}' ({} weekly steps from {} raw SODA rows; stale_fallback={})",
                out.dataset.vintage_id,
                out.dataset.points.len(),
                out.raw_row_count,
                out.used_stale_cache
            );
        }
        "fetch-nwss" => {
            let state_filter = args.get(1).map(String::as_str).unwrap_or("Oregon");
            let out = CdcSodaClient::default().fetch_with_lazy_cache(
                &mut store,
                SodaDatasetKind::NwssWastewater,
                "",
                state_filter,
                800,
                3600 * 6,
            )?;
            println!(
                "Committed NWSS wastewater vintage '{}' ({} dates across {} states)",
                out.dataset.vintage_id,
                out.dataset.points.len(),
                out.state_values_by_fips.len()
            );
        }
        "fetch-faers" => {
            let drug = args.get(1).map(String::as_str).unwrap_or("oseltamivir");
            let out = OpenFdaClient::default().fetch_faers_with_cache(&mut store, drug, 3600 * 12)?;
            println!(
                "Committed openFDA FAERS vintage '{}' ({} weekly points, {} top MedDRA reactions)",
                out.weekly_series.vintage_id,
                out.weekly_series.points.len(),
                out.top_reactions.len()
            );
        }
        "fetch-who" => {
            let indicator = args.get(1).map(String::as_str).unwrap_or("WHS3_62");
            let iso3 = args.get(2).map(String::as_str).unwrap_or("USA");
            let ds =
                WhoGhoClient::default().fetch_indicator_series(&mut store, indicator, iso3, 86_400)?;
            println!(
                "Committed WHO GHO vintage '{}' ({} periods)",
                ds.vintage_id,
                ds.points.len()
            );
        }
        "import-fhir" => {
            let path = args.get(1).ok_or_else(|| {
                EpiError::Schema("Usage: epiledger import-fhir <bundle.json> [code_filter]".to_string())
            })?;
            let flt = args.get(2).map(String::as_str);
            let sum = import_fhir_r4_bundle_file(path, flt)?;
            store.commit_vintage(&sum.dataset)?;
            println!(
                "Imported HL7 FHIR R4 Bundle into vintage '{}' ({} weeks from {} FHIR entries)",
                sum.dataset.vintage_id,
                sum.dataset.points.len(),
                sum.total_resources_scanned
            );
        }
        "export-fhir" => {
            let vid = args.get(1).ok_or_else(|| {
                EpiError::Schema("Usage: epiledger export-fhir <vintage_id> <out_bundle.json>".to_string())
            })?;
            let out_path = args.get(2).ok_or_else(|| {
                EpiError::Schema("Usage: epiledger export-fhir <vintage_id> <out_bundle.json>".to_string())
            })?;
            let ds = store
                .load_vintage(vid)?
                .ok_or_else(|| EpiError::Schema(format!("Vintage '{vid}' not found")))?;
            let bundle = export_timeseries_as_fhir_r4_bundle(&ds, "14189004")?;
            std::fs::write(out_path, bundle)?;
            println!("Exported HL7 FHIR R4 Bundle JSON to {out_path}");
        }
        "survey-xpt" => {
            let path = args.get(1).ok_or_else(|| {
                EpiError::Schema(
                    "Usage: epiledger survey-xpt <file.XPT> <stratum_col> <psu_col> <weight_col> <outcome_col>"
                        .to_string(),
                )
            })?;
            let s_col = args.get(2).map(String::as_str).unwrap_or("SDMVSTRA");
            let p_col = args.get(3).map(String::as_str).unwrap_or("SDMVPSU");
            let w_col = args.get(4).map(String::as_str).unwrap_or("WTMEC2YR");
            let y_col = args.get(5).ok_or_else(|| {
                EpiError::Schema("Specify outcome variable name in 5th argument".to_string())
            })?;
            let bytes = std::fs::read(path)?;
            let xpt = parse_sas_xport_bytes(&bytes, 100_000)?;
            let obs = extract_survey_obs_from_xpt(&xpt, s_col, p_col, w_col, y_col)?;
            let est = compute_taylor_linearized_mean(&obs, LonelyPsuStrategy::GrandMeanCenter)
                .ok_or_else(|| EpiError::Schema("Insufficient valid weighted survey rows".to_string()))?;
            println!(
                "Taylor-Series Linearized Survey Estimate ({} on {}): Mean={:.4}, SE={:.4}, 95% CI=[{:.4}, {:.4}], df={}, DEFF={:.2}",
                y_col, path, est.weighted_mean_or_prop, est.standard_error, est.ci_low_95, est.ci_high_95, est.degrees_of_freedom, est.design_effect_deff
            );
        }
        "faers-prr" => {
            let a: u64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(180);
            let b: u64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4_820);
            let c: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1_200);
            let d: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(993_800);
            let sig = compute_disproportionality_2x2("TargetDrug", "TargetReaction", a, b, c, d);
            println!(
                "2x2 FAERS Disproportionality (a={a}, b={b}, c={c}, d={d}): PRR={:.2} (chi2={:.1}), ROR={:.2} [{:.2}, {:.2}], BCPNN IC={:.2} (IC025={:.2}), EvansSignal={}, BcpnnSignal={}",
                sig.prr, sig.prr_chi2_yates, sig.ror, sig.ror_ci_low_95, sig.ror_ci_high_95, sig.bcpnn_ic, sig.bcpnn_ic_025, sig.flagged_evans_prr, sig.flagged_bcpnn_ic
            );
        }
        "import-csv" => {
            let path = args.get(1).ok_or_else(|| {
                EpiError::Schema("Usage: epiledger import-csv <path.csv>".to_string())
            })?;
            let ds = import_delimited_timeseries(path, None, &ColumnMapping::default())?;
            store.commit_vintage(&ds)?;
            println!(
                "Imported CSV into SQLite vintage '{}' ({} time steps, checksum {})",
                ds.vintage_id,
                ds.points.len(),
                ds.checksum_hex
            );
        }
        "export-csv" => {
            let vid = args.get(1).ok_or_else(|| {
                EpiError::Schema("Usage: epiledger export-csv <vintage_id> <out.csv>".to_string())
            })?;
            let out_path = args.get(2).ok_or_else(|| {
                EpiError::Schema("Usage: epiledger export-csv <vintage_id> <out.csv>".to_string())
            })?;
            let ds = store
                .load_vintage(vid)?
                .ok_or_else(|| EpiError::Schema(format!("Vintage '{vid}' not found")))?;
            export_dataset_csv(&ds, out_path, SuppressionPolicy::Unmasked)?;
            println!("Exported {} rows to {out_path}", ds.points.len());
        }
        "offspring-k" => {
            let src_arg = args.get(1).map(String::as_str).unwrap_or("builtin");
            let mode_arg = args.get(2).map(String::as_str).unwrap_or("all");
            let mode = if mode_arg.eq_ignore_ascii_case("ztnb")
                || mode_arg.eq_ignore_ascii_case("truncated")
            {
                ZeroTruncationMode::ZeroTruncatedClustersOnly
            } else {
                ZeroTruncationMode::StandardAllCases
            };
            let counts: Vec<u32> = if src_arg.eq_ignore_ascii_case("builtin") {
                let recs = builtin_superspreading_outbreak_records();
                let net = analyze_contact_tracing_network(&recs, &OffspringFitConfig::default())?;
                net.case_metrics.iter().map(|c| c.secondary_offspring).collect()
            } else if Path::new(src_arg).exists() {
                let raw = std::fs::read_to_string(src_arg)?;
                raw.split(|c: char| c == ',' || c.is_whitespace() || c == ';')
                    .filter(|s| !s.trim().is_empty())
                    .filter_map(|s| s.trim().parse::<u32>().ok())
                    .collect()
            } else {
                src_arg
                    .split(|c: char| c == ',' || c.is_whitespace() || c == ';')
                    .filter(|s| !s.trim().is_empty())
                    .map(|s| {
                        s.trim().parse::<u32>().map_err(|_| {
                            EpiError::Schema(format!("Invalid integer offspring count '{s}'"))
                        })
                    })
                    .collect::<Result<Vec<u32>, EpiError>>()?
            };

            let cfg = OffspringFitConfig {
                mode,
                ..Default::default()
            };
            let res = fit_negative_binomial_offspring(&counts, &cfg)?;
            let hi_str = if res.k_ci_high_95.is_finite() {
                format!("{:.3}", res.k_ci_high_95)
            } else {
                "inf".to_string()
            };
            println!(
                "Negative Binomial Offspring Fit (n={}, total_z={}, mode={:?}):",
                res.n_cases_evaluated, res.total_secondary_cases, res.mode
            );
            println!(
                "  R0 MLE = {:.3} [95% CI: {:.3}, {:.3}] | Dispersion k MLE = {:.3} [95% Profile CI: {:.3}, {}]",
                res.r0_mle, res.r0_ci_low_95, res.r0_ci_high_95, res.k_mle, res.k_ci_low_95, hi_str
            );
            println!(
                "  LRT vs Poisson = {:.2} (p = {:.3e}) | Top 20% cases -> {:.1}% tx (emp {:.1}%) | 80% tx by {:.1}% cases",
                res.lrt_statistic,
                res.lrt_p_value,
                res.prop_transmission_top_20_theoretical * 100.0,
                res.prop_transmission_top_20_empirical * 100.0,
                res.cases_driving_80_pct_transmission * 100.0
            );
            println!(
                "  Single-Intro Extinction q_ext = {:.1}% | 99% Poisson SSE Threshold Z > {} ({} observed SSEs)",
                res.single_intro_extinction_prob * 100.0,
                res.sse_99_threshold,
                res.observed_sse_count
            );
            println!(
                "{:<8} {:>10} {:>12} {:>14} {:>14}",
                "Z_COUNT", "OBS_CASES", "OBS_PROB_%", "NEGBIN_PMF_%", "POISSON_PMF_%"
            );
            for b in res.histogram.iter().take(16) {
                println!(
                    "{:<8} {:>10} {:>11.2}% {:>13.2}% {:>13.2}%",
                    b.offspring_count,
                    b.observed_cases,
                    b.observed_prob * 100.0,
                    b.fitted_nb_prob * 100.0,
                    b.fitted_poisson_prob * 100.0
                );
            }
        }
        "contact-network" => {
            let src_arg = args.get(1).map(String::as_str).unwrap_or("builtin");
            let dot_out = args.get(2).map(String::as_str);
            let recs = if src_arg.eq_ignore_ascii_case("builtin") {
                builtin_superspreading_outbreak_records()
            } else {
                let raw = std::fs::read_to_string(src_arg)?;
                parse_contact_tracing_csv(&raw)?
            };

            let net = analyze_contact_tracing_network(&recs, &OffspringFitConfig::default())?;
            println!(
                "Bipartite Contact-Tracing Network: {} cases, {} venues, {} directed links, max gen={}, max k-core={}",
                net.total_cases,
                net.total_venues,
                net.total_transmission_links,
                net.max_generation_depth,
                net.max_bipartite_k_core
            );
            println!(
                "  Offspring R0 = {:.2} | Dispersion k = {:.3} [{:.3}, {}] | Top-20% Share = {:.1}% | Empirical SI = {:.2}d +/- {:.2}d (n={})",
                net.offspring_fit.r0_mle,
                net.offspring_fit.k_mle,
                net.offspring_fit.k_ci_low_95,
                if net.offspring_fit.k_ci_high_95.is_finite() {
                    format!("{:.3}", net.offspring_fit.k_ci_high_95)
                } else {
                    "inf".to_string()
                },
                net.offspring_fit.prop_transmission_top_20_theoretical * 100.0,
                net.serial_interval.mean_days,
                net.serial_interval.sd_days,
                net.serial_interval.pair_count
            );
            println!(
                "{:<20} {:<18} {:>8} {:>9} {:>8} {:>8} {:>12}",
                "VENUE_ID", "SETTING", "SEC_TX", "SHARE_%", "SAR_%", "K_CORE", "BETWEENNESS"
            );
            for v in &net.venue_metrics {
                println!(
                    "{:<20} {:<18} {:>8} {:>8.1}% {:>8} {:>8} {:>12.4}",
                    v.venue_id,
                    v.venue_category,
                    v.secondary_transmissions,
                    v.share_of_outbreak_transmissions * 100.0,
                    v.secondary_attack_rate
                        .map(|s| format!("{:.1}%", s * 100.0))
                        .unwrap_or_else(|| "--".to_string()),
                    v.bipartite_k_core,
                    v.bipartite_betweenness
                );
            }
            if let Some(out_p) = dot_out {
                let dot = export_bipartite_network_dot(&net, &recs);
                std::fs::write(out_p, dot)?;
                println!("Exported Graphviz bipartite network .dot to {out_p}");
            }
        }
        "wastewater-demix" => {
            let src_arg = args.get(1).map(String::as_str).unwrap_or("builtin");
            let dw_arg = args.get(2).map(String::as_str).unwrap_or("log2");
            let boot_reps: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(80);
            let depth_weighting = if dw_arg.eq_ignore_ascii_case("sqrt") {
                DepthWeighting::SqrtDepth
            } else if dw_arg.eq_ignore_ascii_case("uniform") {
                DepthWeighting::Uniform
            } else {
                DepthWeighting::Log2DepthPlusOne
            };

            let (lib, smps) = if src_arg.eq_ignore_ascii_case("builtin") {
                builtin_wastewater_barcode_and_samples()
            } else {
                let raw = std::fs::read_to_string(src_arg)?;
                parse_wastewater_snv_csv(&raw)?
            };

            let cfg = WastewaterDemixConfig {
                depth_weighting,
                bootstrap_replicates: boot_reps,
                ..Default::default()
            };
            let res = deconvolve_wastewater_timeseries(&lib, &smps, &cfg)?;
            println!(
                "Freyja L1 Simplex Wastewater Deconvolution ({} samples, {} lineages, weighting={:?}, boot={}):",
                res.samples.len(),
                res.lineage_names.len(),
                cfg.depth_weighting,
                cfg.bootstrap_replicates
            );
            print!("{:<16} {:<10} {:>9} {:>11}", "SAMPLE_ID", "PERIOD", "L1_LOSS", "VIRAL_GC_L");
            for name in &res.lineage_names {
                print!("  {:>18}", name);
            }
            println!("  {:>8}", "CRYPTIC");
            for s in &res.samples {
                print!(
                    "{:<16} {:<10} {:>9.4} {:>11.2e}",
                    s.sample_id, s.date_label, s.weighted_l1_loss, s.viral_copies_per_liter
                );
                for l in &s.lineages {
                    print!(
                        "  {:>6.1}% [{:>4.0}-{:>3.0}]",
                        l.abundance * 100.0,
                        l.ci_low_95 * 100.0,
                        l.ci_high_95 * 100.0
                    );
                }
                println!("  {:>8}", s.cryptic_sites.len());
            }
        }
        "antibiogram" => {
            let src_arg = args.get(1).map(String::as_str).unwrap_or("builtin");
            let dedup_arg = args.get(2).map(String::as_str).unwrap_or("dedup");
            let unit_flt = args.get(3).cloned().unwrap_or_default();
            let deduplication = if dedup_arg.eq_ignore_ascii_case("all_isolates")
                || dedup_arg.eq_ignore_ascii_case("raw")
            {
                IsolateDeduplicationPolicy::AllIsolatesUnfiltered
            } else {
                IsolateDeduplicationPolicy::ClsiM39FirstIsolatePerPatient
            };

            let recs = if src_arg.eq_ignore_ascii_case("builtin") {
                builtin_hospital_amr_isolates()
            } else {
                let raw = std::fs::read_to_string(src_arg)?;
                parse_amr_isolates_csv(&raw)?
            };

            let cfg = AntibiogramConfig {
                deduplication,
                facility_unit_filter: unit_flt,
                ..Default::default()
            };
            let ab = compute_clsi_m39_antibiogram(&recs, &cfg)?;
            println!(
                "CLSI M39 Cumulative Antibiogram (unit={}, policy={:?}, included={} isolates, excluded_duplicates={}):",
                ab.facility_unit_filter,
                ab.policy,
                ab.total_included_isolates,
                ab.total_duplicates_excluded
            );
            println!(
                "{:<24} {:>8} {:>9} {:>8} {:>8}  {}",
                "ORGANISM", "N_DEDUP", "M39_>=30", "MDR_%", "XDR_%", "SENTINEL_RESISTANCE"
            );
            for o in &ab.organism_summaries {
                println!(
                    "{:<24} {:>8} {:>9} {:>7.1}% {:>7.1}%  {}: {:.1}% ({})",
                    o.organism,
                    o.deduplicated_isolates,
                    if o.meets_clsi_m39_n30 { "PASS" } else { "LOW_N<30" },
                    o.pct_mdr,
                    o.pct_xdr,
                    o.sentinel_phenotype_label,
                    o.pct_sentinel_resistant,
                    o.sentinel_resistant_count
                );
            }
            println!(
                "
{:<24} {:<6} {:>7} {:>8} {:>18} {:>14}  {}",
                "ORGANISM", "DRUG", "N_TEST", "%_SUSC", "WILSON_95%_CI", "MIC50/MIC90", "STEWARDSHIP"
            );
            for c in &ab.cells {
                let mic_s = match (c.mic50_ug_ml, c.mic90_ug_ml) {
                    (Some(m50), Some(m90)) => format!("{m50:.2}/{m90:.2}"),
                    _ => "--".to_string(),
                };
                let flag = if !c.meets_clsi_m39_n30 {
                    "CAUTION_N<30"
                } else if c.empiric_caution_below_threshold {
                    "AVOID_MONO_EMPIRIC"
                } else {
                    "EMPIRIC_OK"
                };
                println!(
                    "{:<24} {:<6} {:>7} {:>7.1}% {:>18} {:>14}  {}",
                    c.organism,
                    c.drug_code,
                    c.isolates_tested,
                    c.pct_susceptible,
                    format!("[{:.1}%, {:.1}%]", c.wilson_ci_low_95, c.wilson_ci_high_95),
                    mic_s,
                    flag
                );
            }
        }
        "capture-recapture" => {
            let src_arg = args.get(1).map(String::as_str).unwrap_or("builtin");
            if src_arg.eq_ignore_ascii_case("chapman") {
                let n10: u64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(87);
                let n01: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(27);
                let n11: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(57);
                let ch = compute_two_source_chapman("Source_A", "Source_B", n10, n01, n11);
                println!(
                    "Two-Source Chapman Capture-Recapture (n10={n10}, n01={n01}, n11={n11}, obs={}):",
                    ch.total_observed_unique
                );
                println!(
                    "  Chapman N_hat = {:.2} [95% CI: {:.2}, {:.2}] (Seber SE = {:.2}, unobs n00 = {:.2})",
                    ch.chapman_total_n_hat, ch.ci_low_95, ch.ci_high_95, ch.seber_se, ch.chapman_unobserved_n00
                );
                println!(
                    "  Completeness = {:.1}% | Sens(A) = {:.1}% | Sens(B) = {:.1}%",
                    ch.overall_completeness_pct, ch.source_a_sensitivity_pct, ch.source_b_sensitivity_pct
                );
            } else {
                let counts: ThreeSourceOverlapCounts = if src_arg.eq_ignore_ascii_case("builtin") {
                    builtin_three_source_tb_pertussis_crc()
                } else if Path::new(src_arg).exists() {
                    let raw = std::fs::read_to_string(src_arg)?;
                    parse_crc_case_lists_csv(&raw, "Lab_S1", "Hospital_S2", "Vital_S3")?
                } else {
                    let cells: Vec<u64> = src_arg
                        .split(',')
                        .filter_map(|s| s.trim().parse::<u64>().ok())
                        .collect();
                    if cells.len() != 7 {
                        return Err(EpiError::Schema(
                            "Provide 'builtin', a 4-col CSV file, 'chapman n10 n01 n11', or 7 comma-separated Venn counts n100,n010,n001,n110,n101,n011,n111".to_string(),
                        ));
                    }
                    ThreeSourceOverlapCounts {
                        s1_name: "Source_1".to_string(),
                        s2_name: "Source_2".to_string(),
                        s3_name: "Source_3".to_string(),
                        n_100: cells[0],
                        n_010: cells[1],
                        n_001: cells[2],
                        n_110: cells[3],
                        n_101: cells[4],
                        n_011: cells[5],
                        n_111: cells[6],
                    }
                };

                let crc = fit_three_source_loglinear_crc(&counts)?;
                let best = &crc.models[crc.bic_selected_index];
                println!(
                    "Three-Source Log-Linear Capture-Recapture (Observed = {}, BIC Optimal = {} {}):",
                    crc.total_observed, best.model_id, best.formula_label
                );
                println!(
                    "  BIC Optimal N_hat = {:.1} [95% CI: {:.1}, {:.1}] | Unobserved n000 = {:.1} | Completeness = {:.1}%",
                    best.total_n_hat, best.n_ci_low_95, best.n_ci_high_95, best.unobserved_n000_hat, best.overall_completeness_pct
                );
                println!(
                    "{:<12} {:<30} {:>3} {:>9} {:>9} {:>9} {:>10} {:>22} {:>8}",
                    "MODEL", "INTERACTIONS", "DF", "DEV_G2", "AIC", "BIC", "UNOBS_000", "TOTAL_N_HAT_[95%_CI]", "COMP_%"
                );
                for (i, m) in crc.models.iter().enumerate() {
                    let star = if i == crc.bic_selected_index { "*" } else { " " };
                    println!(
                        "{}{:<11} {:<30} {:>3} {:>9.2} {:>9.2} {:>9.2} {:>10.1} {:>22} {:>7.1}%",
                        star,
                        m.model_id,
                        m.formula_label,
                        m.residual_df,
                        m.deviance_g2,
                        m.aic,
                        m.bic,
                        m.unobserved_n000_hat,
                        format!("{:.1} [{:.1},{:.1}]", m.total_n_hat, m.n_ci_low_95, m.n_ci_high_95),
                        m.overall_completeness_pct
                    );
                }
            }
        }
        "syndromic-classify" => {
            let src_arg = args.get(1).map(String::as_str).unwrap_or("builtin");
            let query_arg = args
                .get(2)
                .map(String::as_str)
                .unwrap_or("(^fever^ OR ^cough^ OR ^j10^) ANDNOT (^immunization^ OR ^heat pack^)");
            let commit_flag = args.get(3).map(String::as_str).unwrap_or("");

            let encs = if src_arg.eq_ignore_ascii_case("builtin") {
                builtin_ed_chief_complaint_encounters()
            } else {
                let raw = std::fs::read_to_string(src_arg)?;
                parse_ed_encounters_csv(&raw)?
            };

            let batch = classify_ed_encounters_essence(&encs, query_arg, 5);
            println!(
                "CDC NSSP ESSENCE CCDD Classification ({} ED visits, {} matched >=1 syndrome, {} NegEx/Exclusions blocked, custom query hits={}):",
                batch.total_encounters,
                batch.encounters_matching_any_syndrome,
                batch.total_negations_blocked,
                batch.custom_query_matches
            );
            println!(
                "{:<20} {:<42} {:>8} {:>9} {:>14}",
                "CODE", "SYNDROME_DEFINITION", "MATCHES", "ED_SHARE", "NEGEX_BLOCKED"
            );
            for r in &batch.syndrome_summaries {
                println!(
                    "{:<20} {:<42} {:>8} {:>8.1}% {:>14}",
                    r.code, r.label, r.matched_visits, r.pct_of_ed_visits, r.negated_false_positives_blocked
                );
            }
            if commit_flag.eq_ignore_ascii_case("commit-ili") {
                let ds = build_syndrome_weekly_dataset(&batch, EssenceSyndromeId::IliRespiratory);
                store.commit_vintage(&ds)?;
                println!(
                    "Committed weekly NSSP ILI time series into SQLite vintage '{}' ({} MMWR weeks).",
                    ds.vintage_id,
                    ds.points.len()
                );
            }
        }
        "cfr-adjust" => {
            let src_arg = args.get(1).map(String::as_str).unwrap_or("builtin");
            let delay_mean: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(13.5);
            let cutoff_arg: Option<usize> = args.get(3).and_then(|s| s.parse().ok());

            let full_series = if src_arg.eq_ignore_ascii_case("builtin") {
                builtin_emerging_wave_cases_and_deaths()
            } else {
                let raw = std::fs::read_to_string(src_arg)?;
                parse_cfr_incidence_deaths_csv(&raw)?
            };
            let end = cutoff_arg
                .unwrap_or_else(|| if src_arg.eq_ignore_ascii_case("builtin") { 26 } else { full_series.len() })
                .clamp(3, full_series.len());

            let cfg = DelayCfrConfig {
                delay_mean,
                ..Default::default()
            };
            let res = compute_delay_adjusted_cfr(&full_series[..end], &cfg)?;
            println!(
                "Nishiura Right-Truncated Delay-Adjusted cCFR (steps=1..{end}, delay_mean={:.1}, delay_sd={:.1}):",
                cfg.delay_mean, cfg.delay_sd
            );
            println!(
                "  Naive Crude CFR = {:.2}% ({:.1} deaths / {:.0} raw cases) | Nishiura cCFR = {:.2}% [95% CI: {:.2}%, {:.2}%] (u_t = {:.1}, peak bias = {:.2}x)",
                res.final_crude_cfr_pct,
                res.final_cumulative_deaths,
                res.final_cumulative_cases,
                res.final_adjusted_ccfr_pct,
                res.final_adjusted_ci_low_95,
                res.final_adjusted_ci_high_95,
                res.final_resolved_denominator,
                res.peak_truncation_bias_factor
            );
            println!(
                "{:<10} {:>10} {:>12} {:>10} {:>11} {:>24} {:>8}",
                "PERIOD", "CUM_CASES", "RESOLVED_UT", "CUM_DEATH", "NAIVE_CFR", "ADJUSTED_CCFR_[95%_CI]", "BIAS"
            );
            let tail_start = res.trajectory.len().saturating_sub(12);
            for p in &res.trajectory[tail_start..] {
                println!(
                    "{:<10} {:>10.0} {:>12.1} {:>10.1} {:>10.2}% {:>24} {:>7.2}x",
                    p.period_label,
                    p.cumulative_cases,
                    p.resolved_denominator_ut,
                    p.cumulative_deaths,
                    p.crude_cfr_pct,
                    format!("{:.2}% [{:.2}%, {:.2}%]", p.delay_adjusted_ccfr_pct, p.adjusted_ci_low_95, p.adjusted_ci_high_95),
                    p.truncation_bias_factor
                );
            }
        }
        "hospital-surge" => {
            let mode_arg = args.get(1).map(String::as_str).unwrap_or("builtin");
            let icu_cap: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(55.0);
            let vent_cap: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(32.0);

            let cfg = HospitalSurgeConfig {
                staffed_icu_beds: icu_cap,
                available_ventilators: vent_cap,
                ..Default::default()
            };

            let surge = if mode_arg.eq_ignore_ascii_case("seir") {
                let vid = args
                    .get(4)
                    .cloned()
                    .unwrap_or_else(|| "cdc_measles_2022_2026_v1".to_string());
                let ds = store
                    .load_vintage(&vid)?
                    .ok_or_else(|| EpiError::Schema(format!("Vintage '{vid}' not found")))?;
                let counts: Vec<f64> = ds.points.iter().map(|p| p.count).collect();
                let fc = calibrate_and_project_seir(&counts, true, None, &SeirConfig::default());
                project_hospital_surge_from_seir(&fc, &cfg)?
            } else {
                let wave = builtin_emerging_wave_cases_and_deaths();
                let cases_d: Vec<f64> = wave.iter().map(|r| r.incident_cases * 6.0).collect();
                let low_d: Vec<f64> = cases_d.iter().map(|&v| v * 0.78).collect();
                let high_d: Vec<f64> = cases_d.iter().map(|&v| v * 1.24).collect();
                project_competing_risk_hospital_surge(&cases_d, &low_d, &high_d, &cfg)?
            };

            println!(
                "Competing-Risk Hospital / ICU / Ventilator Surge Projection (HFR={:.2}%, IFR={:.3}%, Total Admissions={:.1}, Deaths={:.1}):",
                surge.implied_hospital_fatality_pct,
                surge.implied_infection_fatality_pct,
                surge.total_projected_admissions,
                surge.total_projected_hospital_deaths
            );
            println!(
                "  Peak Acute Floor = {:.1} beds (Day {}, {:.0}% util) | Peak ICU = {:.1} beds (Day {}, {:.0}% util, deficit={:.1} bed-days) | Peak Ventilators = {:.1} (Day {}, {:.0}% util)",
                surge.peak_ward_census,
                surge.peak_ward_day,
                surge.peak_ward_utilization_pct,
                surge.peak_icu_census,
                surge.peak_icu_day,
                surge.peak_icu_utilization_pct,
                surge.icu_bed_days_deficit,
                surge.peak_vent_census,
                surge.peak_vent_day,
                surge.peak_vent_utilization_pct
            );
            println!(
                "{:<8} {:>12} {:>14} {:>22} {:>16} {:>12}",
                "DAY", "ADMISSIONS", "ACUTE_FLOOR", "ICU_CENSUS_[95%_PI]", "VENTILATORS", "CUM_DEATHS"
            );
            for p in surge.daily_projections.iter().step_by(5) {
                println!(
                    "Day {:<4} {:>12.1} {:>14.1} {:>22} {:>16.1} {:>12.1}",
                    p.day_index,
                    p.new_admissions_mean,
                    p.ward_census_mean,
                    format!("{:.1} [{:.1},{:.1}]", p.icu_census_mean, p.icu_census_low_95, p.icu_census_high_95),
                    p.vent_census_mean,
                    p.cumulative_deaths_mean
                );
            }
        }
        "inspect-xpt" => {
            let path = args.get(1).ok_or_else(|| {
                EpiError::Schema("Usage: epiledger inspect-xpt <nhanes_or_brfss.XPT>".to_string())
            })?;
            let bytes = std::fs::read(path)?;
            let xpt = parse_sas_xport_bytes(&bytes, 25)?;
            println!(
                "SAS XPORT Member: {} ({} variables, {} preview rows decoded)",
                xpt.member_name,
                xpt.variables.len(),
                xpt.rows.len()
            );
            for v in &xpt.variables {
                println!(
                    "  - {:<10} [{} len={}] {}",
                    v.name,
                    if v.is_char { "CHAR" } else { "IBM370_NUM" },
                    v.length,
                    v.label
                );
            }
        }
        other => {
            println!("EpiLedger Native Rust CLI — Unrecognized subcommand '{other}'");
            println!("Available commands:");
            println!("  epiledger                                                Launch native 16-tab wgpu/egui desktop workstation");
            println!("  epiledger tui                                            Launch headless 8-view ratatui SSH terminal dashboard");
            println!("  epiledger vintages                                       List bitemporal vintages in SQLite");
            println!("  epiledger analyze [vintage_id]                           Run Cori Rt + Farrington GLM on a vintage");
            println!("  epiledger nowcast [vintage_id] [delay_wk]                Run right-truncation reporting-delay nowcast");
            println!("  epiledger seir-forecast [vid] [weeks] [trajectories]     Run scoped-parallel RK4 + stochastic SEIR forecast");
            println!("  epiledger backtest-wis [vid] [horizon] [origins]         Run K=11 WIS backtest & Inverse-WIS ensemble stacking");
            println!("  epiledger export-sitrep [vid] [out.html] [mask]          Export standalone zero-JS HTML5/SVG Situation Report");
            println!("  epiledger dp-release [vid] [eps] [cap] [out.csv]         Release Discrete Laplace DP series & log SQLite budget");
            println!("  epiledger audit-verify [limit]                           Verify genesis-to-tip SHA-256 tamper-evident audit chain");
            println!("  epiledger pprl-bench [num_patients] [blocked|cartesian]  Benchmark LSH + Birth-Year hash-blocked PPRL scaling");
            println!("  epiledger offspring-k [z_list|csv|builtin] [all|ztnb]    Fit Negative Binomial R0 & dispersion k MLE + 80/20 share");
            println!("  epiledger contact-network [pairs.csv|builtin] [out.dot]  Analyze bipartite case-venue k-core, betweenness & SI");
            println!("  epiledger wastewater-demix [snv.csv|builtin] [log2|sqrt] Deconvolve sewage SNV barcodes via Freyja L1 simplex");
            println!("  epiledger antibiogram [amr.csv|builtin] [dedup|raw] [u]  Compute CLSI M39 cumulative antibiogram & MDR/XDR");
            println!("  epiledger capture-recapture [csv|chapman|builtin]        Fit 8-model 3-source log-linear & Chapman completeness");
            println!("  epiledger syndromic-classify [ed.csv|builtin] [query]    Classify NSSP ESSENCE chief complaints with NegEx");
            println!("  epiledger cfr-adjust [csv|builtin] [delay_d] [cutoff]    Compute Nishiura right-truncated delay-adjusted cCFR");
            println!("  epiledger hospital-surge [builtin|seir] [icu] [vent]     Project competing-risk acute, ICU & ventilator surge");
            println!("  epiledger lineage-fit [csv|builtin] [linear|spline]      Fit multinomial variant logistics & variant Rt");
            println!("  epiledger phylo-tree [fasta|builtin] [meta] [nwk]        Reconstruct joint SNP + onset outbreak tree");
            println!("  epiledger dlnm-fit [csv|builtin] [max_lag] [spline|thr]  Fit Gasparrini Cross-Basis DLNM lag-exposure model");
            println!("  epiledger excess-mortality [csv|builtin] [euromomo|wmd]  Compute Box-Cox Z-score excess mortality & YPLL-75");
            println!("  epiledger import-hl7v2 [file.hl7|builtin] [dedup|raw]    Parse HL7 v2.5.1 ELR/NSSP & deduplicate via Bloom PPRL");
            println!("  epiledger pprl-link [dice_cutoff] [episode_days]         Run 512-bit HMAC-SHA256 Bloom + Fellegi-Sunter EM linkage");
            println!("  epiledger sero-foi [csv|builtin] [rev|irrev]             Fit age-structured catalytic / reverse-catalytic FOI");
            println!("  epiledger tnd-ve                                         Estimate confounder-adjusted Test-Negative Design VE");
            println!("  epiledger spatial-scan                                   Rank single US states by Poisson excess LLR");
            println!("  epiledger zone-scan [geojson|builtin] [mc_reps]          Scoped-parallel Kulldorff multi-zone scan (up to M=9999)");
            println!("  epiledger fetch-nndss [disease]                          Live CDC SODA NNDSS fetch + lazy cache");
            println!("  epiledger fetch-nwss [state]                             Live CDC SODA NWSS wastewater fetch");
            println!("  epiledger fetch-delphi [region] [epiweeks]               Live CMU Delphi FluView ILINet fetch");
            println!("  epiledger fetch-faers [drug]                             Live openFDA FAERS weekly series fetch");
            println!("  epiledger fetch-who [ind] [iso3]                         Live WHO GHO OData indicator fetch");
            println!("  epiledger import-fhir <bundle.json> [code]               Import HL7 FHIR R4 Bundle into SQLite");
            println!("  epiledger export-fhir <vid> <bundle.json>                Export vintage as HL7 FHIR R4 Observation Bundle");
            println!("  epiledger survey-xpt <XPT> <s> <p> <w> <y>               Compute Taylor-linearized survey mean, SE & DEFF");
            println!("  epiledger faers-prr <a> <b> <c> <d>                      Compute 2x2 PRR, ROR, and BCPNN IC / IC025");
            println!("  epiledger import-csv <file.csv>                          Import CSV/TSV into SQLite vintage ledger");
            println!("  epiledger export-csv <vid> <out>                         Export vintage to CSV");
            println!("  epiledger inspect-xpt <file.XPT>                         Decode CDC NHANES/BRFSS SAS XPORT v5 file");
        }
    }

    Ok(())
}
