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

use crate::cartogram::US_CARTOGRAM_CELLS;
use crate::cli::extract_survey_obs_from_xpt;
use crate::crc_syndromic_tab::CrcSyndromicWorkbenchState;
use crate::dp_audit_tab::DpAuditWorkbenchState;
use crate::elr_pprl_sero_tab::ElrPprlSeroWorkbenchState;
use crate::env_mortality_tab::EnvMortalityWorkbenchState;
use crate::genomic_tab::GenomicWorkbenchState;
use crate::network_tab::SuperspreadingWorkbenchState;
use crate::scoring_sitrep_tab::ScoringSitrepWorkbenchState;
use crate::severity_surge_tab::SeveritySurgeWorkbenchState;
use crate::wastewater_amr_tab::WastewaterAmrWorkbenchState;
use eframe::egui::{self, Color32, RichText, Stroke, Vec2};
use egui_plot::{Bar, BarChart, Legend, Line, Plot, PlotPoints, Points};
use epi_analytics::{
    build_us_state_scan_zones, calibrate_and_project_seir, compute_aberrations, compute_cori_rt,
    compute_disproportionality_2x2, compute_kulldorff_spatial_clusters,
    compute_taylor_linearized_mean, nowcast_from_reporting_triangle, nowcast_parametric_delay,
    run_kulldorff_circular_scan, seir_pathogen_presets, AberrationConfig, AberrationMethod,
    AberrationSeries, ComplexSurveyEstimate, CoriRtConfig, CoriRtSeries, DisproportionalitySignal,
    LonelyPsuStrategy, MultiZoneCluster, NowcastSeries, ParametricNowcastConfig, ReportingTriangle,
    SeirConfig, SeirForecastResult, SpatialClusterScore, SpatialScanInputZone,
};
use epi_core::seed_data::builtin_fda_seed_recalls;
use epi_core::{
    builtin_measles_weekly_dataset, builtin_pertussis_spatial_records, pathogen_si_presets,
    CellValue, JurisdictionRecord, ReactionFacet, RecallRecord, SuppressionPolicy,
    TimeSeriesDataset, VintageMetadata,
};
use epi_ingest::{
    export_dataset_csv, export_dataset_json, export_timeseries_as_fhir_r4_bundle,
    import_delimited_timeseries, import_fhir_r4_bundle_file, import_json_timeseries,
    parse_geojson_centroids_file, parse_sas_xport_bytes, CdcSodaClient, CmuDelphiClient,
    ColumnMapping, OpenFdaClient, SasXportDataset, SodaDatasetKind, WhoGhoClient,
};
use epi_store::{BackfillRevisionDelta, EpiStore};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::mpsc::{channel, Receiver, Sender};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ActiveTab {
    Surveillance,
    SeirNowcast,
    SpatialScan,
    LiveConnectors,
    LocalImportXpt,
    SurveyAndPrr,
    GenomicPhylo,
    EnvAndMortality,
    ElrPprlSero,
    ForecastScoringSitrep,
    DpAuditAndHpc,
    SuperspreadingNet,
    WastewaterAndAmr,
    CrcAndSyndromic,
    SeverityAndSurge,
    VintageLedger,
}

enum WorkerMessage {
    SodaLoaded {
        dataset: TimeSeriesDataset,
        state_map: HashMap<String, f64>,
        raw_rows: usize,
        stale: bool,
    },
    DelphiLoaded(TimeSeriesDataset),
    FdaFaersLoaded {
        dataset: TimeSeriesDataset,
        reactions: Vec<ReactionFacet>,
        stale: bool,
    },
    FdaRecallsLoaded(Vec<RecallRecord>),
    WhoLoaded(TimeSeriesDataset),
    Error(String),
}

pub struct EpiLedgerGuiApp {
    db_path: PathBuf,
    store: EpiStore,
    active_tab: ActiveTab,
    status_banner: String,
    suppression_policy: SuppressionPolicy,

    // Active time-series & analytics cache
    active_dataset: TimeSeriesDataset,
    vintages: Vec<VintageMetadata>,
    window_weeks: usize,
    rt_cfg: CoriRtConfig,
    selected_si_preset: String,
    ab_cfg: AberrationConfig,
    cached_rt: CoriRtSeries,
    cached_ab: AberrationSeries,

    // Reporting-Delay Nowcasting & SEIR Compartmental Forecasting state
    nowcast_cfg: ParametricNowcastConfig,
    use_triangle_when_available: bool,
    feed_nowcast_into_seir: bool,
    seir_preset_id: String,
    seir_cfg: SeirConfig,
    cached_nowcast: NowcastSeries,
    cached_seir: SeirForecastResult,

    // Spatial state & Multi-Zone Kulldorff Scan
    jurisdictions: Vec<JurisdictionRecord>,
    spatial_use_2025: bool,
    selected_fips: String,
    cached_clusters: Vec<SpatialClusterScore>,
    max_zone_pop_frac: f64,
    cached_zone_clusters: Vec<MultiZoneCluster>,
    highlighted_cluster_members: HashSet<String>,

    // Live connector inputs & results
    soda_is_nwss: bool,
    soda_condition: String,
    soda_jurisdiction: String,
    soda_limit: usize,
    delphi_region: String,
    delphi_epiweeks: String,
    fda_drug: String,
    fda_reactions: Vec<ReactionFacet>,
    fda_recalls: Vec<RecallRecord>,
    who_indicator: String,
    who_iso3: String,
    worker_busy: bool,
    worker_tx: Sender<WorkerMessage>,
    worker_rx: Receiver<WorkerMessage>,

    // Local import / FHIR / GeoJSON / SAS XPT inspector
    import_path_input: String,
    fhir_code_filter: String,
    xpt_preview: Option<SasXportDataset>,

    // Complex Survey Linearization & Pharmacovigilance 2x2 PRR/ROR/IC state
    survey_stratum_col: String,
    survey_psu_col: String,
    survey_weight_col: String,
    survey_outcome_col: String,
    survey_estimate: Option<ComplexSurveyEstimate>,
    pv_drug_name: String,
    pv_reaction_name: String,
    pv_a: u64,
    pv_b: u64,
    pv_c: u64,
    pv_d: u64,
    pv_signal: DisproportionalitySignal,

    // Genomic Surveillance & Outbreak Phylodynamics state
    genomic: GenomicWorkbenchState,
    // Environmental DLNM & EuroMOMO Excess Mortality state
    env_mortality: EnvMortalityWorkbenchState,
    // HL7 v2.5.1 ELR/NSSP, 512-bit Bloom PPRL, & Sero/TND VE state
    elr_pprl_sero: ElrPprlSeroWorkbenchState,
    // K=11 WIS Backtesting, Inverse-WIS Ensembles, & HTML/SVG SitRep state
    scoring_sitrep: ScoringSitrepWorkbenchState,
    // Discrete Laplace DP Ledger, SHA-256 Audit Chain, & Blocked PPRL state
    dp_audit: DpAuditWorkbenchState,
    // Superspreading Offspring Dispersion (k MLE) & Bipartite Contact Network state
    superspreading_net: SuperspreadingWorkbenchState,
    // Freyja L1 Wastewater SNV Demixing & CLSI M39 Antibiograms state
    wastewater_amr: WastewaterAmrWorkbenchState,
    // 3-Source Log-Linear Capture-Recapture & NSSP ESSENCE Syndromic Classifier state
    crc_syndromic: CrcSyndromicWorkbenchState,
    // Right-Truncated Delay-Adjusted cCFR & Competing-Risk Hospital/ICU Surge state
    severity_surge: SeveritySurgeWorkbenchState,

    // Vintage diff
    diff_older_id: String,
    diff_newer_id: String,
    diff_results: Vec<BackfillRevisionDelta>,
}

impl EpiLedgerGuiApp {
    pub fn new(cc: &eframe::CreationContext<'_>, db_path: PathBuf) -> Self {
        Self::new_headless(&cc.egui_ctx, db_path)
    }

    pub fn new_headless(ctx: &egui::Context, db_path: PathBuf) -> Self {
        ctx.set_visuals(egui::Visuals::dark());

        let mut store = EpiStore::open_file(&db_path)
            .or_else(|_| EpiStore::open_in_memory())
            .expect("SQLite initialization");

        let seed_ds = builtin_measles_weekly_dataset();
        if store.list_vintages().map(|v| v.is_empty()).unwrap_or(true) {
            let _ = store.commit_vintage(&seed_ds);
        }

        let vintages = store.list_vintages().unwrap_or_default();
        let active_dataset = store
            .load_vintage(&seed_ds.vintage_id)
            .ok()
            .flatten()
            .unwrap_or(seed_ds);

        let rt_cfg = CoriRtConfig::default();
        let ab_cfg = AberrationConfig::default();
        let nowcast_cfg = ParametricNowcastConfig::default();
        let seir_cfg = SeirConfig::default();

        let counts: Vec<f64> = active_dataset.points.iter().map(|p| p.count).collect();
        let cached_rt = compute_cori_rt(&counts, &rt_cfg);
        let cached_ab = compute_aberrations(&counts, &ab_cfg);
        let cached_nowcast = nowcast_parametric_delay(&counts, &nowcast_cfg);
        let cached_seir = calibrate_and_project_seir(&counts, true, Some(&nowcast_cfg), &seir_cfg);

        let jurisdictions = builtin_pertussis_spatial_records();
        let cached_clusters = compute_kulldorff_spatial_clusters(&jurisdictions, true);
        let state_zones = build_us_state_scan_zones(&jurisdictions, true);
        let cached_zone_clusters = run_kulldorff_circular_scan(&state_zones, 0.25, 499, 6);
        let highlighted_cluster_members = cached_zone_clusters
            .first()
            .map(|c| c.member_ids.iter().cloned().collect())
            .unwrap_or_default();

        let pv_signal = compute_disproportionality_2x2(
            "Oseltamivir",
            "NEUROPSYCHIATRIC EVENT",
            180,
            4_820,
            1_200,
            993_800,
        );

        let (worker_tx, worker_rx) = channel();
        let scoring_sitrep = ScoringSitrepWorkbenchState::new(&active_dataset, &rt_cfg, &seir_cfg);
        let severity_surge = SeveritySurgeWorkbenchState::new(&cached_seir);
        let dp_audit = DpAuditWorkbenchState::new(&store, &active_dataset.vintage_id);
        let first_vid = vintages
            .first()
            .map(|v| v.vintage_id.clone())
            .unwrap_or_default();

        Self {
            db_path,
            store,
            active_tab: ActiveTab::Surveillance,
            status_banner: "Loaded 16-tab workbench with Right-Truncated cCFR, Competing-Risk ICU Surge, Capture-Recapture, ESSENCE & Freyja Demixing.".to_string(),
            suppression_policy: SuppressionPolicy::Unmasked,
            active_dataset,
            vintages,
            window_weeks: 104,
            rt_cfg,
            selected_si_preset: "measles".to_string(),
            ab_cfg,
            cached_rt,
            cached_ab,
            nowcast_cfg,
            use_triangle_when_available: false,
            feed_nowcast_into_seir: true,
            seir_preset_id: "measles".to_string(),
            seir_cfg,
            cached_nowcast,
            cached_seir,
            jurisdictions,
            spatial_use_2025: true,
            selected_fips: "41".to_string(),
            cached_clusters,
            max_zone_pop_frac: 0.25,
            cached_zone_clusters,
            highlighted_cluster_members,
            soda_is_nwss: false,
            soda_condition: "Pertussis".to_string(),
            soda_jurisdiction: String::new(),
            soda_limit: 600,
            delphi_region: "nat".to_string(),
            delphi_epiweeks: "202340-202645".to_string(),
            fda_drug: "oseltamivir".to_string(),
            fda_reactions: Vec::new(),
            fda_recalls: builtin_fda_seed_recalls(),
            who_indicator: "WHS3_62".to_string(),
            who_iso3: "USA".to_string(),
            worker_busy: false,
            worker_tx,
            worker_rx,
            import_path_input: String::new(),
            fhir_code_filter: String::new(),
            xpt_preview: None,
            survey_stratum_col: "SDMVSTRA".to_string(),
            survey_psu_col: "SDMVPSU".to_string(),
            survey_weight_col: "WTMEC2YR".to_string(),
            survey_outcome_col: "BMXBMI".to_string(),
            survey_estimate: None,
            pv_drug_name: "Oseltamivir".to_string(),
            pv_reaction_name: "NEUROPSYCHIATRIC EVENT".to_string(),
            pv_a: 180,
            pv_b: 4_820,
            pv_c: 1_200,
            pv_d: 993_800,
            pv_signal,
            genomic: GenomicWorkbenchState::new(),
            env_mortality: EnvMortalityWorkbenchState::new(),
            elr_pprl_sero: ElrPprlSeroWorkbenchState::new(),
            scoring_sitrep,
            dp_audit,
            superspreading_net: SuperspreadingWorkbenchState::new(),
            wastewater_amr: WastewaterAmrWorkbenchState::new(),
            crc_syndromic: CrcSyndromicWorkbenchState::new(),
            severity_surge,
            diff_older_id: first_vid.clone(),
            diff_newer_id: first_vid,
            diff_results: Vec::new(),
        }
    }

    fn recompute_analytics(&mut self) {
        let counts: Vec<f64> = self.active_dataset.points.iter().map(|p| p.count).collect();
        self.cached_rt = compute_cori_rt(&counts, &self.rt_cfg);
        self.cached_ab = compute_aberrations(&counts, &self.ab_cfg);
        self.recompute_seir_and_nowcast();
    }

    fn recompute_spatial(&mut self) {
        self.cached_clusters =
            compute_kulldorff_spatial_clusters(&self.jurisdictions, self.spatial_use_2025);
        let zones = build_us_state_scan_zones(&self.jurisdictions, self.spatial_use_2025);
        self.cached_zone_clusters =
            run_kulldorff_circular_scan(&zones, self.max_zone_pop_frac, 499, 6);
        self.highlighted_cluster_members = self
            .cached_zone_clusters
            .first()
            .map(|c| c.member_ids.iter().cloned().collect())
            .unwrap_or_default();
    }

    fn recompute_seir_and_nowcast(&mut self) {
        let counts: Vec<f64> = self.active_dataset.points.iter().map(|p| p.count).collect();

        if self.use_triangle_when_available && self.vintages.len() >= 2 {
            let loaded: Vec<TimeSeriesDataset> = self
                .vintages
                .iter()
                .rev()
                .filter_map(|m| self.store.load_vintage(&m.vintage_id).ok().flatten())
                .collect();
            if loaded.len() >= 2 {
                if let Some(tri) =
                    ReportingTriangle::from_vintages(&loaded, self.nowcast_cfg.max_delay_weeks)
                {
                    self.cached_nowcast =
                        nowcast_from_reporting_triangle(&tri, self.nowcast_cfg.overdispersion_phi);
                    let nc_counts: Vec<f64> = self
                        .cached_nowcast
                        .points
                        .iter()
                        .map(|p| {
                            if self.feed_nowcast_into_seir {
                                p.nowcast_mean
                            } else {
                                p.observed_count
                            }
                        })
                        .collect();
                    self.cached_seir =
                        calibrate_and_project_seir(&nc_counts, false, None, &self.seir_cfg);
                    return;
                }
            }
        }

        self.cached_nowcast = nowcast_parametric_delay(&counts, &self.nowcast_cfg);
        self.cached_seir = calibrate_and_project_seir(
            &counts,
            self.feed_nowcast_into_seir,
            Some(&self.nowcast_cfg),
            &self.seir_cfg,
        );
    }

    fn refresh_vintages(&mut self) {
        if let Ok(list) = self.store.list_vintages() {
            self.vintages = list;
        }
    }

    fn poll_worker(&mut self) {
        while let Ok(msg) = self.worker_rx.try_recv() {
            self.worker_busy = false;
            match msg {
                WorkerMessage::SodaLoaded {
                    dataset,
                    state_map,
                    raw_rows,
                    stale,
                } => {
                    let _ = self.store.commit_vintage(&dataset);
                    self.refresh_vintages();
                    if !dataset.points.is_empty() {
                        self.active_dataset = dataset;
                        self.recompute_analytics();
                    }
                    if !state_map.is_empty() {
                        for rec in &mut self.jurisdictions {
                            rec.active_custom_value = state_map.get(&rec.fips).copied();
                        }
                        self.recompute_spatial();
                    }
                    self.status_banner = format!(
                        "CDC SODA sync complete: {} rows parsed (stale_fallback={stale}).",
                        raw_rows
                    );
                }
                WorkerMessage::DelphiLoaded(dataset) => {
                    let _ = self.store.commit_vintage(&dataset);
                    self.refresh_vintages();
                    self.active_dataset = dataset;
                    self.rt_cfg.si_mean_days = 3.6;
                    self.rt_cfg.si_sd_days = 1.6;
                    self.recompute_analytics();
                    self.status_banner = format!(
                        "Loaded CMU Delphi FluView ILINet series '{}' ({} weeks).",
                        self.active_dataset.title,
                        self.active_dataset.points.len()
                    );
                }
                WorkerMessage::FdaFaersLoaded {
                    dataset,
                    reactions,
                    stale,
                } => {
                    let _ = self.store.commit_vintage(&dataset);
                    self.refresh_vintages();
                    self.active_dataset = dataset;
                    if let Some(top_rx) = reactions.first() {
                        self.pv_reaction_name = top_rx.term.clone();
                        self.pv_a = top_rx.count;
                        let sum_other: u64 =
                            reactions.iter().skip(1).map(|r| r.count).sum::<u64>().max(100);
                        self.pv_b = sum_other;
                        self.pv_signal = compute_disproportionality_2x2(
                            &self.fda_drug,
                            &self.pv_reaction_name,
                            self.pv_a,
                            self.pv_b,
                            self.pv_c,
                            self.pv_d,
                        );
                    }
                    self.fda_reactions = reactions;
                    self.rt_cfg.si_mean_days = 3.6;
                    self.rt_cfg.si_sd_days = 1.6;
                    self.recompute_analytics();
                    self.status_banner = format!(
                        "openFDA FAERS weekly series loaded ({} weeks, stale_fallback={stale}).",
                        self.active_dataset.points.len()
                    );
                }
                WorkerMessage::FdaRecallsLoaded(recalls) => {
                    let n = recalls.len();
                    self.fda_recalls = recalls;
                    self.status_banner = format!("Loaded {n} live openFDA drug recall records.");
                }
                WorkerMessage::WhoLoaded(dataset) => {
                    let _ = self.store.commit_vintage(&dataset);
                    self.refresh_vintages();
                    self.active_dataset = dataset;
                    self.recompute_analytics();
                    self.status_banner =
                        format!("Loaded WHO GHO series '{}'.", self.active_dataset.title);
                }
                WorkerMessage::Error(err) => {
                    self.status_banner = format!("Live connector warning: {err}");
                }
            }
        }
    }

    fn render_top_bar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.heading(
                RichText::new("EpiLedger")
                    .strong()
                    .color(Color32::from_rgb(96, 165, 250)),
            );
            ui.label(
                RichText::new(format!("[{}]", self.active_dataset.source_kind.as_str()))
                    .monospace()
                    .color(Color32::from_rgb(52, 211, 153)),
            );
            ui.separator();

            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::Surveillance,
                "1. Rt & Aberration",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::SeirNowcast,
                "2. SEIR & Nowcast",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::SpatialScan,
                "3. Multi-Zone Scan",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::LiveConnectors,
                "4. Live APIs",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::LocalImportXpt,
                "5. CSV / .XPT / FHIR",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::SurveyAndPrr,
                "6. Survey & PRR",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::GenomicPhylo,
                "7. Genomic & Phylo",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::EnvAndMortality,
                "8. Env & Mortality",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::ElrPprlSero,
                "9. ELR, PPRL & VE",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::ForecastScoringSitrep,
                "10. WIS & SitRep",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::DpAuditAndHpc,
                "11. DP, Audit & HPC",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::SuperspreadingNet,
                "12. Superspreading & Net",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::WastewaterAndAmr,
                "13. Wastewater & AMR",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::CrcAndSyndromic,
                "14. CRC & ESSENCE",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::SeverityAndSurge,
                "15. Severity & Surge",
            );
            ui.selectable_value(
                &mut self.active_tab,
                ActiveTab::VintageLedger,
                "16. Vintage Ledger",
            );

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let mut mask_bool = self.suppression_policy == SuppressionPolicy::NchsWonderLe9;
                if ui.checkbox(&mut mask_bool, "NCHS <=9 Mask").changed() {
                    self.suppression_policy = if mask_bool {
                        SuppressionPolicy::NchsWonderLe9
                    } else {
                        SuppressionPolicy::Unmasked
                    };
                }
            });
        });
    }

    fn render_surveillance_tab(&mut self, ui: &mut egui::Ui) {
        let total_pts = self.active_dataset.points.len();
        let start_idx = if self.window_weeks == 0 || self.window_weeks >= total_pts {
            0
        } else {
            total_pts - self.window_weeks
        };

        ui.columns(4, |cols| {
            let last_pt = self.active_dataset.points.last();
            cols[0].group(|ui| {
                ui.label(RichText::new("LATEST INCIDENCE").small().weak());
                if let Some(p) = last_pt {
                    let cv =
                        CellValue::from_count(Some(p.count), p.suppressed, self.suppression_policy);
                    ui.heading(RichText::new(cv.display_string(0)).monospace());
                    ui.label(RichText::new(&p.date_label).small().monospace());
                }
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("CORI INSTANTANEOUS Rt").small().weak());
                if let Some(Some(est)) = self.cached_rt.estimates.last() {
                    let col = if est.mean > 1.0 {
                        Color32::from_rgb(244, 63, 94)
                    } else {
                        Color32::from_rgb(16, 185, 129)
                    };
                    ui.heading(RichText::new(format!("{:.2}", est.mean)).monospace().color(col));
                    ui.label(
                        RichText::new(format!(
                            "95% CrI [{:.2}, {:.2}]",
                            est.cri_low_95, est.cri_high_95
                        ))
                        .small()
                        .monospace(),
                    );
                } else {
                    ui.heading("--");
                }
            });

            cols[2].group(|ui| {
                ui.label(RichText::new(self.ab_cfg.method.label()).small().weak());
                let win_alerts = self.cached_ab.points[start_idx..]
                    .iter()
                    .filter(|a| a.is_alert)
                    .count();
                ui.heading(
                    RichText::new(format!("{win_alerts} alerts"))
                        .monospace()
                        .color(Color32::from_rgb(245, 158, 11)),
                );
                ui.label(RichText::new(format!("z-alpha = {:.2}", self.ab_cfg.z_alpha)).small());
            });

            cols[3].group(|ui| {
                ui.label(RichText::new("ACTIVE VINTAGE").small().weak());
                ui.label(RichText::new(&self.active_dataset.title).strong());
                ui.label(
                    RichText::new(format!(
                        "Hash {} | {} total pts",
                        self.active_dataset.checksum_hex, total_pts
                    ))
                    .small()
                    .monospace(),
                );
            });
        });

        ui.add_space(6.0);

        let mut dirty = false;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("Window:");
                egui::ComboBox::from_id_source("win_combo")
                    .selected_text(match self.window_weeks {
                        0 => "All Points",
                        26 => "Last 26 Steps",
                        52 => "Last 52 Steps",
                        _ => "Last 104 Steps",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.window_weeks, 0, "All Points");
                        ui.selectable_value(&mut self.window_weeks, 104, "Last 104 Steps");
                        ui.selectable_value(&mut self.window_weeks, 52, "Last 52 Steps");
                        ui.selectable_value(&mut self.window_weeks, 26, "Last 26 Steps");
                    });

                ui.separator();
                ui.label("Aberration:");
                for m in [
                    AberrationMethod::FarringtonQuasiPoisson,
                    AberrationMethod::EarsC1,
                    AberrationMethod::EarsC2,
                    AberrationMethod::EarsC3Cusum,
                ] {
                    if ui
                        .selectable_label(self.ab_cfg.method == m, m.label())
                        .clicked()
                    {
                        self.ab_cfg.method = m;
                        dirty = true;
                    }
                }

                ui.separator();
                ui.label("SI Preset:");
                for preset in pathogen_si_presets() {
                    if ui
                        .selectable_label(self.selected_si_preset == preset.id, preset.id)
                        .clicked()
                    {
                        self.selected_si_preset = preset.id.to_string();
                        self.rt_cfg.si_mean_days = preset.mean_days;
                        self.rt_cfg.si_sd_days = preset.sd_days;
                        dirty = true;
                    }
                }
            });

            ui.horizontal_wrapped(|ui| {
                if ui
                    .add(
                        egui::Slider::new(&mut self.rt_cfg.si_mean_days, 1.5..=21.0)
                            .text("SI Mean (d)"),
                    )
                    .changed()
                {
                    self.selected_si_preset = "custom".to_string();
                    dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.rt_cfg.si_sd_days, 0.5..=10.0)
                            .text("SI SD (d)"),
                    )
                    .changed()
                {
                    self.selected_si_preset = "custom".to_string();
                    dirty = true;
                }
                if ui
                    .add(egui::Slider::new(&mut self.rt_cfg.tau_weeks, 1..=8).text("Cori tau (wk)"))
                    .changed()
                {
                    dirty = true;
                }
                if ui
                    .add(egui::Slider::new(&mut self.ab_cfg.z_alpha, 1.65..=3.50).text("z-alpha"))
                    .changed()
                {
                    dirty = true;
                }
            });
        });

        if dirty {
            self.recompute_analytics();
        }

        ui.add_space(6.0);

        let mut bars = Vec::new();
        let mut upper_pts = Vec::new();
        let mut exp_pts = Vec::new();
        let mut alert_pts = Vec::new();

        for (offset, idx) in (start_idx..total_pts).enumerate() {
            let x = offset as f64;
            let raw_c = self.active_dataset.points[idx].count;
            let disp_c = if self.suppression_policy == SuppressionPolicy::NchsWonderLe9
                && raw_c > 0.0
                && raw_c <= 9.0
            {
                0.0
            } else {
                raw_c
            };
            let ab = &self.cached_ab.points[idx];
            let fill = if ab.is_alert {
                Color32::from_rgb(244, 63, 94)
            } else {
                Color32::from_rgb(59, 130, 246)
            };
            bars.push(Bar::new(x, disp_c).width(0.74).fill(fill));

            if let Some(u) = ab.upper_threshold {
                upper_pts.push([x, u]);
            }
            if let Some(e) = ab.expected_baseline {
                exp_pts.push([x, e]);
            }
            if ab.is_alert {
                alert_pts.push([x, disp_c]);
            }
        }

        Plot::new("epi_curve_plot")
            .height(230.0)
            .legend(Legend::default())
            .show(ui, |plot_ui| {
                plot_ui.bar_chart(BarChart::new(bars).name("Observed Incidence"));
                plot_ui.line(
                    Line::new(PlotPoints::from(exp_pts))
                        .name("Expected Baseline")
                        .color(Color32::from_rgb(148, 163, 184)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(upper_pts))
                        .name("Upper Prediction Threshold")
                        .color(Color32::from_rgb(245, 158, 11)),
                );
                plot_ui.points(
                    Points::new(PlotPoints::from(alert_pts))
                        .radius(4.5)
                        .name("Aberration Alert")
                        .color(Color32::from_rgb(251, 113, 133)),
                );
            });

        ui.add_space(6.0);

        let mut rt_mean_pts = Vec::new();
        let mut rt_low_pts = Vec::new();
        let mut rt_high_pts = Vec::new();
        for (offset, idx) in (start_idx..total_pts).enumerate() {
            if let Some(est) = self.cached_rt.estimates[idx] {
                let x = offset as f64;
                rt_mean_pts.push([x, est.mean.min(8.0)]);
                rt_low_pts.push([x, est.cri_low_95.min(8.0)]);
                rt_high_pts.push([x, est.cri_high_95.min(8.0)]);
            }
        }

        let max_x = (total_pts.saturating_sub(start_idx).max(1)) as f64;
        let eq_pts = vec![[0.0, 1.0], [max_x, 1.0]];
        Plot::new("cori_rt_plot")
            .height(195.0)
            .legend(Legend::default())
            .show(ui, |plot_ui| {
                plot_ui.line(
                    Line::new(PlotPoints::from(eq_pts))
                        .name("Equilibrium Rt = 1.0")
                        .color(Color32::from_rgb(16, 185, 129)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(rt_mean_pts))
                        .name("Posterior Mean Rt")
                        .color(Color32::from_rgb(96, 165, 250)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(rt_high_pts))
                        .name("95% CrI Upper")
                        .color(Color32::from_rgba_unmultiplied(96, 165, 250, 110)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(rt_low_pts))
                        .name("95% CrI Lower")
                        .color(Color32::from_rgba_unmultiplied(96, 165, 250, 110)),
                );
            });
    }

    fn render_seir_nowcast_tab(&mut self, ui: &mut egui::Ui) {
        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("LATEST NOWCAST vs RAW").small().weak());
                if let Some(last_nc) = self.cached_nowcast.points.last() {
                    ui.heading(
                        RichText::new(format!("{:.1}", last_nc.nowcast_mean))
                            .monospace()
                            .color(Color32::from_rgb(245, 158, 11)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Raw {:.0} | {:.0}% accrued | +{:.1} unrep",
                            last_nc.observed_count,
                            last_nc.completeness_fraction * 100.0,
                            last_nc.expected_unreported
                        ))
                        .small()
                        .monospace(),
                    );
                } else {
                    ui.heading("--");
                }
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("SEIR INITIAL R_eff & R0").small().weak());
                let col = if self.cached_seir.initial_r_eff > 1.0 {
                    Color32::from_rgb(244, 63, 94)
                } else {
                    Color32::from_rgb(16, 185, 129)
                };
                ui.heading(
                    RichText::new(format!("R_eff = {:.2}", self.cached_seir.initial_r_eff))
                        .monospace()
                        .color(col),
                );
                ui.label(
                    RichText::new(format!("Calibrated R0 = {:.2}", self.cached_seir.calibrated_r0))
                        .small()
                        .monospace(),
                );
            });

            cols[2].group(|ui| {
                ui.label(RichText::new("PROJECTED PEAK INCIDENCE").small().weak());
                ui.heading(
                    RichText::new(format!(
                        "{:.0} / wk (Wk +{})",
                        self.cached_seir.peak_reported_weekly, self.cached_seir.peak_week_offset
                    ))
                    .monospace(),
                );
                ui.label(
                    RichText::new(format!("Horizon: {} weeks", self.seir_cfg.horizon_weeks))
                        .small()
                        .monospace(),
                );
            });

            cols[3].group(|ui| {
                ui.label(RichText::new("HORIZON CUMULATIVE CASES").small().weak());
                ui.heading(
                    RichText::new(format!(
                        "{:.0}",
                        self.cached_seir.cumulative_projected_reported
                    ))
                    .monospace()
                    .color(Color32::from_rgb(96, 165, 250)),
                );
                ui.label(
                    RichText::new(format!(
                        "95% PI [{:.0}, {:.0}]",
                        self.cached_seir.cumulative_pi_low_95,
                        self.cached_seir.cumulative_pi_high_95
                    ))
                    .small()
                    .monospace(),
                );
            });
        });

        ui.add_space(6.0);

        let mut dirty = false;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Nowcast Coupling:").strong());
                if ui
                    .checkbox(
                        &mut self.feed_nowcast_into_seir,
                        "Initialize SEIR from Nowcast-Corrected Tail",
                    )
                    .changed()
                {
                    dirty = true;
                }
                if self.vintages.len() >= 2
                    && ui
                        .checkbox(
                            &mut self.use_triangle_when_available,
                            "Use Multi-Vintage SQLite Reporting Triangle",
                        )
                        .changed()
                {
                    dirty = true;
                }

                ui.separator();
                ui.label(RichText::new("Pathogen SEIR Preset:").strong());
                for p in seir_pathogen_presets() {
                    if ui
                        .selectable_label(self.seir_preset_id == p.id, p.id)
                        .clicked()
                    {
                        self.seir_preset_id = p.id.to_string();
                        self.seir_cfg.r0 = p.r0_default;
                        self.seir_cfg.latent_days = p.latent_days;
                        self.seir_cfg.infectious_days = p.infectious_days;
                        self.seir_cfg.waning_days = p.waning_days;
                        self.seir_cfg.initial_immune_fraction = p.default_immune_fraction;
                        self.seir_cfg.reporting_fraction = p.reporting_fraction;
                        dirty = true;
                    }
                }
            });

            ui.separator();

            ui.horizontal_wrapped(|ui| {
                if ui
                    .add(
                        egui::Slider::new(&mut self.nowcast_cfg.delay_mean_weeks, 0.2..=6.0)
                            .text("Delay Mean (wk)"),
                    )
                    .changed()
                {
                    dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.nowcast_cfg.delay_sd_weeks, 0.2..=4.0)
                            .text("Delay SD (wk)"),
                    )
                    .changed()
                {
                    dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.nowcast_cfg.max_delay_weeks, 1..=12)
                            .text("Max Delay D (wk)"),
                    )
                    .changed()
                {
                    dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.nowcast_cfg.overdispersion_phi, 1.0..=10.0)
                            .text("Nowcast phi"),
                    )
                    .changed()
                {
                    dirty = true;
                }
            });

            ui.horizontal_wrapped(|ui| {
                if ui
                    .checkbox(&mut self.seir_cfg.auto_calibrate_r0_from_rt, "Auto-fit R0 from Rt")
                    .changed()
                {
                    dirty = true;
                }
                if !self.seir_cfg.auto_calibrate_r0_from_rt
                    && ui
                        .add(egui::Slider::new(&mut self.seir_cfg.r0, 0.5..=20.0).text("Manual R0"))
                        .changed()
                {
                    dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.seir_cfg.initial_immune_fraction, 0.0..=0.98)
                            .text("Initial Immune"),
                    )
                    .changed()
                {
                    dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.seir_cfg.latent_days, 0.5..=21.0)
                            .text("Latent 1/sigma (d)"),
                    )
                    .changed()
                {
                    dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.seir_cfg.infectious_days, 1.0..=28.0)
                            .text("Infectious 1/gamma (d)"),
                    )
                    .changed()
                {
                    dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.seir_cfg.contact_reduction, 0.0..=0.80)
                            .text("Contact Reduction"),
                    )
                    .changed()
                {
                    dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.seir_cfg.horizon_weeks, 4..=52)
                            .text("Horizon (wk)"),
                    )
                    .changed()
                {
                    dirty = true;
                }
            });
        });

        if dirty {
            self.recompute_seir_and_nowcast();
        }

        ui.add_space(6.0);

        let n_hist = self.cached_nowcast.points.len();
        let hist_tail = 26usize.min(n_hist);
        let h_start = n_hist.saturating_sub(hist_tail);

        let mut obs_bars = Vec::new();
        let mut nc_line = Vec::new();
        let mut nc_hi = Vec::new();
        let mut nc_lo = Vec::new();

        for (k, idx) in (h_start..n_hist).enumerate() {
            let x = (k as f64) - ((hist_tail.saturating_sub(1)) as f64);
            let pt = &self.cached_nowcast.points[idx];
            let bar_col = if pt.is_right_truncated {
                Color32::from_rgba_unmultiplied(59, 130, 246, 150)
            } else {
                Color32::from_rgb(59, 130, 246)
            };
            obs_bars.push(Bar::new(x, pt.observed_count).width(0.70).fill(bar_col));
            if pt.is_right_truncated || idx + 1 == n_hist {
                nc_line.push([x, pt.nowcast_mean]);
                nc_hi.push([x, pt.pi_high_95]);
                nc_lo.push([x, pt.pi_low_95]);
            }
        }

        let anchor_y = self
            .cached_nowcast
            .points
            .last()
            .map(|p| {
                if self.feed_nowcast_into_seir {
                    p.nowcast_mean
                } else {
                    p.observed_count
                }
            })
            .unwrap_or(0.0);

        let mut seir_mean_pts = vec![[0.0, anchor_y]];
        let mut seir_low_pts = vec![[0.0, anchor_y]];
        let mut seir_high_pts = vec![[0.0, anchor_y]];

        for step in &self.cached_seir.steps {
            let x = step.week_offset as f64;
            seir_mean_pts.push([x, step.reported_weekly_mean]);
            seir_low_pts.push([x, step.reported_pi_low_95]);
            seir_high_pts.push([x, step.reported_pi_high_95]);
        }

        Plot::new("seir_nowcast_trajectory_plot")
            .height(235.0)
            .legend(Legend::default())
            .show(ui, |plot_ui| {
                plot_ui.bar_chart(BarChart::new(obs_bars).name("Observed Provisional (t <= 0)"));
                plot_ui.line(
                    Line::new(PlotPoints::from(nc_line))
                        .name("Nowcast-Adjusted Mean")
                        .color(Color32::from_rgb(245, 158, 11)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(nc_hi))
                        .name("Nowcast 95% PI Upper")
                        .color(Color32::from_rgba_unmultiplied(245, 158, 11, 110)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(nc_lo))
                        .name("Nowcast 95% PI Lower")
                        .color(Color32::from_rgba_unmultiplied(245, 158, 11, 110)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(seir_mean_pts))
                        .name("SEIR RK4 Forecast Mean (t > 0)")
                        .color(Color32::from_rgb(52, 211, 153)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(seir_high_pts))
                        .name("SEIR Stochastic 95% PI Upper")
                        .color(Color32::from_rgba_unmultiplied(52, 211, 153, 115)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(seir_low_pts))
                        .name("SEIR Stochastic 95% PI Lower")
                        .color(Color32::from_rgba_unmultiplied(52, 211, 153, 115)),
                );
            });

        ui.add_space(6.0);

        let mut e_pts = Vec::new();
        let mut i_pts = Vec::new();
        for s in &self.cached_seir.steps {
            let x = s.week_offset as f64;
            e_pts.push([x, s.exposed]);
            i_pts.push([x, s.infectious]);
        }

        Plot::new("seir_compartments_plot")
            .height(175.0)
            .legend(Legend::default())
            .show(ui, |plot_ui| {
                plot_ui.line(
                    Line::new(PlotPoints::from(e_pts))
                        .name("Latent Exposed E(t)")
                        .color(Color32::from_rgb(192, 132, 252)),
                );
                plot_ui.line(
                    Line::new(PlotPoints::from(i_pts))
                        .name("Active Infectious I(t)")
                        .color(Color32::from_rgb(251, 113, 133)),
                );
            });
    }

    fn render_spatial_tab(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.label("Baseline Year:");
            if ui
                .selectable_label(self.spatial_use_2025, "2025 Provisional")
                .clicked()
            {
                self.spatial_use_2025 = true;
                self.recompute_spatial();
            }
            if ui
                .selectable_label(!self.spatial_use_2025, "2024 Provisional")
                .clicked()
            {
                self.spatial_use_2025 = false;
                self.recompute_spatial();
            }
            if ui
                .add(
                    egui::Slider::new(&mut self.max_zone_pop_frac, 0.08..=0.50)
                        .text("Max Cluster Pop Cap"),
                )
                .changed()
            {
                self.recompute_spatial();
            }
            if ui.button("Import County/State GeoJSON Scan...").clicked() {
                if let Some(p) = rfd::FileDialog::new()
                    .add_filter("GeoJSON", &["geojson", "json"])
                    .pick_file()
                {
                    match parse_geojson_centroids_file(&p) {
                        Ok(feats) => {
                            let zones: Vec<SpatialScanInputZone> = feats
                                .into_iter()
                                .map(|f| SpatialScanInputZone {
                                    id: f.geoid,
                                    name: f.name,
                                    lat_deg: f.lat_deg,
                                    lon_deg: f.lon_deg,
                                    observed_cases: f.cases.unwrap_or(0.0),
                                    population: f.population.unwrap_or(100_000.0),
                                })
                                .collect();
                            self.cached_zone_clusters = run_kulldorff_circular_scan(
                                &zones,
                                self.max_zone_pop_frac,
                                999,
                                6,
                            );
                            self.status_banner = format!(
                                "Completed M=999 Kulldorff circular scan on {} GeoJSON polygons.",
                                zones.len()
                            );
                        }
                        Err(e) => self.status_banner = format!("GeoJSON error: {e}"),
                    }
                }
            }
            if ui.button("Clear Live Overlay").clicked() {
                for r in &mut self.jurisdictions {
                    r.active_custom_value = None;
                }
                self.recompute_spatial();
            }
        });

        ui.add_space(6.0);

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.label(
                    RichText::new("US 51-Jurisdiction Cartogram (Emerald Border = Selected Multi-State Cluster)")
                        .strong(),
                );
                let cell_w = 40.0;
                let cell_h = 34.0;
                let gap = 4.0;
                let desired = Vec2::new(12.0 * (cell_w + gap), 8.0 * (cell_h + gap) + 8.0);
                let (rect, response) = ui.allocate_exact_size(desired, egui::Sense::click());
                let painter = ui.painter_at(rect);

                let mut by_fips: HashMap<&str, &JurisdictionRecord> = HashMap::new();
                let mut max_rate = 1.0_f64;
                for j in &self.jurisdictions {
                    let r = j.active_custom_value.unwrap_or(if self.spatial_use_2025 {
                        j.rate_2025
                    } else {
                        j.rate_2024
                    });
                    max_rate = max_rate.max(r);
                    by_fips.insert(j.fips.as_str(), j);
                }

                for &(col, row, fips, abbr) in US_CARTOGRAM_CELLS {
                    let x0 = rect.min.x + (col as f32) * (cell_w + gap);
                    let y0 = rect.min.y + 4.0 + (row as f32) * (cell_h + gap);
                    let cell_rect =
                        egui::Rect::from_min_size(egui::pos2(x0, y0), Vec2::new(cell_w, cell_h));

                    if response.clicked() {
                        if let Some(pos) = response.interact_pointer_pos() {
                            if cell_rect.contains(pos) {
                                self.selected_fips = fips.to_string();
                            }
                        }
                    }

                    let fill = if let Some(j) = by_fips.get(fips) {
                        let raw_cases = if self.spatial_use_2025 {
                            j.cases_2025
                        } else {
                            j.cases_2024
                        };
                        if self.suppression_policy == SuppressionPolicy::NchsWonderLe9
                            && raw_cases > 0
                            && raw_cases <= 9
                        {
                            Color32::from_rgb(51, 65, 85)
                        } else {
                            let val = j.active_custom_value.unwrap_or(if self.spatial_use_2025 {
                                j.rate_2025
                            } else {
                                j.rate_2024
                            });
                            let t = (val / max_rate).clamp(0.0, 1.0) as f32;
                            Color32::from_rgb(
                                (30.0 + 30.0 * (1.0 - t)) as u8,
                                (64.0 + 110.0 * t) as u8,
                                (120.0 + 130.0 * t) as u8,
                            )
                        }
                    } else {
                        Color32::from_rgb(30, 41, 59)
                    };

                    painter.rect_filled(cell_rect, 4.0, fill);
                    if self.selected_fips == fips {
                        painter.rect_stroke(
                            cell_rect,
                            4.0,
                            Stroke::new(2.2, Color32::from_rgb(245, 158, 11)),
                        );
                    } else if self.highlighted_cluster_members.contains(fips) {
                        painter.rect_stroke(
                            cell_rect,
                            4.0,
                            Stroke::new(1.8, Color32::from_rgb(52, 211, 153)),
                        );
                    }
                    painter.text(
                        cell_rect.center(),
                        egui::Align2::CENTER_CENTER,
                        abbr,
                        egui::FontId::monospace(11.5),
                        Color32::WHITE,
                    );
                }

                if let Some(sel) = self
                    .jurisdictions
                    .iter()
                    .find(|j| j.fips == self.selected_fips)
                {
                    ui.separator();
                    ui.label(
                        RichText::new(format!(
                            "{} (FIPS {}) — 2025 Census Pop: {}",
                            sel.name, sel.fips, sel.census_pop_2025
                        ))
                        .strong(),
                    );
                    ui.label(format!(
                        "2025 Cases: {} ({:.2}/100k) | 2024 Cases: {} ({:.2}/100k)",
                        sel.cases_2025, sel.rate_2025, sel.cases_2024, sel.rate_2024
                    ));
                }
            });

            cols[1].group(|ui| {
                ui.label(
                    RichText::new("Kulldorff Circular Multi-Zone Clusters (Monte Carlo p-values)")
                        .strong(),
                );
                let mut select_members: Option<HashSet<String>> = None;
                egui::ScrollArea::vertical().max_height(190.0).show(ui, |ui| {
                    for (idx, zc) in self.cached_zone_clusters.iter().enumerate() {
                        ui.group(|ui| {
                            ui.horizontal_wrapped(|ui| {
                                if ui
                                    .button(format!("#{} {}", idx + 1, zc.center_name))
                                    .clicked()
                                {
                                    select_members = Some(zc.member_ids.iter().cloned().collect());
                                }
                                ui.monospace(format!(
                                    "R={:.0}km | RR={:.2}x | LLR={:.1} | p={:.3}",
                                    zc.radius_km,
                                    zc.relative_risk,
                                    zc.kulldorff_llr,
                                    zc.monte_carlo_p_value
                                ));
                            });
                            ui.label(
                                RichText::new(format!(
                                    "Zones ({}): {}",
                                    zc.member_names.len(),
                                    zc.member_names.join(", ")
                                ))
                                .small()
                                .weak(),
                            );
                        });
                    }
                });
                if let Some(m) = select_members {
                    self.highlighted_cluster_members = m;
                }

                ui.separator();
                ui.label(RichText::new("Single-Jurisdiction Excess LLR Table").strong());
                egui::ScrollArea::vertical().max_height(160.0).show(ui, |ui| {
                    egui::Grid::new("kulldorff_single_grid")
                        .striped(true)
                        .show(ui, |ui| {
                            ui.label(RichText::new("State").strong());
                            ui.label(RichText::new("Obs").strong());
                            ui.label(RichText::new("Exp").strong());
                            ui.label(RichText::new("RR").strong());
                            ui.label(RichText::new("LLR").strong());
                            ui.end_row();
                            for c in &self.cached_clusters {
                                if ui.link(&c.name).clicked() {
                                    self.selected_fips = c.fips.clone();
                                }
                                ui.monospace(format!("{:.0}", c.observed_cases));
                                ui.monospace(format!("{:.1}", c.expected_cases));
                                ui.monospace(format!("{:.2}x", c.relative_risk));
                                ui.monospace(format!("{:.1}", c.kulldorff_llr));
                                ui.end_row();
                            }
                        });
                });
            });
        });
    }

    fn render_live_connectors_tab(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.heading("CDC Socrata SODA & CMU Delphi Epidata");
                ui.horizontal(|ui| {
                    ui.radio_value(&mut self.soda_is_nwss, false, "NNDSS Weekly (x9gk-5huc)");
                    ui.radio_value(&mut self.soda_is_nwss, true, "NWSS Wastewater (2ew6-ywp6)");
                });
                ui.horizontal(|ui| {
                    ui.label("Disease / Query ($q):");
                    ui.text_edit_singleline(&mut self.soda_condition);
                });
                ui.horizontal(|ui| {
                    ui.label("Jurisdiction:");
                    ui.text_edit_singleline(&mut self.soda_jurisdiction);
                });
                ui.add(egui::Slider::new(&mut self.soda_limit, 100..=2500).text("Batch Limit"));

                if ui
                    .add_enabled(!self.worker_busy, egui::Button::new("Fetch Live CDC SODA"))
                    .clicked()
                {
                    self.worker_busy = true;
                    let tx = self.worker_tx.clone();
                    let db = self.db_path.clone();
                    let kind = if self.soda_is_nwss {
                        SodaDatasetKind::NwssWastewater
                    } else {
                        SodaDatasetKind::NndssWeekly
                    };
                    let cond = self.soda_condition.clone();
                    let jur = self.soda_jurisdiction.clone();
                    let lim = self.soda_limit;
                    let ctx_clone = ctx.clone();

                    std::thread::spawn(move || {
                        let res = (|| {
                            let mut st = EpiStore::open_file(&db)
                                .or_else(|_| EpiStore::open_in_memory())?;
                            CdcSodaClient::default()
                                .fetch_with_lazy_cache(&mut st, kind, &cond, &jur, lim, 21_600)
                        })();
                        match res {
                            Ok(out) => {
                                let _ = tx.send(WorkerMessage::SodaLoaded {
                                    dataset: out.dataset,
                                    state_map: out.state_values_by_fips,
                                    raw_rows: out.raw_row_count,
                                    stale: out.used_stale_cache,
                                });
                            }
                            Err(e) => {
                                let _ = tx.send(WorkerMessage::Error(e.to_string()));
                            }
                        }
                        ctx_clone.request_repaint();
                    });
                }

                ui.separator();
                ui.label(RichText::new("CMU Delphi FluView ILINet Connector").strong());
                ui.horizontal(|ui| {
                    ui.label("Region (nat/or/ca):");
                    ui.text_edit_singleline(&mut self.delphi_region);
                    ui.label("Epiweeks:");
                    ui.text_edit_singleline(&mut self.delphi_epiweeks);
                });
                if ui
                    .add_enabled(
                        !self.worker_busy,
                        egui::Button::new("Fetch CMU Delphi FluView ILINet"),
                    )
                    .clicked()
                {
                    self.worker_busy = true;
                    let tx = self.worker_tx.clone();
                    let db = self.db_path.clone();
                    let reg = self.delphi_region.clone();
                    let wks = self.delphi_epiweeks.clone();
                    let ctx_clone = ctx.clone();
                    std::thread::spawn(move || {
                        let res = (|| {
                            let mut st = EpiStore::open_file(&db)
                                .or_else(|_| EpiStore::open_in_memory())?;
                            CmuDelphiClient::default()
                                .fetch_fluview_ilinet(&mut st, &reg, &wks, true, 21_600)
                        })();
                        match res {
                            Ok(ds) => {
                                let _ = tx.send(WorkerMessage::DelphiLoaded(ds));
                            }
                            Err(e) => {
                                let _ = tx.send(WorkerMessage::Error(e.to_string()));
                            }
                        }
                        ctx_clone.request_repaint();
                    });
                }
            });

            cols[1].group(|ui| {
                ui.heading("openFDA FAERS & WHO GHO OData");
                ui.horizontal(|ui| {
                    ui.label("Generic Drug (FAERS):");
                    ui.text_edit_singleline(&mut self.fda_drug);
                });
                ui.horizontal(|ui| {
                    if ui
                        .add_enabled(
                            !self.worker_busy,
                            egui::Button::new("Fetch FAERS Weekly + Reactions"),
                        )
                        .clicked()
                    {
                        self.worker_busy = true;
                        let tx = self.worker_tx.clone();
                        let db = self.db_path.clone();
                        let drug = self.fda_drug.clone();
                        let ctx_clone = ctx.clone();
                        std::thread::spawn(move || {
                            let res = (|| {
                                let mut st = EpiStore::open_file(&db)
                                    .or_else(|_| EpiStore::open_in_memory())?;
                                OpenFdaClient::default().fetch_faers_with_cache(
                                    &mut st, &drug, 43_200,
                                )
                            })();
                            match res {
                                Ok(out) => {
                                    let _ = tx.send(WorkerMessage::FdaFaersLoaded {
                                        dataset: out.weekly_series,
                                        reactions: out.top_reactions,
                                        stale: out.used_stale_cache,
                                    });
                                }
                                Err(e) => {
                                    let _ = tx.send(WorkerMessage::Error(e.to_string()));
                                }
                            }
                            ctx_clone.request_repaint();
                        });
                    }

                    if ui
                        .add_enabled(!self.worker_busy, egui::Button::new("Fetch Drug Recalls"))
                        .clicked()
                    {
                        self.worker_busy = true;
                        let tx = self.worker_tx.clone();
                        let ctx_clone = ctx.clone();
                        std::thread::spawn(move || {
                            match OpenFdaClient::default().fetch_drug_recalls(25) {
                                Ok(list) => {
                                    let _ = tx.send(WorkerMessage::FdaRecallsLoaded(list));
                                }
                                Err(e) => {
                                    let _ = tx.send(WorkerMessage::Error(e.to_string()));
                                }
                            }
                            ctx_clone.request_repaint();
                        });
                    }
                });

                ui.separator();
                ui.horizontal(|ui| {
                    ui.label("WHO Indicator:");
                    ui.text_edit_singleline(&mut self.who_indicator);
                    ui.label("ISO3:");
                    ui.text_edit_singleline(&mut self.who_iso3);
                });
                if ui
                    .add_enabled(!self.worker_busy, egui::Button::new("Fetch WHO GHO Series"))
                    .clicked()
                {
                    self.worker_busy = true;
                    let tx = self.worker_tx.clone();
                    let db = self.db_path.clone();
                    let ind = self.who_indicator.clone();
                    let iso = self.who_iso3.clone();
                    let ctx_clone = ctx.clone();
                    std::thread::spawn(move || {
                        let res = (|| {
                            let mut st = EpiStore::open_file(&db)
                                .or_else(|_| EpiStore::open_in_memory())?;
                            WhoGhoClient::default().fetch_indicator_series(&mut st, &ind, &iso, 86_400)
                        })();
                        match res {
                            Ok(ds) => {
                                let _ = tx.send(WorkerMessage::WhoLoaded(ds));
                            }
                            Err(e) => {
                                let _ = tx.send(WorkerMessage::Error(e.to_string()));
                            }
                        }
                        ctx_clone.request_repaint();
                    });
                }
            });
        });

        ui.add_space(8.0);
        ui.group(|ui| {
            ui.label(RichText::new("openFDA Drug Recall Enforcement Stream").strong());
            egui::ScrollArea::vertical().max_height(190.0).show(ui, |ui| {
                egui::Grid::new("fda_recalls_grid").striped(true).show(ui, |ui| {
                    ui.label(RichText::new("Recall ID").strong());
                    ui.label(RichText::new("Class").strong());
                    ui.label(RichText::new("Date").strong());
                    ui.label(RichText::new("Firm").strong());
                    ui.label(RichText::new("Reason").strong());
                    ui.end_row();

                    for r in &self.fda_recalls {
                        ui.monospace(&r.recall_number);
                        ui.label(&r.classification);
                        ui.monospace(&r.recall_initiation_date);
                        ui.label(format!("{} ({})", r.recalling_firm, r.state));
                        ui.label(&r.reason_for_recall);
                        ui.end_row();
                    }
                });
            });
        });
    }

    fn render_local_import_tab(&mut self, ui: &mut egui::Ui) {
        ui.group(|ui| {
            ui.heading("Offline Public Dataset & HL7 FHIR R4 Bundle Ingestion");
            ui.horizontal(|ui| {
                ui.label("File Path:");
                ui.text_edit_singleline(&mut self.import_path_input);
                if ui.button("Browse File...").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter(
                            "Public Health Files",
                            &["csv", "tsv", "txt", "json", "xpt", "XPT", "geojson"],
                        )
                        .pick_file()
                    {
                        self.import_path_input = p.display().to_string();
                    }
                }
                ui.label("FHIR Code Filter:");
                ui.text_edit_singleline(&mut self.fhir_code_filter);
            });

            ui.horizontal_wrapped(|ui| {
                if ui.button("Import CSV / TSV / JSON Series").clicked() {
                    let path = PathBuf::from(self.import_path_input.trim());
                    let is_json = path
                        .extension()
                        .and_then(|s| s.to_str())
                        .map(|e| e.eq_ignore_ascii_case("json"))
                        .unwrap_or(false);
                    let res = if is_json {
                        import_json_timeseries(&path)
                    } else {
                        import_delimited_timeseries(&path, None, &ColumnMapping::default())
                    };
                    match res {
                        Ok(ds) => {
                            let _ = self.store.commit_vintage(&ds);
                            self.refresh_vintages();
                            self.active_dataset = ds;
                            self.recompute_analytics();
                            self.active_tab = ActiveTab::Surveillance;
                            self.status_banner =
                                "Imported local dataset and committed to SQLite vintage ledger."
                                    .to_string();
                        }
                        Err(e) => self.status_banner = format!("Import failed: {e}"),
                    }
                }

                if ui.button("Import HL7 FHIR R4 Bundle (.json)").clicked() {
                    let path = PathBuf::from(self.import_path_input.trim());
                    let flt = if self.fhir_code_filter.trim().is_empty() {
                        None
                    } else {
                        Some(self.fhir_code_filter.trim())
                    };
                    match import_fhir_r4_bundle_file(&path, flt) {
                        Ok(sum) => {
                            let _ = self.store.commit_vintage(&sum.dataset);
                            self.refresh_vintages();
                            if !sum.state_counts_by_fips.is_empty() {
                                for r in &mut self.jurisdictions {
                                    r.active_custom_value =
                                        sum.state_counts_by_fips.get(&r.fips).copied();
                                }
                                self.recompute_spatial();
                            }
                            self.active_dataset = sum.dataset;
                            self.recompute_analytics();
                            self.active_tab = ActiveTab::Surveillance;
                            self.status_banner = format!(
                                "Imported HL7 FHIR R4 Bundle ({} resources scanned).",
                                sum.total_resources_scanned
                            );
                        }
                        Err(e) => self.status_banner = format!("FHIR R4 import error: {e}"),
                    }
                }

                if ui.button("Decode CDC NHANES / BRFSS SAS .XPT").clicked() {
                    let path = PathBuf::from(self.import_path_input.trim());
                    match std::fs::read(&path)
                        .map_err(epi_core::EpiError::from)
                        .and_then(|bytes| parse_sas_xport_bytes(&bytes, 100_000))
                    {
                        Ok(xpt) => {
                            self.status_banner = format!(
                                "Decoded SAS XPORT '{}' ({} variables, {} rows loaded for Taylor Survey Engine).",
                                xpt.member_name,
                                xpt.variables.len(),
                                xpt.rows.len()
                            );
                            self.xpt_preview = Some(xpt);
                        }
                        Err(e) => self.status_banner = format!("SAS .XPT decode error: {e}"),
                    }
                }

                if ui.button("Export Active Series CSV...").clicked() {
                    if let Some(out) = rfd::FileDialog::new()
                        .set_file_name("epiledger_export.csv")
                        .save_file()
                    {
                        match export_dataset_csv(
                            &self.active_dataset,
                            &out,
                            self.suppression_policy,
                        ) {
                            Ok(()) => {
                                self.status_banner = format!("Exported CSV to {}", out.display())
                            }
                            Err(e) => self.status_banner = format!("CSV export failed: {e}"),
                        }
                    }
                }

                if ui.button("Export Active Vintage as HL7 FHIR R4 Bundle...").clicked() {
                    if let Some(out) = rfd::FileDialog::new()
                        .set_file_name("epiledger_fhir_r4_bundle.json")
                        .save_file()
                    {
                        match export_timeseries_as_fhir_r4_bundle(&self.active_dataset, "14189004")
                            .and_then(|txt| std::fs::write(&out, txt).map_err(Into::into))
                        {
                            Ok(()) => {
                                self.status_banner =
                                    format!("Exported HL7 FHIR R4 Bundle to {}", out.display())
                            }
                            Err(e) => self.status_banner = format!("FHIR export failed: {e}"),
                        }
                    }
                }

                if ui.button("Export Active Vintage JSON...").clicked() {
                    if let Some(out) = rfd::FileDialog::new()
                        .set_file_name("epiledger_vintage.json")
                        .save_file()
                    {
                        match export_dataset_json(&self.active_dataset, &out) {
                            Ok(()) => {
                                self.status_banner = format!("Exported JSON to {}", out.display())
                            }
                            Err(e) => self.status_banner = format!("JSON export failed: {e}"),
                        }
                    }
                }
            });
        });

        if let Some(xpt) = &self.xpt_preview {
            ui.add_space(8.0);
            ui.group(|ui| {
                ui.label(
                    RichText::new(format!(
                        "SAS XPORT v5 Binary Inspector — Member: {} ({} Variables, {} Rows)",
                        xpt.member_name,
                        xpt.variables.len(),
                        xpt.rows.len()
                    ))
                    .strong(),
                );
                egui::ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                    egui::Grid::new("xpt_vars_grid").striped(true).show(ui, |ui| {
                        ui.label(RichText::new("Variable").strong());
                        ui.label(RichText::new("Type").strong());
                        ui.label(RichText::new("Bytes").strong());
                        ui.label(RichText::new("SAS Label").strong());
                        ui.end_row();
                        for v in &xpt.variables {
                            ui.monospace(&v.name);
                            ui.monospace(if v.is_char { "CHAR" } else { "IBM370_FLOAT" });
                            ui.monospace(v.length.to_string());
                            ui.label(&v.label);
                            ui.end_row();
                        }
                    });
                });
            });
        }
    }

    fn render_survey_and_prr_tab(&mut self, ui: &mut egui::Ui) {
        ui.columns(2, |cols| {
            // Left: Complex Survey Taylor-Series Linearization (NHANES / BRFSS .XPT)
            cols[0].group(|ui| {
                ui.heading("Complex Survey Taylor-Series Linearization (NHANES / BRFSS)");
                ui.label(
                    RichText::new(
                        "Computes stratified cluster-adjusted mean/prevalence, linearized SE, Korn-Graubard df, and Design Effect (DEFF) with NOMCAR preservation.",
                    )
                    .small()
                    .weak(),
                );
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.label("Stratum Col:");
                    ui.text_edit_singleline(&mut self.survey_stratum_col);
                    ui.label("PSU Col:");
                    ui.text_edit_singleline(&mut self.survey_psu_col);
                });
                ui.horizontal(|ui| {
                    ui.label("Weight Col:");
                    ui.text_edit_singleline(&mut self.survey_weight_col);
                    ui.label("Outcome Col:");
                    ui.text_edit_singleline(&mut self.survey_outcome_col);
                });

                ui.horizontal(|ui| {
                    if ui.button("NHANES Preset").clicked() {
                        self.survey_stratum_col = "SDMVSTRA".to_string();
                        self.survey_psu_col = "SDMVPSU".to_string();
                        self.survey_weight_col = "WTMEC2YR".to_string();
                    }
                    if ui.button("BRFSS Preset").clicked() {
                        self.survey_stratum_col = "_STSTR".to_string();
                        self.survey_psu_col = "_PSU".to_string();
                        self.survey_weight_col = "_LLCPWT".to_string();
                    }
                    if ui
                        .button("Compute Linearized Estimate on Loaded .XPT")
                        .clicked()
                    {
                        if let Some(xpt) = &self.xpt_preview {
                            match extract_survey_obs_from_xpt(
                                xpt,
                                &self.survey_stratum_col,
                                &self.survey_psu_col,
                                &self.survey_weight_col,
                                &self.survey_outcome_col,
                            ) {
                                Ok(obs) => {
                                    self.survey_estimate = compute_taylor_linearized_mean(
                                        &obs,
                                        LonelyPsuStrategy::GrandMeanCenter,
                                    );
                                    if self.survey_estimate.is_some() {
                                        self.status_banner = format!(
                                            "Computed Taylor-linearized survey estimate across {} weighted rows.",
                                            obs.len()
                                        );
                                    } else {
                                        self.status_banner =
                                            "Insufficient valid survey observations.".to_string();
                                    }
                                }
                                Err(e) => self.status_banner = format!("Survey column error: {e}"),
                            }
                        } else {
                            self.status_banner =
                                "Decode a SAS .XPT file in Tab 5 first to run Taylor linearization."
                                    .to_string();
                        }
                    }
                });

                if let Some(est) = &self.survey_estimate {
                    ui.separator();
                    ui.monospace(format!(
                        "Weighted Mean / Prop : {:.4} (95% CI [{:.4}, {:.4}])",
                        est.weighted_mean_or_prop, est.ci_low_95, est.ci_high_95
                    ));
                    ui.monospace(format!(
                        "Taylor Linearized SE : {:.4} | Design Effect (DEFF) = {:.2}",
                        est.standard_error, est.design_effect_deff
                    ));
                    ui.monospace(format!(
                        "Design Structure     : {} strata, {} PSUs (df = {}, lonely strata = {})",
                        est.n_strata, est.n_psus, est.degrees_of_freedom, est.lonely_strata_count
                    ));
                    ui.monospace(format!(
                        "Subpop Sample n      : {} | Sum of Weights = {:.0}",
                        est.subpop_sample_size, est.sum_of_weights
                    ));
                }
            });

            // Right: Pharmacovigilance 2x2 Disproportionality (PRR, ROR, BCPNN IC / IC025)
            cols[1].group(|ui| {
                ui.heading("FAERS Disproportionality Signal Detector (PRR, ROR & BCPNN IC)");
                ui.label(
                    RichText::new(
                        "Evaluates 2x2 drug-reaction contingency tables using Evans PRR criteria, Woolf ROR 95% CI, and Bayesian shrinkage Information Component (IC025).",
                    )
                    .small()
                    .weak(),
                );
                ui.add_space(4.0);

                let mut pv_dirty = false;
                ui.horizontal(|ui| {
                    ui.label("Drug:");
                    pv_dirty |= ui.text_edit_singleline(&mut self.pv_drug_name).changed();
                    ui.label("Reaction:");
                    pv_dirty |= ui.text_edit_singleline(&mut self.pv_reaction_name).changed();
                });
                ui.horizontal(|ui| {
                    pv_dirty |= ui
                        .add(egui::DragValue::new(&mut self.pv_a).prefix("a (Drug+Rx): "))
                        .changed();
                    pv_dirty |= ui
                        .add(egui::DragValue::new(&mut self.pv_b).prefix("b (Drug+Other): "))
                        .changed();
                });
                ui.horizontal(|ui| {
                    pv_dirty |= ui
                        .add(egui::DragValue::new(&mut self.pv_c).prefix("c (Other+Rx): "))
                        .changed();
                    pv_dirty |= ui
                        .add(egui::DragValue::new(&mut self.pv_d).prefix("d (Background): "))
                        .changed();
                });

                if pv_dirty {
                    self.pv_signal = compute_disproportionality_2x2(
                        &self.pv_drug_name,
                        &self.pv_reaction_name,
                        self.pv_a,
                        self.pv_b,
                        self.pv_c,
                        self.pv_d,
                    );
                }

                ui.separator();
                let s = &self.pv_signal;
                ui.monospace(format!(
                    "Proportional Reporting Ratio (PRR) : {:.2}  (Yates chi2 = {:.1})",
                    s.prr, s.prr_chi2_yates
                ));
                ui.monospace(format!(
                    "Reporting Odds Ratio (ROR)         : {:.2}  (95% CI [{:.2}, {:.2}])",
                    s.ror, s.ror_ci_low_95, s.ror_ci_high_95
                ));
                ui.monospace(format!(
                    "Bayesian Information Component     : IC = {:+.2} bits  (IC025 = {:+.2})",
                    s.bcpnn_ic, s.bcpnn_ic_025
                ));
                ui.horizontal(|ui| {
                    let ev_col = if s.flagged_evans_prr {
                        Color32::from_rgb(244, 63, 94)
                    } else {
                        Color32::from_rgb(16, 185, 129)
                    };
                    let ic_col = if s.flagged_bcpnn_ic {
                        Color32::from_rgb(244, 63, 94)
                    } else {
                        Color32::from_rgb(16, 185, 129)
                    };
                    ui.label(
                        RichText::new(format!(
                            "Evans PRR Signal: {}",
                            if s.flagged_evans_prr { "FLAGGED" } else { "Below Threshold" }
                        ))
                        .monospace()
                        .color(ev_col),
                    );
                    ui.separator();
                    ui.label(
                        RichText::new(format!(
                            "BCPNN IC025 > 0 Signal: {}",
                            if s.flagged_bcpnn_ic { "FLAGGED" } else { "Below Threshold" }
                        ))
                        .monospace()
                        .color(ic_col),
                    );
                });
            });
        });
    }

    fn render_vintage_ledger_tab(&mut self, ui: &mut egui::Ui) {
        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.heading("SQLite WAL Bitemporal Vintage Ledger");
                let mut load_id: Option<String> = None;
                egui::ScrollArea::vertical().max_height(340.0).show(ui, |ui| {
                    for v in &self.vintages {
                        let is_cur = v.vintage_id == self.active_dataset.vintage_id;
                        ui.group(|ui| {
                            ui.horizontal(|ui| {
                                ui.strong(&v.title);
                                ui.monospace(format!("[{}]", v.source_kind));
                            });
                            ui.monospace(format!(
                                "ID: {} | {} pts | SHA256-32: {}",
                                v.vintage_id, v.row_count, v.checksum_hex
                            ));
                            if !is_cur && ui.button("Activate Vintage").clicked() {
                                load_id = Some(v.vintage_id.clone());
                            } else if is_cur {
                                ui.label(
                                    RichText::new("Active in Workbench")
                                        .color(Color32::from_rgb(52, 211, 153)),
                                );
                            }
                        });
                    }
                });
                if let Some(vid) = load_id {
                    if let Ok(Some(ds)) = self.store.load_vintage(&vid) {
                        self.active_dataset = ds;
                        self.recompute_analytics();
                        self.status_banner = format!("Switched active vintage to '{vid}'.");
                    }
                }
            });

            cols[1].group(|ui| {
                ui.heading("Retroactive Backfill Revision Diff");
                ui.label(
                    RichText::new(
                        "Compare two vintages over overlapping dates to audit state backfill revisions.",
                    )
                    .small()
                    .weak(),
                );
                ui.horizontal(|ui| {
                    ui.label("Older Vintage ID:");
                    ui.text_edit_singleline(&mut self.diff_older_id);
                });
                ui.horizontal(|ui| {
                    ui.label("Newer Vintage ID:");
                    ui.text_edit_singleline(&mut self.diff_newer_id);
                });
                if ui.button("Compute Backfill Revisions").clicked() {
                    match self.store.diff_vintages(
                        self.diff_older_id.trim(),
                        self.diff_newer_id.trim(),
                    ) {
                        Ok(diffs) => {
                            self.status_banner =
                                format!("Found {} revised periods between vintages.", diffs.len());
                            self.diff_results = diffs;
                        }
                        Err(e) => self.status_banner = format!("Diff error: {e}"),
                    }
                }

                egui::ScrollArea::vertical().max_height(240.0).show(ui, |ui| {
                    egui::Grid::new("diff_grid").striped(true).show(ui, |ui| {
                        ui.label(RichText::new("Period").strong());
                        ui.label(RichText::new("Older").strong());
                        ui.label(RichText::new("Newer").strong());
                        ui.label(RichText::new("Delta").strong());
                        ui.end_row();
                        for d in &self.diff_results {
                            ui.monospace(&d.date_label);
                            ui.monospace(format!("{:.1}", d.older_count));
                            ui.monospace(format!("{:.1}", d.newer_count));
                            ui.monospace(format!("{:+.1}", d.delta));
                            ui.end_row();
                        }
                    });
                });
            });
        });
    }
}

impl eframe::App for EpiLedgerGuiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.update_ui(ctx);
    }
}

impl EpiLedgerGuiApp {
    pub fn configure_screenshot_state(&mut self) {
        self.configure_screenshot_view(1);
    }

    pub fn configure_screenshot_view(&mut self, view: u8) {
        match view {
            2 => {
                self.active_tab = ActiveTab::WastewaterAndAmr;
            }
            3 => {
                self.active_tab = ActiveTab::DpAuditAndHpc;
                self.dp_audit.trigger_dp_release_demo(
                    &self.store,
                    &self.active_dataset,
                    &mut self.status_banner,
                );
            }
            _ => {
                self.active_tab = ActiveTab::Surveillance;
                self.window_weeks = 0;
                self.recompute_analytics();
            }
        }
    }

    pub fn update_ui(&mut self, ctx: &egui::Context) {
        self.poll_worker();

        egui::TopBottomPanel::top("top_nav").show(ctx, |ui| {
            self.render_top_bar(ui);
        });

        egui::TopBottomPanel::bottom("bottom_status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.worker_busy {
                    ui.spinner();
                }
                ui.label(RichText::new(&self.status_banner).small().monospace());
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| match self.active_tab {
                ActiveTab::Surveillance => self.render_surveillance_tab(ui),
                ActiveTab::SeirNowcast => self.render_seir_nowcast_tab(ui),
                ActiveTab::SpatialScan => self.render_spatial_tab(ui),
                ActiveTab::LiveConnectors => self.render_live_connectors_tab(ui, ctx),
                ActiveTab::LocalImportXpt => self.render_local_import_tab(ui),
                ActiveTab::SurveyAndPrr => self.render_survey_and_prr_tab(ui),
                ActiveTab::GenomicPhylo => self.genomic.render(ui, &mut self.status_banner),
                ActiveTab::EnvAndMortality => self.env_mortality.render(ui, &mut self.status_banner),
                ActiveTab::ElrPprlSero => {
                    let act = self.elr_pprl_sero.render(ui, &mut self.status_banner);
                    if let Some(ds) = act.commit_hl7_dataset {
                        let vid = ds.vintage_id.clone();
                        let _ = self.store.commit_vintage(&ds);
                        self.refresh_vintages();
                        self.active_dataset = ds;
                        self.recompute_analytics();
                        self.status_banner = format!("Committed deduplicated HL7 v2.5.1 vintage '{vid}' into SQLite.");
                    }
                    if let Some((immune_frac, waning_days)) = act.apply_sero_to_seir {
                        self.seir_cfg.initial_immune_fraction = immune_frac.clamp(0.0, 0.99);
                        self.seir_cfg.waning_days = waning_days.max(0.0);
                        self.recompute_seir_and_nowcast();
                        self.active_tab = ActiveTab::SeirNowcast;
                        self.status_banner = format!(
                            "Calibrated SEIR initial immunity ({:.1}%) & waning ({:.0}d) from Sero-Catalytic FOI fit.",
                            immune_frac * 100.0,
                            waning_days
                        );
                    }
                }
                ActiveTab::ForecastScoringSitrep => self.scoring_sitrep.render(
                    ui,
                    &self.active_dataset,
                    self.suppression_policy,
                    &self.rt_cfg,
                    &self.cached_rt,
                    &self.cached_ab,
                    &self.cached_nowcast,
                    &self.seir_cfg,
                    &self.cached_zone_clusters,
                    &mut self.status_banner,
                ),
                ActiveTab::DpAuditAndHpc => {
                    let act = self.dp_audit.render(
                        ui,
                        &self.store,
                        &self.active_dataset,
                        &mut self.status_banner,
                    );
                    if let Some(dp_ds) = act.commit_dp_vintage {
                        let vid = dp_ds.vintage_id.clone();
                        let _ = self.store.commit_vintage(&dp_ds);
                        self.refresh_vintages();
                        self.active_dataset = dp_ds;
                        self.recompute_analytics();
                        self.dp_audit.refresh_from_store(&self.store, &vid);
                    }
                }
                ActiveTab::SuperspreadingNet => {
                    let act = self.superspreading_net.render(ui, &mut self.status_banner);
                    if let Some((si_mean, si_sd)) = act.apply_si_to_cori_rt {
                        self.rt_cfg.si_mean_days = si_mean.clamp(0.5, 40.0);
                        self.rt_cfg.si_sd_days = si_sd.clamp(0.25, 25.0);
                        self.selected_si_preset = "custom_contact_network".to_string();
                        self.recompute_analytics();
                        self.active_tab = ActiveTab::Surveillance;
                        self.status_banner = format!(
                            "Applied empirical contact-tracing Serial Interval (mean={:.2}d, SD={:.2}d) to Cori Rt.",
                            self.rt_cfg.si_mean_days, self.rt_cfg.si_sd_days
                        );
                    }
                }
                ActiveTab::WastewaterAndAmr => {
                    let act = self.wastewater_amr.render(ui, &mut self.status_banner);
                    if let Some(lin_tbl) = act.push_wastewater_to_genomic {
                        let n_wks = lin_tbl.rows.len();
                        let n_clades = lin_tbl.clade_names.len();
                        self.genomic.lineage_table = lin_tbl;
                        let _ = self.genomic.recompute_lineage();
                        self.active_tab = ActiveTab::GenomicPhylo;
                        self.status_banner = format!(
                            "Loaded Freyja-deconvolved sewage trajectories ({n_wks} weeks, {n_clades} lineages) into Multinomial Genomic Fitness."
                        );
                    }
                }
                ActiveTab::CrcAndSyndromic => {
                    let act = self.crc_syndromic.render(ui, &mut self.status_banner);
                    if let Some(syn_ds) = act.commit_syndrome_vintage {
                        let vid = syn_ds.vintage_id.clone();
                        let _ = self.store.commit_vintage(&syn_ds);
                        self.refresh_vintages();
                        self.active_dataset = syn_ds;
                        self.recompute_analytics();
                        self.active_tab = ActiveTab::Surveillance;
                        self.status_banner = format!(
                            "Committed NSSP ESSENCE syndrome series '{vid}' into SQLite and updated Cori Rt / Farrington."
                        );
                    }
                }
                ActiveTab::SeverityAndSurge => {
                    self.severity_surge
                        .render(ui, &self.cached_seir, &mut self.status_banner);
                }
                ActiveTab::VintageLedger => self.render_vintage_ledger_tab(ui),
            });
        });
    }
}
