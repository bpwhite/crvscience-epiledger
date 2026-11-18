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
use egui_plot::{Bar, BarChart, Legend, Line, Plot, PlotPoints, Points, VLine};
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
    builtin_measles_weekly_dataset, builtin_pertussis_spatial_records,
    parse_flexible_date_to_mmwr, pathogen_si_presets, CellValue, JurisdictionRecord,
    ReactionFacet, RecallRecord, SuppressionPolicy, SurveillanceCampaignPreset, TimeSeriesDataset,
    VintageMetadata,
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
use std::sync::Arc;

fn major_step_for_span(span: f64) -> usize {
    if span > 180.0 {
        26
    } else if span > 80.0 {
        13
    } else if span > 32.0 {
        4
    } else if span > 14.0 {
        2
    } else {
        1
    }
}

pub(crate) fn mmwr_week_grid_spacer(
    input: egui_plot::GridInput,
    n_pts: usize,
) -> Vec<egui_plot::GridMark> {
    let (min_x, max_x) = input.bounds;
    let span = (max_x - min_x).abs().max(1.0);
    let major_step = major_step_for_span(span);
    let minor_step: usize = match major_step {
        26 => 13,
        13 => 4,
        4 => 1,
        _ => 1,
    };
    let start = (min_x.floor() as isize).max(0) as usize;
    let end = ((max_x.ceil() as isize).max(0) as usize).min(n_pts.saturating_sub(1));
    let mut marks = Vec::new();
    for idx in start..=end {
        if idx % major_step == 0 {
            marks.push(egui_plot::GridMark {
                value: idx as f64,
                step_size: major_step as f64,
            });
        } else if minor_step < major_step && idx % minor_step == 0 {
            marks.push(egui_plot::GridMark {
                value: idx as f64,
                step_size: minor_step as f64,
            });
        }
    }
    marks
}

pub(crate) fn format_mmwr_axis_mark(
    mark: egui_plot::GridMark,
    range: &std::ops::RangeInclusive<f64>,
    labels: &[String],
) -> String {
    let span = (range.end() - range.start()).abs().max(1.0);
    let major_step = major_step_for_span(span) as f64;
    if mark.step_size + 1e-6 < major_step {
        return String::new();
    }
    let rounded = mark.value.round();
    if (mark.value - rounded).abs() > 0.05 || rounded < 0.0 {
        return String::new();
    }
    labels
        .get(rounded as usize)
        .cloned()
        .unwrap_or_default()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DomainPillar {
    OutbreakAndRt,
    SpatialAndSyndromic,
    GenomicWastewaterLab,
    ForecastingAuditLedger,
    ConnectorsAndSurveys,
}

impl DomainPillar {
    const ALL: [DomainPillar; 5] = [
        Self::OutbreakAndRt,
        Self::SpatialAndSyndromic,
        Self::GenomicWastewaterLab,
        Self::ForecastingAuditLedger,
        Self::ConnectorsAndSurveys,
    ];

    fn label(self) -> &'static str {
        match self {
            Self::OutbreakAndRt => "1. Outbreak & Rt",
            Self::SpatialAndSyndromic => "2. Spatial & Syndromic",
            Self::GenomicWastewaterLab => "3. Genomic, Wastewater & Lab",
            Self::ForecastingAuditLedger => "4. Forecasting, Audit & Ledger",
            Self::ConnectorsAndSurveys => "5. Connectors & Surveys",
        }
    }

    fn sub_tabs(self) -> &'static [(ActiveTab, &'static str)] {
        match self {
            Self::OutbreakAndRt => &[
                (ActiveTab::Surveillance, "Rt & Aberration"),
                (ActiveTab::SeirNowcast, "SEIR & Delay Nowcast"),
                (ActiveTab::SuperspreadingNet, "Superspreading & Network"),
            ],
            Self::SpatialAndSyndromic => &[
                (ActiveTab::SpatialScan, "Kulldorff Spatial Scan"),
                (ActiveTab::CrcAndSyndromic, "Capture-Recapture & ESSENCE"),
                (ActiveTab::EnvAndMortality, "Env DLNM & EuroMOMO"),
            ],
            Self::GenomicWastewaterLab => &[
                (ActiveTab::GenomicPhylo, "Genomic & Phylodynamics"),
                (ActiveTab::WastewaterAndAmr, "Wastewater & CLSI AMR"),
                (ActiveTab::ElrPprlSero, "ELR, PPRL & Sero/VE"),
            ],
            Self::ForecastingAuditLedger => &[
                (ActiveTab::ForecastScoringSitrep, "WIS Scoring & SitRep"),
                (ActiveTab::SeverityAndSurge, "cCFR & ICU Surge"),
                (ActiveTab::DpAuditAndHpc, "DP Privacy & Audit"),
                (ActiveTab::VintageLedger, "Vintage Ledger & Diff"),
            ],
            Self::ConnectorsAndSurveys => &[
                (ActiveTab::LiveConnectors, "Live CDC / Delphi / WHO APIs"),
                (ActiveTab::LocalImportXpt, "Local CSV / .XPT / FHIR"),
                (ActiveTab::SurveyAndPrr, "Survey (NHANES) & PRR (FAERS)"),
            ],
        }
    }
}

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

impl ActiveTab {
    const ALL_SEARCHABLE: [(ActiveTab, &'static str, &'static str); 16] = [
        (Self::Surveillance, "Rt & Aberration", "cori reproduction aberration farrington ears measles pertussis influenza"),
        (Self::SeirNowcast, "SEIR & Delay Nowcast", "seir compartment reporting delay nowcast forecast"),
        (Self::SuperspreadingNet, "Superspreading & Network", "offspring dispersion k mle contact tracing network"),
        (Self::SpatialScan, "Kulldorff Spatial Scan", "spatial cluster states map fips rr llr"),
        (Self::CrcAndSyndromic, "Capture-Recapture & ESSENCE", "crc chapman essence syndromic chief complaint"),
        (Self::EnvAndMortality, "Env DLNM & EuroMOMO", "dlnm wildfire pm2.5 euromomo excess mortality"),
        (Self::GenomicPhylo, "Genomic & Phylodynamics", "clades variants lineage growth advantage outbreak tree fasta"),
        (Self::WastewaterAndAmr, "Wastewater & CLSI AMR", "wastewater freyja snv antibiogram clsi m39 mdr xdr"),
        (Self::ElrPprlSero, "ELR, PPRL & Sero/VE", "hl7 elr pprl bloom filter serosurvey vaccine effectiveness"),
        (Self::ForecastScoringSitrep, "WIS Scoring & SitRep", "wis interval score ensemble sitrep report"),
        (Self::SeverityAndSurge, "cCFR & ICU Surge", "cfr case fatality hospital icu ventilator surge"),
        (Self::DpAuditAndHpc, "DP Privacy & Audit", "differential privacy laplace epsilon budget audit chain"),
        (Self::VintageLedger, "Vintage Ledger & Diff", "sqlite vintage bitemporal backfill diff"),
        (Self::LiveConnectors, "Live CDC / Delphi / WHO APIs", "soda nndss nwss delphi openfda who gho"),
        (Self::LocalImportXpt, "Local CSV / .XPT / FHIR", "csv tsv sas xpt fhir import geojson"),
        (Self::SurveyAndPrr, "Survey (NHANES) & PRR (FAERS)", "nhanes complex survey taylor faers prr ror ic025"),
    ];

    fn pillar(self) -> DomainPillar {
        match self {
            Self::Surveillance | Self::SeirNowcast | Self::SuperspreadingNet => {
                DomainPillar::OutbreakAndRt
            }
            Self::SpatialScan | Self::CrcAndSyndromic | Self::EnvAndMortality => {
                DomainPillar::SpatialAndSyndromic
            }
            Self::GenomicPhylo | Self::WastewaterAndAmr | Self::ElrPprlSero => {
                DomainPillar::GenomicWastewaterLab
            }
            Self::ForecastScoringSitrep
            | Self::SeverityAndSurge
            | Self::DpAuditAndHpc
            | Self::VintageLedger => DomainPillar::ForecastingAuditLedger,
            Self::LiveConnectors | Self::LocalImportXpt | Self::SurveyAndPrr => {
                DomainPillar::ConnectorsAndSurveys
            }
        }
    }

    fn binding_badge(self) -> (bool, &'static str) {
        match self {
            Self::SurveyAndPrr => (false, "Auxiliary Sub-Dataset: NHANES / FAERS"),
            Self::EnvAndMortality => (false, "Auxiliary Sub-Dataset: Env PM2.5 / EuroMOMO"),
            Self::ElrPprlSero => (false, "Auxiliary Sub-Dataset: ELR / Sero Archive"),
            Self::LiveConnectors => (false, "Auxiliary Sub-Dataset: Live External Feeds"),
            Self::LocalImportXpt => (false, "Auxiliary Sub-Dataset: Local File Import"),
            _ => (true, "Bound to Active Cohort"),
        }
    }
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
    active_campaign: SurveillanceCampaignPreset,
    active_window_cutoff_idx: usize,
    global_search_query: String,
    focus_global_search: bool,
    show_about_modal: bool,
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

        let mut app = Self {
            db_path,
            store,
            active_tab: ActiveTab::Surveillance,
            active_campaign: SurveillanceCampaignPreset::UsMeasles2022_2026,
            active_window_cutoff_idx: SurveillanceCampaignPreset::UsMeasles2022_2026
                .default_active_cutoff_idx(),
            global_search_query: String::new(),
            focus_global_search: false,
            show_about_modal: false,
            status_banner: "Synchronized 5-Domain Surveillance Workbench ready (Active Cohort: US Measles 2022–2026).".to_string(),
            suppression_policy: SuppressionPolicy::Unmasked,
            active_dataset,
            vintages,
            window_weeks: 0,
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
            selected_fips: "48".to_string(),
            cached_clusters,
            max_zone_pop_frac: 0.25,
            cached_zone_clusters,
            highlighted_cluster_members,
            soda_is_nwss: false,
            soda_condition: "Measles".to_string(),
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
        };
        app.apply_campaign_preset(SurveillanceCampaignPreset::UsMeasles2022_2026);
        app
    }

    fn apply_campaign_preset(&mut self, campaign: SurveillanceCampaignPreset) {
        self.active_campaign = campaign;
        let ds = campaign.build_weekly_dataset();
        let vid = ds.vintage_id.clone();
        let _ = self.store.commit_vintage(&ds);
        self.refresh_vintages();
        self.active_dataset = ds;
        self.active_window_cutoff_idx = campaign
            .default_active_cutoff_idx()
            .min(self.active_dataset.points.len().saturating_sub(1));

        self.selected_si_preset = campaign.si_preset_id().to_string();
        if let Some(sp) = pathogen_si_presets()
            .into_iter()
            .find(|p| p.id == campaign.si_preset_id())
        {
            self.rt_cfg.si_mean_days = sp.mean_days;
            self.rt_cfg.si_sd_days = sp.sd_days;
        }
        self.seir_preset_id = campaign.si_preset_id().to_string();
        if let Some(seir_p) = seir_pathogen_presets()
            .into_iter()
            .find(|p| p.id == campaign.si_preset_id())
        {
            self.seir_cfg.r0 = seir_p.r0_default;
            self.seir_cfg.latent_days = seir_p.latent_days;
            self.seir_cfg.infectious_days = seir_p.infectious_days;
            self.seir_cfg.waning_days = seir_p.waning_days;
            self.seir_cfg.initial_immune_fraction = seir_p.default_immune_fraction;
            self.seir_cfg.reporting_fraction = seir_p.reporting_fraction;
        }
        self.soda_condition = match campaign {
            SurveillanceCampaignPreset::UsMeasles2022_2026 => "Measles",
            SurveillanceCampaignPreset::UsPertussis2022_2026 => "Pertussis",
            SurveillanceCampaignPreset::InfluenzaIliWastewater2022_2026 => "Influenza",
            SurveillanceCampaignPreset::SarsCoV2VariantCohort2022_2026 => "COVID-19",
        }
        .to_string();

        self.jurisdictions = campaign.build_jurisdictions();
        self.recompute_analytics();
        self.recompute_spatial();
        if let Some(top) = self.cached_clusters.first() {
            self.selected_fips = top.fips.clone();
        }

        self.genomic.sync_to_campaign(campaign, &self.active_dataset);
        self.wastewater_amr
            .sync_to_campaign(campaign, &self.active_dataset);
        self.severity_surge
            .sync_to_campaign(campaign, &self.active_dataset, &self.cached_seir);
        self.crc_syndromic
            .sync_to_campaign(campaign, &self.active_dataset);
        self.superspreading_net.sync_to_campaign(campaign);
        self.dp_audit.refresh_from_store(&self.store, &vid);
        self.status_banner = format!(
            "Synchronized Active Cohort '{}': pathogen {}, {} weeks, {} US jurisdictions & all sub-modules.",
            campaign.label(),
            campaign.pathogen(),
            self.active_dataset.points.len(),
            self.jurisdictions.len()
        );
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

    fn render_top_bar(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        if ctx.input(|i| i.modifiers.command && i.key_pressed(egui::Key::K)) {
            self.focus_global_search = true;
        }

        // Row 1: 5 Clean Domain Tabs + Global Search & Filter Bar (Ctrl+K)
        let cur_pillar = self.active_tab.pillar();
        ui.horizontal_wrapped(|ui| {
            ui.heading(
                RichText::new("EpiLedger")
                    .strong()
                    .color(Color32::from_rgb(96, 165, 250)),
            );
            ui.separator();

            for pillar in DomainPillar::ALL {
                let selected = cur_pillar == pillar;
                if ui
                    .selectable_label(selected, RichText::new(pillar.label()).strong())
                    .clicked()
                {
                    if cur_pillar != pillar {
                        if let Some(&(first_tab, _)) = pillar.sub_tabs().first() {
                            self.active_tab = first_tab;
                        }
                    }
                }
            }

            ui.separator();
            ui.label(RichText::new("Search (Ctrl+K):").small().weak());
            let search_resp = ui.add(
                egui::TextEdit::singleline(&mut self.global_search_query)
                    .desired_width(175.0)
                    .hint_text("Module, state, clade, pathogen..."),
            );
            if self.focus_global_search {
                search_resp.request_focus();
                self.focus_global_search = false;
            }
            if !self.global_search_query.is_empty() && ui.small_button("Clear").clicked() {
                self.global_search_query.clear();
            }

            let mut mask_bool = self.suppression_policy == SuppressionPolicy::NchsWonderLe9;
            if ui.checkbox(&mut mask_bool, "NCHS <=9 Mask").changed() {
                self.suppression_policy = if mask_bool {
                    SuppressionPolicy::NchsWonderLe9
                } else {
                    SuppressionPolicy::Unmasked
                };
            }
        });

        // Quick matching-module jump strip whenever Ctrl+K Global Search is active
        let q_low = self.global_search_query.trim().to_ascii_lowercase();
        if !q_low.is_empty() {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Matching Modules:").small().color(Color32::from_rgb(56, 189, 248)));
                let mut hits = 0usize;
                for &(tab, title, tags) in &ActiveTab::ALL_SEARCHABLE {
                    let hay = format!("{} {} {}", tab.pillar().label(), title, tags).to_ascii_lowercase();
                    if hay.contains(&q_low) {
                        hits += 1;
                        if ui.small_button(format!("-> {title}")).clicked() {
                            self.active_tab = tab;
                        }
                    }
                }
                if hits == 0 {
                    ui.weak("Filtering active table rows (no module title match).");
                }
            });
        }

        ui.separator();

        // Row 2: Pinned ActiveSurveillanceHeader (Cohort Selector + Metadata + Binding Badge)
        let mut next_campaign: Option<SurveillanceCampaignPreset> = None;
        let active_cutoff_label = self
            .active_dataset
            .points
            .get(self.active_window_cutoff_idx)
            .map(|p| {
                let wk = p
                    .epiweek
                    .map(|w| w.compact_label())
                    .unwrap_or_else(|| "w--".to_string());
                format!("{} ({})", p.date_label, wk)
            })
            .unwrap_or_else(|| "Full Range".to_string());

        ui.horizontal_wrapped(|ui| {
            ui.label(RichText::new("Active Cohort:").small().strong().color(Color32::from_rgb(148, 163, 184)));
            egui::ComboBox::from_id_source("active_campaign_combo")
                .selected_text(self.active_campaign.label())
                .width(245.0)
                .show_ui(ui, |ui| {
                    for preset in SurveillanceCampaignPreset::ALL {
                        if ui
                            .selectable_label(self.active_campaign == preset, preset.label())
                            .clicked()
                        {
                            next_campaign = Some(preset);
                        }
                    }
                });

            ui.separator();
            ui.label(
                RichText::new(format!(
                    "Pathogen: {}  |  Region: {}  |  {} · {}",
                    self.active_campaign.pathogen(),
                    self.active_campaign.region(),
                    self.active_campaign.date_span_iso(),
                    self.active_campaign.mmwr_span_label()
                ))
                .small()
                .monospace(),
            );

            ui.separator();
            let (is_bound, badge_txt) = self.active_tab.binding_badge();
            let (badge_col, badge_bg, prefix) = if is_bound {
                (
                    Color32::from_rgb(52, 211, 153),
                    Color32::from_rgb(6, 78, 59),
                    "[*] ",
                )
            } else {
                (
                    Color32::from_rgb(251, 191, 36),
                    Color32::from_rgb(120, 53, 15),
                    "[o] ",
                )
            };
            egui::Frame::none()
                .fill(badge_bg)
                .rounding(4.0)
                .inner_margin(egui::Margin::symmetric(7.0, 2.0))
                .show(ui, |ui| {
                    ui.label(
                        RichText::new(format!("{prefix}{badge_txt}"))
                            .small()
                            .strong()
                            .color(badge_col),
                    );
                });
        });

        if let Some(c) = next_campaign {
            if c != self.active_campaign {
                self.apply_campaign_preset(c);
            }
        }

        // Row 3: Sub-View Selector for Current Domain + Dual Time-Window Split Controls
        let n_pts = self.active_dataset.points.len().max(2);
        ui.horizontal_wrapped(|ui| {
            ui.label(RichText::new("View:").small().weak());
            for &(sub_tab, sub_label) in cur_pillar.sub_tabs() {
                ui.selectable_value(&mut self.active_tab, sub_tab, sub_label);
            }

            ui.separator();
            ui.label(RichText::new("Baseline / Active Cutoff:").small().weak());
            ui.add(
                egui::Slider::new(&mut self.active_window_cutoff_idx, 12..=(n_pts - 4))
                    .show_value(false),
            );
            ui.monospace(RichText::new(&active_cutoff_label).small().color(Color32::from_rgb(34, 211, 238)));
            if ui.small_button("52w Active").clicked() {
                self.active_window_cutoff_idx = n_pts.saturating_sub(52);
            }
            if ui.small_button("26w Active").clicked() {
                self.active_window_cutoff_idx = n_pts.saturating_sub(26);
            }
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
                    if self.cached_rt.estimates.last().is_some() {
                        ui.label(
                            RichText::new("Gated (n < 12 cases in window)")
                                .small()
                                .weak(),
                        );
                    }
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

        let mut hist_bars = Vec::new();
        let mut active_bars = Vec::new();
        let mut alert_bars = Vec::new();
        let mut upper_pts = Vec::new();
        let mut exp_pts = Vec::new();
        let mut alert_pts = Vec::new();

        let cutoff_idx = self
            .active_window_cutoff_idx
            .clamp(1, total_pts.saturating_sub(1));
        let cutoff_x_opt = if cutoff_idx >= start_idx && cutoff_idx < total_pts {
            Some((cutoff_idx - start_idx) as f64 - 0.5)
        } else {
            None
        };
        let cutoff_label = self
            .active_dataset
            .points
            .get(cutoff_idx)
            .map(|p| format!("Active Outbreak Cutoff ({})", p.date_label))
            .unwrap_or_else(|| "Active Outbreak Cutoff".to_string());

        let x_labels: Arc<Vec<String>> = Arc::new(
            self.active_dataset.points[start_idx..total_pts]
                .iter()
                .map(|p| {
                    p.epiweek
                        .map(|w| w.compact_label())
                        .or_else(|| {
                            parse_flexible_date_to_mmwr(&p.date_label).map(|w| w.compact_label())
                        })
                        .unwrap_or_else(|| p.date_label.clone())
                })
                .collect(),
        );
        let n_vis = total_pts.saturating_sub(start_idx).max(1);
        let left_end = (n_vis / 3).max(1);
        let right_start = n_vis.saturating_sub((n_vis / 3).max(1));
        let mut epi_y_max = 0.0_f64;
        let mut epi_left_max = 0.0_f64;
        let mut epi_right_max = 0.0_f64;

        for (offset, idx) in (start_idx..total_pts).enumerate() {
            let x = offset as f64;
            let raw_pt = &self.active_dataset.points[idx];
            let raw_c = raw_pt.count;
            let disp_c = if self.suppression_policy == SuppressionPolicy::NchsWonderLe9
                && (raw_pt.suppressed || (raw_c > 0.0 && raw_c <= 9.0))
            {
                0.0
            } else {
                raw_c
            };
            let ab = &self.cached_ab.points[idx];
            let local_peak = disp_c
                .max(ab.upper_threshold.unwrap_or(0.0))
                .max(ab.expected_baseline.unwrap_or(0.0));
            epi_y_max = epi_y_max.max(local_peak);
            if offset < left_end {
                epi_left_max = epi_left_max.max(local_peak);
            }
            if offset >= right_start {
                epi_right_max = epi_right_max.max(local_peak);
            }

            if idx < cutoff_idx {
                let col = if ab.is_alert {
                    Color32::from_rgba_unmultiplied(148, 163, 184, 210)
                } else {
                    Color32::from_rgba_unmultiplied(100, 116, 139, 155)
                };
                hist_bars.push(Bar::new(x, disp_c).width(0.74).fill(col));
            } else if ab.is_alert {
                alert_bars.push(
                    Bar::new(x, disp_c)
                        .width(0.74)
                        .fill(Color32::from_rgb(244, 63, 94)),
                );
            } else {
                active_bars.push(
                    Bar::new(x, disp_c)
                        .width(0.74)
                        .fill(Color32::from_rgb(34, 211, 238)),
                );
            }

            if let Some(u) = ab.upper_threshold {
                upper_pts.push([x, u]);
            }
            if let Some(e) = ab.expected_baseline {
                exp_pts.push([x, e]);
            }
            if ab.is_alert && idx >= cutoff_idx {
                alert_pts.push([x, disp_c]);
            }
        }

        let avail_h = (ui.clip_rect().bottom() - ui.cursor().top().max(ui.clip_rect().top()))
            .max(ui.available_height());
        let usable_h = (avail_h - 10.0).max(425.0);
        let epi_h = (usable_h * 0.56).max(230.0);
        let rt_h = (usable_h - epi_h).max(195.0);

        let epi_legend_corner = if epi_left_max <= epi_right_max {
            egui_plot::Corner::LeftTop
        } else {
            egui_plot::Corner::RightTop
        };
        let x_labels_epi = Arc::clone(&x_labels);
        let x_link_id = egui::Id::new("surveillance_x_link");

        let mut epi_plot = Plot::new("epi_curve_plot")
            .height(epi_h)
            .legend(
                Legend::default()
                    .position(epi_legend_corner)
                    .background_alpha(0.85),
            )
            .link_axis(x_link_id, true, false)
            .x_grid_spacer(move |input| mmwr_week_grid_spacer(input, n_vis))
            .x_axis_formatter(move |mark, range| format_mmwr_axis_mark(mark, range, &x_labels_epi));
        if epi_y_max > 0.0 {
            epi_plot = epi_plot.include_y(epi_y_max * 1.18);
        }
        epi_plot.show(ui, |plot_ui| {
            plot_ui.bar_chart(BarChart::new(hist_bars).name("Historical Baseline Archive"));
            plot_ui.bar_chart(BarChart::new(active_bars).name("Current / Active Outbreak Window"));
            if !alert_bars.is_empty() {
                plot_ui.bar_chart(BarChart::new(alert_bars).name("Active Outbreak Alert"));
            }
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
            if !alert_pts.is_empty() {
                plot_ui.points(
                    Points::new(PlotPoints::from(alert_pts))
                        .radius(4.5)
                        .name("Aberration Alert Peak")
                        .color(Color32::from_rgb(251, 113, 133)),
                );
            }
            if let Some(cx) = cutoff_x_opt {
                plot_ui.vline(
                    VLine::new(cx)
                        .name(cutoff_label.clone())
                        .color(Color32::from_rgb(251, 191, 36))
                        .width(2.0),
                );
            }
        });

        ui.add_space(6.0);

        let mut rt_hist_segs: Vec<Vec<[f64; 2]>> = Vec::new();
        let mut rt_mean_segs: Vec<Vec<[f64; 2]>> = Vec::new();
        let mut rt_low_segs: Vec<Vec<[f64; 2]>> = Vec::new();
        let mut rt_high_segs: Vec<Vec<[f64; 2]>> = Vec::new();
        let mut cur_hist: Vec<[f64; 2]> = Vec::new();
        let mut cur_mean: Vec<[f64; 2]> = Vec::new();
        let mut cur_low: Vec<[f64; 2]> = Vec::new();
        let mut cur_high: Vec<[f64; 2]> = Vec::new();
        let mut prev_hist_pt: Option<[f64; 2]> = None;

        let mut rt_y_max = 1.5_f64;
        let mut rt_left_max = 1.0_f64;
        let mut rt_right_max = 1.0_f64;

        for (offset, idx) in (start_idx..total_pts).enumerate() {
            if let Some(est) = self.cached_rt.estimates[idx] {
                let x = offset as f64;
                let mean_y = est.mean.min(8.0);
                let low_y = est.cri_low_95.min(8.0);
                let high_y = est.cri_high_95.min(8.0);
                let local_peak = mean_y.max(high_y);
                rt_y_max = rt_y_max.max(local_peak);
                if offset < left_end {
                    rt_left_max = rt_left_max.max(local_peak);
                }
                if offset >= right_start {
                    rt_right_max = rt_right_max.max(local_peak);
                }

                if idx < cutoff_idx {
                    let pt = [x, mean_y];
                    cur_hist.push(pt);
                    prev_hist_pt = Some(pt);
                } else {
                    if cur_mean.is_empty() {
                        if let Some(prev) = prev_hist_pt.take() {
                            cur_mean.push(prev);
                        }
                    }
                    cur_mean.push([x, mean_y]);
                }
                cur_low.push([x, low_y]);
                cur_high.push([x, high_y]);
            } else {
                prev_hist_pt = None;
                if !cur_hist.is_empty() {
                    rt_hist_segs.push(std::mem::take(&mut cur_hist));
                }
                if !cur_mean.is_empty() {
                    rt_mean_segs.push(std::mem::take(&mut cur_mean));
                }
                if !cur_low.is_empty() {
                    rt_low_segs.push(std::mem::take(&mut cur_low));
                }
                if !cur_high.is_empty() {
                    rt_high_segs.push(std::mem::take(&mut cur_high));
                }
            }
        }
        if !cur_hist.is_empty() {
            rt_hist_segs.push(cur_hist);
        }
        if !cur_mean.is_empty() {
            rt_mean_segs.push(cur_mean);
        }
        if !cur_low.is_empty() {
            rt_low_segs.push(cur_low);
        }
        if !cur_high.is_empty() {
            rt_high_segs.push(cur_high);
        }

        let rt_legend_corner = if rt_left_max <= rt_right_max {
            egui_plot::Corner::LeftTop
        } else {
            egui_plot::Corner::RightTop
        };
        let x_labels_rt = Arc::clone(&x_labels);
        let max_x = (total_pts.saturating_sub(start_idx).max(1)) as f64;
        let eq_pts = vec![[0.0, 1.0], [max_x, 1.0]];
        Plot::new("cori_rt_plot")
            .height(rt_h)
            .legend(
                Legend::default()
                    .position(rt_legend_corner)
                    .background_alpha(0.85),
            )
            .include_y((rt_y_max * 1.18).max(2.0))
            .link_axis(x_link_id, true, false)
            .x_grid_spacer(move |input| mmwr_week_grid_spacer(input, n_vis))
            .x_axis_formatter(move |mark, range| format_mmwr_axis_mark(mark, range, &x_labels_rt))
            .x_axis_label("MMWR Week (YYYY-Www)")
            .show(ui, |plot_ui| {
                plot_ui.line(
                    Line::new(PlotPoints::from(eq_pts))
                        .name("Equilibrium Rt = 1.0")
                        .color(Color32::from_rgb(16, 185, 129)),
                );
                for seg in rt_hist_segs {
                    plot_ui.line(
                        Line::new(PlotPoints::from(seg))
                            .name("Baseline Archive Rt")
                            .color(Color32::from_rgb(148, 163, 184)),
                    );
                }
                for seg in rt_mean_segs {
                    plot_ui.line(
                        Line::new(PlotPoints::from(seg))
                            .name("Active Outbreak Mean Rt")
                            .color(Color32::from_rgb(34, 211, 238))
                            .width(2.0),
                    );
                }
                for seg in rt_high_segs {
                    plot_ui.line(
                        Line::new(PlotPoints::from(seg))
                            .name("95% CrI Upper")
                            .color(Color32::from_rgba_unmultiplied(96, 165, 250, 110)),
                    );
                }
                for seg in rt_low_segs {
                    plot_ui.line(
                        Line::new(PlotPoints::from(seg))
                            .name("95% CrI Lower")
                            .color(Color32::from_rgba_unmultiplied(96, 165, 250, 110)),
                    );
                }
                if let Some(cx) = cutoff_x_opt {
                    plot_ui.vline(
                        VLine::new(cx)
                            .name(cutoff_label.clone())
                            .color(Color32::from_rgb(251, 191, 36))
                            .width(2.0),
                    );
                }
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
                            let row_q = self.global_search_query.trim().to_ascii_lowercase();
                            for c in &self.cached_clusters {
                                if !row_q.is_empty()
                                    && !c.name.to_ascii_lowercase().contains(&row_q)
                                    && !c.fips.to_ascii_lowercase().contains(&row_q)
                                {
                                    continue;
                                }
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
                    let vq = self.global_search_query.trim().to_ascii_lowercase();
                    for v in &self.vintages {
                        if !vq.is_empty()
                            && !v.title.to_ascii_lowercase().contains(&vq)
                            && !v.vintage_id.to_ascii_lowercase().contains(&vq)
                            && !v.source_kind.to_ascii_lowercase().contains(&vq)
                        {
                            continue;
                        }
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
    pub fn configure_screenshot_view(&mut self, view: u8) {
        match view {
            2 => {
                self.apply_campaign_preset(
                    SurveillanceCampaignPreset::SarsCoV2VariantCohort2022_2026,
                );
                self.active_tab = ActiveTab::WastewaterAndAmr;
            }
            3 => {
                self.apply_campaign_preset(SurveillanceCampaignPreset::UsPertussis2022_2026);
                self.active_tab = ActiveTab::DpAuditAndHpc;
                self.dp_audit.trigger_dp_release_demo(
                    &self.store,
                    &self.active_dataset,
                    &mut self.status_banner,
                );
            }
            _ => {
                self.apply_campaign_preset(SurveillanceCampaignPreset::UsMeasles2022_2026);
                self.active_tab = ActiveTab::Surveillance;
                self.window_weeks = 0;
                self.recompute_analytics();
            }
        }
    }

    pub fn update_ui(&mut self, ctx: &egui::Context) {
        self.poll_worker();

        egui::TopBottomPanel::top("top_nav").show(ctx, |ui| {
            self.render_top_bar(ui, ctx);
        });

        egui::TopBottomPanel::bottom("bottom_status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.worker_busy {
                    ui.spinner();
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.small_button("About").clicked() {
                        self.show_about_modal = true;
                    }
                    ui.separator();
                    ui.label(
                        RichText::new("Copyright CRV Science Network")
                            .small()
                            .weak(),
                    );
                    ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
                        ui.add(
                            egui::Label::new(
                                RichText::new(&self.status_banner).small().monospace(),
                            )
                            .truncate(),
                        );
                    });
                });
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
                ActiveTab::GenomicPhylo => {
                    self.genomic
                        .render_with_filter(ui, &mut self.status_banner, &self.global_search_query);
                }
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
                        self.suppression_policy,
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
                    let act = self.wastewater_amr.render_with_filter(
                        ui,
                        &mut self.status_banner,
                        &self.global_search_query,
                    );
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

        if self.show_about_modal {
            let mut open = self.show_about_modal;
            egui::Window::new("About EpiLedger")
                .open(&mut open)
                .collapsible(false)
                .resizable(false)
                .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
                .show(ctx, |ui| {
                    ui.heading("EpiLedger");
                    ui.label(format!("Version {}", env!("CARGO_PKG_VERSION")));
                    ui.separator();
                    ui.label(
                        "Native Rust public health informatics workbench for outbreak surveillance, Cori Rt, Farrington/EARS aberrations, reporting-delay nowcasting, wastewater deconvolution, and differential privacy auditing.",
                    );
                    ui.add_space(6.0);
                    ui.label(RichText::new("Copyright CRV Science Network").strong());
                    ui.label("License: GNU General Public License v3.0 or later (GPL-3.0-or-later)");
                    ui.add_space(4.0);
                    ui.hyperlink_to(
                        "CRV Science Network — https://crvscience.com",
                        "https://crvscience.com",
                    );
                    ui.hyperlink_to(
                        "CRV Science Apps — https://crvscience.com/apps",
                        "https://crvscience.com/apps",
                    );
                    ui.add_space(8.0);
                    if ui.button("Close").clicked() {
                        self.show_about_modal = false;
                    }
                });
            if !open {
                self.show_about_modal = false;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn campaign_preset_synchronizes_all_modules_and_all_16_tabs_render() {
        let ctx = egui::Context::default();
        let mut app = EpiLedgerGuiApp::new_headless(&ctx, PathBuf::from(":memory:"));

        for preset in SurveillanceCampaignPreset::ALL {
            app.apply_campaign_preset(preset);
            assert_eq!(app.active_campaign, preset);
            assert_eq!(app.active_dataset.default_si_preset, preset.si_preset_id());
            assert_eq!(app.jurisdictions.len(), 51);
            assert!(!app.wastewater_amr.ww_result.samples.is_empty());
            assert!(app.genomic.cached_lineage.is_some());
        }

        app.global_search_query = "TX".to_string();
        for &(tab, _, _) in &ActiveTab::ALL_SEARCHABLE {
            app.active_tab = tab;
            let _ = ctx.run(egui::RawInput::default(), |c| {
                app.update_ui(c);
            });
        }
    }

    #[test]
    fn test_mmwr_x_axis_formatter_and_about_modal() {
        let labels: Vec<String> = (1..=253)
            .map(|i| format!("{}-W{:02}", 2022 + (i - 1) / 52, ((i - 1) % 52) + 1))
            .collect();
        let marks = mmwr_week_grid_spacer(
            egui_plot::GridInput {
                bounds: (0.0, 252.0),
                base_step_size: 10.0,
            },
            labels.len(),
        );
        assert!(marks.iter().any(|m| m.value == 0.0 && m.step_size == 26.0));
        assert!(marks.iter().any(|m| m.value == 26.0 && m.step_size == 26.0));
        let major_lbl = format_mmwr_axis_mark(
            egui_plot::GridMark {
                value: 26.0,
                step_size: 26.0,
            },
            &(0.0..=252.0),
            &labels,
        );
        assert_eq!(major_lbl, "2022-W27");
        let minor_lbl = format_mmwr_axis_mark(
            egui_plot::GridMark {
                value: 13.0,
                step_size: 13.0,
            },
            &(0.0..=252.0),
            &labels,
        );
        assert!(minor_lbl.is_empty(), "minor ticks must not emit overlapping labels");

        let ctx = egui::Context::default();
        let mut app = EpiLedgerGuiApp::new_headless(&ctx, PathBuf::from(":memory:"));
        assert!(!app.show_about_modal);
        app.show_about_modal = true;
        let out = ctx.run(egui::RawInput::default(), |c| {
            app.update_ui(c);
        });
        assert!(!out.shapes.is_empty());
    }
}
