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

#![forbid(unsafe_code)]

pub mod aberration;
pub mod amr_antibiogram;
pub mod capture_recapture;
pub mod cfr_and_surge;
pub mod dlnm_env;
pub mod dp_risk_digest;
pub mod euromomo_mortality;
pub mod forecast_scoring;
pub mod lineage_logistics;
pub mod math;
pub mod nowcast;
pub mod pharmacovigilance;
pub mod phylo_outbreak;
pub mod rates_and_lag;
pub mod renewal_rt;
pub mod seir;
pub mod sero_tnd_ve;
pub mod spatial_scan_zones;
pub mod superspreading_network;
pub mod survey_weights;
pub mod wastewater_demix;

pub use aberration::{
    compute_aberrations, AberrationConfig, AberrationMethod, AberrationPoint, AberrationSeries,
};
pub use amr_antibiogram::{
    builtin_hospital_amr_isolates, compute_clsi_m39_antibiogram, parse_amr_isolates_csv,
    wilson_95_ci_pct, AmrIsolateRecord, AntibiogramCell, AntibiogramConfig, AntimicrobialTestCall,
    CumulativeAntibiogramResult, IsolateDeduplicationPolicy, OrganismResistanceSummary,
    SirInterpretation,
};
pub use capture_recapture::{
    builtin_three_source_tb_pertussis_crc, compute_two_source_chapman,
    fit_three_source_loglinear_crc, parse_crc_case_lists_csv, LogLinearCrcModelFit,
    ThreeSourceCrcAnalysis, ThreeSourceOverlapCounts, TwoSourceChapmanResult,
};
pub use cfr_and_surge::{
    build_competing_risk_admission_kernels, builtin_emerging_wave_cases_and_deaths,
    compute_delay_adjusted_cfr, discrete_gamma_cdf, discrete_gamma_pmf,
    parse_cfr_incidence_deaths_csv, project_competing_risk_hospital_surge,
    project_hospital_surge_from_seir, CaseDeathObservation, CfrTrajectoryPoint,
    CompetingRiskCohortPoint, DelayAdjustedCfrResult, DelayCfrConfig, HospitalSurgeConfig,
    HospitalSurgeResult, SurgeDayProjection,
};
pub use dlnm_env::{
    builtin_wildfire_pm25_dlnm_dataset, fit_dlnm_quasi_poisson, parse_dlnm_exposure_csv,
    CumulativeExposurePoint, DlnmConfig, DlnmDataset, DlnmFitResult, DlnmObservation,
    ExposureBasisKind, LagBasisKind, LagSlicePoint,
};
pub use dp_risk_digest::{
    builtin_multi_jurisdiction_dp_inputs, builtin_multi_jurisdiction_dp_inputs_for_campaign,
    compute_multi_jurisdiction_dp_risk_digest, export_dp_risk_digest_csv,
    export_dp_risk_digest_json, DpRiskDigestConfig, DpRiskTier, JurisdictionDpRiskSummary,
    JurisdictionSeriesInput, MultiJurisdictionDpRiskDigest,
};
pub use euromomo_mortality::{
    builtin_all_cause_weekly_mortality_dataset, compute_euromomo_excess_mortality,
    parse_weekly_mortality_csv, ExcessMortalityConfig, ExcessMortalityPoint,
    ExcessMortalityResult, MortalityBaselineMode, WeeklyMortalityDataset, WeeklyMortalityRow,
};
pub use forecast_scoring::{
    combine_median_ensemble, combine_weighted_ensemble, compute_wis_11,
    run_rolling_wis_backtest_and_ensemble, BacktestConfig, EnsembleBacktestResult, ModelScorecard,
    QuantileForecast23, StackedHorizonStep, FORECAST_HUB_QUANTILES, WIS_ALPHAS_11,
};
pub use lineage_logistics::{
    builtin_respiratory_lineage_table, fit_multinomial_lineage_logistics, parse_lineage_counts_csv,
    CladeFitnessSummary, LineageBasisMode, LineageCountRow, LineageCountTable,
    LineageLogisticsConfig, LineageLogisticsResult, LineageTimePoint,
};
pub use nowcast::{
    build_reporting_delay_distribution, nowcast_from_reporting_triangle, nowcast_parametric_delay,
    NowcastPoint, NowcastSeries, ParametricNowcastConfig, ReportingTriangle,
};
pub use pharmacovigilance::{compute_disproportionality_2x2, DisproportionalitySignal};
pub use phylo_outbreak::{
    builtin_outbreak_fasta_and_metadata, export_outbreak_graphviz_dot, parse_fasta_and_metadata,
    reconstruct_outbreak_tree, CaseOffspringMetric, OutbreakIsolateInput, OutbreakTreeConfig,
    OutbreakTreeResult, PairwiseGenomicDistance, TransmissionEdge,
};
pub use rates_and_lag::{
    compute_direct_age_standardized_rate, compute_kulldorff_spatial_clusters,
    compute_wastewater_lead_lag, AgeStandardizedResult, LagCorrelationPoint, SpatialClusterScore,
};
pub use renewal_rt::{
    build_discrete_gamma_si, compute_cori_rt, CoriRtConfig, CoriRtEstimate, CoriRtSeries,
};
pub use seir::{
    calibrate_and_project_seir, seir_pathogen_presets, SeirConfig, SeirForecastResult,
    SeirPathogenPreset, SeirWeeklyStep,
};
pub use sero_tnd_ve::{
    builtin_age_serosurvey_dataset, builtin_tnd_vaccine_study_dataset,
    estimate_tnd_vaccine_effectiveness, fit_catalytic_seroprevalence, parse_serosurvey_csv,
    AgeSeroprevalenceStratum, CatalyticCurvePoint, CatalyticFitResult, CatalyticModelKind,
    SerosurveyDataset, TndCellRecord, TndDataset, TndVeResult, TndWaningStratumVe,
};
pub use spatial_scan_zones::{
    build_us_state_scan_zones, haversine_km, run_kulldorff_circular_scan, MultiZoneCluster,
    SpatialScanInputZone,
};
pub use superspreading_network::{
    analyze_contact_tracing_network, branching_extinction_probability,
    builtin_superspreading_outbreak_records, compute_brandes_betweenness,
    compute_undirected_k_cores, export_bipartite_network_dot, fit_negative_binomial_offspring,
    lloyd_smith_cases_for_transmission_share, lloyd_smith_top_share, negbin_pmf,
    parse_contact_tracing_csv, poisson_99_sse_cutoff, poisson_pmf, BipartiteNetworkAnalysis,
    CaseNetworkMetric, ContactTracingRecord, OffspringDispersionResult, OffspringFitConfig,
    OffspringHistogramBin, SerialIntervalSummary, VenueAmplificationMetric, ZeroTruncationMode,
};
pub use survey_weights::{
    compute_taylor_linearized_mean, ComplexSurveyEstimate, LonelyPsuStrategy, SurveyObservation,
};
pub use wastewater_demix::{
    builtin_wastewater_barcode_and_samples, deconvolve_wastewater_sample,
    deconvolve_wastewater_timeseries, parse_wastewater_snv_csv, project_onto_probability_simplex,
    solve_weighted_l1_simplex, BarcodeLibrary, CrypticResidualSite, DepthWeighting,
    LineageAbundanceEstimate, LongitudinalWastewaterResult, SampleDemixResult, SnvObservation,
    WastewaterDemixConfig, WastewaterSampleInput,
};
