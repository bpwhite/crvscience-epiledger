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

use crate::amr_antibiogram::wilson_95_ci_pct;
use crate::math::regularized_gamma_p;
use crate::seir::SeirForecastResult;
use epi_core::EpiError;
use serde::{Deserialize, Serialize};

/// Right-Truncated Delay-Adjusted Case Fatality Risk (`cCFR` via Nishiura Convolution)
/// and Competing-Risk Hospital Acute Ward, `ICU`, and Invasive Mechanical Ventilator (`IMV`)
/// Sojourn & Surge Census Projections.
///
/// References:
/// - Delay-adjusted age- and sex-specific case fatality rates under right truncation: https://pmc.ncbi.nlm.nih.gov/articles/PMC7529598/
/// - Real-time outbreak cCFR adjustment for onset-to-death lag: https://pmc.ncbi.nlm.nih.gov/articles/PMC7258483/
/// - Nowcast-coupled stochastic SEIR incidence to hospital capacity projection: https://pmc.ncbi.nlm.nih.gov/articles/PMC9200328/

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DelayCfrConfig {
    /// Mean onset-to-death (or diagnosis-to-death) delay in time-series steps (days or weeks).
    pub delay_mean: f64,
    /// Standard deviation of onset-to-death delay in time-series steps.
    pub delay_sd: f64,
    /// Maximum support horizon for discretizing the Gamma delay distribution.
    pub max_delay_steps: usize,
}

impl Default for DelayCfrConfig {
    fn default() -> Self {
        Self {
            delay_mean: 13.5,
            delay_sd: 6.8,
            max_delay_steps: 60,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CaseDeathObservation {
    pub period_label: String,
    pub incident_cases: f64,
    pub incident_deaths: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CfrTrajectoryPoint {
    pub step_index: usize,
    pub period_label: String,
    pub cumulative_cases: f64,
    pub cumulative_deaths: f64,
    /// Nishiura resolved-outcome effective denominator u_t = sum_{tau<=t} c_tau * F_{t - tau}.
    pub resolved_denominator_ut: f64,
    pub crude_cfr_pct: f64,
    pub crude_ci_low_95: f64,
    pub crude_ci_high_95: f64,
    pub delay_adjusted_ccfr_pct: f64,
    pub adjusted_ci_low_95: f64,
    pub adjusted_ci_high_95: f64,
    /// Ratio cCFR / crude_CFR = C_cum / u_t quantifying downward right-truncation bias in naive CFR.
    pub truncation_bias_factor: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DelayAdjustedCfrResult {
    pub config: DelayCfrConfig,
    pub final_cumulative_cases: f64,
    pub final_cumulative_deaths: f64,
    pub final_resolved_denominator: f64,
    pub final_crude_cfr_pct: f64,
    pub final_adjusted_ccfr_pct: f64,
    pub final_adjusted_ci_low_95: f64,
    pub final_adjusted_ci_high_95: f64,
    pub peak_truncation_bias_factor: f64,
    pub delay_cdf: Vec<f64>,
    pub trajectory: Vec<CfrTrajectoryPoint>,
}

/// Competing-risk clinical progression & regional hospital bed/ventilator capacity parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HospitalSurgeConfig {
    /// Probability an incident infection (or reported case) requires acute hospital admission.
    pub prob_hospital_admission: f64,
    /// Mean lag in days from symptom onset / case report to hospital admission.
    pub onset_to_admit_mean_days: f64,
    pub onset_to_admit_sd_days: f64,
    /// Competing-risk probability that an admitted patient escalates to ICU.
    pub prob_icu_given_hosp: f64,
    /// Probability an ICU patient requires invasive mechanical ventilation (IMV).
    pub prob_vent_given_icu: f64,
    /// Competing-risk in-hospital fatality probabilities across clinical acuity tiers.
    pub prob_death_ward_only: f64,
    pub prob_death_icu_nonvent: f64,
    pub prob_death_vent: f64,
    /// Mean length-of-stay (LOS) in days for each competing pathway compartment.
    pub ward_only_los_days: f64,
    pub pre_icu_ward_days: f64,
    pub icu_nonvent_los_days: f64,
    pub vent_los_days: f64,
    pub post_icu_stepdown_ward_days: f64,
    /// Regional staffed capacity thresholds for surge breach evaluation.
    pub staffed_ward_beds: f64,
    pub staffed_icu_beds: f64,
    pub available_ventilators: f64,
    /// Whether to drive surge from `true_weekly_infections` (using Infection Hospitalization Ratio)
    /// rather than `reported_weekly_mean` (Case Hospitalization Ratio) when coupled to SEIR.
    pub use_true_seir_infections: bool,
}

impl Default for HospitalSurgeConfig {
    fn default() -> Self {
        Self {
            prob_hospital_admission: 0.045,
            onset_to_admit_mean_days: 4.5,
            onset_to_admit_sd_days: 2.4,
            prob_icu_given_hosp: 0.22,
            prob_vent_given_icu: 0.48,
            prob_death_ward_only: 0.025,
            prob_death_icu_nonvent: 0.14,
            prob_death_vent: 0.42,
            ward_only_los_days: 6.2,
            pre_icu_ward_days: 1.8,
            icu_nonvent_los_days: 7.5,
            vent_los_days: 11.2,
            post_icu_stepdown_ward_days: 4.0,
            staffed_ward_beds: 320.0,
            staffed_icu_beds: 55.0,
            available_ventilators: 32.0,
            use_true_seir_infections: false,
        }
    }
}

/// Aalen-Johansen style competing-risk Cumulative Incidence Functions (CIF) per admitted cohort over days since admission.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompetingRiskCohortPoint {
    pub day_since_admission: usize,
    pub prop_in_ward: f64,
    pub prop_in_icu_nonvent: f64,
    pub prop_on_ventilator: f64,
    pub cif_discharged_alive: f64,
    pub cif_in_hospital_death: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurgeDayProjection {
    pub day_index: usize,
    pub new_admissions_mean: f64,
    pub ward_census_mean: f64,
    pub ward_census_low_95: f64,
    pub ward_census_high_95: f64,
    pub icu_census_mean: f64,
    pub icu_census_low_95: f64,
    pub icu_census_high_95: f64,
    pub vent_census_mean: f64,
    pub vent_census_low_95: f64,
    pub vent_census_high_95: f64,
    pub cumulative_deaths_mean: f64,
    pub ward_breach: bool,
    pub icu_breach: bool,
    pub vent_breach: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HospitalSurgeResult {
    pub implied_hospital_fatality_pct: f64,
    pub implied_infection_fatality_pct: f64,
    pub peak_ward_census: f64,
    pub peak_ward_day: usize,
    pub peak_ward_utilization_pct: f64,
    pub first_ward_breach_day: Option<usize>,
    pub peak_icu_census: f64,
    pub peak_icu_day: usize,
    pub peak_icu_utilization_pct: f64,
    pub first_icu_breach_day: Option<usize>,
    pub peak_vent_census: f64,
    pub peak_vent_day: usize,
    pub peak_vent_utilization_pct: f64,
    pub first_vent_breach_day: Option<usize>,
    pub total_projected_admissions: f64,
    pub total_projected_hospital_deaths: f64,
    pub icu_bed_days_deficit: f64,
    pub cohort_competing_cif: Vec<CompetingRiskCohortPoint>,
    pub daily_projections: Vec<SurgeDayProjection>,
}

/// Build discrete Gamma CDF F_k = P(T <= k + 0.5) for k = 0..=max_steps.
pub fn discrete_gamma_cdf(mean: f64, sd: f64, max_steps: usize) -> Vec<f64> {
    let m = mean.max(0.25);
    let s = sd.max(0.15);
    let shape = (m / s).powi(2).max(0.2);
    let rate = shape / m;
    (0..=max_steps)
        .map(|k| regularized_gamma_p(shape, rate * ((k as f64) + 0.5)).clamp(0.0, 1.0))
        .collect()
}

/// Discrete Gamma PMF f_k = F_k - F_{k-1} normalized to sum to 1 over 0..=max_steps.
pub fn discrete_gamma_pmf(mean: f64, sd: f64, max_steps: usize) -> Vec<f64> {
    let cdf = discrete_gamma_cdf(mean, sd, max_steps);
    let mut pmf = Vec::with_capacity(cdf.len());
    let mut prev = 0.0;
    for &c in &cdf {
        pmf.push((c - prev).max(0.0));
        prev = c;
    }
    let sum: f64 = pmf.iter().sum();
    if sum > 1e-12 {
        for v in &mut pmf {
            *v /= sum;
        }
    } else if !pmf.is_empty() {
        pmf[0] = 1.0;
    }
    pmf
}

/// Compute Nishiura right-truncation delay-adjusted Case Fatality Risk (`cCFR`) trajectory and 95% CIs.
pub fn compute_delay_adjusted_cfr(
    series: &[CaseDeathObservation],
    cfg: &DelayCfrConfig,
) -> Result<DelayAdjustedCfrResult, EpiError> {
    if series.is_empty() {
        return Err(EpiError::Schema(
            "Case-death series cannot be empty for delay-adjusted CFR".to_string(),
        ));
    }

    let max_lag = cfg.max_delay_steps.max(7).min(240);
    let delay_cdf = discrete_gamma_cdf(cfg.delay_mean, cfg.delay_sd, max_lag);

    let n = series.len();
    let mut cum_c = 0.0_f64;
    let mut cum_d = 0.0_f64;
    let mut peak_bias = 1.0_f64;
    let mut trajectory = Vec::with_capacity(n);

    for t in 0..n {
        let c_t = series[t].incident_cases.max(0.0);
        let d_t = series[t].incident_deaths.max(0.0);
        cum_c += c_t;
        cum_d += d_t;

        let mut u_t = 0.0_f64;
        for tau in 0..=t {
            let lag = (t - tau).min(max_lag);
            u_t += series[tau].incident_cases.max(0.0) * delay_cdf[lag];
        }

        let crude_pct = if cum_c > 0.0 {
            ((cum_d / cum_c) * 100.0).clamp(0.0, 100.0)
        } else {
            0.0
        };
        let (crude_lo, crude_hi) = wilson_95_ci_pct(
            (cum_d.round() as usize).min(cum_c.round() as usize),
            cum_c.round().max(1.0) as usize,
        );

        let adj_pct = if u_t > 0.25 {
            ((cum_d / u_t) * 100.0).clamp(0.0, 100.0)
        } else {
            0.0
        };
        let u_int = u_t.round().max(1.0) as usize;
        let d_int = (cum_d.round() as usize).min(u_int);
        let (adj_lo, adj_hi) = if u_t >= 1.0 {
            wilson_95_ci_pct(d_int, u_int)
        } else {
            (0.0, 100.0)
        };

        let bias = if u_t > 1e-6 && cum_c > 0.0 {
            (cum_c / u_t).max(1.0)
        } else {
            1.0
        };
        if cum_d >= 2.0 && bias > peak_bias {
            peak_bias = bias;
        }

        trajectory.push(CfrTrajectoryPoint {
            step_index: t,
            period_label: series[t].period_label.clone(),
            cumulative_cases: cum_c,
            cumulative_deaths: cum_d,
            resolved_denominator_ut: u_t,
            crude_cfr_pct: crude_pct,
            crude_ci_low_95: crude_lo,
            crude_ci_high_95: crude_hi,
            delay_adjusted_ccfr_pct: adj_pct,
            adjusted_ci_low_95: adj_lo,
            adjusted_ci_high_95: adj_hi,
            truncation_bias_factor: bias,
        });
    }

    let last = trajectory.last().unwrap().clone();
    Ok(DelayAdjustedCfrResult {
        config: cfg.clone(),
        final_cumulative_cases: last.cumulative_cases,
        final_cumulative_deaths: last.cumulative_deaths,
        final_resolved_denominator: last.resolved_denominator_ut,
        final_crude_cfr_pct: last.crude_cfr_pct,
        final_adjusted_ccfr_pct: last.delay_adjusted_ccfr_pct,
        final_adjusted_ci_low_95: last.adjusted_ci_low_95,
        final_adjusted_ci_high_95: last.adjusted_ci_high_95,
        peak_truncation_bias_factor: peak_bias,
        delay_cdf,
        trajectory,
    })
}

/// Build competing-risk daily occupancy kernels (tau = 0..=max_tau days post-admission) for Ward, non-vent ICU, Vent, Alive Discharge, and Death.
pub fn build_competing_risk_admission_kernels(
    cfg: &HospitalSurgeConfig,
    max_tau: usize,
) -> Vec<CompetingRiskCohortPoint> {
    let p_icu = cfg.prob_icu_given_hosp.clamp(0.0, 1.0);
    let p_ward_only = 1.0 - p_icu;
    let p_vent_given_icu = cfg.prob_vent_given_icu.clamp(0.0, 1.0);
    let p_icu_nv = p_icu * (1.0 - p_vent_given_icu);
    let p_icu_v = p_icu * p_vent_given_icu;

    let d_ward = cfg.prob_death_ward_only.clamp(0.0, 1.0);
    let d_icu_nv = cfg.prob_death_icu_nonvent.clamp(0.0, 1.0);
    let d_vent = cfg.prob_death_vent.clamp(0.0, 1.0);

    // Gamma sojourn survival functions 1 - F(tau) for each clinical stage (CV ~ 0.45 realistic right-skewed LOS)
    let ward_only_cdf = discrete_gamma_cdf(
        cfg.ward_only_los_days.max(1.0),
        (cfg.ward_only_los_days * 0.45).max(0.8),
        max_tau,
    );
    let pre_icu_pmf = discrete_gamma_pmf(
        cfg.pre_icu_ward_days.max(0.5),
        (cfg.pre_icu_ward_days * 0.5).max(0.4),
        max_tau,
    );
    let pre_icu_cdf = discrete_gamma_cdf(
        cfg.pre_icu_ward_days.max(0.5),
        (cfg.pre_icu_ward_days * 0.5).max(0.4),
        max_tau,
    );
    let icu_nv_pmf = discrete_gamma_pmf(
        cfg.icu_nonvent_los_days.max(1.0),
        (cfg.icu_nonvent_los_days * 0.45).max(0.8),
        max_tau,
    );
    let icu_nv_cdf = discrete_gamma_cdf(
        cfg.icu_nonvent_los_days.max(1.0),
        (cfg.icu_nonvent_los_days * 0.45).max(0.8),
        max_tau,
    );
    let vent_pmf = discrete_gamma_pmf(
        cfg.vent_los_days.max(1.0),
        (cfg.vent_los_days * 0.45).max(0.8),
        max_tau,
    );
    let vent_cdf = discrete_gamma_cdf(
        cfg.vent_los_days.max(1.0),
        (cfg.vent_los_days * 0.45).max(0.8),
        max_tau,
    );
    let stepdown_cdf = discrete_gamma_cdf(
        cfg.post_icu_stepdown_ward_days.max(0.5),
        (cfg.post_icu_stepdown_ward_days * 0.45).max(0.5),
        max_tau,
    );

    // Convolve pre-ICU ward arrival PMF with ICU sojourn survival (1 - F_icu) and step-down ward recovery
    let mut out = Vec::with_capacity(max_tau + 1);

    for tau in 0..=max_tau {
        // 1. Ward-only branch: still in ward with probability (1 - F_ward(tau))
        let w_only_in = p_ward_only * (1.0 - ward_only_cdf[tau]);
        let w_only_exited = p_ward_only * ward_only_cdf[tau];
        let dead_w = w_only_exited * d_ward;
        let disch_w = w_only_exited * (1.0 - d_ward);

        // 2. Pre-ICU ward occupancy before transfer to ICU
        let w_pre_icu = p_icu * (1.0 - pre_icu_cdf[tau]);

        // 3. Non-ventilated ICU occupancy + exit convolution
        let mut in_icu_nv = 0.0;
        let mut exited_icu_nv_at_or_before = 0.0;
        let mut in_stepdown_nv = 0.0;
        for s1 in 0..=tau {
            let p_enter_icu = pre_icu_pmf[s1];
            let rem = tau - s1;
            in_icu_nv += p_icu_nv * p_enter_icu * (1.0 - icu_nv_cdf[rem]);
            exited_icu_nv_at_or_before += p_icu_nv * p_enter_icu * icu_nv_cdf[rem];

            // Survivors exiting ICU at exact lag s2 enter post-ICU step-down ward
            for s2 in 0..=rem {
                let p_exit_icu = icu_nv_pmf[s2];
                let step_rem = rem - s2;
                in_stepdown_nv += p_icu_nv
                    * (1.0 - d_icu_nv)
                    * p_enter_icu
                    * p_exit_icu
                    * (1.0 - stepdown_cdf[step_rem]);
            }
        }
        let dead_icu_nv = exited_icu_nv_at_or_before * d_icu_nv;
        let disch_icu_nv = (exited_icu_nv_at_or_before * (1.0 - d_icu_nv) - in_stepdown_nv).max(0.0);

        // 4. Ventilated ICU occupancy + exit convolution
        let mut in_vent = 0.0;
        let mut exited_vent_at_or_before = 0.0;
        let mut in_stepdown_v = 0.0;
        for s1 in 0..=tau {
            let p_enter_icu = pre_icu_pmf[s1];
            let rem = tau - s1;
            in_vent += p_icu_v * p_enter_icu * (1.0 - vent_cdf[rem]);
            exited_vent_at_or_before += p_icu_v * p_enter_icu * vent_cdf[rem];

            for s2 in 0..=rem {
                let p_exit_v = vent_pmf[s2];
                let step_rem = rem - s2;
                in_stepdown_v += p_icu_v
                    * (1.0 - d_vent)
                    * p_enter_icu
                    * p_exit_v
                    * (1.0 - stepdown_cdf[step_rem]);
            }
        }
        let dead_vent = exited_vent_at_or_before * d_vent;
        let disch_vent = (exited_vent_at_or_before * (1.0 - d_vent) - in_stepdown_v).max(0.0);

        let prop_ward = (w_only_in + w_pre_icu + in_stepdown_nv + in_stepdown_v).clamp(0.0, 1.0);
        let prop_icu_nv = in_icu_nv.clamp(0.0, 1.0);
        let prop_vent = in_vent.clamp(0.0, 1.0);
        let cif_death = (dead_w + dead_icu_nv + dead_vent).clamp(0.0, 1.0);
        let cif_disch = (disch_w + disch_icu_nv + disch_vent).clamp(0.0, 1.0);

        out.push(CompetingRiskCohortPoint {
            day_since_admission: tau,
            prop_in_ward: prop_ward,
            prop_in_icu_nonvent: prop_icu_nv,
            prop_on_ventilator: prop_vent,
            cif_discharged_alive: cif_disch,
            cif_in_hospital_death: cif_death,
        });
    }

    out
}

/// Project daily competing-risk hospital Ward, Total ICU (`non-vent + vent`), and Ventilator bed occupancy
/// from daily incident case trajectories (`mean`, `low_95`, `high_95`).
pub fn project_competing_risk_hospital_surge(
    daily_cases_mean: &[f64],
    daily_cases_low_95: &[f64],
    daily_cases_high_95: &[f64],
    cfg: &HospitalSurgeConfig,
) -> Result<HospitalSurgeResult, EpiError> {
    let n_days = daily_cases_mean.len();
    if n_days == 0 {
        return Err(EpiError::Schema(
            "Daily incidence vector cannot be empty for hospital surge projection".to_string(),
        ));
    }

    let max_kernel_days = 65_usize;
    let admit_pmf = discrete_gamma_pmf(
        cfg.onset_to_admit_mean_days,
        cfg.onset_to_admit_sd_days,
        28,
    );
    let cohort_cif = build_competing_risk_admission_kernels(cfg, max_kernel_days);

    let p_hosp = cfg.prob_hospital_admission.clamp(0.0, 1.0);
    let convolve_admissions = |inc: &[f64]| -> Vec<f64> {
        let mut adm = vec![0.0_f64; n_days];
        for t in 0..n_days {
            let mut s = 0.0;
            for (lag, &w_lag) in admit_pmf.iter().enumerate() {
                if lag <= t {
                    s += inc[t - lag].max(0.0) * p_hosp * w_lag;
                }
            }
            adm[t] = s;
        }
        adm
    };

    let adm_mean = convolve_admissions(daily_cases_mean);
    let low_fallback: Vec<f64> = if daily_cases_low_95.len() == n_days {
        daily_cases_low_95.to_vec()
    } else {
        daily_cases_mean.iter().map(|&v| 0.75 * v).collect()
    };
    let high_fallback: Vec<f64> = if daily_cases_high_95.len() == n_days {
        daily_cases_high_95.to_vec()
    } else {
        daily_cases_mean.iter().map(|&v| 1.25 * v).collect()
    };
    let adm_low = convolve_admissions(&low_fallback);
    let adm_high = convolve_admissions(&high_fallback);

    let convolve_census = |adm: &[f64]| -> Vec<(f64, f64, f64, f64)> {
        let mut series = Vec::with_capacity(n_days);
        for t in 0..n_days {
            let mut w_occ = 0.0;
            let mut icu_occ = 0.0;
            let mut vent_occ = 0.0;
            let mut cum_dead = 0.0;
            for tau in 0..=t {
                let a_src = adm[t - tau];
                if a_src <= 0.0 {
                    continue;
                }
                if tau <= max_kernel_days {
                    let kpt = &cohort_cif[tau];
                    w_occ += a_src * kpt.prop_in_ward;
                    icu_occ += a_src * (kpt.prop_in_icu_nonvent + kpt.prop_on_ventilator);
                    vent_occ += a_src * kpt.prop_on_ventilator;
                    cum_dead += a_src * kpt.cif_in_hospital_death;
                } else {
                    cum_dead += a_src * cohort_cif[max_kernel_days].cif_in_hospital_death;
                }
            }
            series.push((w_occ, icu_occ, vent_occ, cum_dead));
        }
        series
    };

    let cen_mean = convolve_census(&adm_mean);
    let cen_low = convolve_census(&adm_low);
    let cen_high = convolve_census(&adm_high);

    let mut peak_ward = 0.0_f64;
    let mut peak_ward_day = 0_usize;
    let mut first_ward_breach = None;

    let mut peak_icu = 0.0_f64;
    let mut peak_icu_day = 0_usize;
    let mut first_icu_breach = None;

    let mut peak_vent = 0.0_f64;
    let mut peak_vent_day = 0_usize;
    let mut first_vent_breach = None;

    let mut icu_deficit = 0.0_f64;
    let mut daily_projections = Vec::with_capacity(n_days);

    for d in 0..n_days {
        let (wm, im, vm, dm) = cen_mean[d];
        let (wl, il, vl, _) = cen_low[d];
        let (wh, ih, vh, _) = cen_high[d];

        if wm > peak_ward {
            peak_ward = wm;
            peak_ward_day = d;
        }
        if im > peak_icu {
            peak_icu = im;
            peak_icu_day = d;
        }
        if vm > peak_vent {
            peak_vent = vm;
            peak_vent_day = d;
        }

        let w_br = wm > cfg.staffed_ward_beds;
        let i_br = im > cfg.staffed_icu_beds;
        let v_br = vm > cfg.available_ventilators;
        if w_br && first_ward_breach.is_none() {
            first_ward_breach = Some(d);
        }
        if i_br && first_icu_breach.is_none() {
            first_icu_breach = Some(d);
        }
        if v_br && first_vent_breach.is_none() {
            first_vent_breach = Some(d);
        }
        if im > cfg.staffed_icu_beds {
            icu_deficit += im - cfg.staffed_icu_beds;
        }

        daily_projections.push(SurgeDayProjection {
            day_index: d,
            new_admissions_mean: adm_mean[d],
            ward_census_mean: wm,
            ward_census_low_95: wl.min(wm),
            ward_census_high_95: wh.max(wm),
            icu_census_mean: im,
            icu_census_low_95: il.min(im),
            icu_census_high_95: ih.max(im),
            vent_census_mean: vm,
            vent_census_low_95: vl.min(vm),
            vent_census_high_95: vh.max(vm),
            cumulative_deaths_mean: dm,
            ward_breach: w_br,
            icu_breach: i_br,
            vent_breach: v_br,
        });
    }

    let p_icu = cfg.prob_icu_given_hosp.clamp(0.0, 1.0);
    let pv = cfg.prob_vent_given_icu.clamp(0.0, 1.0);
    let exact_hfr = (1.0 - p_icu) * cfg.prob_death_ward_only
        + p_icu * (1.0 - pv) * cfg.prob_death_icu_nonvent
        + p_icu * pv * cfg.prob_death_vent;
    let exact_ifr = p_hosp * exact_hfr;

    let total_adm: f64 = adm_mean.iter().sum();
    let total_deaths = daily_projections
        .last()
        .map(|p| p.cumulative_deaths_mean)
        .unwrap_or(0.0);

    Ok(HospitalSurgeResult {
        implied_hospital_fatality_pct: exact_hfr * 100.0,
        implied_infection_fatality_pct: exact_ifr * 100.0,
        peak_ward_census: peak_ward,
        peak_ward_day,
        peak_ward_utilization_pct: (peak_ward * 100.0) / cfg.staffed_ward_beds.max(1.0),
        first_ward_breach_day: first_ward_breach,
        peak_icu_census: peak_icu,
        peak_icu_day,
        peak_icu_utilization_pct: (peak_icu * 100.0) / cfg.staffed_icu_beds.max(1.0),
        first_icu_breach_day: first_icu_breach,
        peak_vent_census: peak_vent,
        peak_vent_day,
        peak_vent_utilization_pct: (peak_vent * 100.0) / cfg.available_ventilators.max(1.0),
        first_vent_breach_day: first_vent_breach,
        total_projected_admissions: total_adm,
        total_projected_hospital_deaths: total_deaths,
        icu_bed_days_deficit: icu_deficit,
        cohort_competing_cif: cohort_cif,
        daily_projections,
    })
}

/// Couple an existing `SeirForecastResult` directly into the competing-risk hospital/ICU/ventilator surge engine.
pub fn project_hospital_surge_from_seir(
    seir: &SeirForecastResult,
    cfg: &HospitalSurgeConfig,
) -> Result<HospitalSurgeResult, EpiError> {
    if seir.steps.is_empty() {
        return Err(EpiError::Schema(
            "SEIR forecast contains zero steps".to_string(),
        ));
    }
    let mut daily_mean = Vec::with_capacity(seir.steps.len() * 7);
    let mut daily_low = Vec::with_capacity(seir.steps.len() * 7);
    let mut daily_high = Vec::with_capacity(seir.steps.len() * 7);

    for step in &seir.steps {
        let rep_m = step.reported_weekly_mean.max(1e-6);
        let scale_lo = (step.reported_pi_low_95 / rep_m).clamp(0.0, 1.0);
        let scale_hi = (step.reported_pi_high_95 / rep_m).max(1.0);
        let wk_base = if cfg.use_true_seir_infections {
            step.true_weekly_infections
        } else {
            step.reported_weekly_mean
        };
        let d_m = (wk_base / 7.0).max(0.0);
        for _ in 0..7 {
            daily_mean.push(d_m);
            daily_low.push(d_m * scale_lo);
            daily_high.push(d_m * scale_hi);
        }
    }

    project_competing_risk_hospital_surge(&daily_mean, &daily_low, &daily_high, cfg)
}

/// Parse daily or weekly case-fatality CSV (`period_label,incident_cases,incident_deaths`).
pub fn parse_cfr_incidence_deaths_csv(
    csv_text: &str,
) -> Result<Vec<CaseDeathObservation>, EpiError> {
    let mut out = Vec::new();
    for (idx, raw_line) in csv_text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        if idx == 0
            && (cols[0].eq_ignore_ascii_case("period_label")
                || cols[0].eq_ignore_ascii_case("date"))
        {
            continue;
        }
        if cols.len() < 3 {
            return Err(EpiError::Schema(format!(
                "CFR CSV line {} requires 3 columns (period_label,incident_cases,incident_deaths)",
                idx + 1
            )));
        }
        let cases: f64 = cols[1]
            .parse()
            .map_err(|_| EpiError::Schema(format!("Invalid incident_cases '{}'", cols[1])))?;
        let deaths: f64 = cols[2]
            .parse()
            .map_err(|_| EpiError::Schema(format!("Invalid incident_deaths '{}'", cols[2])))?;
        out.push(CaseDeathObservation {
            period_label: cols[0].to_string(),
            incident_cases: cases.max(0.0),
            incident_deaths: deaths.max(0.0),
        });
    }
    if out.is_empty() {
        return Err(EpiError::Schema(
            "No valid case-death observations found in CSV".to_string(),
        ));
    }
    Ok(out)
}

/// Built-in 60-day acute outbreak wave with true infection-to-death cCFR = 2.40% and mean onset-to-death lag = 13.5 days,
/// demonstrating severe early-wave downward bias in naive CFR (~0.9% on Day 20 vs ~2.4% Nishiura cCFR).
pub fn builtin_emerging_wave_cases_and_deaths() -> Vec<CaseDeathObservation> {
    let n_days = 60_usize;
    let true_cfr = 0.024_f64;
    let pmf = discrete_gamma_pmf(13.5, 6.8, 60);

    let mut cases = Vec::with_capacity(n_days);
    for d in 0..n_days {
        // Logistic-bell epidemic wave peaking around Day 28 at ~620 cases/day
        let x = ((d as f64) - 28.0) / 7.2;
        let inc = (620.0 * (-0.5 * x * x).exp() + 12.0).round();
        cases.push(inc);
    }

    let mut series = Vec::with_capacity(n_days);
    for t in 0..n_days {
        let mut exp_d = 0.0;
        for tau in 0..=t {
            exp_d += cases[tau] * true_cfr * pmf[t - tau];
        }
        series.push(CaseDeathObservation {
            period_label: format!("Day-{:02}", t + 1),
            incident_cases: cases[t],
            incident_deaths: (exp_d * 10.0).round() / 10.0,
        });
    }
    series
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nishiura_delay_adjusted_cfr_corrects_early_wave_right_truncation() {
        let wave = builtin_emerging_wave_cases_and_deaths();
        // Truncate right during active exponential growth at Day 22!
        let early_slice = &wave[..22];
        let res = compute_delay_adjusted_cfr(early_slice, &DelayCfrConfig::default())
            .expect("cCFR succeeds");

        // On Day 22, naive crude CFR is severely biased downward (< 1.1%), whereas Nishiura cCFR recovers ~2.4%!
        assert!(
            res.final_crude_cfr_pct < 1.25,
            "Expected right-truncated naive CFR < 1.25%, got {}",
            res.final_crude_cfr_pct
        );
        assert!(
            (res.final_adjusted_ccfr_pct - 2.40).abs() < 0.20,
            "Expected Nishiura cCFR ~ 2.40%, got {}",
            res.final_adjusted_ccfr_pct
        );
        assert!(res.peak_truncation_bias_factor > 2.0);
    }

    #[test]
    fn test_competing_risk_hospital_surge_conserves_cohort_probabilities() {
        let cfg = HospitalSurgeConfig::default();
        let cif = build_competing_risk_admission_kernels(&cfg, 65);
        let tail = cif.last().unwrap();
        let total_mass = tail.prop_in_ward
            + tail.prop_in_icu_nonvent
            + tail.prop_on_ventilator
            + tail.cif_discharged_alive
            + tail.cif_in_hospital_death;
        assert!(
            (total_mass - 1.0).abs() < 0.03,
            "Expected competing-risk cohort probabilities to sum ~ 1.0, got {total_mass}"
        );

        let wave = builtin_emerging_wave_cases_and_deaths();
        let cases_d: Vec<f64> = wave.iter().map(|r| r.incident_cases * 6.0).collect();
        let surge = project_competing_risk_hospital_surge(&cases_d, &[], &[], &cfg)
            .expect("surge succeeds");
        assert!(surge.peak_ward_census > surge.peak_icu_census);
        assert!(surge.peak_icu_census > surge.peak_vent_census);
    }
}
