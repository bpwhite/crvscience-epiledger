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

use crate::math::{gamma_quantile, normal_quantile};
use crate::renewal_rt::{build_discrete_gamma_si, compute_cori_rt, CoriRtConfig};
use crate::seir::{calibrate_and_project_seir, SeirConfig};
use epi_core::EpiError;
use serde::{Deserialize, Serialize};

/// Proper Scoring Rules for Probabilistic Epidemic Forecasts (`K = 11` Weighted Interval Score `WIS`,
/// `log(y + 1)` Transformed `WIS`, Sharpness / Overprediction / Underprediction Decomposition,
/// Empirical PI Coverage) & Out-of-Sample Multi-Model Probabilistic Ensemble Stacking.
///
/// References:
/// - Bracher J, Ray EL, Gneiting T, Reich NG (2021), Evaluating epidemic forecasts in an interval format: https://pmc.ncbi.nlm.nih.gov/articles/PMC7880475/
/// - Bosse NI et al. (2023), Scoring epidemiological forecasts on transformed scales: https://pmc.ncbi.nlm.nih.gov/articles/PMC10495027/
/// - Taylor JW, Taylor KS (2022), Inverse-score & median interval forecast combining: https://pmc.ncbi.nlm.nih.gov/articles/PMC8963571/
/// - Ray EL et al. (2023), Comparing trained and untrained probabilistic ensemble forecasts: https://pmc.ncbi.nlm.nih.gov/articles/PMC9247236/

/// Standard CDC Forecast Hub 23 probability levels:
/// 11 lower bounds, median `0.50` (index 11), and 11 symmetric upper bounds,
/// corresponding to `K = 11` central prediction intervals at `alpha = [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]`.
pub const FORECAST_HUB_QUANTILES: [f64; 23] = [
    0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550,
    0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 0.975, 0.990,
];

pub const WIS_ALPHAS_11: [f64; 11] = [
    0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
];

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuantileForecast23 {
    /// Non-decreasing quantile predictions at `FORECAST_HUB_QUANTILES` (`[0.01 .. 0.99]`).
    pub quantiles: [f64; 23],
}

impl QuantileForecast23 {
    pub fn median(&self) -> f64 {
        self.quantiles[11]
    }
    pub fn pi50_low(&self) -> f64 {
        self.quantiles[6] // 0.25
    }
    pub fn pi50_high(&self) -> f64 {
        self.quantiles[16] // 0.75
    }
    pub fn pi95_low(&self) -> f64 {
        self.quantiles[1] // 0.025
    }
    pub fn pi95_high(&self) -> f64 {
        self.quantiles[21] // 0.975
    }

    /// Construct a monotonic 23-quantile predictive distribution from mean `mu >= 0` and overdispersed variance `var`.
    pub fn from_mean_and_variance(mean: f64, variance: f64) -> Self {
        let mu = mean.max(0.05);
        let var = variance.max(mu * 1.05);
        // Gamma match on (mean = shape/rate, var = shape/rate^2) -> shape = mu^2 / var, rate = mu / var
        let shape = (mu * mu / var).clamp(0.15, 5_000.0);
        let rate = (mu / var).max(1e-6);

        let mut q = [0.0_f64; 23];
        if shape >= 80.0 {
            let sd = var.sqrt();
            for (idx, &p) in FORECAST_HUB_QUANTILES.iter().enumerate() {
                q[idx] = (mu + normal_quantile(p) * sd).max(0.0);
            }
        } else {
            for (idx, &p) in FORECAST_HUB_QUANTILES.iter().enumerate() {
                q[idx] = gamma_quantile(p, shape, rate).max(0.0);
            }
        }
        // Enforce strict monotonicity across quantile levels
        for i in 1..23 {
            if q[i] < q[i - 1] {
                q[i] = q[i - 1];
            }
        }
        Self { quantiles: q }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct WisBreakdown {
    pub total_wis: f64,
    pub sharpness: f64,
    pub underprediction: f64,
    pub overprediction: f64,
    pub log_wis: f64,
    pub abs_error_median: f64,
    pub covered_50: bool,
    pub covered_95: bool,
}

/// Evaluate exact `K = 11` Weighted Interval Score (`WIS_{0,K}`) and `log(y + 1)` transformed `WIS`
/// per Bracher et al. (2021, PMC7880475) and Bosse et al. (2023, PMC10495027).
pub fn compute_wis_11(fc: &QuantileForecast23, observed_y: f64) -> WisBreakdown {
    let y = observed_y.max(0.0);
    let denom = (WIS_ALPHAS_11.len() as f64) + 0.5; // 11.5

    let eval_scale = |q_arr: &[f64; 23], target: f64| -> (f64, f64, f64) {
        let med = q_arr[11];
        // Median absolute error term w_0 * |y - m| (w_0 = 0.5) splits cleanly into under/over-prediction
        let mut sharp_sum = 0.0;
        let mut under_sum = 0.5 * (target - med).max(0.0);
        let mut over_sum = 0.5 * (med - target).max(0.0);

        for k in 0..11 {
            let alpha = WIS_ALPHAS_11[k];
            let w_k = 0.5 * alpha;
            let l_k = q_arr[k];
            let u_k = q_arr[22 - k];

            // Note: w_k * IS_alpha = (alpha / 2) * (u - l) + max(0, y - u) + max(0, l - y)
            sharp_sum += w_k * (u_k - l_k).max(0.0);
            under_sum += (target - u_k).max(0.0);
            over_sum += (l_k - target).max(0.0);
        }

        (sharp_sum / denom, under_sum / denom, over_sum / denom)
    };

    let (sharp, under, over) = eval_scale(&fc.quantiles, y);

    let mut log_q = [0.0_f64; 23];
    for i in 0..23 {
        log_q[i] = (fc.quantiles[i].max(0.0) + 1.0).ln();
    }
    let (l_sharp, l_under, l_over) = eval_scale(&log_q, (y + 1.0).ln());

    WisBreakdown {
        total_wis: sharp + under + over,
        sharpness: sharp,
        underprediction: under,
        overprediction: over,
        log_wis: l_sharp + l_under + l_over,
        abs_error_median: (y - fc.median()).abs(),
        covered_50: y >= fc.pi50_low() && y <= fc.pi50_high(),
        covered_95: y >= fc.pi95_low() && y <= fc.pi95_high(),
    }
}

/// Combine constituent 23-quantile forecasts via quantile-wise weighted linear pool (`sum_m w_m * q_m(tau)`).
pub fn combine_weighted_ensemble(
    constituents: &[QuantileForecast23],
    weights: &[f64],
) -> QuantileForecast23 {
    if constituents.is_empty() {
        return QuantileForecast23 {
            quantiles: [0.0; 23],
        };
    }
    let w_sum: f64 = weights.iter().copied().sum();
    let mut out = [0.0_f64; 23];
    for tau_i in 0..23 {
        let mut acc = 0.0;
        for (m, fc) in constituents.iter().enumerate() {
            let wm = if w_sum > 1e-12 {
                weights.get(m).copied().unwrap_or(0.0) / w_sum
            } else {
                1.0 / (constituents.len() as f64)
            };
            acc += wm * fc.quantiles[tau_i];
        }
        out[tau_i] = acc.max(0.0);
    }
    for i in 1..23 {
        if out[i] < out[i - 1] {
            out[i] = out[i - 1];
        }
    }
    QuantileForecast23 { quantiles: out }
}

/// Combine constituent 23-quantile forecasts via quantile-wise median (`Untrained Median Ensemble`, PMC8963571).
pub fn combine_median_ensemble(constituents: &[QuantileForecast23]) -> QuantileForecast23 {
    if constituents.is_empty() {
        return QuantileForecast23 {
            quantiles: [0.0; 23],
        };
    }
    let mut out = [0.0_f64; 23];
    let m_len = constituents.len();
    for tau_i in 0..23 {
        let mut vals: Vec<f64> = constituents.iter().map(|c| c.quantiles[tau_i]).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        out[tau_i] = if m_len % 2 == 1 {
            vals[m_len / 2]
        } else {
            0.5 * (vals[m_len / 2 - 1] + vals[m_len / 2])
        };
    }
    for i in 1..23 {
        if out[i] < out[i - 1] {
            out[i] = out[i - 1];
        }
    }
    QuantileForecast23 { quantiles: out }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BacktestConfig {
    pub horizon_weeks: usize,
    pub num_origins: usize,
    pub si_mean_days: f64,
    pub si_sd_days: f64,
}

impl Default for BacktestConfig {
    fn default() -> Self {
        Self {
            horizon_weeks: 4,
            num_origins: 16,
            si_mean_days: 11.7,
            si_sd_days: 2.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModelScorecard {
    pub model_id: String,
    pub model_label: String,
    pub evaluated_forecasts: usize,
    pub mean_wis: f64,
    pub mean_sharpness: f64,
    pub mean_underprediction: f64,
    pub mean_overprediction: f64,
    pub mean_log_wis: f64,
    pub mean_mae: f64,
    /// Relative WIS (`mean_wis / naive_mean_wis`; `< 1.0` beats naive persistence benchmark).
    pub relative_wis_vs_naive: f64,
    pub coverage_50_pct: f64,
    pub coverage_95_pct: f64,
    /// Final learned Inverse-WIS stacking weight (`0.0` for benchmarks/unweighted).
    pub final_ensemble_weight: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StackedHorizonStep {
    pub horizon_week: usize,
    pub inverse_wis_ensemble: QuantileForecast23,
    pub median_ensemble: QuantileForecast23,
    pub seir_median: f64,
    pub renewal_median: f64,
    pub trend_median: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnsembleBacktestResult {
    pub horizon_weeks: usize,
    pub origins_evaluated: usize,
    pub best_model_label: String,
    pub scorecards: Vec<ModelScorecard>,
    /// Final trained weights `[w_seir, w_renewal, w_trend]` summing to `1.0`.
    pub trained_constituent_weights: [f64; 3],
    /// Live forward `1..=horizon_weeks` multi-model forecast from the end of the series.
    pub forward_steps: Vec<StackedHorizonStep>,
}

/// Generate `1..=horizon` 23-quantile forecasts from Naive Persistence + 3 constituent models
/// (`SEIR`, `Cori Renewal`, `Damped Log-Linear Trend`) given history `hist[0..=origin]`.
fn forecast_constituents_at_origin(
    hist: &[f64],
    horizon: usize,
    si_mean_days: f64,
    si_sd_days: f64,
    seir_base: &SeirConfig,
) -> (
    Vec<QuantileForecast23>, // Naive Persistence
    Vec<QuantileForecast23>, // 0: Mechanistic SEIR
    Vec<QuantileForecast23>, // 1: Cori Renewal
    Vec<QuantileForecast23>, // 2: Damped Log-Linear Trend
) {
    let n = hist.len();
    let last_y = hist.last().copied().unwrap_or(1.0).max(0.1);

    // Empirical trailing 8-step first-difference variance for random-walk naive baseline
    let win = n.min(8);
    let mut diff_sq = 0.0;
    for i in (n - win + 1)..n {
        let d = hist[i] - hist[i - 1];
        diff_sq += d * d;
    }
    let rw_step_var = (diff_sq / (win.saturating_sub(1).max(1) as f64)).max(last_y * 1.25);

    let mut naive_fcs = Vec::with_capacity(horizon);
    for h in 1..=horizon {
        naive_fcs.push(QuantileForecast23::from_mean_and_variance(
            last_y,
            rw_step_var * (h as f64),
        ));
    }

    // Constituent 0: Mechanistic SEIR (reduced ensemble runs inside rolling loop for crisp latency)
    let mut seir_cfg = *seir_base;
    seir_cfg.horizon_weeks = horizon;
    seir_cfg.ensemble_runs = 48;
    let seir_res = calibrate_and_project_seir(hist, false, None, &seir_cfg);
    let mut seir_fcs = Vec::with_capacity(horizon);
    for h in 0..horizon {
        if let Some(st) = seir_res.steps.get(h) {
            let mu = st.reported_weekly_mean.max(0.1);
            let approx_sd = ((st.reported_pi_high_95 - st.reported_pi_low_95) / 3.92).max(1.0);
            seir_fcs.push(QuantileForecast23::from_mean_and_variance(
                mu,
                (approx_sd * approx_sd).max(mu * 1.15),
            ));
        } else {
            seir_fcs.push(naive_fcs[h].clone());
        }
    }

    // Constituent 1: Cori Renewal Autoregressive Projection
    let rt_cfg = CoriRtConfig {
        si_mean_days,
        si_sd_days,
        tau_weeks: 3,
        prior_a: 1.0,
        prior_b: 5.0,
    };
    let rt_series = compute_cori_rt(hist, &rt_cfg);
    let (rt_mean, rt_sd) = rt_series
        .estimates
        .iter()
        .rev()
        .flatten()
        .next()
        .map(|e| {
            // Mild damping toward 1.0 over multi-week horizon to avoid explosive divergence
            let m = 0.82 * e.mean.clamp(0.25, 2.8) + 0.18;
            let s = ((e.cri_high_95 - e.cri_low_95) / 3.92).clamp(0.08, 0.85);
            (m, s)
        })
        .unwrap_or((1.0, 0.25));

    let max_wk_lag = 8usize;
    let w_wk = build_discrete_gamma_si(si_mean_days.max(3.5), si_sd_days.max(1.0), 7.0, max_wk_lag);
    let mut ext_counts = hist.to_vec();
    let mut renewal_fcs = Vec::with_capacity(horizon);
    for h in 1..=horizon {
        let cur_t = ext_counts.len();
        let mut inf_potential = 0.0;
        for s in 1..=max_wk_lag {
            if cur_t >= s {
                inf_potential += ext_counts[cur_t - s].max(0.0) * w_wk[s];
            }
        }
        if inf_potential < 0.05 {
            inf_potential = last_y;
        }
        let damp_rt = 1.0 + (rt_mean - 1.0) * 0.88_f64.powi((h - 1) as i32);
        let mu_h = (damp_rt * inf_potential).max(0.1);
        let rel_var = (rt_sd * (h as f64).sqrt()).powi(2);
        let var_h = mu_h * 1.25 + mu_h * mu_h * rel_var;
        ext_counts.push(mu_h);
        renewal_fcs.push(QuantileForecast23::from_mean_and_variance(mu_h, var_h));
    }

    // Constituent 2: Damped Local Quasi-Poisson Log-Linear Trend (trailing 6 weeks)
    let tr_win = n.min(6);
    let slice = &hist[n - tr_win..];
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for (i, &v) in slice.iter().enumerate() {
        let x = i as f64;
        let ly = (v.max(0.0) + 1.0).ln();
        sx += x;
        sy += ly;
        sxx += x * x;
        sxy += x * ly;
    }
    let wf = tr_win as f64;
    let denom = (wf * sxx - sx * sx).max(1e-6);
    let slope = ((wf * sxy - sx * sy) / denom).clamp(-0.32, 0.32);
    let intercept = (sy - slope * sx) / wf;
    let mut resid_sq = 0.0;
    for (i, &v) in slice.iter().enumerate() {
        let pred = (intercept + slope * (i as f64)).exp() - 1.0;
        resid_sq += (v - pred).powi(2) / pred.max(1.0);
    }
    let phi_tr = (resid_sq / (tr_win.saturating_sub(2).max(1) as f64)).clamp(1.1, 12.0);

    let mut trend_fcs = Vec::with_capacity(horizon);
    for h in 1..=horizon {
        // Damped trend factor 0.85^h prevents runaway exponential overshoot
        let eff_steps = (1..=h).map(|k| 0.85_f64.powi(k as i32)).sum::<f64>();
        let x_last = (tr_win - 1) as f64;
        let base_last = (intercept + slope * x_last).exp();
        let mu_h = (base_last * (slope * eff_steps).exp() - 1.0).max(0.1);
        let var_h = phi_tr * mu_h * (1.0 + 0.35 * (h as f64));
        trend_fcs.push(QuantileForecast23::from_mean_and_variance(mu_h, var_h));
    }

    (naive_fcs, seir_fcs, renewal_fcs, trend_fcs)
}

/// Run rolling out-of-sample backtesting across `cfg.num_origins` historical cutpoints,
/// training causal Inverse-WIS stacking weights (`w_m proportional to 1 / WIS_m`) without future lookahead,
/// and emit live forward `1..=horizon_weeks` stacked forecasts from the full series.
pub fn run_rolling_wis_backtest_and_ensemble(
    counts: &[f64],
    cfg: &BacktestConfig,
    seir_base: &SeirConfig,
) -> Result<EnsembleBacktestResult, EpiError> {
    let n = counts.len();
    let horizon = cfg.horizon_weeks.clamp(1, 12);
    if n < horizon + 10 {
        return Err(EpiError::Schema(format!(
            "Backtesting requires at least {} time steps (got {})",
            horizon + 10,
            n
        )));
    }

    let max_last_origin = n - horizon - 1;
    let min_first_origin = 8usize;
    let avail_origins = max_last_origin.saturating_sub(min_first_origin) + 1;
    let num_origins = cfg.num_origins.clamp(3, 40).min(avail_origins);
    let start_origin = max_last_origin + 1 - num_origins;

    // 6 tracked models:
    // 0: Naive Persistence, 1: SEIR, 2: Cori Renewal, 3: Damped Trend, 4: Median Ensemble, 5: Inverse-WIS Ensemble
    let model_meta = [
        ("naive", "Naive Persistence Benchmark"),
        ("seir", "Mechanistic SEIR"),
        ("renewal", "Cori Renewal AR"),
        ("trend", "Damped Quasi-Poisson Trend"),
        ("median_ens", "Untrained Median Ensemble"),
        ("inv_wis_ens", "Trained Inverse-WIS Ensemble"),
    ];

    let mut acc_wis = [0.0_f64; 6];
    let mut acc_sharp = [0.0_f64; 6];
    let mut acc_under = [0.0_f64; 6];
    let mut acc_over = [0.0_f64; 6];
    let mut acc_log_wis = [0.0_f64; 6];
    let mut acc_mae = [0.0_f64; 6];
    let mut acc_cov50 = [0usize; 6];
    let mut acc_cov95 = [0usize; 6];
    let mut total_evals = 0usize;

    // Store historical constituent forecasts `(origin_idx, h_idx, [seir, renewal, trend], realized_y)`
    // so that at origin `t`, Inverse-WIS weights use ONLY targets `origin + h <= t` (strictly causal!).
    struct PastConstituentRecord {
        target_index: usize,
        constituent_wis: [f64; 3],
    }
    let mut completed_records: Vec<PastConstituentRecord> = Vec::new();

    for origin in start_origin..=max_last_origin {
        // Compute strictly causal Inverse-WIS weights using only past targets `<= origin`
        let mut causal_sum_wis = [0.0_f64; 3];
        let mut causal_cnt = 0usize;
        for rec in &completed_records {
            if rec.target_index <= origin {
                for m in 0..3 {
                    causal_sum_wis[m] += rec.constituent_wis[m];
                }
                causal_cnt += 1;
            }
        }
        let weights: [f64; 3] = if causal_cnt >= 2 {
            let inv0 = 1.0 / (causal_sum_wis[0] / (causal_cnt as f64) + 0.25);
            let inv1 = 1.0 / (causal_sum_wis[1] / (causal_cnt as f64) + 0.25);
            let inv2 = 1.0 / (causal_sum_wis[2] / (causal_cnt as f64) + 0.25);
            let s = inv0 + inv1 + inv2;
            [inv0 / s, inv1 / s, inv2 / s]
        } else {
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]
        };

        let hist = &counts[..=origin];
        let (naive_h, seir_h, ren_h, tr_h) = forecast_constituents_at_origin(
            hist,
            horizon,
            cfg.si_mean_days,
            cfg.si_sd_days,
            seir_base,
        );

        for h in 0..horizon {
            let target_idx = origin + 1 + h;
            let y_true = counts[target_idx].max(0.0);
            let const_slice = [seir_h[h].clone(), ren_h[h].clone(), tr_h[h].clone()];
            let med_ens = combine_median_ensemble(&const_slice);
            let inv_ens = combine_weighted_ensemble(&const_slice, &weights);

            let all_six = [
                &naive_h[h],
                &const_slice[0],
                &const_slice[1],
                &const_slice[2],
                &med_ens,
                &inv_ens,
            ];

            let mut c_wis = [0.0_f64; 3];
            for (m_idx, fc) in all_six.iter().enumerate() {
                let sc = compute_wis_11(fc, y_true);
                acc_wis[m_idx] += sc.total_wis;
                acc_sharp[m_idx] += sc.sharpness;
                acc_under[m_idx] += sc.underprediction;
                acc_over[m_idx] += sc.overprediction;
                acc_log_wis[m_idx] += sc.log_wis;
                acc_mae[m_idx] += sc.abs_error_median;
                if sc.covered_50 {
                    acc_cov50[m_idx] += 1;
                }
                if sc.covered_95 {
                    acc_cov95[m_idx] += 1;
                }
                if (1..=3).contains(&m_idx) {
                    c_wis[m_idx - 1] = sc.total_wis;
                }
            }
            completed_records.push(PastConstituentRecord {
                target_index: target_idx,
                constituent_wis: c_wis,
            });
            total_evals += 1;
        }
    }

    // Final trained Inverse-WIS weights using all completed backtest evaluations for live forward projection
    let denom_evals = (total_evals.max(1)) as f64;
    let final_inv = [
        1.0 / (acc_wis[1] / denom_evals + 0.25),
        1.0 / (acc_wis[2] / denom_evals + 0.25),
        1.0 / (acc_wis[3] / denom_evals + 0.25),
    ];
    let inv_sum = final_inv[0] + final_inv[1] + final_inv[2];
    let trained_weights = [
        final_inv[0] / inv_sum,
        final_inv[1] / inv_sum,
        final_inv[2] / inv_sum,
    ];

    let naive_mean_wis = (acc_wis[0] / denom_evals).max(1e-6);
    let mut scorecards = Vec::with_capacity(6);
    for m_idx in 0..6 {
        let mw = acc_wis[m_idx] / denom_evals;
        let w_assigned = if (1..=3).contains(&m_idx) {
            trained_weights[m_idx - 1]
        } else {
            0.0
        };
        scorecards.push(ModelScorecard {
            model_id: model_meta[m_idx].0.to_string(),
            model_label: model_meta[m_idx].1.to_string(),
            evaluated_forecasts: total_evals,
            mean_wis: mw,
            mean_sharpness: acc_sharp[m_idx] / denom_evals,
            mean_underprediction: acc_under[m_idx] / denom_evals,
            mean_overprediction: acc_over[m_idx] / denom_evals,
            mean_log_wis: acc_log_wis[m_idx] / denom_evals,
            mean_mae: acc_mae[m_idx] / denom_evals,
            relative_wis_vs_naive: mw / naive_mean_wis,
            coverage_50_pct: (acc_cov50[m_idx] as f64 / denom_evals) * 100.0,
            coverage_95_pct: (acc_cov95[m_idx] as f64 / denom_evals) * 100.0,
            final_ensemble_weight: w_assigned,
        });
    }

    let best_model_label = scorecards
        .iter()
        .min_by(|a, b| {
            a.mean_wis
                .partial_cmp(&b.mean_wis)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|s| s.model_label.clone())
        .unwrap_or_else(|| "Trained Inverse-WIS Ensemble".to_string());

    // Live forward stacked forecast from full series `counts`
    let (_, f_seir, f_ren, f_tr) = forecast_constituents_at_origin(
        counts,
        horizon,
        cfg.si_mean_days,
        cfg.si_sd_days,
        seir_base,
    );
    let mut forward_steps = Vec::with_capacity(horizon);
    for h in 0..horizon {
        let c_slice = [f_seir[h].clone(), f_ren[h].clone(), f_tr[h].clone()];
        let inv_fc = combine_weighted_ensemble(&c_slice, &trained_weights);
        let med_fc = combine_median_ensemble(&c_slice);
        forward_steps.push(StackedHorizonStep {
            horizon_week: h + 1,
            inverse_wis_ensemble: inv_fc,
            median_ensemble: med_fc,
            seir_median: c_slice[0].median(),
            renewal_median: c_slice[1].median(),
            trend_median: c_slice[2].median(),
        });
    }

    Ok(EnsembleBacktestResult {
        horizon_weeks: horizon,
        origins_evaluated: num_origins,
        best_model_label,
        scorecards,
        trained_constituent_weights: trained_weights,
        forward_steps,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use epi_core::builtin_measles_weekly_dataset;

    #[test]
    fn test_wis_11_exact_additive_decomposition() {
        let fc = QuantileForecast23::from_mean_and_variance(85.0, 220.0);
        for &y in &[20.0, 85.0, 165.0] {
            let w = compute_wis_11(&fc, y);
            let sum_parts = w.sharpness + w.underprediction + w.overprediction;
            assert!((w.total_wis - sum_parts).abs() < 1e-11);
            assert!(w.log_wis > 0.0);
        }
        // Score must be minimized near the predictive median
        let w_center = compute_wis_11(&fc, fc.median());
        let w_tail = compute_wis_11(&fc, fc.median() + 120.0);
        assert!(w_center.total_wis < w_tail.total_wis);
    }

    #[test]
    fn test_rolling_inverse_wis_backtest_on_measles() {
        let ds = builtin_measles_weekly_dataset();
        let counts: Vec<f64> = ds.points.iter().map(|p| p.count).collect();
        let cfg = BacktestConfig {
            horizon_weeks: 4,
            num_origins: 12,
            si_mean_days: 11.7,
            si_sd_days: 2.0,
        };
        let res =
            run_rolling_wis_backtest_and_ensemble(&counts, &cfg, &SeirConfig::default()).unwrap();
        assert_eq!(res.scorecards.len(), 6);
        assert_eq!(res.forward_steps.len(), 4);
        let w_sum: f64 = res.trained_constituent_weights.iter().sum();
        assert!((w_sum - 1.0).abs() < 1e-9);
    }
}
