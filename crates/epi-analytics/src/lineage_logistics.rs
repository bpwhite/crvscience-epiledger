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

use crate::math::normal_quantile;
use crate::renewal_rt::{compute_cori_rt, CoriRtConfig, CoriRtEstimate};
use epi_core::EpiError;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// Multinomial logistic variant competition, natural cubic spline transition dynamics,
/// short-term clade share forecasting, and joint variant-specific reproduction numbers (`R_{t,k}`).
///
/// References:
/// - Multinomial logistic & cubic-spline variant transition dynamics: https://pmc.ncbi.nlm.nih.gov/articles/PMC10065418/
/// - Short-term evolutionary forecasting via multinomial logistic fitness models: https://pmc.ncbi.nlm.nih.gov/articles/PMC10705624/
/// - Joint inference of variant-specific Rt and growth advantage: https://pmc.ncbi.nlm.nih.gov/articles/PMC13349379/

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LineageBasisMode {
    /// Log-linear relative fitness model `eta_{t,k} = alpha_k + s_k * (t - t_0)`.
    LogLinearFitness,
    /// 3-basis natural cubic spline over normalized time for non-monotonic wave transitions.
    NaturalCubicSpline,
}

impl LineageBasisMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::LogLinearFitness => "Multinomial Log-Linear Fitness",
            Self::NaturalCubicSpline => "Multinomial Natural Cubic Spline (df=3)",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageCountRow {
    pub period_label: String,
    pub total_cases: f64,
    pub clade_counts: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageCountTable {
    pub title: String,
    pub clade_names: Vec<String>,
    pub rows: Vec<LineageCountRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageLogisticsConfig {
    pub basis_mode: LineageBasisMode,
    pub reference_clade_index: usize,
    pub forecast_weeks: usize,
    pub generation_time_days: f64,
    pub generation_sd_days: f64,
    pub ridge_l2_penalty: f64,
    pub max_irls_iters: usize,
}

impl Default for LineageLogisticsConfig {
    fn default() -> Self {
        Self {
            basis_mode: LineageBasisMode::LogLinearFitness,
            reference_clade_index: 0,
            forecast_weeks: 6,
            generation_time_days: 5.9,
            generation_sd_days: 4.8,
            ridge_l2_penalty: 1e-3,
            max_irls_iters: 60,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CladeFitnessSummary {
    pub clade_name: String,
    pub is_reference: bool,
    pub total_sequences: f64,
    pub latest_observed_share: f64,
    pub latest_fitted_share: f64,
    pub forecast_end_share: f64,
    /// Weekly logistic growth advantage `s_k` relative to reference clade (`week^-1`).
    pub growth_advantage_per_week: f64,
    pub growth_se_per_week: f64,
    pub growth_ci_low_95: f64,
    pub growth_ci_high_95: f64,
    /// Multiplicative selection coefficient per generation `exp(s_k * T_g / 7) - 1`.
    pub multiplicative_advantage_pct: f64,
    /// Estimated week index (0-based relative to start of table) where clade crosses 50% vs reference.
    pub crossover_week_vs_ref: Option<f64>,
    pub peak_fitted_share: f64,
    pub peak_period_label: String,
    pub latest_variant_rt: Option<CoriRtEstimate>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageTimePoint {
    pub step_index: usize,
    pub period_label: String,
    pub is_forecast: bool,
    pub total_sequences: f64,
    pub aggregate_cases: f64,
    pub observed_shares: Vec<Option<f64>>,
    pub fitted_shares: Vec<f64>,
    pub share_ci_low_95: Vec<f64>,
    pub share_ci_high_95: Vec<f64>,
    pub decomposed_cases: Vec<f64>,
    pub variant_rt: Vec<Option<CoriRtEstimate>>,
    /// Shannon diversity `H = -sum_k p_k log2(p_k)` in bits.
    pub shannon_bits: f64,
    /// Effective number of circulating clades `2^H`.
    pub effective_clades: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageLogisticsResult {
    pub clade_names: Vec<String>,
    pub reference_clade: String,
    pub basis_mode: LineageBasisMode,
    pub log_likelihood: f64,
    pub aic: f64,
    pub converged: bool,
    pub iterations_used: usize,
    pub clade_summaries: Vec<CladeFitnessSummary>,
    pub timeline: Vec<LineageTimePoint>,
}

/// Evaluate design basis vector `x(t)` of length `P` (`P = 2` for log-linear, `P = 4` for 3-df natural cubic spline).
fn basis_vector(t: f64, t_min: f64, t_max: f64, mode: LineageBasisMode) -> Vec<f64> {
    let span = (t_max - t_min).max(1.0);
    let u = (t - t_min) / span;
    match mode {
        LineageBasisMode::LogLinearFitness => {
            let t_mid = 0.5 * (t_min + t_max);
            vec![1.0, t - t_mid]
        }
        LineageBasisMode::NaturalCubicSpline => {
            let k1 = 0.0_f64;
            let k2 = 1.0 / 3.0;
            let k3 = 2.0 / 3.0;
            let k4 = 1.0_f64;
            let pos3 = |v: f64| if v > 0.0 { v * v * v } else { 0.0 };
            let d = |v: f64, knot: f64| -> f64 {
                (pos3(v - knot) - pos3(v - k4)) / (k4 - knot).max(1e-9)
            };
            let d1 = d(u, k1);
            let b2 = d(u, k2) - d1;
            let b3 = d(u, k3) - d1;
            vec![1.0, u, b2, b3]
        }
    }
}

/// Solve symmetric positive-definite linear system `A x = b` via Gauss-Jordan elimination with partial pivoting,
/// returning `(solution, A^-1)` for asymptotic covariance estimation.
fn solve_spd_and_invert(a_in: &[Vec<f64>], b_in: &[f64]) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = b_in.len();
    if n == 0 || a_in.len() != n {
        return None;
    }
    let mut aug = vec![vec![0.0; 2 * n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a_in[i][j];
        }
        aug[i][n + i] = 1.0;
        aug[i][2 * n] = b_in[i];
    }

    for col in 0..n {
        let mut pivot = col;
        let mut max_abs = aug[col][col].abs();
        for row in (col + 1)..n {
            if aug[row][col].abs() > max_abs {
                max_abs = aug[row][col].abs();
                pivot = row;
            }
        }
        if max_abs < 1e-13 {
            return None;
        }
        if pivot != col {
            aug.swap(col, pivot);
        }
        let diag = aug[col][col];
        for j in col..=(2 * n) {
            aug[col][j] /= diag;
        }
        for row in 0..n {
            if row != col {
                let factor = aug[row][col];
                if factor.abs() > 0.0 {
                    for j in col..=(2 * n) {
                        aug[row][j] -= factor * aug[col][j];
                    }
                }
            }
        }
    }

    let mut sol = vec![0.0; n];
    let mut inv = vec![vec![0.0; n]; n];
    for i in 0..n {
        sol[i] = aug[i][2 * n];
        for j in 0..n {
            inv[i][j] = aug[i][n + j];
        }
    }
    Some((sol, inv))
}

fn softmax_with_ref(etas: &[f64], ref_idx: usize) -> Vec<f64> {
    let k_len = etas.len() + 1;
    let mut full = Vec::with_capacity(k_len);
    let mut ptr = 0;
    for k in 0..k_len {
        if k == ref_idx {
            full.push(0.0);
        } else {
            full.push(etas[ptr]);
            ptr += 1;
        }
    }
    let max_e = full.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut exps: Vec<f64> = full.iter().map(|&e| (e - max_e).exp()).collect();
    let sum_e: f64 = exps.iter().sum();
    if sum_e > 0.0 {
        for v in &mut exps {
            *v /= sum_e;
        }
    }
    exps
}

/// Fit a multinomial logistic variant dynamics model, project short-term clade shares,
/// and jointly estimate variant-specific effective reproduction numbers `R_{t,k}`.
pub fn fit_multinomial_lineage_logistics(
    table: &LineageCountTable,
    cfg: &LineageLogisticsConfig,
) -> Result<LineageLogisticsResult, EpiError> {
    let k_total = table.clade_names.len();
    let t_obs = table.rows.len();
    if k_total < 2 {
        return Err(EpiError::Schema(
            "Multinomial lineage model requires at least 2 clades".to_string(),
        ));
    }
    if t_obs < 3 {
        return Err(EpiError::Schema(
            "Multinomial lineage model requires at least 3 time periods".to_string(),
        ));
    }
    for r in &table.rows {
        if r.clade_counts.len() != k_total {
            return Err(EpiError::Schema(format!(
                "Row '{}' has {} clade counts, expected {}",
                r.period_label,
                r.clade_counts.len(),
                k_total
            )));
        }
    }

    let ref_idx = cfg.reference_clade_index.min(k_total - 1);
    let non_ref_indices: Vec<usize> = (0..k_total).filter(|&k| k != ref_idx).collect();
    let m_free = non_ref_indices.len();

    let t_min = 0.0_f64;
    let t_max = (t_obs.saturating_sub(1)) as f64;
    let basis_rows: Vec<Vec<f64>> = (0..t_obs)
        .map(|t| basis_vector(t as f64, t_min, t_max, cfg.basis_mode))
        .collect();
    let p_dim = basis_rows[0].len();
    let total_params = m_free * p_dim;

    let mut beta = vec![0.0_f64; total_params];
    let ridge = cfg.ridge_l2_penalty.max(1e-6);

    let eval_log_lik = |b_vec: &[f64]| -> f64 {
        let mut ll = 0.0;
        for (t, row) in table.rows.iter().enumerate() {
            let x = &basis_rows[t];
            let mut etas = vec![0.0; m_free];
            for c_pos in 0..m_free {
                let mut dot = 0.0;
                for p in 0..p_dim {
                    dot += b_vec[c_pos * p_dim + p] * x[p];
                }
                etas[c_pos] = dot;
            }
            let probs = softmax_with_ref(&etas, ref_idx);
            for k in 0..k_total {
                let y = row.clade_counts[k].max(0.0);
                if y > 0.0 {
                    ll += y * probs[k].max(1e-15).ln();
                }
            }
        }
        for c_pos in 0..m_free {
            for p in 1..p_dim {
                let v = b_vec[c_pos * p_dim + p];
                ll -= 0.5 * ridge * v * v;
            }
        }
        ll
    };

    let mut converged = false;
    let mut iters_used = 0usize;
    let mut cov_matrix = vec![vec![0.0; total_params]; total_params];

    for iter in 1..=cfg.max_irls_iters.max(10) {
        iters_used = iter;
        let mut grad = vec![0.0_f64; total_params];
        let mut info = vec![vec![0.0_f64; total_params]; total_params];

        for (t, row) in table.rows.iter().enumerate() {
            let n_t: f64 = row.clade_counts.iter().map(|v| v.max(0.0)).sum();
            if n_t <= 0.0 {
                continue;
            }
            let x = &basis_rows[t];
            let mut etas = vec![0.0; m_free];
            for c_pos in 0..m_free {
                let mut dot = 0.0;
                for p in 0..p_dim {
                    dot += beta[c_pos * p_dim + p] * x[p];
                }
                etas[c_pos] = dot;
            }
            let probs = softmax_with_ref(&etas, ref_idx);

            for a_pos in 0..m_free {
                let k_a = non_ref_indices[a_pos];
                let diff = row.clade_counts[k_a].max(0.0) - n_t * probs[k_a];
                for p in 0..p_dim {
                    grad[a_pos * p_dim + p] += diff * x[p];
                }

                for b_pos in 0..m_free {
                    let k_b = non_ref_indices[b_pos];
                    let w_ab = if a_pos == b_pos {
                        n_t * probs[k_a] * (1.0 - probs[k_a])
                    } else {
                        -n_t * probs[k_a] * probs[k_b]
                    };
                    for p in 0..p_dim {
                        for q in 0..p_dim {
                            info[a_pos * p_dim + p][b_pos * p_dim + q] += w_ab * x[p] * x[q];
                        }
                    }
                }
            }
        }

        for c_pos in 0..m_free {
            info[c_pos * p_dim][c_pos * p_dim] += 1e-6;
            for p in 1..p_dim {
                let idx = c_pos * p_dim + p;
                grad[idx] -= ridge * beta[idx];
                info[idx][idx] += ridge;
            }
        }

        let Some((step, inv)) = solve_spd_and_invert(&info, &grad) else {
            break;
        };
        cov_matrix = inv;

        let cur_ll = eval_log_lik(&beta);
        let mut alpha_step = 1.0;
        let mut accepted = false;
        for _ in 0..10 {
            let cand: Vec<f64> = beta
                .iter()
                .zip(&step)
                .map(|(b, d)| b + alpha_step * d)
                .collect();
            if eval_log_lik(&cand) >= cur_ll - 1e-9 {
                beta = cand;
                accepted = true;
                break;
            }
            alpha_step *= 0.5;
        }
        if !accepted {
            break;
        }

        let max_delta = step
            .iter()
            .map(|d| (d * alpha_step).abs())
            .fold(0.0_f64, f64::max);
        if max_delta < 1e-6 {
            converged = true;
            break;
        }
    }

    let mut raw_ll = 0.0;
    for (t, row) in table.rows.iter().enumerate() {
        let x = &basis_rows[t];
        let mut etas = vec![0.0; m_free];
        for c_pos in 0..m_free {
            etas[c_pos] = (0..p_dim).map(|p| beta[c_pos * p_dim + p] * x[p]).sum();
        }
        let probs = softmax_with_ref(&etas, ref_idx);
        for k in 0..k_total {
            let y = row.clade_counts[k].max(0.0);
            if y > 0.0 {
                raw_ll += y * probs[k].max(1e-15).ln();
            }
        }
    }
    let aic = 2.0 * (total_params as f64) - 2.0 * raw_ll;

    let recent_win = t_obs.min(5);
    let last_cases = table
        .rows
        .last()
        .map(|r| r.total_cases.max(0.0))
        .unwrap_or(0.0);
    let first_recent = table.rows[t_obs - recent_win].total_cases.max(1.0);
    let case_weekly_ratio = if recent_win > 1 && last_cases > 0.0 {
        ((last_cases.max(1.0) / first_recent).ln() / ((recent_win - 1) as f64))
            .clamp(-0.35, 0.35)
            .exp()
    } else {
        1.0
    };

    let total_steps = t_obs + cfg.forecast_weeks;
    let z95 = normal_quantile(0.975);

    let mut fitted_shares_all: Vec<Vec<f64>> = Vec::with_capacity(total_steps);
    let mut ci_low_all: Vec<Vec<f64>> = Vec::with_capacity(total_steps);
    let mut ci_high_all: Vec<Vec<f64>> = Vec::with_capacity(total_steps);
    let mut agg_cases_all: Vec<f64> = Vec::with_capacity(total_steps);

    for step in 0..total_steps {
        let x = basis_vector(step as f64, t_min, t_max, cfg.basis_mode);
        let mut etas = vec![0.0; m_free];
        for c_pos in 0..m_free {
            etas[c_pos] = (0..p_dim).map(|p| beta[c_pos * p_dim + p] * x[p]).sum();
        }
        let probs = softmax_with_ref(&etas, ref_idx);

        let mut lows = vec![0.0; k_total];
        let mut highs = vec![0.0; k_total];
        for k in 0..k_total {
            let pk = probs[k].clamp(1e-8, 1.0 - 1e-8);
            let mut g_logit = vec![0.0; total_params];
            for a_pos in 0..m_free {
                let k_a = non_ref_indices[a_pos];
                let d_pk_d_eta = if k == k_a {
                    pk * (1.0 - pk)
                } else {
                    -pk * probs[k_a]
                };
                let scale = d_pk_d_eta / (pk * (1.0 - pk));
                for p in 0..p_dim {
                    g_logit[a_pos * p_dim + p] = scale * x[p];
                }
            }
            let mut var_logit = 0.0;
            for i in 0..total_params {
                if g_logit[i].abs() > 0.0 {
                    for j in 0..total_params {
                        var_logit += g_logit[i] * cov_matrix[i][j] * g_logit[j];
                    }
                }
            }
            let se_logit = var_logit.max(0.0).sqrt();
            let logit_k = (pk / (1.0 - pk)).ln();
            lows[k] = 1.0 / (1.0 + (-(logit_k - z95 * se_logit)).exp());
            highs[k] = 1.0 / (1.0 + (-(logit_k + z95 * se_logit)).exp());
        }

        let cases_t = if step < t_obs {
            let raw = table.rows[step].total_cases;
            if raw > 0.0 {
                raw
            } else {
                table.rows[step].clade_counts.iter().sum()
            }
        } else {
            let h = (step - t_obs + 1) as i32;
            last_cases * case_weekly_ratio.powi(h)
        };

        fitted_shares_all.push(probs);
        ci_low_all.push(lows);
        ci_high_all.push(highs);
        agg_cases_all.push(cases_t);
    }

    let rt_cfg = CoriRtConfig {
        si_mean_days: cfg.generation_time_days.max(1.5),
        si_sd_days: cfg.generation_sd_days.max(0.5),
        tau_weeks: 3,
        prior_a: 1.0,
        prior_b: 5.0,
    };
    let mut variant_rt_series: Vec<Vec<Option<CoriRtEstimate>>> = Vec::with_capacity(k_total);
    for k in 0..k_total {
        let clade_inc: Vec<f64> = (0..t_obs)
            .map(|t| agg_cases_all[t] * fitted_shares_all[t][k])
            .collect();
        let rt_res = compute_cori_rt(&clade_inc, &rt_cfg);
        variant_rt_series.push(rt_res.estimates);
    }

    let mut timeline = Vec::with_capacity(total_steps);
    for step in 0..total_steps {
        let is_forecast = step >= t_obs;
        let period_label = if !is_forecast {
            table.rows[step].period_label.clone()
        } else {
            format!("Forecast +{}w", step - t_obs + 1)
        };
        let total_seq = if !is_forecast {
            table.rows[step]
                .clade_counts
                .iter()
                .map(|v| v.max(0.0))
                .sum()
        } else {
            0.0
        };
        let obs_shares = if !is_forecast && total_seq > 0.0 {
            table.rows[step]
                .clade_counts
                .iter()
                .map(|v| Some(v.max(0.0) / total_seq))
                .collect()
        } else {
            vec![None; k_total]
        };
        let probs = &fitted_shares_all[step];
        let decomp: Vec<f64> = probs.iter().map(|&p| p * agg_cases_all[step]).collect();
        let vrt: Vec<Option<CoriRtEstimate>> = (0..k_total)
            .map(|k| {
                if step < t_obs {
                    variant_rt_series[k][step]
                } else {
                    None
                }
            })
            .collect();

        let mut shannon_bits = 0.0;
        for &pk in probs {
            if pk > 1e-12 {
                shannon_bits -= pk * pk.log2();
            }
        }
        let effective_clades = 2.0_f64.powf(shannon_bits);

        timeline.push(LineageTimePoint {
            step_index: step,
            period_label,
            is_forecast,
            total_sequences: total_seq,
            aggregate_cases: agg_cases_all[step],
            observed_shares: obs_shares,
            fitted_shares: probs.clone(),
            share_ci_low_95: ci_low_all[step].clone(),
            share_ci_high_95: ci_high_all[step].clone(),
            decomposed_cases: decomp,
            variant_rt: vrt,
            shannon_bits,
            effective_clades,
        });
    }

    let gen_weeks = (cfg.generation_time_days / 7.0).max(0.1);
    let mut clade_summaries = Vec::with_capacity(k_total);
    for k in 0..k_total {
        let total_seq: f64 = table
            .rows
            .iter()
            .map(|r| r.clade_counts[k].max(0.0))
            .sum();
        let latest_obs = timeline[t_obs - 1].observed_shares[k].unwrap_or(0.0);
        let latest_fit = timeline[t_obs - 1].fitted_shares[k];
        let end_fit = timeline.last().map(|pt| pt.fitted_shares[k]).unwrap_or(latest_fit);

        let mut peak_share = -1.0;
        let mut peak_label = String::new();
        for pt in &timeline {
            if pt.fitted_shares[k] > peak_share {
                peak_share = pt.fitted_shares[k];
                peak_label = pt.period_label.clone();
            }
        }

        let latest_vrt = (0..t_obs).rev().find_map(|t| variant_rt_series[k][t]);

        if k == ref_idx {
            clade_summaries.push(CladeFitnessSummary {
                clade_name: table.clade_names[k].clone(),
                is_reference: true,
                total_sequences: total_seq,
                latest_observed_share: latest_obs,
                latest_fitted_share: latest_fit,
                forecast_end_share: end_fit,
                growth_advantage_per_week: 0.0,
                growth_se_per_week: 0.0,
                growth_ci_low_95: 0.0,
                growth_ci_high_95: 0.0,
                multiplicative_advantage_pct: 0.0,
                crossover_week_vs_ref: None,
                peak_fitted_share: peak_share.max(0.0),
                peak_period_label: peak_label,
                latest_variant_rt: latest_vrt,
            });
        } else {
            let c_pos = non_ref_indices.iter().position(|&idx| idx == k).unwrap();
            let (s_wk, se_wk, crossover) = match cfg.basis_mode {
                LineageBasisMode::LogLinearFitness => {
                    let b0 = beta[c_pos * p_dim];
                    let b1 = beta[c_pos * p_dim + 1];
                    let se1 = cov_matrix[c_pos * p_dim + 1][c_pos * p_dim + 1]
                        .max(0.0)
                        .sqrt();
                    let t_mid = 0.5 * (t_min + t_max);
                    let cross = if b1.abs() > 1e-4 {
                        let t_star = t_mid - b0 / b1;
                        if (-12.0..=(total_steps as f64 + 24.0)).contains(&t_star) {
                            Some(t_star)
                        } else {
                            None
                        }
                    } else {
                        None
                    };
                    (b1, se1, cross)
                }
                LineageBasisMode::NaturalCubicSpline => {
                    let t_hi = t_max;
                    let t_lo = (t_max - 1.0).max(0.0);
                    let x_hi = basis_vector(t_hi, t_min, t_max, cfg.basis_mode);
                    let x_lo = basis_vector(t_lo, t_min, t_max, cfg.basis_mode);
                    let dx: Vec<f64> = x_hi.iter().zip(&x_lo).map(|(a, b)| a - b).collect();
                    let slope: f64 = (0..p_dim).map(|p| beta[c_pos * p_dim + p] * dx[p]).sum();
                    let mut var_s = 0.0;
                    for p in 0..p_dim {
                        for q in 0..p_dim {
                            var_s += dx[p]
                                * cov_matrix[c_pos * p_dim + p][c_pos * p_dim + q]
                                * dx[q];
                        }
                    }
                    let mut cross = None;
                    for w in 1..timeline.len() {
                        let prev_diff = timeline[w - 1].fitted_shares[k]
                            - timeline[w - 1].fitted_shares[ref_idx];
                        let cur_diff =
                            timeline[w].fitted_shares[k] - timeline[w].fitted_shares[ref_idx];
                        if prev_diff <= 0.0 && cur_diff > 0.0 {
                            let frac = (-prev_diff) / (cur_diff - prev_diff).max(1e-9);
                            cross = Some((w - 1) as f64 + frac);
                            break;
                        }
                    }
                    (slope, var_s.max(0.0).sqrt(), cross)
                }
            };

            let mult_pct = ((s_wk * gen_weeks).exp() - 1.0) * 100.0;
            clade_summaries.push(CladeFitnessSummary {
                clade_name: table.clade_names[k].clone(),
                is_reference: false,
                total_sequences: total_seq,
                latest_observed_share: latest_obs,
                latest_fitted_share: latest_fit,
                forecast_end_share: end_fit,
                growth_advantage_per_week: s_wk,
                growth_se_per_week: se_wk,
                growth_ci_low_95: s_wk - z95 * se_wk,
                growth_ci_high_95: s_wk + z95 * se_wk,
                multiplicative_advantage_pct: mult_pct,
                crossover_week_vs_ref: crossover,
                peak_fitted_share: peak_share.max(0.0),
                peak_period_label: peak_label,
                latest_variant_rt: latest_vrt,
            });
        }
    }

    Ok(LineageLogisticsResult {
        clade_names: table.clade_names.clone(),
        reference_clade: table.clade_names[ref_idx].clone(),
        basis_mode: cfg.basis_mode,
        log_likelihood: raw_ll,
        aic,
        converged,
        iterations_used: iters_used,
        clade_summaries,
        timeline,
    })
}

/// Parse either a wide clade count table (`week,total_cases,CladeA,CladeB,...`) or a tidy long table
/// (`week,lineage,count[,total_cases]`) from CSV or TSV text.
pub fn parse_lineage_counts_csv(text: &str, title: &str) -> Result<LineageCountTable, EpiError> {
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let header_line = lines
        .next()
        .ok_or_else(|| EpiError::Schema("Empty lineage CSV/TSV file".to_string()))?;
    let delim = if header_line.contains('\t') { '\t' } else { ',' };
    let headers: Vec<String> = header_line
        .split(delim)
        .map(|s| s.trim().trim_matches('"').to_string())
        .collect();
    if headers.len() < 3 {
        return Err(EpiError::Schema(
            "Lineage CSV requires at least 3 columns".to_string(),
        ));
    }

    let h_lower: Vec<String> = headers.iter().map(|h| h.to_ascii_lowercase()).collect();
    let is_tidy_long = h_lower
        .iter()
        .any(|h| h == "lineage" || h == "clade" || h == "variant" || h == "pango_lineage");

    if is_tidy_long {
        let wk_col = h_lower
            .iter()
            .position(|h| h.contains("week") || h.contains("date") || h == "period")
            .unwrap_or(0);
        let lin_col = h_lower
            .iter()
            .position(|h| h == "lineage" || h == "clade" || h == "variant" || h == "pango_lineage")
            .unwrap();
        let cnt_col = h_lower
            .iter()
            .position(|h| h.contains("seq") || h == "count" || h == "n" || h == "specimens")
            .unwrap_or(2);
        let case_col = h_lower
            .iter()
            .position(|h| h == "total_cases" || h == "cases");

        let mut clades_set = BTreeSet::new();
        let mut period_order = Vec::new();
        let mut period_map: BTreeMap<String, (f64, BTreeMap<String, f64>)> = BTreeMap::new();

        for line in lines {
            let cols: Vec<&str> = line.split(delim).map(|s| s.trim().trim_matches('"')).collect();
            if cols.len() <= lin_col.max(cnt_col).max(wk_col) {
                continue;
            }
            let wk = cols[wk_col].to_string();
            let clade = cols[lin_col].to_string();
            let cnt: f64 = cols[cnt_col].parse().unwrap_or(0.0);
            let tot_c: f64 = case_col
                .and_then(|c| cols.get(c))
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0);

            if !period_map.contains_key(&wk) {
                period_order.push(wk.clone());
            }
            clades_set.insert(clade.clone());
            let entry = period_map
                .entry(wk)
                .or_insert_with(|| (tot_c, BTreeMap::new()));
            if tot_c > entry.0 {
                entry.0 = tot_c;
            }
            *entry.1.entry(clade).or_insert(0.0) += cnt.max(0.0);
        }

        let clade_names: Vec<String> = clades_set.into_iter().collect();
        let mut rows = Vec::with_capacity(period_order.len());
        for wk in period_order {
            if let Some((tot_c, cmap)) = period_map.remove(&wk) {
                let counts: Vec<f64> = clade_names
                    .iter()
                    .map(|cn| cmap.get(cn).copied().unwrap_or(0.0))
                    .collect();
                let sum_seq: f64 = counts.iter().sum();
                rows.push(LineageCountRow {
                    period_label: wk,
                    total_cases: if tot_c > 0.0 { tot_c } else { sum_seq },
                    clade_counts: counts,
                });
            }
        }
        Ok(LineageCountTable {
            title: title.to_string(),
            clade_names,
            rows,
        })
    } else {
        let has_total_cases = h_lower[1] == "total_cases" || h_lower[1] == "cases";
        let clade_start = if has_total_cases { 2 } else { 1 };
        let clade_names: Vec<String> = headers[clade_start..].to_vec();
        if clade_names.len() < 2 {
            return Err(EpiError::Schema(
                "Wide lineage CSV needs at least 2 clade columns".to_string(),
            ));
        }
        let mut rows = Vec::new();
        for line in lines {
            let cols: Vec<&str> = line.split(delim).map(|s| s.trim().trim_matches('"')).collect();
            if cols.len() < clade_start + clade_names.len() {
                continue;
            }
            let counts: Vec<f64> = (0..clade_names.len())
                .map(|k| cols[clade_start + k].parse::<f64>().unwrap_or(0.0).max(0.0))
                .collect();
            let sum_seq: f64 = counts.iter().sum();
            let tot_c = if has_total_cases {
                cols[1].parse::<f64>().unwrap_or(sum_seq).max(0.0)
            } else {
                sum_seq
            };
            rows.push(LineageCountRow {
                period_label: cols[0].to_string(),
                total_cases: tot_c,
                clade_counts: counts,
            });
        }
        Ok(LineageCountTable {
            title: title.to_string(),
            clade_names,
            rows,
        })
    }
}

/// Built-in 20-week multi-clade respiratory genomic surveillance benchmark (`JN.1`, `KP.2`, `KP.3.1.1`, `XEC`, `LP.8.1`).
pub fn builtin_respiratory_lineage_table() -> LineageCountTable {
    let raw: &[(&str, f64, [f64; 5])] = &[
        ("2026-W24", 1_420.0, [148.0, 96.0, 42.0, 6.0, 1.0]),
        ("2026-W25", 1_510.0, [132.0, 108.0, 58.0, 9.0, 1.0]),
        ("2026-W26", 1_640.0, [115.0, 118.0, 79.0, 13.0, 2.0]),
        ("2026-W27", 1_790.0, [96.0, 124.0, 105.0, 19.0, 3.0]),
        ("2026-W28", 1_960.0, [79.0, 122.0, 136.0, 27.0, 4.0]),
        ("2026-W29", 2_180.0, [64.0, 114.0, 168.0, 38.0, 6.0]),
        ("2026-W30", 2_390.0, [50.0, 102.0, 198.0, 53.0, 9.0]),
        ("2026-W31", 2_540.0, [39.0, 88.0, 224.0, 74.0, 13.0]),
        ("2026-W32", 2_610.0, [29.0, 74.0, 241.0, 98.0, 19.0]),
        ("2026-W33", 2_520.0, [22.0, 61.0, 246.0, 126.0, 28.0]),
        ("2026-W34", 2_380.0, [16.0, 48.0, 238.0, 156.0, 41.0]),
        ("2026-W35", 2_290.0, [12.0, 38.0, 220.0, 188.0, 58.0]),
        ("2026-W36", 2_340.0, [9.0, 29.0, 196.0, 218.0, 82.0]),
        ("2026-W37", 2_490.0, [7.0, 22.0, 169.0, 244.0, 114.0]),
        ("2026-W38", 2_680.0, [5.0, 17.0, 142.0, 261.0, 154.0]),
        ("2026-W39", 2_890.0, [4.0, 13.0, 116.0, 268.0, 202.0]),
        ("2026-W40", 3_120.0, [3.0, 9.0, 92.0, 262.0, 256.0]),
        ("2026-W41", 3_290.0, [2.0, 7.0, 71.0, 246.0, 312.0]),
        ("2026-W42", 3_410.0, [2.0, 5.0, 54.0, 224.0, 368.0]),
        ("2026-W43", 3_480.0, [1.0, 4.0, 41.0, 198.0, 418.0]),
    ];

    LineageCountTable {
        title: "Multi-Clade Respiratory Surveillance Benchmark (2026-W24..W43)".to_string(),
        clade_names: vec![
            "JN.1".to_string(),
            "KP.2".to_string(),
            "KP.3.1.1".to_string(),
            "XEC".to_string(),
            "LP.8.1".to_string(),
        ],
        rows: raw
            .iter()
            .map(|(wk, cases, c)| LineageCountRow {
                period_label: (*wk).to_string(),
                total_cases: *cases,
                clade_counts: c.to_vec(),
            })
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multinomial_lineage_growth_advantage_and_spline() {
        let tbl = builtin_respiratory_lineage_table();
        let mut cfg = LineageLogisticsConfig::default();
        let res = fit_multinomial_lineage_logistics(&tbl, &cfg).expect("log-linear fit");
        assert!(res.converged);
        assert_eq!(res.timeline.len(), 20 + 6);
        assert!(
            res.clade_summaries[4].growth_advantage_per_week
                > res.clade_summaries[3].growth_advantage_per_week
        );
        assert!(res.clade_summaries[3].growth_advantage_per_week > 0.20);
        for pt in &res.timeline {
            let sum_p: f64 = pt.fitted_shares.iter().sum();
            assert!((sum_p - 1.0).abs() < 1e-6);
        }

        cfg.basis_mode = LineageBasisMode::NaturalCubicSpline;
        let spline_res = fit_multinomial_lineage_logistics(&tbl, &cfg).expect("spline fit");
        assert_eq!(spline_res.timeline.len(), 26);
        assert!(spline_res.log_likelihood > res.log_likelihood - 5.0);
    }
}
