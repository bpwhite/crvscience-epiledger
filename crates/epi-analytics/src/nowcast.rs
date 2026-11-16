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

use crate::math::gamma_quantile;
use crate::renewal_rt::build_discrete_gamma_si;
use epi_core::TimeSeriesDataset;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Configuration for parametric Gamma-delay Bayesian smoothing nowcasting
/// inspired by NobBS (McGough et al. 2020, PMC7162546).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ParametricNowcastConfig {
    /// Mean reporting delay in weeks (e.g. 1.6 weeks for lab-confirmed NNDSS notifications).
    pub delay_mean_weeks: f64,
    /// Standard deviation of reporting delay in weeks.
    pub delay_sd_weeks: f64,
    /// Maximum delay horizon D (weeks) beyond which counts are treated as 100% accrued.
    pub max_delay_weeks: usize,
    /// Negative-Binomial / Gamma-Poisson overdispersion factor (>= 1.0).
    pub overdispersion_phi: f64,
    /// Number of recent pre-truncation weeks used to anchor the local smooth prior rate.
    pub smooth_prior_window: usize,
}

impl Default for ParametricNowcastConfig {
    fn default() -> Self {
        Self {
            delay_mean_weeks: 1.6,
            delay_sd_weeks: 1.1,
            max_delay_weeks: 6,
            overdispersion_phi: 2.0,
            smooth_prior_window: 4,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct NowcastPoint {
    pub observed_count: f64,
    pub completeness_fraction: f64,
    pub expected_unreported: f64,
    pub nowcast_mean: f64,
    pub nowcast_median: f64,
    pub pi_low_95: f64,
    pub pi_high_95: f64,
    pub is_right_truncated: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NowcastSeries {
    pub points: Vec<NowcastPoint>,
    pub delay_pmf: Vec<f64>,
    pub cumulative_completeness: Vec<f64>,
    pub total_imputed_unreported: f64,
}

/// Empirical `N x (D + 1)` incremental reporting triangle (`n[t][d]` = cases at reference step `t`
/// reported after `d` delay steps; `None` in the unobserved bottom-right right-truncated corner).
/// Implements the empirical chain-ladder baseline framework (Wolffram et al. / PMC13005074).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportingTriangle {
    pub reference_labels: Vec<String>,
    pub max_delay: usize,
    /// Row-major `n_rows x (max_delay + 1)` incremental counts (`None` when `t + d >= n_rows`).
    pub matrix: Vec<Vec<Option<f64>>>,
}

impl ReportingTriangle {
    /// Reconstruct an empirical reporting triangle from chronological bitemporal `TimeSeriesDataset`
    /// vintages stored in `epi-store` (ordered oldest-ingested to newest-ingested).
    pub fn from_vintages(vintages: &[TimeSeriesDataset], max_delay: usize) -> Option<Self> {
        if vintages.len() < 2 {
            return None;
        }
        let latest = vintages.last()?;
        if latest.points.len() < 4 {
            return None;
        }

        let reference_labels: Vec<String> =
            latest.points.iter().map(|p| p.date_label.clone()).collect();
        let n_rows = reference_labels.len();
        let d_cap = max_delay.clamp(1, 12);

        // Map each vintage's cumulative reported count by date_label
        let vintage_maps: Vec<BTreeMap<&str, f64>> = vintages
            .iter()
            .map(|ds| {
                ds.points
                    .iter()
                    .map(|p| (p.date_label.as_str(), p.count.max(0.0)))
                    .collect()
            })
            .collect();

        // If we have K ordered vintages representing successive snapshots, map snapshot index to delay
        let mut matrix = vec![vec![None; d_cap + 1]; n_rows];
        for (t, label) in reference_labels.iter().enumerate() {
            let max_obs_d = (n_rows - 1 - t).min(d_cap);
            // Locate earliest vintage containing `label`
            let first_v = vintage_maps
                .iter()
                .position(|m| m.contains_key(label.as_str()));

            if let Some(v0) = first_v {
                let mut prev_cum = 0.0;
                for d in 0..=max_obs_d {
                    let v_idx = (v0 + d).min(vintage_maps.len() - 1);
                    let cum = vintage_maps[v_idx]
                        .get(label.as_str())
                        .copied()
                        .unwrap_or(prev_cum)
                        .max(prev_cum);
                    matrix[t][d] = Some((cum - prev_cum).max(0.0));
                    prev_cum = cum;
                }
            } else {
                let final_c = latest.points[t].count.max(0.0);
                matrix[t][0] = Some(final_c);
                for d in 1..=max_obs_d {
                    matrix[t][d] = Some(0.0);
                }
            }
        }

        Some(Self {
            reference_labels,
            max_delay: d_cap,
            matrix,
        })
    }
}

/// Build discrete delay PMF `p[0..=max_delay]` (allowing non-zero same-week `d = 0` reporting!)
/// and cumulative completeness `F[d] = sum_{u=0..=d} p[u]`.
pub fn build_reporting_delay_distribution(
    mean_weeks: f64,
    sd_weeks: f64,
    max_delay: usize,
) -> (Vec<f64>, Vec<f64>) {
    let d_max = max_delay.max(1);
    // Shift by +0.5 weeks so d=0 covers [0, 1) weeks, d=1 covers [1, 2) weeks, etc.
    let shifted_weights = build_discrete_gamma_si(
        (mean_weeks + 0.5).max(0.35),
        sd_weeks.max(0.20),
        1.0,
        d_max + 1,
    );
    let mut pmf: Vec<f64> = shifted_weights.into_iter().skip(1).collect();
    let sum: f64 = pmf.iter().sum();
    if sum > 0.0 {
        for p in &mut pmf {
            *p /= sum;
        }
    } else {
        pmf[0] = 1.0;
    }

    let mut cdf = Vec::with_capacity(pmf.len());
    let mut acc = 0.0;
    for &p in &pmf {
        acc = (acc + p).min(1.0);
        cdf.push(acc.max(1e-4));
    }
    if let Some(last) = cdf.last_mut() {
        *last = 1.0;
    }
    (pmf, cdf)
}

/// Correct right-truncated surveillance counts using a Bayesian Gamma-Poisson smoothed delay model.
///
/// For each reference time `t` within `max_delay_weeks` of the end `T - 1`:
/// - `lag = (T - 1) - t` has cumulative reporting completeness `F = cdf[lag]`.
/// - Combines a local smooth epidemic rate prior `Gamma(a0, b0)` anchored on recent mature weeks
///   with observed partial count `y_t ~ Poisson(lambda_t * F)` -> posterior `lambda_t | y_t ~ Gamma(a0 + y_t, b0 + F)`.
/// - Unreported count `N_unrep = lambda_t * (1 - F)` is inflated by `overdispersion_phi` for 95% predictive bounds.
pub fn nowcast_parametric_delay(counts: &[f64], cfg: &ParametricNowcastConfig) -> NowcastSeries {
    let n = counts.len();
    let d_max = cfg.max_delay_weeks.clamp(1, 16);
    let (pmf, cdf) =
        build_reporting_delay_distribution(cfg.delay_mean_weeks, cfg.delay_sd_weeks, d_max);

    if n == 0 {
        return NowcastSeries {
            points: Vec::new(),
            delay_pmf: pmf,
            cumulative_completeness: cdf,
            total_imputed_unreported: 0.0,
        };
    }

    // Estimate local baseline incidence prior from mature/near-mature steps just before the truncated tail
    let mature_end = n.saturating_sub(d_max);
    let prior_win = cfg.smooth_prior_window.max(2);
    let prior_start = mature_end.saturating_sub(prior_win);
    let prior_slice = if mature_end > prior_start {
        &counts[prior_start..mature_end]
    } else {
        &counts[..n]
    };
    let mut prior_mean = if prior_slice.is_empty() {
        1.0
    } else {
        (prior_slice.iter().map(|v| v.max(0.0)).sum::<f64>() / (prior_slice.len() as f64)).max(0.5)
    };

    let phi = cfg.overdispersion_phi.max(1.0);
    let mut points = Vec::with_capacity(n);
    let mut total_unreported = 0.0;

    for (t, &raw_y) in counts.iter().enumerate() {
        let y = raw_y.max(0.0);
        let lag = (n - 1) - t;

        if lag >= d_max {
            // Fully accrued historical observation
            prior_mean = 0.7 * prior_mean + 0.3 * y.max(0.5);
            points.push(NowcastPoint {
                observed_count: y,
                completeness_fraction: 1.0,
                expected_unreported: 0.0,
                nowcast_mean: y,
                nowcast_median: y,
                pi_low_95: y,
                pi_high_95: y,
                is_right_truncated: false,
            });
        } else {
            let f_comp = cdf[lag].clamp(0.02, 1.0);
            // Weakly informative time-correlated Gamma prior with prior weight = 1.0 / phi
            let b0 = (1.0 / phi).max(0.05);
            let a0 = (prior_mean * b0).max(0.10);

            // Conjugate Gamma posterior for eventual total rate lambda_t
            let a_post = a0 + y;
            let b_post = b0 + f_comp;
            let lambda_mean = a_post / b_post;

            let expected_unrep = (lambda_mean * (1.0 - f_comp)).max(0.0);
            let nowcast_mean = y + expected_unrep;

            // Predictive distribution for unreported tail incorporating overdispersion phi
            let (pi_low, pi_high, median) = if expected_unrep > 0.05 {
                let var_unrep = expected_unrep * phi + (expected_unrep * expected_unrep) / a_post;
                let k_u = ((expected_unrep * expected_unrep) / var_unrep).max(0.25);
                let r_u = (expected_unrep / var_unrep).max(1e-4);
                let u_low = gamma_quantile(0.025, k_u, r_u);
                let u_med = gamma_quantile(0.50, k_u, r_u);
                let u_high = gamma_quantile(0.975, k_u, r_u);
                (y + u_low, y + u_high, y + u_med)
            } else {
                (y, y, y)
            };

            // Update rolling prior smoothly with posterior estimate (NobBS random-walk coupling)
            prior_mean = 0.65 * prior_mean + 0.35 * nowcast_mean.max(0.5);
            total_unreported += expected_unrep;

            points.push(NowcastPoint {
                observed_count: y,
                completeness_fraction: f_comp,
                expected_unreported: expected_unrep,
                nowcast_mean,
                nowcast_median: median,
                pi_low_95: pi_low,
                pi_high_95: pi_high,
                is_right_truncated: f_comp < 0.999,
            });
        }
    }

    NowcastSeries {
        points,
        delay_pmf: pmf,
        cumulative_completeness: cdf,
        total_imputed_unreported: total_unreported,
    }
}

/// Empirical chain-ladder reporting-triangle nowcaster (Wolffram et al. / PMC13005074).
/// Estimates step-to-step development factors from observed upper-left cells and imputes the missing bottom-right corner.
pub fn nowcast_from_reporting_triangle(
    triangle: &ReportingTriangle,
    overdispersion_phi: f64,
) -> NowcastSeries {
    let n = triangle.reference_labels.len();
    let d_max = triangle.max_delay.max(1);
    if n == 0 {
        return NowcastSeries {
            points: Vec::new(),
            delay_pmf: vec![1.0],
            cumulative_completeness: vec![1.0],
            total_imputed_unreported: 0.0,
        };
    }

    // Impute full matrix left-to-right via empirical cumulative link ratios
    let mut full = vec![vec![0.0; d_max + 1]; n];
    for t in 0..n {
        for d in 0..=d_max {
            if let Some(v) = triangle.matrix[t][d] {
                full[t][d] = v.max(0.0);
            }
        }
    }

    for d in 1..=d_max {
        let mut sum_prev_cum = 0.0;
        let mut sum_curr_inc = 0.0;
        for t in 0..n {
            if triangle.matrix[t][d].is_some() {
                let prev_c: f64 = full[t][0..d].iter().sum();
                sum_prev_cum += prev_c;
                sum_curr_inc += full[t][d];
            }
        }
        let hazard = if sum_prev_cum > 1e-6 {
            (sum_curr_inc / sum_prev_cum).clamp(0.0, 5.0)
        } else {
            0.0
        };

        for t in 0..n {
            if triangle.matrix[t][d].is_none() {
                let prev_c: f64 = full[t][0..d].iter().sum();
                full[t][d] = prev_c * hazard;
            }
        }
    }

    // Derive empirical delay PMF and CDF across all rows
    let mut col_sums = vec![0.0; d_max + 1];
    for row in &full {
        for (d, &val) in row.iter().enumerate() {
            col_sums[d] += val;
        }
    }
    let grand_total: f64 = col_sums.iter().sum();
    let pmf: Vec<f64> = if grand_total > 1e-6 {
        col_sums.iter().map(|&c| c / grand_total).collect()
    } else {
        let mut p = vec![0.0; d_max + 1];
        p[0] = 1.0;
        p
    };

    let mut cdf = Vec::with_capacity(d_max + 1);
    let mut acc = 0.0;
    for &p in &pmf {
        acc = (acc + p).min(1.0);
        cdf.push(acc.max(1e-4));
    }
    if let Some(last) = cdf.last_mut() {
        *last = 1.0;
    }

    let phi = overdispersion_phi.max(1.0);
    let mut points = Vec::with_capacity(n);
    let mut total_unreported = 0.0;

    for t in 0..n {
        let mut obs = 0.0;
        let mut unrep = 0.0;
        let mut max_obs_d = 0usize;
        for d in 0..=d_max {
            if let Some(v) = triangle.matrix[t][d] {
                obs += v.max(0.0);
                max_obs_d = d;
            } else {
                unrep += full[t][d].max(0.0);
            }
        }
        let f_comp = cdf[max_obs_d.min(d_max)];
        let is_trunc = max_obs_d < d_max;
        let mean_total = obs + unrep;
        total_unreported += unrep;

        let (low, high, med) = if is_trunc && unrep > 0.05 {
            let var = unrep * phi + 0.15 * unrep * unrep;
            let k = ((unrep * unrep) / var).max(0.25);
            let r = (unrep / var).max(1e-4);
            (
                obs + gamma_quantile(0.025, k, r),
                obs + gamma_quantile(0.975, k, r),
                obs + gamma_quantile(0.50, k, r),
            )
        } else {
            (obs, obs, obs)
        };

        points.push(NowcastPoint {
            observed_count: obs,
            completeness_fraction: f_comp,
            expected_unreported: unrep,
            nowcast_mean: mean_total,
            nowcast_median: med,
            pi_low_95: low,
            pi_high_95: high,
            is_right_truncated: is_trunc,
        });
    }

    NowcastSeries {
        points,
        delay_pmf: pmf,
        cumulative_completeness: cdf,
        total_imputed_unreported: total_unreported,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parametric_nowcast_inflates_truncated_tail() {
        let mut series = vec![50.0; 20];
        // Simulate right-truncated drop on the final 3 weeks
        series.extend_from_slice(&[42.0, 31.0, 18.0]);
        let out = nowcast_parametric_delay(&series, &ParametricNowcastConfig::default());
        let last = out.points.last().unwrap();
        assert!(last.is_right_truncated);
        assert!(last.completeness_fraction < 0.70);
        assert!(last.nowcast_mean > last.observed_count + 10.0);
        assert!(last.pi_high_95 > last.nowcast_mean);
    }
}
