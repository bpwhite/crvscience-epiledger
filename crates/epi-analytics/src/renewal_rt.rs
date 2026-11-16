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

use crate::math::{gamma_pdf, gamma_quantile};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CoriRtConfig {
    pub si_mean_days: f64,
    pub si_sd_days: f64,
    pub tau_weeks: usize,
    pub prior_a: f64,
    pub prior_b: f64,
}

impl Default for CoriRtConfig {
    fn default() -> Self {
        Self {
            si_mean_days: 11.7,
            si_sd_days: 2.0,
            tau_weeks: 3,
            prior_a: 1.0,
            prior_b: 5.0,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CoriRtEstimate {
    pub mean: f64,
    pub median: f64,
    pub cri_low_95: f64,
    pub cri_high_95: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoriRtSeries {
    pub estimates: Vec<Option<CoriRtEstimate>>,
    pub discrete_kernel: Vec<(String, f64)>,
    pub used_daily_em_disaggregation: bool,
}

/// Build normalized discrete Gamma serial interval weights `w(1..=max_steps)` (`w(0) = 0`).
pub fn build_discrete_gamma_si(
    mean_days: f64,
    sd_days: f64,
    step_days: f64,
    max_steps: usize,
) -> Vec<f64> {
    let mu = (mean_days / step_days).max(0.20);
    let sigma = (sd_days / step_days).max(0.10);
    let shape = (mu * mu) / (sigma * sigma);
    let rate = mu / (sigma * sigma);

    let mut weights = vec![0.0; max_steps + 1];
    let sub = 16usize;
    let mut sum = 0.0;

    for s in 1..=max_steps {
        let low = ((s as f64) - 0.5).max(1e-4);
        let high = (s as f64) + 0.5;
        let h = (high - low) / (sub as f64);
        let mut area = 0.5 * (gamma_pdf(low, shape, rate) + gamma_pdf(high, shape, rate));
        for k in 1..sub {
            area += gamma_pdf(low + (k as f64) * h, shape, rate);
        }
        let w = area * h;
        weights[s] = w;
        sum += w;
    }

    if sum > 0.0 {
        for w in weights.iter_mut().skip(1) {
            *w /= sum;
        }
    } else if max_steps >= 1 {
        weights[1] = 1.0;
    }
    weights
}

/// Estimate time-varying reproduction number `Rt` via Cori et al. (`EpiEstim`, PMC3816335).
/// When `si_mean_days < 7.0` on weekly series, reconstructs daily incidence via mass-conserving
/// EM-inspired smoothing before evaluating the daily renewal equation (Nash et al., PMC10491397).
pub fn compute_cori_rt(weekly_counts: &[f64], cfg: &CoriRtConfig) -> CoriRtSeries {
    let n = weekly_counts.len();
    if n == 0 {
        return CoriRtSeries {
            estimates: Vec::new(),
            discrete_kernel: Vec::new(),
            used_daily_em_disaggregation: false,
        };
    }

    let tau = cfg.tau_weeks.max(1);
    if cfg.si_mean_days >= 7.0 {
        let max_lag = 10usize;
        let w = build_discrete_gamma_si(cfg.si_mean_days, cfg.si_sd_days, 7.0, max_lag);
        let mut lambda = vec![0.0; n];
        for t in 0..n {
            let mut inf = 0.0;
            for s in 1..=max_lag {
                if t >= s {
                    inf += weekly_counts[t - s].max(0.0) * w[s];
                }
            }
            lambda[t] = inf;
        }

        let mut estimates = vec![None; n];
        for t in 1..n {
            let start = t.saturating_sub(tau - 1).max(1);
            let mut sum_i = 0.0;
            let mut sum_l = 0.0;
            for u in start..=t {
                sum_i += weekly_counts[u].max(0.0);
                sum_l += lambda[u];
            }
            if sum_l > 0.05 {
                let a_post = cfg.prior_a + sum_i;
                let b_post = (1.0 / cfg.prior_b) + sum_l;
                estimates[t] = Some(CoriRtEstimate {
                    mean: a_post / b_post,
                    median: gamma_quantile(0.50, a_post, b_post),
                    cri_low_95: gamma_quantile(0.025, a_post, b_post),
                    cri_high_95: gamma_quantile(0.975, a_post, b_post),
                });
            }
        }

        let discrete_kernel = w
            .iter()
            .enumerate()
            .skip(1)
            .map(|(s, &val)| (format!("Wk +{s}"), val))
            .collect();

        CoriRtSeries {
            estimates,
            discrete_kernel,
            used_daily_em_disaggregation: false,
        }
    } else {
        // Mass-conserving daily reconstruction preserving exact weekly totals
        let daily_len = n * 7;
        let mut daily = vec![0.0; daily_len];
        for (wk, &c) in weekly_counts.iter().enumerate() {
            let per_day = c.max(0.0) / 7.0;
            for d in 0..7 {
                daily[wk * 7 + d] = per_day;
            }
        }

        for _ in 0..3 {
            let mut smoothed = daily.clone();
            for i in 3..daily_len.saturating_sub(3) {
                let win_sum: f64 = daily[i - 3..=i + 3].iter().sum();
                smoothed[i] = win_sum / 7.0;
            }
            for (wk, &c) in weekly_counts.iter().enumerate() {
                let wk_sum: f64 = smoothed[wk * 7..wk * 7 + 7].iter().sum();
                let scale = if wk_sum > 1e-9 { c.max(0.0) / wk_sum } else { 0.0 };
                for d in 0..7 {
                    daily[wk * 7 + d] = smoothed[wk * 7 + d] * scale;
                }
            }
        }

        let max_day_lag = 28usize;
        let w_day = build_discrete_gamma_si(cfg.si_mean_days, cfg.si_sd_days, 1.0, max_day_lag);
        let mut lambda_day = vec![0.0; daily_len];
        for d in 0..daily_len {
            let mut inf = 0.0;
            for s in 1..=max_day_lag {
                if d >= s {
                    inf += daily[d - s] * w_day[s];
                }
            }
            lambda_day[d] = inf;
        }

        let tau_days = tau * 7;
        let mut estimates = vec![None; n];
        for wk in 1..n {
            let end_d = wk * 7 + 6;
            let start_d = end_d.saturating_sub(tau_days - 1).max(1);
            let sum_i: f64 = daily[start_d..=end_d].iter().sum();
            let sum_l: f64 = lambda_day[start_d..=end_d].iter().sum();
            if sum_l > 0.05 {
                let a_post = cfg.prior_a + sum_i;
                let b_post = (1.0 / cfg.prior_b) + sum_l;
                estimates[wk] = Some(CoriRtEstimate {
                    mean: a_post / b_post,
                    median: gamma_quantile(0.50, a_post, b_post),
                    cri_low_95: gamma_quantile(0.025, a_post, b_post),
                    cri_high_95: gamma_quantile(0.975, a_post, b_post),
                });
            }
        }

        let discrete_kernel = w_day
            .iter()
            .enumerate()
            .skip(1)
            .take(14)
            .map(|(s, &val)| (format!("d+{s}"), val))
            .collect();

        CoriRtSeries {
            estimates,
            discrete_kernel,
            used_daily_em_disaggregation: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cori_rt_expanding_vs_declining() {
        let growing = vec![2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0];
        let cfg = CoriRtConfig::default();
        let res = compute_cori_rt(&growing, &cfg);
        let last = res.estimates.last().unwrap().unwrap();
        assert!(last.mean > 1.3);

        let shrinking = vec![128.0, 64.0, 32.0, 16.0, 8.0, 4.0, 2.0];
        let res_down = compute_cori_rt(&shrinking, &cfg);
        let last_down = res_down.estimates.last().unwrap().unwrap();
        assert!(last_down.mean < 0.8);
    }
}
