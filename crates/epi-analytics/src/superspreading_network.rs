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

use crate::math::{gamma_quantile, ln_gamma, normal_quantile, regularized_gamma_p};
use epi_core::EpiError;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// Observation regime for individual secondary offspring counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ZeroTruncationMode {
    /// Complete line-list including terminal / non-transmitting cases (z_i >= 0).
    StandardAllCases,
    /// Cluster-ascertained observations conditioning on z_i >= 1 (ZTNB).
    ZeroTruncatedClustersOnly,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OffspringFitConfig {
    pub mode: ZeroTruncationMode,
    /// Minimum grid search bound for dispersion parameter k.
    pub min_k: f64,
    /// Maximum grid search bound before treating offspring distribution as Poisson.
    pub max_k: f64,
}

impl Default for OffspringFitConfig {
    fn default() -> Self {
        Self {
            mode: ZeroTruncationMode::StandardAllCases,
            min_k: 0.01,
            max_k: 250.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OffspringHistogramBin {
    pub offspring_count: u32,
    pub observed_cases: usize,
    pub observed_prob: f64,
    pub fitted_nb_prob: f64,
    pub fitted_poisson_prob: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OffspringDispersionResult {
    pub n_cases_evaluated: usize,
    pub total_secondary_cases: u64,
    pub mode: ZeroTruncationMode,
    pub r0_mle: f64,
    pub r0_ci_low_95: f64,
    pub r0_ci_high_95: f64,
    pub k_mle: f64,
    pub k_ci_low_95: f64,
    pub k_ci_high_95: f64,
    pub log_likelihood_nb: f64,
    pub log_likelihood_poisson: f64,
    pub lrt_statistic: f64,
    pub lrt_p_value: f64,
    /// Lloyd-Smith theoretical share of all transmissions caused by the top 20% most infectious cases.
    pub prop_transmission_top_20_theoretical: f64,
    /// Empirical share of transmissions caused by the top 20% of observed index cases.
    pub prop_transmission_top_20_empirical: f64,
    /// Theoretical fraction of most infectious cases responsible for 80% of all transmissions (t_80).
    pub cases_driving_80_pct_transmission: f64,
    /// 99th-percentile homogeneous Poisson superspreading event (SSE) cutoff Z^(99).
    pub sse_99_threshold: u32,
    pub observed_sse_count: usize,
    /// Single-introduction Galton-Watson branching-process extinction probability q_ext in [0, 1].
    pub single_intro_extinction_prob: f64,
    /// Empirical Lorenz curve points [cum_fraction_of_cases_desc, cum_fraction_of_transmissions].
    pub empirical_lorenz_curve: Vec<[f64; 2]>,
    /// Theoretical Gamma-offspring Lorenz curve points [u, 1 - P(k+1, x_u)].
    pub theoretical_lorenz_curve: Vec<[f64; 2]>,
    /// Profile log-likelihood curve [k, log_likelihood - max_log_likelihood].
    pub profile_ll_curve: Vec<[f64; 2]>,
    pub histogram: Vec<OffspringHistogramBin>,
}

/// Directed case-to-case contact-tracing record with shared exposure venue metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactTracingRecord {
    /// Source case ID, or empty / "INDEX" / "PRIMARY" for primary introductions.
    pub infector_id: String,
    /// Newly infected case ID.
    pub infectee_id: String,
    /// Exposure setting / venue identifier where transmission occurred.
    pub venue_id: String,
    /// Venue setting category (e.g., "Choir Rehearsal", "Food Processing", "Household").
    pub venue_category: String,
    pub infector_onset_day: Option<f64>,
    pub infectee_onset_day: Option<f64>,
    /// Traced susceptible contacts exposed at this venue event (optional denominator for SAR).
    pub venue_exposed_contacts: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerialIntervalSummary {
    pub pair_count: usize,
    pub mean_days: f64,
    pub sd_days: f64,
    pub median_days: f64,
    pub iqr_low_days: f64,
    pub iqr_high_days: f64,
    pub prop_non_positive: f64,
    /// Positive-shifted Gamma shape/rate suitable for Cori Rt renewal equations.
    pub positive_gamma_mean_days: f64,
    pub positive_gamma_sd_days: f64,
    pub empirical_deltas: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CaseNetworkMetric {
    pub case_id: String,
    pub is_primary_index: bool,
    pub generation_depth: usize,
    pub secondary_offspring: u32,
    pub primary_venue_id: String,
    pub bipartite_degree: usize,
    pub bipartite_k_core: usize,
    pub bipartite_betweenness: f64,
    pub is_superspreader_99: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenueAmplificationMetric {
    pub venue_id: String,
    pub venue_category: String,
    pub secondary_transmissions: usize,
    pub distinct_cases_linked: usize,
    pub exposed_contacts: Option<u32>,
    pub secondary_attack_rate: Option<f64>,
    pub share_of_outbreak_transmissions: f64,
    pub bipartite_degree: usize,
    pub bipartite_k_core: usize,
    pub bipartite_betweenness: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BipartiteNetworkAnalysis {
    pub total_cases: usize,
    pub total_venues: usize,
    pub total_transmission_links: usize,
    pub primary_introductions: usize,
    pub max_generation_depth: usize,
    pub max_bipartite_k_core: usize,
    pub bipartite_density: f64,
    pub offspring_fit: OffspringDispersionResult,
    pub serial_interval: SerialIntervalSummary,
    pub case_metrics: Vec<CaseNetworkMetric>,
    pub venue_metrics: Vec<VenueAmplificationMetric>,
}

/// Negative binomial PMF P(Z = z | r0, k).
pub fn negbin_pmf(z: u32, r0: f64, k: f64) -> f64 {
    if r0 <= 0.0 {
        return if z == 0 { 1.0 } else { 0.0 };
    }
    if k <= 0.0 || !k.is_finite() {
        return poisson_pmf(z, r0);
    }
    let zf = f64::from(z);
    let ln_p = ln_gamma(zf + k)
        - ln_gamma(k)
        - ln_gamma(zf + 1.0)
        + k * (k / (k + r0)).ln()
        + zf * (r0 / (k + r0)).ln();
    ln_p.exp().clamp(0.0, 1.0)
}

/// Homogeneous Poisson PMF P(Z = z | r0).
pub fn poisson_pmf(z: u32, r0: f64) -> f64 {
    if r0 <= 0.0 {
        return if z == 0 { 1.0 } else { 0.0 };
    }
    let zf = f64::from(z);
    (zf * r0.ln() - r0 - ln_gamma(zf + 1.0)).exp().clamp(0.0, 1.0)
}

fn nb_log_likelihood_all(counts: &[u32], r0: f64, k: f64) -> f64 {
    if r0 <= 0.0 || k <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let ln_p0 = k * (k / (k + r0)).ln();
    let ln_q = (r0 / (k + r0)).ln();
    let mut sum = 0.0;
    for &z in counts {
        let zf = f64::from(z);
        sum += ln_gamma(zf + k) - ln_gamma(k) - ln_gamma(zf + 1.0) + ln_p0 + zf * ln_q;
    }
    sum
}

fn poisson_log_likelihood_all(counts: &[u32], r0: f64) -> f64 {
    if r0 <= 0.0 {
        return if counts.iter().all(|&z| z == 0) {
            0.0
        } else {
            f64::NEG_INFINITY
        };
    }
    let ln_r0 = r0.ln();
    counts
        .iter()
        .map(|&z| {
            let zf = f64::from(z);
            zf * ln_r0 - r0 - ln_gamma(zf + 1.0)
        })
        .sum()
}

/// Zero-truncated negative binomial profile log-likelihood for fixed k (optimizing r0 > 0).
fn ztnb_profile_r0_and_ll(pos_counts: &[u32], k: f64) -> (f64, f64) {
    let mean_pos =
        pos_counts.iter().map(|&z| f64::from(z)).sum::<f64>() / (pos_counts.len().max(1) as f64);
    let eval_r0 = |r0: f64| -> f64 {
        if r0 <= 1e-6 {
            return f64::NEG_INFINITY;
        }
        let p0 = (k / (k + r0)).powf(k);
        let one_minus_p0 = (1.0 - p0).max(1e-15);
        nb_log_likelihood_all(pos_counts, r0, k) - (pos_counts.len() as f64) * one_minus_p0.ln()
    };

    // Golden-section maximize over log(r0) in [1e-3, mean_pos * 2.5]
    let mut a = 1e-3_f64.ln();
    let mut b = (mean_pos * 2.5).max(1.0).ln();
    let gr = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    let mut fc = eval_r0(c.exp());
    let mut fd = eval_r0(d.exp());
    for _ in 0..55 {
        if fc < fd {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = eval_r0(d.exp());
        } else {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = eval_r0(c.exp());
        }
    }
    let r0_best = ((a + b) * 0.5).exp();
    let ll_best = eval_r0(r0_best);
    (r0_best, ll_best)
}

fn zt_poisson_r0_and_ll(pos_counts: &[u32]) -> (f64, f64) {
    let mean_pos =
        pos_counts.iter().map(|&z| f64::from(z)).sum::<f64>() / (pos_counts.len().max(1) as f64);
    // Solve r0 / (1 - exp(-r0)) = mean_pos via Newton-Raphson
    let mut r0 = (mean_pos - 0.5).max(0.05);
    for _ in 0..30 {
        let e = (-r0).exp();
        let denom = (1.0 - e).max(1e-15);
        let f = r0 / denom - mean_pos;
        let df = (denom - r0 * e) / (denom * denom);
        if df.abs() > 1e-12 {
            r0 = (r0 - f / df).max(1e-5);
        }
    }
    let one_minus_p0 = (1.0 - (-r0).exp()).max(1e-15);
    let ll = poisson_log_likelihood_all(pos_counts, r0)
        - (pos_counts.len() as f64) * one_minus_p0.ln();
    (r0, ll)
}

/// Lloyd-Smith theoretical share of total transmission generated by top fraction `u` (e.g. 0.20) of cases.
pub fn lloyd_smith_top_share(u: f64, r0: f64, k: f64) -> f64 {
    let u_clamped = u.clamp(1e-6, 1.0 - 1e-6);
    if r0 <= 0.0 || k <= 0.0 || !k.is_finite() {
        return u_clamped;
    }
    let rate = k / r0;
    let nu_thresh = gamma_quantile(1.0 - u_clamped, k, rate);
    (1.0 - regularized_gamma_p(k + 1.0, rate * nu_thresh)).clamp(0.0, 1.0)
}

/// Invert Lloyd-Smith curve to find fraction of cases responsible for target share (e.g., 0.80) of transmission.
pub fn lloyd_smith_cases_for_transmission_share(target_tx_share: f64, r0: f64, k: f64) -> f64 {
    let target = target_tx_share.clamp(1e-6, 1.0 - 1e-6);
    if r0 <= 0.0 || k <= 0.0 || !k.is_finite() {
        return target;
    }
    let mut lo = 1e-6;
    let mut hi = 1.0 - 1e-6;
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        if lloyd_smith_top_share(mid, r0, k) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Single-introduction Galton-Watson branching-process extinction probability solving q = (1 + (R0/k)(1 - q))^(-k).
pub fn branching_extinction_probability(r0: f64, k: f64) -> f64 {
    if r0 <= 1.0 {
        return 1.0;
    }
    if !k.is_finite() || k > 1e5 {
        // Homogeneous Poisson limit: q = exp(R0 * (q - 1))
        let mut q = 0.2;
        for _ in 0..80 {
            q = (r0 * (q - 1.0)).exp();
        }
        return q.clamp(0.0, 1.0);
    }
    let mut q = 0.5_f64;
    for _ in 0..120 {
        let next = (1.0 + (r0 / k) * (1.0 - q)).powf(-k);
        if (next - q).abs() < 1e-12 {
            q = next;
            break;
        }
        q = next;
    }
    q.clamp(0.0, 1.0)
}

/// 99th-percentile homogeneous Poisson superspreading cutoff Z^(99).
pub fn poisson_99_sse_cutoff(r0: f64) -> u32 {
    if r0 <= 0.0 {
        return 1;
    }
    let mut cdf = 0.0;
    for z in 0..500_u32 {
        cdf += poisson_pmf(z, r0);
        if cdf >= 0.99 {
            return z.max(1);
        }
    }
    500
}

/// Fit negative binomial offspring mean R0 and dispersion parameter k via exact MLE + profile-likelihood 95% CIs.
pub fn fit_negative_binomial_offspring(
    raw_counts: &[u32],
    cfg: &OffspringFitConfig,
) -> Result<OffspringDispersionResult, EpiError> {
    let counts: Vec<u32> = match cfg.mode {
        ZeroTruncationMode::StandardAllCases => raw_counts.to_vec(),
        ZeroTruncationMode::ZeroTruncatedClustersOnly => {
            raw_counts.iter().copied().filter(|&z| z >= 1).collect()
        }
    };

    if counts.len() < 3 {
        return Err(EpiError::Schema(
            "At least 3 index case offspring observations are required to estimate R0 and k"
                .to_string(),
        ));
    }

    let n = counts.len();
    let total_sec: u64 = counts.iter().map(|&z| u64::from(z)).sum();
    if total_sec == 0 {
        return Err(EpiError::Schema(
            "All offspring counts are zero (R0 = 0); dispersion parameter k is undefined"
                .to_string(),
        ));
    }

    let min_k = cfg.min_k.max(1e-3);
    let max_k = cfg.max_k.max(10.0);

    let (r0_poisson, ll_poisson) = match cfg.mode {
        ZeroTruncationMode::StandardAllCases => {
            let r0_m = (total_sec as f64) / (n as f64);
            (r0_m, poisson_log_likelihood_all(&counts, r0_m))
        }
        ZeroTruncationMode::ZeroTruncatedClustersOnly => zt_poisson_r0_and_ll(&counts),
    };

    let eval_profile_k = |k: f64| -> (f64, f64) {
        match cfg.mode {
            ZeroTruncationMode::StandardAllCases => {
                let r0_m = (total_sec as f64) / (n as f64);
                (r0_m, nb_log_likelihood_all(&counts, r0_m, k))
            }
            ZeroTruncationMode::ZeroTruncatedClustersOnly => ztnb_profile_r0_and_ll(&counts, k),
        }
    };

    // 1. Coarse log-spaced scan across [min_k, max_k] to bracket global maximum
    let grid_n = 120_usize;
    let ln_min = min_k.ln();
    let ln_max = max_k.ln();
    let mut best_idx = 0_usize;
    let mut grid_pts: Vec<(f64, f64, f64)> = Vec::with_capacity(grid_n);

    for i in 0..grid_n {
        let frac = (i as f64) / ((grid_n - 1) as f64);
        let k = (ln_min + frac * (ln_max - ln_min)).exp();
        let (r0_k, ll_k) = eval_profile_k(k);
        if i == 0 || ll_k > grid_pts[best_idx].2 {
            best_idx = i;
        }
        grid_pts.push((k, r0_k, ll_k));
    }

    // 2. Golden-section refinement around best grid bracket
    let lo_idx = best_idx.saturating_sub(2);
    let hi_idx = (best_idx + 2).min(grid_n - 1);
    let mut a = grid_pts[lo_idx].0.ln();
    let mut b = grid_pts[hi_idx].0.ln();
    let gr = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    let mut fc = eval_profile_k(c.exp()).1;
    let mut fd = eval_profile_k(d.exp()).1;
    for _ in 0..65 {
        if fc < fd {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = eval_profile_k(d.exp()).1;
        } else {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = eval_profile_k(c.exp()).1;
        }
    }

    let mut k_mle = ((a + b) * 0.5).exp();
    let (mut r0_mle, mut ll_nb) = eval_profile_k(k_mle);

    // Check if homogeneous Poisson limit exceeds max_k boundary
    if ll_poisson >= ll_nb {
        k_mle = max_k;
        r0_mle = r0_poisson;
        ll_nb = ll_poisson;
    }

    // 3. Likelihood-ratio test against homogeneous Poisson (chi-bar-square 0.5 * P(chi^2_1 > LRT))
    let lrt = (2.0 * (ll_nb - ll_poisson)).max(0.0);
    let p_chi2_1 = 1.0 - regularized_gamma_p(0.5, 0.5 * lrt);
    let lrt_p_value = if lrt <= 1e-9 {
        1.0
    } else {
        (0.5 * p_chi2_1).clamp(0.0, 1.0)
    };

    // 4. Profile-likelihood 95% confidence interval [k_low_95, k_high_95] at threshold ll_nb - 1.920729
    let target_ll = ll_nb - 0.5 * 3.8414588;

    let k_ci_low_95 = {
        let (_, ll_at_min) = eval_profile_k(min_k);
        if ll_at_min >= target_ll {
            min_k
        } else {
            let mut lo = min_k.ln();
            let mut hi = k_mle.ln();
            for _ in 0..55 {
                let mid = 0.5 * (lo + hi);
                if eval_profile_k(mid.exp()).1 < target_ll {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            (0.5 * (lo + hi)).exp()
        }
    };

    let k_ci_high_95 = if ll_poisson >= target_ll {
        f64::INFINITY
    } else {
        let (_, ll_at_max) = eval_profile_k(max_k);
        if ll_at_max >= target_ll {
            f64::INFINITY
        } else {
            let mut lo = k_mle.ln();
            let mut hi = max_k.ln();
            for _ in 0..55 {
                let mid = 0.5 * (lo + hi);
                if eval_profile_k(mid.exp()).1 >= target_ll {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            (0.5 * (lo + hi)).exp()
        }
    };

    // 5. Asymptotic log-scale 95% CI for R0 incorporating overdispersion Var(R0_hat) = (R0 + R0^2 / k) / n
    let z975 = normal_quantile(0.975);
    let var_r0 = (r0_mle + (r0_mle * r0_mle) / k_mle.max(1e-3)) / (n as f64);
    let se_ln_r0 = (var_r0.sqrt() / r0_mle.max(1e-6)).min(4.0);
    let r0_ci_low_95 = (r0_mle.ln() - z975 * se_ln_r0).exp();
    let r0_ci_high_95 = (r0_mle.ln() + z975 * se_ln_r0).exp();

    // 6. Empirical & theoretical Lorenz curves + 80/20 superspreading metrics
    let mut sorted_desc = counts.clone();
    sorted_desc.sort_unstable_by(|a, b| b.cmp(a));
    let mut empirical_lorenz = Vec::with_capacity(n + 1);
    empirical_lorenz.push([0.0, 0.0]);
    let mut running_tx = 0_u64;
    for (idx, &z) in sorted_desc.iter().enumerate() {
        running_tx += u64::from(z);
        empirical_lorenz.push([
            ((idx + 1) as f64) / (n as f64),
            (running_tx as f64) / (total_sec as f64),
        ]);
    }

    // Interpolate empirical Lorenz curve at x = 0.20
    let mut prop_emp_20 = 0.20;
    for w in empirical_lorenz.windows(2) {
        if w[0][0] <= 0.20 && w[1][0] >= 0.20 {
            let span = (w[1][0] - w[0][0]).max(1e-12);
            let t = (0.20 - w[0][0]) / span;
            prop_emp_20 = w[0][1] + t * (w[1][1] - w[0][1]);
            break;
        }
    }

    let prop_theo_20 = lloyd_smith_top_share(0.20, r0_mle, k_mle);
    let t_80 = lloyd_smith_cases_for_transmission_share(0.80, r0_mle, k_mle);

    let mut theoretical_lorenz = Vec::with_capacity(51);
    theoretical_lorenz.push([0.0, 0.0]);
    for step in 1..50 {
        let u = (step as f64) / 50.0;
        theoretical_lorenz.push([u, lloyd_smith_top_share(u, r0_mle, k_mle)]);
    }
    theoretical_lorenz.push([1.0, 1.0]);

    let sse_99 = poisson_99_sse_cutoff(r0_mle);
    let observed_sse_count = counts.iter().filter(|&&z| z > sse_99).count();
    let q_ext = branching_extinction_probability(r0_mle, k_mle);

    let profile_ll_curve: Vec<[f64; 2]> = grid_pts
        .into_iter()
        .map(|(k, _, ll)| [k, (ll - ll_nb).max(-25.0)])
        .collect();

    // 7. Offspring histogram comparing Observed vs Fitted NB vs Homogeneous Poisson
    let max_z = sorted_desc.first().copied().unwrap_or(0).min(40);
    let mut freq_map = vec![0_usize; (max_z as usize) + 1];
    for &z in &counts {
        let bucket = (z as usize).min(max_z as usize);
        freq_map[bucket] += 1;
    }
    let p0_nb = (k_mle / (k_mle + r0_mle)).powf(k_mle);
    let p0_pois = (-r0_mle).exp();
    let mut histogram = Vec::with_capacity((max_z as usize) + 1);
    for z in 0..=max_z {
        let obs_c = freq_map[z as usize];
        let (nb_p, pois_p) = match cfg.mode {
            ZeroTruncationMode::StandardAllCases => {
                (negbin_pmf(z, r0_mle, k_mle), poisson_pmf(z, r0_mle))
            }
            ZeroTruncationMode::ZeroTruncatedClustersOnly => {
                if z == 0 {
                    (0.0, 0.0)
                } else {
                    (
                        negbin_pmf(z, r0_mle, k_mle) / (1.0 - p0_nb).max(1e-12),
                        poisson_pmf(z, r0_mle) / (1.0 - p0_pois).max(1e-12),
                    )
                }
            }
        };
        histogram.push(OffspringHistogramBin {
            offspring_count: z,
            observed_cases: obs_c,
            observed_prob: (obs_c as f64) / (n as f64),
            fitted_nb_prob: nb_p,
            fitted_poisson_prob: pois_p,
        });
    }

    Ok(OffspringDispersionResult {
        n_cases_evaluated: n,
        total_secondary_cases: total_sec,
        mode: cfg.mode,
        r0_mle,
        r0_ci_low_95,
        r0_ci_high_95,
        k_mle,
        k_ci_low_95,
        k_ci_high_95,
        log_likelihood_nb: ll_nb,
        log_likelihood_poisson: ll_poisson,
        lrt_statistic: lrt,
        lrt_p_value,
        prop_transmission_top_20_theoretical: prop_theo_20,
        prop_transmission_top_20_empirical: prop_emp_20,
        cases_driving_80_pct_transmission: t_80,
        sse_99_threshold: sse_99,
        observed_sse_count,
        single_intro_extinction_prob: q_ext,
        empirical_lorenz_curve: empirical_lorenz,
        theoretical_lorenz_curve: theoretical_lorenz,
        profile_ll_curve,
        histogram,
    })
}

/// Exact Matula-Beck linear-time k-core decomposition on an undirected adjacency list.
pub fn compute_undirected_k_cores(adj: &[Vec<usize>]) -> Vec<usize> {
    let n = adj.len();
    if n == 0 {
        return Vec::new();
    }
    let mut deg: Vec<usize> = adj.iter().map(Vec::len).collect();
    let max_d = deg.iter().copied().max().unwrap_or(0);

    let mut bin = vec![0_usize; max_d + 1];
    for &d in &deg {
        bin[d] += 1;
    }
    let mut start = 0_usize;
    for d in 0..=max_d {
        let num = bin[d];
        bin[d] = start;
        start += num;
    }

    let mut pos = vec![0_usize; n];
    let mut vert = vec![0_usize; n];
    for v in 0..n {
        pos[v] = bin[deg[v]];
        vert[pos[v]] = v;
        bin[deg[v]] += 1;
    }
    for d in (1..=max_d).rev() {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;

    for i in 0..n {
        let v = vert[i];
        for &u in &adj[v] {
            if deg[u] > deg[v] {
                let du = deg[u];
                let pu = pos[u];
                let pw = bin[du];
                let w = vert[pw];
                if u != w {
                    pos[u] = pw;
                    vert[pu] = w;
                    pos[w] = pu;
                    vert[pw] = u;
                }
                bin[du] += 1;
                deg[u] -= 1;
            }
        }
    }
    deg
}

/// Exact Brandes O(|V| * |E|) shortest-path betweenness centrality normalized to [0, 1].
pub fn compute_brandes_betweenness(adj: &[Vec<usize>]) -> Vec<f64> {
    let n = adj.len();
    let mut cb = vec![0.0_f64; n];
    if n < 3 {
        return cb;
    }

    let mut stack = Vec::with_capacity(n);
    let mut preds: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut sigma = vec![0.0_f64; n];
    let mut dist = vec![-1_i32; n];
    let mut delta = vec![0.0_f64; n];
    let mut queue = VecDeque::with_capacity(n);

    for s in 0..n {
        stack.clear();
        for p in &mut preds {
            p.clear();
        }
        sigma.fill(0.0);
        dist.fill(-1);
        sigma[s] = 1.0;
        dist[s] = 0;
        queue.push_back(s);

        while let Some(v) = queue.pop_front() {
            stack.push(v);
            let dv = dist[v];
            for &w in &adj[v] {
                if dist[w] < 0 {
                    queue.push_back(w);
                    dist[w] = dv + 1;
                }
                if dist[w] == dv + 1 {
                    sigma[w] += sigma[v];
                    preds[w].push(v);
                }
            }
        }

        delta.fill(0.0);
        while let Some(w) = stack.pop() {
            for &v in &preds[w] {
                if sigma[w] > 0.0 {
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                }
            }
            if w != s {
                cb[w] += delta[w];
            }
        }
    }

    // Undirected normalization: divide by (n - 1)(n - 2) because each pair was counted in both directions
    let norm = ((n - 1) * (n - 2)) as f64;
    if norm > 0.0 {
        for val in &mut cb {
            *val /= norm;
        }
    }
    cb
}

fn is_primary_marker(infector_id: &str) -> bool {
    let s = infector_id.trim();
    s.is_empty()
        || s.eq_ignore_ascii_case("INDEX")
        || s.eq_ignore_ascii_case("PRIMARY")
        || s.eq_ignore_ascii_case("ROOT")
        || s.eq_ignore_ascii_case("NA")
        || s == "-"
}

/// Full epidemiological analysis of a contact-tracing line list and shared exposure venues.
pub fn analyze_contact_tracing_network(
    records: &[ContactTracingRecord],
    offspring_cfg: &OffspringFitConfig,
) -> Result<BipartiteNetworkAnalysis, EpiError> {
    if records.is_empty() {
        return Err(EpiError::Schema(
            "Contact-tracing line list cannot be empty".to_string(),
        ));
    }

    let mut case_set: BTreeSet<String> = BTreeSet::new();
    let mut venue_meta: BTreeMap<String, (String, Option<u32>, usize)> = BTreeMap::new();
    let mut out_offspring: BTreeMap<String, u32> = BTreeMap::new();
    let mut parent_map: BTreeMap<String, String> = BTreeMap::new();
    let mut case_primary_venue: BTreeMap<String, String> = BTreeMap::new();
    let mut serial_deltas: Vec<f64> = Vec::new();
    let mut total_tx_links = 0_usize;

    for rec in records {
        let infectee = rec.infectee_id.trim().to_string();
        if infectee.is_empty() {
            continue;
        }
        case_set.insert(infectee.clone());
        let venue = if rec.venue_id.trim().is_empty() {
            "VEN-UNSPECIFIED".to_string()
        } else {
            rec.venue_id.trim().to_string()
        };
        let cat = if rec.venue_category.trim().is_empty() {
            "Community".to_string()
        } else {
            rec.venue_category.trim().to_string()
        };
        case_primary_venue
            .entry(infectee.clone())
            .or_insert_with(|| venue.clone());

        let v_entry = venue_meta
            .entry(venue.clone())
            .or_insert_with(|| (cat, rec.venue_exposed_contacts, 0));
        if let Some(exp_c) = rec.venue_exposed_contacts {
            v_entry.1 = Some(v_entry.1.map_or(exp_c, |prev| prev.max(exp_c)));
        }

        if !is_primary_marker(&rec.infector_id) {
            let infector = rec.infector_id.trim().to_string();
            case_set.insert(infector.clone());
            *out_offspring.entry(infector.clone()).or_insert(0) += 1;
            parent_map.entry(infectee.clone()).or_insert(infector.clone());
            case_primary_venue
                .entry(infector)
                .or_insert_with(|| venue.clone());
            v_entry.2 += 1;
            total_tx_links += 1;

            if let (Some(t1), Some(t2)) = (rec.infector_onset_day, rec.infectee_onset_day) {
                if t1.is_finite() && t2.is_finite() {
                    serial_deltas.push(t2 - t1);
                }
            }
        }
    }

    let case_ids: Vec<String> = case_set.into_iter().collect();
    let venue_ids: Vec<String> = venue_meta.keys().cloned().collect();
    let n_cases = case_ids.len();
    let n_venues = venue_ids.len();

    let case_idx: BTreeMap<String, usize> = case_ids
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i))
        .collect();
    let venue_idx: BTreeMap<String, usize> = venue_ids
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), n_cases + i))
        .collect();

    // Build undirected bipartite graph V = V_cases (0..n_cases) U V_venues (n_cases..n_cases+n_venues)
    let total_nodes = n_cases + n_venues;
    let mut adj_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); total_nodes];
    for rec in records {
        let infectee = rec.infectee_id.trim();
        if infectee.is_empty() {
            continue;
        }
        let venue = if rec.venue_id.trim().is_empty() {
            "VEN-UNSPECIFIED"
        } else {
            rec.venue_id.trim()
        };
        if let (Some(&c_u), Some(&v_u)) = (case_idx.get(infectee), venue_idx.get(venue)) {
            adj_sets[c_u].insert(v_u);
            adj_sets[v_u].insert(c_u);
        }
        if !is_primary_marker(&rec.infector_id) {
            let infector = rec.infector_id.trim();
            if let (Some(&c_src), Some(&v_u)) = (case_idx.get(infector), venue_idx.get(venue)) {
                adj_sets[c_src].insert(v_u);
                adj_sets[v_u].insert(c_src);
            }
        }
    }

    let adj: Vec<Vec<usize>> = adj_sets
        .iter()
        .map(|s| s.iter().copied().collect())
        .collect();
    let bipartite_edges: usize = adj.iter().map(Vec::len).sum::<usize>() / 2;
    let bipartite_density = if n_cases > 0 && n_venues > 0 {
        (bipartite_edges as f64) / ((n_cases * n_venues) as f64)
    } else {
        0.0
    };

    let k_cores = compute_undirected_k_cores(&adj);
    let max_bipartite_k_core = k_cores.iter().copied().max().unwrap_or(0);
    let betweenness = compute_brandes_betweenness(&adj);

    // Offspring vector over all cases in the outbreak (leaf infectees correctly carry z_i = 0)
    let offspring_vec: Vec<u32> = case_ids
        .iter()
        .map(|cid| out_offspring.get(cid).copied().unwrap_or(0))
        .collect();
    let offspring_fit = fit_negative_binomial_offspring(&offspring_vec, offspring_cfg)?;

    // Compute generation depth in transmission forest (guarding against cycles)
    let mut max_gen_depth = 0_usize;
    let mut primary_intros = 0_usize;
    let mut case_metrics = Vec::with_capacity(n_cases);

    for (i, cid) in case_ids.iter().enumerate() {
        let is_primary = !parent_map.contains_key(cid);
        if is_primary {
            primary_intros += 1;
        }
        let mut depth = 0_usize;
        let mut cursor = cid.as_str();
        let mut seen = BTreeSet::new();
        seen.insert(cursor);
        while let Some(par) = parent_map.get(cursor) {
            if !seen.insert(par.as_str()) {
                break;
            }
            depth += 1;
            cursor = par.as_str();
        }
        if depth > max_gen_depth {
            max_gen_depth = depth;
        }
        let sec = out_offspring.get(cid).copied().unwrap_or(0);
        case_metrics.push(CaseNetworkMetric {
            case_id: cid.clone(),
            is_primary_index: is_primary,
            generation_depth: depth,
            secondary_offspring: sec,
            primary_venue_id: case_primary_venue
                .get(cid)
                .cloned()
                .unwrap_or_else(|| "VEN-UNSPECIFIED".to_string()),
            bipartite_degree: adj[i].len(),
            bipartite_k_core: k_cores[i],
            bipartite_betweenness: betweenness[i],
            is_superspreader_99: sec > offspring_fit.sse_99_threshold,
        });
    }

    case_metrics.sort_by(|a, b| {
        b.secondary_offspring
            .cmp(&a.secondary_offspring)
            .then_with(|| {
                b.bipartite_betweenness
                    .partial_cmp(&a.bipartite_betweenness)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.case_id.cmp(&b.case_id))
    });

    let mut venue_metrics = Vec::with_capacity(n_venues);
    for (j, vid) in venue_ids.iter().enumerate() {
        let node_u = n_cases + j;
        let (cat, exposed_opt, sec_tx) = venue_meta
            .get(vid)
            .cloned()
            .unwrap_or_else(|| ("Community".to_string(), None, 0));
        let sar = exposed_opt.and_then(|exp_n| {
            if exp_n > 0 {
                Some(((sec_tx as f64) / f64::from(exp_n)).clamp(0.0, 1.0))
            } else {
                None
            }
        });
        venue_metrics.push(VenueAmplificationMetric {
            venue_id: vid.clone(),
            venue_category: cat,
            secondary_transmissions: sec_tx,
            distinct_cases_linked: adj[node_u].len(),
            exposed_contacts: exposed_opt,
            secondary_attack_rate: sar,
            share_of_outbreak_transmissions: if total_tx_links > 0 {
                (sec_tx as f64) / (total_tx_links as f64)
            } else {
                0.0
            },
            bipartite_degree: adj[node_u].len(),
            bipartite_k_core: k_cores[node_u],
            bipartite_betweenness: betweenness[node_u],
        });
    }

    venue_metrics.sort_by(|a, b| {
        b.secondary_transmissions
            .cmp(&a.secondary_transmissions)
            .then_with(|| {
                b.bipartite_betweenness
                    .partial_cmp(&a.bipartite_betweenness)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.venue_id.cmp(&b.venue_id))
    });

    let serial_interval = summarize_serial_intervals(&serial_deltas);

    Ok(BipartiteNetworkAnalysis {
        total_cases: n_cases,
        total_venues: n_venues,
        total_transmission_links: total_tx_links,
        primary_introductions: primary_intros,
        max_generation_depth: max_gen_depth,
        max_bipartite_k_core,
        bipartite_density,
        offspring_fit,
        serial_interval,
        case_metrics,
        venue_metrics,
    })
}

fn summarize_serial_intervals(deltas: &[f64]) -> SerialIntervalSummary {
    if deltas.is_empty() {
        return SerialIntervalSummary {
            pair_count: 0,
            mean_days: 4.5,
            sd_days: 2.2,
            median_days: 4.2,
            iqr_low_days: 3.0,
            iqr_high_days: 5.8,
            prop_non_positive: 0.0,
            positive_gamma_mean_days: 4.5,
            positive_gamma_sd_days: 2.2,
            empirical_deltas: Vec::new(),
        };
    }

    let mut sorted = deltas.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len() as f64;
    let mean = sorted.iter().sum::<f64>() / n;
    let var = if sorted.len() > 1 {
        sorted.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)
    } else {
        2.25
    };
    let sd = var.sqrt().max(0.25);

    let quantile_sorted = |q: f64| -> f64 {
        let pos = q.clamp(0.0, 1.0) * ((sorted.len() - 1) as f64);
        let lo = pos.floor() as usize;
        let hi = pos.ceil() as usize;
        if lo == hi {
            sorted[lo]
        } else {
            sorted[lo] + (pos - (lo as f64)) * (sorted[hi] - sorted[lo])
        }
    };

    let non_pos = sorted.iter().filter(|&&d| d <= 0.0).count();
    let pos_vals: Vec<f64> = sorted.iter().copied().filter(|&d| d > 0.25).collect();
    let (pos_mean, pos_sd) = if pos_vals.len() >= 2 {
        let pm = pos_vals.iter().sum::<f64>() / (pos_vals.len() as f64);
        let pv = pos_vals.iter().map(|&x| (x - pm).powi(2)).sum::<f64>()
            / ((pos_vals.len() - 1) as f64);
        (pm, pv.sqrt().max(0.5))
    } else {
        (mean.max(1.5), sd.max(0.8))
    };

    SerialIntervalSummary {
        pair_count: sorted.len(),
        mean_days: mean,
        sd_days: sd,
        median_days: quantile_sorted(0.50),
        iqr_low_days: quantile_sorted(0.25),
        iqr_high_days: quantile_sorted(0.75),
        prop_non_positive: (non_pos as f64) / n,
        positive_gamma_mean_days: pos_mean,
        positive_gamma_sd_days: pos_sd,
        empirical_deltas: sorted,
    }
}

/// Parse contact-tracing CSV (`infector_id,infectee_id,venue_id,venue_category,infector_onset_day,infectee_onset_day,venue_exposed_contacts`).
pub fn parse_contact_tracing_csv(csv_text: &str) -> Result<Vec<ContactTracingRecord>, EpiError> {
    let mut rows = Vec::new();
    for (line_idx, raw_line) in csv_text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        if line_idx == 0
            && (cols[0].eq_ignore_ascii_case("infector_id")
                || cols[0].eq_ignore_ascii_case("source_id"))
        {
            continue;
        }
        if cols.len() < 3 {
            return Err(EpiError::Schema(format!(
                "Contact-tracing CSV line {} requires at least 3 columns (infector_id,infectee_id,venue_id)",
                line_idx + 1
            )));
        }
        let infector_onset_day = cols.get(4).and_then(|s| s.parse::<f64>().ok());
        let infectee_onset_day = cols.get(5).and_then(|s| s.parse::<f64>().ok());
        let venue_exposed_contacts = cols.get(6).and_then(|s| s.parse::<u32>().ok());
        rows.push(ContactTracingRecord {
            infector_id: cols[0].to_string(),
            infectee_id: cols[1].to_string(),
            venue_id: cols[2].to_string(),
            venue_category: cols.get(3).copied().unwrap_or("Community").to_string(),
            infector_onset_day,
            infectee_onset_day,
            venue_exposed_contacts,
        });
    }
    if rows.is_empty() {
        return Err(EpiError::Schema(
            "No valid contact-tracing rows found in CSV".to_string(),
        ));
    }
    Ok(rows)
}

/// Export styled Graphviz `.dot` representation of the bipartite case-venue and directed transmission network.
pub fn export_bipartite_network_dot(analysis: &BipartiteNetworkAnalysis, records: &[ContactTracingRecord]) -> String {
    let mut out = String::from(
        "digraph EpiLedgerSuperspreadingBipartite {\n  rankdir=LR;\n  node [fontname=\"Helvetica\", fontsize=10];\n",
    );
    for vm in &analysis.venue_metrics {
        out.push_str(&format!(
            "  \"{}\" [shape=box, style=\"rounded,filled\", fillcolor=\"#1e293b\", fontcolor=\"#f8fafc\", label=\"{}\\n({} | tx={} | core={})\"];\n",
            vm.venue_id, vm.venue_id, vm.venue_category, vm.secondary_transmissions, vm.bipartite_k_core
        ));
    }
    for cm in &analysis.case_metrics {
        let fill = if cm.is_superspreader_99 {
            "#f43f5e"
        } else if cm.secondary_offspring >= 2 {
            "#f59e0b"
        } else if cm.is_primary_index {
            "#38bdf8"
        } else {
            "#334155"
        };
        out.push_str(&format!(
            "  \"{}\" [shape=ellipse, style=filled, fillcolor=\"{}\", fontcolor=\"#ffffff\", label=\"{} (z={})\"];\n",
            cm.case_id, fill, cm.case_id, cm.secondary_offspring
        ));
    }
    for rec in records {
        if !is_primary_marker(&rec.infector_id) {
            out.push_str(&format!(
                "  \"{}\" -> \"{}\" [color=\"#64748b\", label=\"{}\"];\n",
                rec.infector_id.trim(),
                rec.infectee_id.trim(),
                rec.venue_id.trim()
            ));
        }
    }
    out.push_str("}\n");
    out
}

/// Realistic multi-cluster respiratory superspreading line list (choir rehearsal, cold-pack shift, clinical unit, and households; k ~ 0.25).
pub fn builtin_superspreading_outbreak_records() -> Vec<ContactTracingRecord> {
    let mut recs = Vec::new();

    // Primary index cases
    recs.push(ContactTracingRecord {
        infector_id: "INDEX".to_string(),
        infectee_id: "CASE-001".to_string(),
        venue_id: "VEN-IMPORT-AIRPORT".to_string(),
        venue_category: "Travel Import".to_string(),
        infector_onset_day: None,
        infectee_onset_day: Some(0.0),
        venue_exposed_contacts: Some(18),
    });
    recs.push(ContactTracingRecord {
        infector_id: "INDEX".to_string(),
        infectee_id: "CASE-040".to_string(),
        venue_id: "VEN-IMPORT-CONF".to_string(),
        venue_category: "Regional Seminar".to_string(),
        infector_onset_day: None,
        infectee_onset_day: Some(2.0),
        venue_exposed_contacts: Some(30),
    });

    // Superspreading Event 1: CASE-001 attends Indoor Choir Rehearsal (14 secondary cases out of 24 attendees)
    for k in 2..=15 {
        let onset = 3.0 + ((k * 3) % 7) as f64 * 0.7;
        recs.push(ContactTracingRecord {
            infector_id: "CASE-001".to_string(),
            infectee_id: format!("CASE-{k:03}"),
            venue_id: "VEN-CHOIR-HALL".to_string(),
            venue_category: "Indoor Choir".to_string(),
            infector_onset_day: Some(0.0),
            infectee_onset_day: Some(onset),
            venue_exposed_contacts: Some(24),
        });
    }

    // Bridge case CASE-004 works at Cold-Storage Packing Shift B -> triggers SSE 2 (9 secondary cases out of 36 workers)
    for (idx, cid) in (16..=24).enumerate() {
        let onset = 7.4 + (idx as f64) * 0.55;
        recs.push(ContactTracingRecord {
            infector_id: "CASE-004".to_string(),
            infectee_id: format!("CASE-{cid:03}"),
            venue_id: "VEN-COLDPACK-B".to_string(),
            venue_category: "Shift Workplace".to_string(),
            infector_onset_day: Some(4.4),
            infectee_onset_day: Some(onset),
            venue_exposed_contacts: Some(36),
        });
    }

    // Secondary household & dining transmissions (mostly small chains of 1-3, while most leaf cases produce 0)
    let secondary_pairs: &[(&str, &str, &str, &str, f64, f64, u32)] = &[
        ("CASE-002", "CASE-025", "VEN-HH-NORTH-1", "Household", 3.7, 7.2, 4),
        ("CASE-002", "CASE-026", "VEN-HH-NORTH-1", "Household", 3.7, 8.1, 4),
        ("CASE-007", "CASE-027", "VEN-HH-WEST-3", "Household", 5.1, 8.9, 3),
        ("CASE-016", "CASE-028", "VEN-BISTRO-OAK", "Indoor Dining", 7.4, 11.0, 14),
        ("CASE-016", "CASE-029", "VEN-BISTRO-OAK", "Indoor Dining", 7.4, 11.8, 14),
        ("CASE-016", "CASE-030", "VEN-BISTRO-OAK", "Indoor Dining", 7.4, 12.5, 14),
        // Presymptomatic transmission example (delta_t = -0.4d)
        ("CASE-019", "CASE-031", "VEN-HH-SOUTH-2", "Household", 9.05, 8.65, 3),
        // CASE-040 seeds a 6-case Clinical Unit cluster (SSE 3)
        ("CASE-040", "CASE-041", "VEN-CLINIC-3B", "Clinical Unit", 2.0, 5.5, 22),
        ("CASE-040", "CASE-042", "VEN-CLINIC-3B", "Clinical Unit", 2.0, 6.1, 22),
        ("CASE-040", "CASE-043", "VEN-CLINIC-3B", "Clinical Unit", 2.0, 6.8, 22),
        ("CASE-040", "CASE-044", "VEN-CLINIC-3B", "Clinical Unit", 2.0, 7.4, 22),
        ("CASE-040", "CASE-045", "VEN-CLINIC-3B", "Clinical Unit", 2.0, 8.0, 22),
        ("CASE-040", "CASE-046", "VEN-CLINIC-3B", "Clinical Unit", 2.0, 8.7, 22),
        // Bridge between Clinic 3B and Bistro Oak via CASE-042 dining at VEN-BISTRO-OAK
        ("CASE-042", "CASE-047", "VEN-BISTRO-OAK", "Indoor Dining", 6.1, 10.2, 14),
        ("CASE-044", "CASE-048", "VEN-HH-EAST-7", "Household", 7.4, 11.4, 2),
    ];

    for &(src, dst, vid, cat, t1, t2, exp_n) in secondary_pairs {
        recs.push(ContactTracingRecord {
            infector_id: src.to_string(),
            infectee_id: dst.to_string(),
            venue_id: vid.to_string(),
            venue_category: cat.to_string(),
            infector_onset_day: Some(t1),
            infectee_onset_day: Some(t2),
            venue_exposed_contacts: Some(exp_n),
        });
    }

    recs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_overdispersed_offspring_mle_recovers_low_k_and_80_20_rule() {
        // 48 non-transmitters (0), 8 small transmitters (1..2), and 4 superspreaders (8, 11, 14, 19)
        let mut z = vec![0_u32; 48];
        z.extend([1, 1, 1, 1, 1, 2, 2, 2, 8, 11, 14, 19]);

        let res = fit_negative_binomial_offspring(&z, &OffspringFitConfig::default())
            .expect("valid NB fit");
        assert_eq!(res.n_cases_evaluated, 60);
        assert!((res.r0_mle - 1.05).abs() < 1e-6);
        assert!(
            res.k_mle > 0.05 && res.k_mle < 0.45,
            "Expected strong overdispersion k in (0.05, 0.45), got {}",
            res.k_mle
        );
        assert!(res.k_ci_low_95 <= res.k_mle && res.k_mle <= res.k_ci_high_95);
        assert!(res.lrt_p_value < 1e-4, "Should strongly reject homogeneous Poisson");
        assert!(
            res.prop_transmission_top_20_theoretical > 0.75,
            "Top 20% should drive >75% of transmission under k ~ 0.18"
        );
    }

    #[test]
    fn test_bipartite_contact_network_k_core_and_betweenness() {
        let recs = builtin_superspreading_outbreak_records();
        let net = analyze_contact_tracing_network(&recs, &OffspringFitConfig::default())
            .expect("network analysis succeeds");

        assert!(net.total_cases >= 40);
        assert!(net.total_venues >= 7);
        assert!(net.offspring_fit.k_mle < 0.50);
        assert_eq!(net.case_metrics[0].case_id, "CASE-001");
        assert_eq!(net.case_metrics[0].secondary_offspring, 14);
        assert_eq!(net.venue_metrics[0].venue_id, "VEN-CHOIR-HALL");
        assert!(net.max_bipartite_k_core >= 1);
        assert!(net.serial_interval.mean_days > 2.5 && net.serial_interval.mean_days < 6.5);
    }
}
