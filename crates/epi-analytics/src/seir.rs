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

use crate::nowcast::{nowcast_parametric_delay, ParametricNowcastConfig};
use crate::renewal_rt::{compute_cori_rt, CoriRtConfig};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeirPathogenPreset {
    pub id: &'static str,
    pub label: &'static str,
    pub r0_default: f64,
    pub latent_days: f64,
    pub infectious_days: f64,
    pub waning_days: f64, // 0.0 = lifelong immunity (SEIR); >0.0 = waning immunity (SEIRS)
    pub default_immune_fraction: f64,
    pub reporting_fraction: f64,
}

pub fn seir_pathogen_presets() -> Vec<SeirPathogenPreset> {
    vec![
        SeirPathogenPreset {
            id: "measles",
            label: "Measles (SEIR: Latent 8d, Inf 5d)",
            r0_default: 14.0,
            latent_days: 8.0,
            infectious_days: 5.0,
            waning_days: 0.0,
            default_immune_fraction: 0.91,
            reporting_fraction: 0.75,
        },
        SeirPathogenPreset {
            id: "pertussis",
            label: "Pertussis (SEIRS: Latent 7d, Inf 14d, Waning 5y)",
            r0_default: 12.0,
            latent_days: 7.0,
            infectious_days: 14.0,
            waning_days: 1825.0,
            default_immune_fraction: 0.88,
            reporting_fraction: 0.25,
        },
        SeirPathogenPreset {
            id: "sarscov2",
            label: "SARS-CoV-2 (SEIRS: Latent 3d, Inf 5d, Waning 240d)",
            r0_default: 3.2,
            latent_days: 3.0,
            infectious_days: 5.0,
            waning_days: 240.0,
            default_immune_fraction: 0.60,
            reporting_fraction: 0.15,
        },
        SeirPathogenPreset {
            id: "influenza",
            label: "Seasonal Influenza A (SEIRS: Latent 1.8d, Inf 3.2d)",
            r0_default: 1.65,
            latent_days: 1.8,
            infectious_days: 3.2,
            waning_days: 365.0,
            default_immune_fraction: 0.35,
            reporting_fraction: 0.12,
        },
    ]
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SeirConfig {
    pub population: f64,
    pub initial_immune_fraction: f64,
    pub r0: f64,
    pub latent_days: f64,
    pub infectious_days: f64,
    /// `0.0` = permanent immunity (`SEIR`); `> 0.0` = mean immunity duration in days (`SEIRS`).
    pub waning_days: f64,
    /// Daily vaccination rate moving `S -> R` (`0.0..=0.05`).
    pub daily_vaccination_rate: f64,
    /// Fraction of true incident infections captured by surveillance (`0.01..=1.0`).
    pub reporting_fraction: f64,
    /// Fractional reduction in transmission rate `beta` from interventions (`0.0..=0.90`).
    pub contact_reduction: f64,
    /// Forward projection horizon in MMWR weeks (`1..=104`).
    pub horizon_weeks: usize,
    /// Number of stochastic tau-leaping ensemble trajectories (`10..=2000`, parallelized via `std::thread::scope`).
    pub ensemble_runs: usize,
    /// Whether to calibrate initial `R_eff` and `(E0, I0)` from the Cori `Rt` estimate of the series.
    pub auto_calibrate_r0_from_rt: bool,
}

impl Default for SeirConfig {
    fn default() -> Self {
        Self {
            population: 341_784_857.0,
            initial_immune_fraction: 0.91,
            r0: 14.0,
            latent_days: 8.0,
            infectious_days: 5.0,
            waning_days: 0.0,
            daily_vaccination_rate: 0.0,
            reporting_fraction: 0.75,
            contact_reduction: 0.0,
            horizon_weeks: 16,
            ensemble_runs: 120,
            auto_calibrate_r0_from_rt: true,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SeirWeeklyStep {
    pub week_offset: usize,
    pub susceptible: f64,
    pub exposed: f64,
    pub infectious: f64,
    pub recovered: f64,
    pub true_weekly_infections: f64,
    pub reported_weekly_mean: f64,
    pub reported_pi_low_95: f64,
    pub reported_pi_median: f64,
    pub reported_pi_high_95: f64,
    pub effective_rt: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeirForecastResult {
    pub steps: Vec<SeirWeeklyStep>,
    pub calibrated_r0: f64,
    pub initial_r_eff: f64,
    pub peak_week_offset: usize,
    pub peak_reported_weekly: f64,
    pub cumulative_projected_reported: f64,
    pub cumulative_pi_low_95: f64,
    pub cumulative_pi_high_95: f64,
}

#[derive(Clone, Copy)]
struct StateVec {
    s: f64,
    e: f64,
    i: f64,
    r: f64,
    cum_inf: f64,
}

impl StateVec {
    fn add_scaled(self, deriv: StateVec, dt: f64) -> Self {
        Self {
            s: (self.s + deriv.s * dt).max(0.0),
            e: (self.e + deriv.e * dt).max(0.0),
            i: (self.i + deriv.i * dt).max(0.0),
            r: (self.r + deriv.r * dt).max(0.0),
            cum_inf: (self.cum_inf + deriv.cum_inf * dt).max(0.0),
        }
    }
}

fn seir_derivatives(
    st: StateVec,
    pop: f64,
    beta_eff: f64,
    sigma: f64,
    gamma: f64,
    omega: f64,
    nu: f64,
) -> StateVec {
    let n_safe = pop.max(1.0);
    let force_of_infection = beta_eff * (st.i / n_safe);
    let new_exposed = force_of_infection * st.s;
    let new_infectious = sigma * st.e;
    let new_recovered = gamma * st.i;
    let waned = omega * st.r;
    let vaccinated = nu * st.s;

    StateVec {
        s: -new_exposed - vaccinated + waned,
        e: new_exposed - new_infectious,
        i: new_infectious - new_recovered,
        r: new_recovered + vaccinated - waned,
        cum_inf: new_infectious,
    }
}

fn splitmix64_traj_seed(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    let out = z ^ (z >> 31);
    if out == 0 {
        0xCAFE_BABE_1234_5678
    } else {
        out
    }
}

/// Calibrate initial `(S0, E0, I0, R0)` and optionally `R0` from a weekly incidence series,
/// then integrate forward using 4th-order Runge-Kutta (`dt = 1/14` day) + multi-threaded `std::thread::scope`
/// stochastic tau-leaping ensemble (deterministic per-trajectory SplitMix64 seeding).
pub fn calibrate_and_project_seir(
    weekly_counts: &[f64],
    use_nowcast_tail: bool,
    nowcast_cfg: Option<&ParametricNowcastConfig>,
    cfg: &SeirConfig,
) -> SeirForecastResult {
    let effective_series: Vec<f64> = if use_nowcast_tail && !weekly_counts.is_empty() {
        let nc = nowcast_parametric_delay(
            weekly_counts,
            nowcast_cfg.unwrap_or(&ParametricNowcastConfig::default()),
        );
        nc.points.into_iter().map(|p| p.nowcast_mean).collect()
    } else {
        weekly_counts.iter().map(|&v| v.max(0.0)).collect()
    };

    let pop = cfg.population.max(100.0);
    let rep_frac = cfg.reporting_fraction.clamp(0.01, 1.0);
    let latent_d = cfg.latent_days.max(0.5);
    let inf_d = cfg.infectious_days.max(0.5);
    let sigma = 1.0 / latent_d;
    let gamma = 1.0 / inf_d;
    let omega = if cfg.waning_days > 1.0 {
        1.0 / cfg.waning_days
    } else {
        0.0
    };
    let nu = cfg.daily_vaccination_rate.clamp(0.0, 0.05);

    let si_mean = (latent_d + 0.5 * inf_d).max(1.5);
    let si_sd = (0.35 * si_mean).max(0.8);
    let rt_out = compute_cori_rt(
        &effective_series,
        &CoriRtConfig {
            si_mean_days: si_mean,
            si_sd_days: si_sd,
            tau_weeks: 3,
            prior_a: 1.0,
            prior_b: 5.0,
        },
    );
    let recent_rt = rt_out
        .estimates
        .iter()
        .rev()
        .flatten()
        .next()
        .map(|e| e.mean.clamp(0.2, 12.0))
        .unwrap_or(1.15);

    let last_weekly_obs = effective_series.last().copied().unwrap_or(10.0).max(1.0);
    let true_daily_inc = (last_weekly_obs / rep_frac) / 7.0;

    let i0 = (true_daily_inc * inf_d).clamp(1.0, pop * 0.25);
    let e0 = (true_daily_inc * latent_d * recent_rt.powf(0.5)).clamp(1.0, pop * 0.25);
    let susc_frac = (1.0 - cfg.initial_immune_fraction.clamp(0.0, 0.995)).max(0.002);
    let s0 = ((pop - e0 - i0) * susc_frac).max(1.0);
    let r0_comp = (pop - s0 - e0 - i0).max(0.0);

    let calibrated_r0 = if cfg.auto_calibrate_r0_from_rt {
        (recent_rt / (s0 / pop).max(0.001)).clamp(0.3, 35.0)
    } else {
        cfg.r0.max(0.1)
    };

    let beta_eff = calibrated_r0 * gamma * (1.0 - cfg.contact_reduction.clamp(0.0, 0.95));
    let horizon = cfg.horizon_weeks.clamp(1, 104);

    // 1. Deterministic RK4 integration (sub_steps = 14 per day -> dt = 1/14 day)
    let sub_steps_per_day = 14usize;
    let dt = 1.0 / (sub_steps_per_day as f64);
    let mut st = StateVec {
        s: s0,
        e: e0,
        i: i0,
        r: r0_comp,
        cum_inf: 0.0,
    };

    let mut det_steps = Vec::with_capacity(horizon);
    let mut prev_cum = 0.0;

    for wk in 1..=horizon {
        for _ in 0..(7 * sub_steps_per_day) {
            let k1 = seir_derivatives(st, pop, beta_eff, sigma, gamma, omega, nu);
            let k2 = seir_derivatives(
                st.add_scaled(k1, 0.5 * dt),
                pop,
                beta_eff,
                sigma,
                gamma,
                omega,
                nu,
            );
            let k3 = seir_derivatives(
                st.add_scaled(k2, 0.5 * dt),
                pop,
                beta_eff,
                sigma,
                gamma,
                omega,
                nu,
            );
            let k4 = seir_derivatives(
                st.add_scaled(k3, dt),
                pop,
                beta_eff,
                sigma,
                gamma,
                omega,
                nu,
            );

            st.s = (st.s + (dt / 6.0) * (k1.s + 2.0 * k2.s + 2.0 * k3.s + k4.s)).max(0.0);
            st.e = (st.e + (dt / 6.0) * (k1.e + 2.0 * k2.e + 2.0 * k3.e + k4.e)).max(0.0);
            st.i = (st.i + (dt / 6.0) * (k1.i + 2.0 * k2.i + 2.0 * k3.i + k4.i)).max(0.0);
            st.r = (st.r + (dt / 6.0) * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r)).max(0.0);
            st.cum_inf +=
                (dt / 6.0) * (k1.cum_inf + 2.0 * k2.cum_inf + 2.0 * k3.cum_inf + k4.cum_inf);
        }

        let true_wk = (st.cum_inf - prev_cum).max(0.0);
        prev_cum = st.cum_inf;
        let rep_wk = true_wk * rep_frac;
        let r_eff_wk = (beta_eff / gamma) * (st.s / pop);

        det_steps.push(SeirWeeklyStep {
            week_offset: wk,
            susceptible: st.s,
            exposed: st.e,
            infectious: st.i,
            recovered: st.r,
            true_weekly_infections: true_wk,
            reported_weekly_mean: rep_wk,
            reported_pi_low_95: rep_wk,
            reported_pi_median: rep_wk,
            reported_pi_high_95: rep_wk,
            effective_rt: r_eff_wk,
        });
    }

    // 2. Parallel Stochastic Tau-Leaping Ensemble via std::thread::scope
    let runs = cfg.ensemble_runs.clamp(10, 2_000);
    let base_seed = 0x9E37_79B9_7F4A_7C15_u64 ^ last_weekly_obs.to_bits();

    let num_workers = if runs < 64 {
        1
    } else {
        std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(4)
            .clamp(1, 16)
            .min(runs)
    };
    let chunk_size = runs.div_ceil(num_workers);

    let traj_outputs: Vec<(Vec<f64>, f64)> = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(num_workers);
        for w_id in 0..num_workers {
            let r_start = w_id * chunk_size;
            let r_end = ((w_id + 1) * chunk_size).min(runs);
            if r_start >= r_end {
                continue;
            }
            handles.push(scope.spawn(move || {
                let mut local_trajs = Vec::with_capacity(r_end - r_start);
                let tau = 0.5_f64;
                for run_idx in r_start..r_end {
                    let mut rng =
                        XorShift64::new(splitmix64_traj_seed(base_seed ^ (run_idx as u64)));
                    let beta_noise = (1.0 + 0.12 * rng.standard_normal()).clamp(0.65, 1.45);
                    let b_run = beta_eff * beta_noise;

                    let mut s = s0;
                    let mut e = e0;
                    let mut i = i0;
                    let mut r = r0_comp;
                    let mut run_cum_rep = 0.0;
                    let mut wk_draws = Vec::with_capacity(horizon);

                    for _wk in 0..horizon {
                        let mut wk_inf = 0.0;
                        for _ in 0..14 {
                            let inf_rate = b_run * (i / pop) * s * tau;
                            let n_se = rng.poisson_approx(inf_rate).min(s);
                            let n_ei = rng.poisson_approx(sigma * e * tau).min(e);
                            let n_ir = rng.poisson_approx(gamma * i * tau).min(i);
                            let n_rs = rng.poisson_approx(omega * r * tau).min(r);
                            let n_vac = (nu * s * tau).min(s - n_se);

                            s = (s - n_se - n_vac + n_rs).max(0.0);
                            e = (e + n_se - n_ei).max(0.0);
                            i = (i + n_ei - n_ir).max(0.0);
                            r = (r + n_ir + n_vac - n_rs).max(0.0);
                            wk_inf += n_ei;
                        }
                        let rep_draw = rng.poisson_approx(wk_inf * rep_frac);
                        wk_draws.push(rep_draw);
                        run_cum_rep += rep_draw;
                    }
                    local_trajs.push((wk_draws, run_cum_rep));
                }
                local_trajs
            }));
        }

        let mut all = Vec::with_capacity(runs);
        for h in handles {
            if let Ok(part) = h.join() {
                all.extend(part);
            }
        }
        all
    });

    let mut weekly_samples: Vec<Vec<f64>> = vec![Vec::with_capacity(runs); horizon];
    let mut cum_samples: Vec<f64> = Vec::with_capacity(runs);
    for (wk_draws, c_rep) in traj_outputs {
        for (wk, v) in wk_draws.into_iter().enumerate() {
            weekly_samples[wk].push(v);
        }
        cum_samples.push(c_rep);
    }

    for (wk, step) in det_steps.iter_mut().enumerate() {
        let slice = &mut weekly_samples[wk];
        slice.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        step.reported_pi_low_95 = percentile_sorted(slice, 0.025);
        step.reported_pi_median = percentile_sorted(slice, 0.50);
        step.reported_pi_high_95 = percentile_sorted(slice, 0.975);
    }

    cum_samples.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let cum_low = percentile_sorted(&cum_samples, 0.025);
    let cum_high = percentile_sorted(&cum_samples, 0.975);

    let (peak_week_offset, peak_reported_weekly) =
        det_steps
            .iter()
            .fold((1usize, 0.0_f64), |(best_w, best_v), s| {
                if s.reported_weekly_mean > best_v {
                    (s.week_offset, s.reported_weekly_mean)
                } else {
                    (best_w, best_v)
                }
            });

    let cumulative_projected_reported: f64 =
        det_steps.iter().map(|s| s.reported_weekly_mean).sum();
    let initial_r_eff = (beta_eff / gamma) * (s0 / pop);

    SeirForecastResult {
        steps: det_steps,
        calibrated_r0,
        initial_r_eff,
        peak_week_offset,
        peak_reported_weekly,
        cumulative_projected_reported,
        cumulative_pi_low_95: cum_low,
        cumulative_pi_high_95: cum_high,
    }
}

fn percentile_sorted(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let pos = q.clamp(0.0, 1.0) * ((sorted.len() - 1) as f64);
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let frac = pos - (lo as f64);
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0xCAFE_BABE_1234_5678 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn next_f64_open01(&mut self) -> f64 {
        let u = (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64);
        u.clamp(1e-12, 1.0 - 1e-12)
    }

    fn standard_normal(&mut self) -> f64 {
        let u1 = self.next_f64_open01();
        let u2 = self.next_f64_open01();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    fn poisson_approx(&mut self, lambda: f64) -> f64 {
        if lambda <= 0.0 || !lambda.is_finite() {
            return 0.0;
        }
        if lambda < 30.0 {
            let l = (-lambda).exp();
            let mut k = 0u32;
            let mut p = 1.0;
            loop {
                k += 1;
                p *= self.next_f64_open01();
                if p <= l || k > 500 {
                    break;
                }
            }
            f64::from(k - 1)
        } else {
            (lambda + lambda.sqrt() * self.standard_normal())
                .round()
                .max(0.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seir_parallel_reproducibility_and_conservation() {
        let history = vec![15.0, 22.0, 34.0, 52.0, 78.0, 110.0];
        let base_cfg = SeirConfig {
            population: 1_000_000.0,
            initial_immune_fraction: 0.80,
            r0: 10.0,
            auto_calibrate_r0_from_rt: false,
            contact_reduction: 0.0,
            horizon_weeks: 12,
            ensemble_runs: 240,
            ..Default::default()
        };
        let r1 = calibrate_and_project_seir(&history, true, None, &base_cfg);
        let r2 = calibrate_and_project_seir(&history, true, None, &base_cfg);
        assert!((r1.cumulative_pi_high_95 - r2.cumulative_pi_high_95).abs() < 1e-12);

        for step in &r1.steps {
            let sum_pop = step.susceptible + step.exposed + step.infectious + step.recovered;
            assert!((sum_pop - 1_000_000.0).abs() < 1.0);
        }
    }
}
