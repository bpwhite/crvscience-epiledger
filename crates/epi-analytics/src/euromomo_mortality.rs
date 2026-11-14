use epi_core::EpiError;
use serde::{Deserialize, Serialize};

/// All-Cause Excess Mortality Surveillance via EuroMOMO / Serfling Harmonic Quasi-Poisson Regression
/// with 2/3-Power Box-Cox Variance-Stabilizing Transformation (`Z > 2` substantial excess threshold),
/// World Mortality Dataset full-year baseline comparison, P-scores, and Years of Potential Life Lost (`YPLL-75`).
///
/// References:
/// - EuroMOMO harmonic Poisson + 2/3-power Box-Cox Z-score excess mortality: https://pmc.ncbi.nlm.nih.gov/articles/PMC12944855/
/// - Karlinsky & Kobak (2021), The World Mortality Dataset baseline modeling: https://pmc.ncbi.nlm.nih.gov/articles/PMC7852240/

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MortalityBaselineMode {
    /// EuroMOMO / Serfling baseline trained on spring/autumn non-epidemic weeks (`W15..W25` & `W36..W45`)
    /// plus iterative exclusion of residual outliers (`Z > 2.0`), preventing prior influenza/heat waves from inflating baseline.
    EuroMomoSpringAutumn,
    /// World Mortality Dataset style baseline fit across all historical weeks (typical seasonal average).
    WorldMortalityFullYear,
}

impl MortalityBaselineMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::EuroMomoSpringAutumn => "EuroMOMO / Serfling Non-Epidemic Baseline",
            Self::WorldMortalityFullYear => "World Mortality Dataset Full-Season Baseline",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WeeklyMortalityRow {
    pub period_label: String,
    pub iso_week_of_year: u32,
    pub observed_deaths: f64,
    pub population: f64,
    /// Mean age at death in this stratum/week (used to compute `YPLL-75` when < 75).
    pub mean_age_years: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WeeklyMortalityDataset {
    pub title: String,
    pub rows: Vec<WeeklyMortalityRow>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ExcessMortalityConfig {
    pub baseline_mode: MortalityBaselineMode,
    /// Number of trailing weeks held out from baseline fitting as the active surveillance window.
    pub evaluation_window_weeks: usize,
    pub seasonal_harmonics: usize,
    /// EuroMOMO substantial excess Z-score threshold (standard = `2.0`).
    pub z_substantial_threshold: f64,
    /// EuroMOMO high/extreme excess Z-score threshold (standard = `4.0`).
    pub z_extreme_threshold: f64,
}

impl Default for ExcessMortalityConfig {
    fn default() -> Self {
        Self {
            baseline_mode: MortalityBaselineMode::EuroMomoSpringAutumn,
            evaluation_window_weeks: 52,
            seasonal_harmonics: 2,
            z_substantial_threshold: 2.0,
            z_extreme_threshold: 4.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExcessMortalityPoint {
    pub step_index: usize,
    pub period_label: String,
    pub iso_week: u32,
    pub used_in_baseline_fit: bool,
    pub observed_deaths: f64,
    pub expected_baseline: f64,
    /// Upper threshold corresponding to `Z = +z_substantial_threshold` (`Z = 2`) inverse Box-Cox transformed.
    pub substantial_upper_2z: f64,
    /// Upper threshold corresponding to `Z = +z_extreme_threshold` (`Z = 4`) inverse Box-Cox transformed.
    pub extreme_upper_4z: f64,
    pub excess_deaths: f64,
    pub p_score_pct: f64,
    /// EuroMOMO 2/3-power Box-Cox variance-stabilized Z-score.
    pub box_cox_z_score: f64,
    pub excess_ypll_75: f64,
    pub is_substantial_excess: bool,
    pub is_extreme_excess: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExcessMortalityResult {
    pub title: String,
    pub baseline_mode: MortalityBaselineMode,
    pub dispersion_phi: f64,
    pub baseline_weeks_used: usize,
    pub evaluation_cumulative_excess: f64,
    pub evaluation_mean_p_score_pct: f64,
    pub evaluation_ypll_75: f64,
    pub substantial_alert_weeks: usize,
    pub extreme_alert_weeks: usize,
    pub peak_z_score: f64,
    pub peak_period_label: String,
    pub timeline: Vec<ExcessMortalityPoint>,
}

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
        let mut best = aug[col][col].abs();
        for row in (col + 1)..n {
            if aug[row][col].abs() > best {
                best = aug[row][col].abs();
                pivot = row;
            }
        }
        if best < 1e-13 {
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
                let f = aug[row][col];
                if f.abs() > 0.0 {
                    for j in col..=(2 * n) {
                        aug[row][j] -= f * aug[col][j];
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

fn build_harmonic_row(step: usize, iso_week: u32, n_total: usize, n_harm: usize) -> Vec<f64> {
    let mut r = Vec::with_capacity(2 + 2 * n_harm);
    r.push(1.0);
    r.push(((step as f64) - 0.5 * (n_total as f64)) / (n_total.max(1) as f64));
    for h in 1..=n_harm {
        let ang = 2.0 * std::f64::consts::PI * (h as f64) * (iso_week as f64) / 52.1775;
        r.push(ang.sin());
        r.push(ang.cos());
    }
    r
}

fn is_euromomo_clean_season_week(iso_week: u32) -> bool {
    // Spring (W15..=W25) and Autumn (W36..=W45) baseline windows outside winter flu & peak summer heat
    (15..=25).contains(&iso_week) || (36..=45).contains(&iso_week)
}

/// Fit harmonic Quasi-Poisson excess mortality model with 2/3-power Box-Cox variance-stabilized Z-scores.
pub fn compute_euromomo_excess_mortality(
    dataset: &WeeklyMortalityDataset,
    cfg: &ExcessMortalityConfig,
) -> Result<ExcessMortalityResult, EpiError> {
    let n = dataset.rows.len();
    if n < 24 {
        return Err(EpiError::Schema(
            "Excess mortality modeling requires at least 24 weekly rows".to_string(),
        ));
    }

    let eval_win = cfg.evaluation_window_weeks.min(n / 2);
    let train_cutoff = n - eval_win;
    let n_harm = cfg.seasonal_harmonics.clamp(1, 3);
    let p_dim = 2 + 2 * n_harm;

    let design_rows: Vec<Vec<f64>> = (0..n)
        .map(|t| build_harmonic_row(t, dataset.rows[t].iso_week_of_year, n, n_harm))
        .collect();
    let log_offsets: Vec<f64> = dataset
        .rows
        .iter()
        .map(|r| (r.population.max(1.0) / 100_000.0).ln())
        .collect();

    let mut mask: Vec<bool> = (0..n)
        .map(|t| {
            if t >= train_cutoff {
                false
            } else {
                match cfg.baseline_mode {
                    MortalityBaselineMode::EuroMomoSpringAutumn => {
                        is_euromomo_clean_season_week(dataset.rows[t].iso_week_of_year)
                    }
                    MortalityBaselineMode::WorldMortalityFullYear => true,
                }
            }
        })
        .collect();

    // Fallback if series only contains non-spring/autumn weeks
    if mask.iter().filter(|&&b| b).count() < p_dim + 4 {
        for slot in mask.iter_mut().take(train_cutoff) {
            *slot = true;
        }
    }

    let mean_rate = (0..train_cutoff)
        .map(|t| {
            dataset.rows[t].observed_deaths.max(0.5)
                / (dataset.rows[t].population.max(1.0) / 100_000.0)
        })
        .sum::<f64>()
        / (train_cutoff.max(1) as f64);

    let mut beta = vec![0.0_f64; p_dim];
    beta[0] = mean_rate.max(1e-3).ln();
    let mut xtwx_inv = vec![vec![0.0; p_dim]; p_dim];
    let mut phi = 1.0_f64;

    // Two-pass IRLS in EuroMOMO mode (second pass downweights/excludes training weeks with residual Z > 2.5)
    let outer_passes = match cfg.baseline_mode {
        MortalityBaselineMode::EuroMomoSpringAutumn => 2,
        MortalityBaselineMode::WorldMortalityFullYear => 1,
    };

    for pass in 0..outer_passes {
        for _iter in 0..40 {
            let mut grad = vec![0.0_f64; p_dim];
            let mut info = vec![vec![0.0_f64; p_dim]; p_dim];

            for t in 0..n {
                if !mask[t] {
                    continue;
                }
                let xr = &design_rows[t];
                let eta: f64 = log_offsets[t]
                    + xr.iter()
                        .zip(&beta)
                        .map(|(a, b)| a * b)
                        .sum::<f64>()
                        .clamp(-15.0, 15.0);
                let mu = eta.exp();
                let y = dataset.rows[t].observed_deaths.max(0.0);
                let resid = y - mu;
                for a in 0..p_dim {
                    grad[a] += resid * xr[a];
                    for b in a..p_dim {
                        info[a][b] += mu * xr[a] * xr[b];
                    }
                }
            }
            for a in 0..p_dim {
                for b in 0..a {
                    info[a][b] = info[b][a];
                }
                let pen = if a == 0 { 1e-8 } else { 1e-4 };
                grad[a] -= pen * beta[a];
                info[a][a] += pen;
            }

            let Some((step, inv)) = solve_spd_and_invert(&info, &grad) else {
                break;
            };
            xtwx_inv = inv;
            let max_s = step.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            let damp = if max_s > 1.5 { 1.5 / max_s } else { 1.0 };
            for a in 0..p_dim {
                beta[a] += damp * step[a];
            }
            if max_s * damp < 1e-6 {
                break;
            }
        }

        let mut chi2 = 0.0;
        let mut n_used = 0usize;
        for t in 0..n {
            if !mask[t] {
                continue;
            }
            let eta = log_offsets[t]
                + design_rows[t]
                    .iter()
                    .zip(&beta)
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
            let mu = eta.exp().max(1e-6);
            let y = dataset.rows[t].observed_deaths.max(0.0);
            chi2 += (y - mu).powi(2) / mu;
            n_used += 1;
        }
        let df = n_used.saturating_sub(p_dim).max(1) as f64;
        phi = (chi2 / df).max(1.0);

        if pass == 0 && outer_passes > 1 {
            // Prune unseasonal spikes (Z > 2.5) from training set if sufficient baseline weeks remain
            let mut candidate_mask = mask.clone();
            for t in 0..train_cutoff {
                if candidate_mask[t] {
                    let eta = log_offsets[t]
                        + design_rows[t]
                            .iter()
                            .zip(&beta)
                            .map(|(a, b)| a * b)
                            .sum::<f64>();
                    let mu = eta.exp().max(1e-6);
                    let y = dataset.rows[t].observed_deaths.max(0.0);
                    let sd_bc = (2.0 / 3.0) * mu.powf(-1.0 / 3.0) * (phi * mu).sqrt();
                    let z = (y.powf(2.0 / 3.0) - mu.powf(2.0 / 3.0)) / sd_bc.max(1e-9);
                    if z > 2.5 {
                        candidate_mask[t] = false;
                    }
                }
            }
            if candidate_mask.iter().filter(|&&b| b).count() >= p_dim + 6 {
                mask = candidate_mask;
            }
        }
    }

    let baseline_weeks_used = mask.iter().filter(|&&b| b).count();
    let mut timeline = Vec::with_capacity(n);
    let mut eval_cum_excess = 0.0;
    let mut eval_p_sum = 0.0;
    let mut eval_count = 0usize;
    let mut eval_ypll = 0.0;
    let mut substantial_weeks = 0usize;
    let mut extreme_weeks = 0usize;
    let mut peak_z = f64::NEG_INFINITY;
    let mut peak_label = String::new();

    for t in 0..n {
        let xr = &design_rows[t];
        let eta = log_offsets[t]
            + xr.iter()
                .zip(&beta)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                .clamp(-15.0, 15.0);
        let mu = eta.exp().max(1e-6);

        // Model prediction variance on log scale: var(eta_t) = phi * x_t^T (X^T W X)^-1 x_t
        let mut var_eta = 0.0;
        for a in 0..p_dim {
            for b in 0..p_dim {
                var_eta += xr[a] * (phi * xtwx_inv[a][b]) * xr[b];
            }
        }
        // Total prediction variance on count scale: Var(O_t - mu_hat_t) = phi * mu + mu^2 * var_eta
        let total_var_y = (phi * mu + mu * mu * var_eta.max(0.0)).max(1e-9);

        // EuroMOMO 2/3-power Box-Cox transformation: g(y) = y^(2/3), g'(mu) = (2/3) * mu^(-1/3)
        let obs = dataset.rows[t].observed_deaths.max(0.0);
        let mu_bc = mu.powf(2.0 / 3.0);
        let obs_bc = obs.powf(2.0 / 3.0);
        let sd_bc = ((2.0 / 3.0) * mu.powf(-1.0 / 3.0) * total_var_y.sqrt()).max(1e-9);
        let z_bc = (obs_bc - mu_bc) / sd_bc;

        let upper_2z = (mu_bc + cfg.z_substantial_threshold * sd_bc)
            .max(0.0)
            .powf(1.5);
        let upper_4z = (mu_bc + cfg.z_extreme_threshold * sd_bc)
            .max(0.0)
            .powf(1.5);

        let excess = obs - mu;
        let p_score = (excess / mu) * 100.0;
        let years_lost_per_death = (75.0 - dataset.rows[t].mean_age_years).max(0.0);
        let ypll_75 = excess.max(0.0) * years_lost_per_death;

        let is_sub = z_bc > cfg.z_substantial_threshold;
        let is_ext = z_bc > cfg.z_extreme_threshold;

        if t >= train_cutoff {
            eval_cum_excess += excess;
            eval_p_sum += p_score;
            eval_count += 1;
            eval_ypll += ypll_75;
            if is_sub {
                substantial_weeks += 1;
            }
            if is_ext {
                extreme_weeks += 1;
            }
        }

        if z_bc > peak_z {
            peak_z = z_bc;
            peak_label = dataset.rows[t].period_label.clone();
        }

        timeline.push(ExcessMortalityPoint {
            step_index: t,
            period_label: dataset.rows[t].period_label.clone(),
            iso_week: dataset.rows[t].iso_week_of_year,
            used_in_baseline_fit: mask[t],
            observed_deaths: obs,
            expected_baseline: mu,
            substantial_upper_2z: upper_2z,
            extreme_upper_4z: upper_4z,
            excess_deaths: excess,
            p_score_pct: p_score,
            box_cox_z_score: z_bc,
            excess_ypll_75: ypll_75,
            is_substantial_excess: is_sub,
            is_extreme_excess: is_ext,
        });
    }

    Ok(ExcessMortalityResult {
        title: dataset.title.clone(),
        baseline_mode: cfg.baseline_mode,
        dispersion_phi: phi,
        baseline_weeks_used,
        evaluation_cumulative_excess: eval_cum_excess,
        evaluation_mean_p_score_pct: if eval_count > 0 {
            eval_p_sum / (eval_count as f64)
        } else {
            0.0
        },
        evaluation_ypll_75: eval_ypll,
        substantial_alert_weeks: substantial_weeks,
        extreme_alert_weeks: extreme_weeks,
        peak_z_score: peak_z,
        peak_period_label: peak_label,
        timeline,
    })
}

/// Parse weekly all-cause mortality CSV/TSV (`period,iso_week,deaths[,population,mean_age]`).
pub fn parse_weekly_mortality_csv(
    text: &str,
    title: &str,
) -> Result<WeeklyMortalityDataset, EpiError> {
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let hdr = lines
        .next()
        .ok_or_else(|| EpiError::Schema("Empty mortality CSV file".to_string()))?;
    let delim = if hdr.contains('\t') { '\t' } else { ',' };

    let mut rows = Vec::new();
    for (idx, line) in lines.enumerate() {
        let cols: Vec<&str> = line.split(delim).map(|s| s.trim().trim_matches('"')).collect();
        if cols.len() < 2 {
            continue;
        }
        let (iso_w, deaths, pop, mean_age) = if cols.len() >= 3 {
            let iw = cols[1].parse::<u32>().unwrap_or(((idx % 52) + 1) as u32);
            let d = cols[2].parse::<f64>().unwrap_or(0.0);
            let p = cols.get(3).and_then(|s| s.parse().ok()).unwrap_or(4_273_586.0);
            let ma = cols.get(4).and_then(|s| s.parse().ok()).unwrap_or(71.5);
            (iw.clamp(1, 53), d.max(0.0), p, ma)
        } else {
            let iw = ((idx % 52) + 1) as u32;
            let d = cols[1].parse::<f64>().unwrap_or(0.0);
            (iw, d.max(0.0), 4_273_586.0, 71.5)
        };
        rows.push(WeeklyMortalityRow {
            period_label: cols[0].to_string(),
            iso_week_of_year: iso_w,
            observed_deaths: deaths,
            population: pop,
            mean_age_years: mean_age,
        });
    }
    if rows.len() < 24 {
        return Err(EpiError::Schema(
            "Mortality CSV requires at least 24 weekly rows".to_string(),
        ));
    }
    Ok(WeeklyMortalityDataset {
        title: title.to_string(),
        rows,
    })
}

/// Built-in 208-week (`2023-W01..2026-W52`) regional all-cause weekly mortality surveillance series
/// featuring annual sinusoidal winter cycles, demographic population growth, a severe 2026 winter respiratory surge,
/// and a 2026 late-summer heatwave mortality spike.
pub fn builtin_all_cause_weekly_mortality_dataset() -> WeeklyMortalityDataset {
    let mut rows = Vec::with_capacity(208);
    let mut seed: u64 = 0xE020_4040_2026;
    let mut next_u01 = || -> f64 {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((seed >> 11) as f64) / ((1u64 << 53) as f64)
    };

    for yr_idx in 0..4usize {
        let year = 2023 + yr_idx;
        for wk in 1..=52u32 {
            let t = yr_idx * 52 + (wk as usize - 1);
            let pop = 4_220_000.0 + (t as f64) * 310.0;
            let ang = 2.0 * std::f64::consts::PI * ((wk as f64) - 3.0) / 52.1775;
            // Clean seasonal mortality rate per 100k (~18.2 baseline + winter harmonic)
            let base_rate = 18.2 + 2.35 * ang.cos() + 0.35 * (2.0 * ang).cos();
            let base_deaths = base_rate * (pop / 100_000.0);

            // Inject realistic winter influenza spikes and a 2026 summer heatwave excess
            let winter_wave = match (year, wk) {
                (2024, 1..=5) => 75.0 * (-(((wk as f64) - 3.0) / 1.8).powi(2)).exp(),
                (2025, 2..=7) => 95.0 * (-(((wk as f64) - 4.5) / 2.0).powi(2)).exp(),
                (2026, 1..=8) => 185.0 * (-(((wk as f64) - 4.0) / 2.3).powi(2)).exp(),
                _ => 0.0,
            };
            let heat_spike = if year == 2026 && (31..=34).contains(&wk) {
                118.0 * (-(((wk as f64) - 32.5) / 1.1).powi(2)).exp()
            } else {
                0.0
            };

            let z = (next_u01() + next_u01() - 1.0) * 1.4;
            let obs = (base_deaths + winter_wave + heat_spike + z * base_deaths.sqrt())
                .round()
                .max(10.0);
            let mean_age = if heat_spike > 20.0 { 68.4 } else { 71.8 };

            rows.push(WeeklyMortalityRow {
                period_label: format!("{year}-W{wk:02}"),
                iso_week_of_year: wk,
                observed_deaths: obs,
                population: pop,
                mean_age_years: mean_age,
            });
        }
    }

    WeeklyMortalityDataset {
        title: "Regional All-Cause Weekly Mortality Surveillance (2023-W01..2026-W52)".to_string(),
        rows,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euromomo_vs_world_mortality_excess_detection() {
        let ds = builtin_all_cause_weekly_mortality_dataset();
        let cfg = ExcessMortalityConfig::default();
        let euro = compute_euromomo_excess_mortality(&ds, &cfg).expect("EuroMOMO fit");
        assert_eq!(euro.timeline.len(), 208);
        assert!(euro.substantial_alert_weeks >= 4);
        assert!(euro.peak_z_score > 4.0);
        assert!(euro.evaluation_ypll_75 > 500.0);

        let wmd_cfg = ExcessMortalityConfig {
            baseline_mode: MortalityBaselineMode::WorldMortalityFullYear,
            ..cfg
        };
        let wmd = compute_euromomo_excess_mortality(&ds, &wmd_cfg).expect("WMD fit");
        // EuroMOMO excludes prior winter flu spikes from baseline, yielding higher cumulative excess than full-year average
        assert!(euro.evaluation_cumulative_excess > wmd.evaluation_cumulative_excess);
    }
}
