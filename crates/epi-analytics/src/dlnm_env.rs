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
use epi_core::EpiError;
use serde::{Deserialize, Serialize};

/// Distributed Lag Non-Linear Models (`DLNM`) via bi-dimensional cross-basis tensor products
/// across environmental exposure intensity (heat index, wildfire PM2.5, extreme precipitation)
/// and lag dimensions (`0..=max_lag`), estimated via overdispersed Quasi-Poisson IRLS.
///
/// References:
/// - Gasparrini A, Armstrong B, Kenward MG (2010), Distributed lag non-linear models: https://pmc.ncbi.nlm.nih.gov/articles/PMC2998707/
/// - Gasparrini A (2011), dlnm cross-basis functional space specification: https://pmc.ncbi.nlm.nih.gov/articles/PMC3191524/
/// - Spatial Bayesian small-area DLNM epidemiological extensions: https://pmc.ncbi.nlm.nih.gov/articles/PMC11031409/

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExposureBasisKind {
    /// 3-df natural cubic spline centered at `reference_exposure`.
    NaturalCubicSpline3Df,
    /// Piecewise linear + quadratic threshold terms below and above `reference_exposure`.
    QuadraticThreshold,
}

impl ExposureBasisKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::NaturalCubicSpline3Df => "Natural Cubic Spline (df=3)",
            Self::QuadraticThreshold => "Piecewise Quadratic Threshold",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LagBasisKind {
    /// 3-df natural cubic spline over `log(lag + 1)` concentrating knots at acute short lags.
    NaturalCubicLogLag3Df,
    /// Quadratic polynomial basis `(1, l/L, (l/L)^2)`.
    QuadraticPolynomial,
}

impl LagBasisKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::NaturalCubicLogLag3Df => "Log-Knots Cubic Spline (df=3)",
            Self::QuadraticPolynomial => "Quadratic Lag Polynomial (df=3)",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DlnmObservation {
    pub date_label: String,
    pub exposure: f64,
    pub outcome_count: f64,
    pub covariate_temp_or_confounder: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DlnmDataset {
    pub title: String,
    pub exposure_unit: String,
    pub outcome_label: String,
    pub observations: Vec<DlnmObservation>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DlnmConfig {
    pub exposure_basis: ExposureBasisKind,
    pub lag_basis: LagBasisKind,
    pub max_lag: usize,
    /// Centering exposure `x_ref` where `RR = 1.0` (if `None`, uses median exposure).
    pub reference_exposure: Option<f64>,
    /// Seasonal cycle period in steps (`365.25` for daily series, `52.18` for weekly series).
    pub seasonal_period_steps: f64,
    pub seasonal_harmonics: usize,
    pub ridge_l2_penalty: f64,
}

impl Default for DlnmConfig {
    fn default() -> Self {
        Self {
            exposure_basis: ExposureBasisKind::NaturalCubicSpline3Df,
            lag_basis: LagBasisKind::NaturalCubicLogLag3Df,
            max_lag: 14,
            reference_exposure: None,
            seasonal_period_steps: 365.25,
            seasonal_harmonics: 2,
            ridge_l2_penalty: 1e-4,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CumulativeExposurePoint {
    pub exposure: f64,
    pub cumulative_rr: f64,
    pub rr_ci_low_95: f64,
    pub rr_ci_high_95: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LagSlicePoint {
    pub lag: usize,
    pub rr: f64,
    pub rr_ci_low_95: f64,
    pub rr_ci_high_95: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DlnmFitResult {
    pub title: String,
    pub exposure_unit: String,
    pub effective_rows_fitted: usize,
    pub max_lag: usize,
    pub reference_exposure: f64,
    pub minimum_risk_exposure: f64,
    pub p95_exposure: f64,
    pub p95_cumulative_rr: f64,
    pub p95_rr_ci_low_95: f64,
    pub p95_rr_ci_high_95: f64,
    pub dispersion_phi: f64,
    pub deviance: f64,
    /// Gasparrini-Leone forward attributable events across observed series.
    pub attributable_events: f64,
    pub attributable_fraction_pct: f64,
    pub cumulative_curve: Vec<CumulativeExposurePoint>,
    /// Lag-response curve `RR(p95_exposure, lag)` for `lag = 0..=max_lag`.
    pub p95_lag_slice: Vec<LagSlicePoint>,
    /// 2D exposure-lag `log(RR)` surface `[exposure_grid_idx][lag]` for contour/heatmap rendering.
    pub surface_exposures: Vec<f64>,
    pub surface_rr_grid: Vec<Vec<f64>>,
}

fn raw_exposure_basis(x: f64, x_min: f64, x_max: f64, x_ref: f64, kind: ExposureBasisKind) -> [f64; 3] {
    let span = (x_max - x_min).max(1e-3);
    match kind {
        ExposureBasisKind::NaturalCubicSpline3Df => {
            let u = (x - x_min) / span;
            let k1 = 0.0_f64;
            let k2 = 0.35_f64;
            let k3 = 0.70_f64;
            let k4 = 1.0_f64;
            let pos3 = |v: f64| if v > 0.0 { v * v * v } else { 0.0 };
            let d = |v: f64, knot: f64| (pos3(v - knot) - pos3(v - k4)) / (k4 - knot).max(1e-9);
            let d1 = d(u, k1);
            [u, d(u, k2) - d1, d(u, k3) - d1]
        }
        ExposureBasisKind::QuadraticThreshold => {
            let diff = (x - x_ref) / span;
            let above = diff.max(0.0);
            let below = (-diff).max(0.0);
            [diff, above * above, below * below]
        }
    }
}

/// Centered exposure basis `b_x(x) - b_x(x_ref)` so `f(x_ref) = 0` identically.
fn centered_exposure_basis(
    x: f64,
    x_min: f64,
    x_max: f64,
    x_ref: f64,
    kind: ExposureBasisKind,
) -> [f64; 3] {
    let bx = raw_exposure_basis(x, x_min, x_max, x_ref, kind);
    let bref = raw_exposure_basis(x_ref, x_min, x_max, x_ref, kind);
    [bx[0] - bref[0], bx[1] - bref[1], bx[2] - bref[2]]
}

fn lag_basis_vec(lag: usize, max_lag: usize, kind: LagBasisKind) -> [f64; 3] {
    let l_max = (max_lag.max(1)) as f64;
    match kind {
        LagBasisKind::NaturalCubicLogLag3Df => {
            let u_max = (l_max + 1.0).ln().max(1e-6);
            let u = ((lag as f64) + 1.0).ln() / u_max;
            let pos3 = |v: f64| if v > 0.0 { v * v * v } else { 0.0 };
            // Intercept + linear log-lag + interior knot at u = 0.50
            [1.0, u, pos3(u - 0.50)]
        }
        LagBasisKind::QuadraticPolynomial => {
            let s = (lag as f64) / l_max;
            [1.0, s, s * s]
        }
    }
}

fn solve_linear_and_invert(a_in: &[Vec<f64>], b_in: &[f64]) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
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

fn sorted_quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let pos = q.clamp(0.0, 1.0) * ((sorted.len() - 1) as f64);
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let frac = pos - (lo as f64);
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

/// Fit a Gasparrini Cross-Basis Distributed Lag Non-Linear Model (`DLNM`) via Quasi-Poisson IRLS,
/// returning cumulative exposure-response curves, acute-vs-harvesting lag slices, 2D RR surface,
/// and forward attributable burden.
pub fn fit_dlnm_quasi_poisson(
    dataset: &DlnmDataset,
    cfg: &DlnmConfig,
) -> Result<DlnmFitResult, EpiError> {
    let n_total = dataset.observations.len();
    let max_lag = cfg.max_lag.clamp(1, 45);
    if n_total <= max_lag + 12 {
        return Err(EpiError::Schema(format!(
            "DLNM requires more than max_lag + 12 observations (got {}, max_lag={})",
            n_total, max_lag
        )));
    }

    let mut sorted_exp: Vec<f64> = dataset.observations.iter().map(|o| o.exposure).collect();
    sorted_exp.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let x_min = *sorted_exp.first().unwrap_or(&0.0);
    let x_max = *sorted_exp.last().unwrap_or(&1.0);
    if (x_max - x_min).abs() < 1e-9 {
        return Err(EpiError::Schema(
            "Exposure series is constant; cannot construct cross-basis".to_string(),
        ));
    }

    let x_median = sorted_quantile(&sorted_exp, 0.50);
    let x_p95 = sorted_quantile(&sorted_exp, 0.95);
    let x_ref = cfg.reference_exposure.unwrap_or(x_median).clamp(x_min, x_max);

    let v_dim = 3usize;
    let l_dim = 3usize;
    let cb_dim = v_dim * l_dim; // 9 cross-basis parameters
    let has_confounder = dataset
        .observations
        .iter()
        .all(|o| o.covariate_temp_or_confounder.is_some());
    let n_harm = cfg.seasonal_harmonics.min(4);
    // Design columns: [0..cb_dim] cross-basis, then intercept, linear trend, 2*n_harm sin/cos, optional confounder
    let extra_cols = 2 + 2 * n_harm + usize::from(has_confounder);
    let p_total = cb_dim + extra_cols;

    let lag_mats: Vec<[f64; 3]> = (0..=max_lag)
        .map(|l| lag_basis_vec(l, max_lag, cfg.lag_basis))
        .collect();

    let n_eff = n_total - max_lag;
    let mut design_x = Vec::with_capacity(n_eff);
    let mut y_vec = Vec::with_capacity(n_eff);
    let period = cfg.seasonal_period_steps.max(4.0);

    for t in max_lag..n_total {
        let mut row = vec![0.0_f64; p_total];
        for lag in 0..=max_lag {
            let x_tl = dataset.observations[t - lag].exposure;
            let bx = centered_exposure_basis(x_tl, x_min, x_max, x_ref, cfg.exposure_basis);
            let bl = &lag_mats[lag];
            for j in 0..v_dim {
                for m in 0..l_dim {
                    row[j * l_dim + m] += bx[j] * bl[m];
                }
            }
        }
        let mut col_ptr = cb_dim;
        row[col_ptr] = 1.0; // Intercept
        col_ptr += 1;
        row[col_ptr] = ((t as f64) - 0.5 * (n_total as f64)) / (n_total as f64); // Centered trend
        col_ptr += 1;

        for h in 1..=n_harm {
            let ang = 2.0 * std::f64::consts::PI * (h as f64) * (t as f64) / period;
            row[col_ptr] = ang.sin();
            row[col_ptr + 1] = ang.cos();
            col_ptr += 2;
        }
        if has_confounder {
            row[col_ptr] = dataset.observations[t]
                .covariate_temp_or_confounder
                .unwrap_or(0.0);
        }

        design_x.push(row);
        y_vec.push(dataset.observations[t].outcome_count.max(0.0));
    }

    let mean_y = (y_vec.iter().sum::<f64>() / (n_eff as f64)).max(0.5);
    let mut beta = vec![0.0_f64; p_total];
    beta[cb_dim] = mean_y.ln();

    let ridge = cfg.ridge_l2_penalty.max(1e-6);
    let mut xtwx_inv = vec![vec![0.0; p_total]; p_total];

    for _iter in 0..45 {
        let mut grad = vec![0.0_f64; p_total];
        let mut info = vec![vec![0.0_f64; p_total]; p_total];

        for i in 0..n_eff {
            let xr = &design_x[i];
            let eta: f64 = xr
                .iter()
                .zip(&beta)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                .clamp(-15.0, 15.0);
            let mu = eta.exp();
            let resid = y_vec[i] - mu;
            for a in 0..p_total {
                grad[a] += resid * xr[a];
                for b in a..p_total {
                    info[a][b] += mu * xr[a] * xr[b];
                }
            }
        }
        for a in 0..p_total {
            for b in 0..a {
                info[a][b] = info[b][a];
            }
            // Do not penalize intercept (col == cb_dim)
            let pen = if a == cb_dim { 1e-8 } else { ridge };
            grad[a] -= pen * beta[a];
            info[a][a] += pen;
        }

        let Some((step, inv)) = solve_linear_and_invert(&info, &grad) else {
            break;
        };
        xtwx_inv = inv;

        let max_step = step.iter().map(|s| s.abs()).fold(0.0_f64, f64::max);
        let damp = if max_step > 2.0 { 2.0 / max_step } else { 1.0 };
        for a in 0..p_total {
            beta[a] += damp * step[a];
        }
        if max_step * damp < 1e-6 {
            break;
        }
    }

    // Compute Pearson overdispersion phi and Poisson deviance
    let mut pearson_chi2 = 0.0;
    let mut deviance = 0.0;
    let mut mu_hat = vec![0.0; n_eff];
    for i in 0..n_eff {
        let eta: f64 = design_x[i]
            .iter()
            .zip(&beta)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            .clamp(-15.0, 15.0);
        let mu = eta.exp().max(1e-9);
        mu_hat[i] = mu;
        let y = y_vec[i];
        pearson_chi2 += (y - mu).powi(2) / mu;
        deviance += if y > 0.0 {
            2.0 * (y * (y / mu).ln() - (y - mu))
        } else {
            2.0 * mu
        };
    }
    let df_resid = (n_eff.saturating_sub(p_total)).max(1) as f64;
    let phi = (pearson_chi2 / df_resid).max(1.0);

    // Extract cross-basis parameter subvector theta [0..cb_dim] and covariance V_theta = phi * inv[0..cb_dim][0..cb_dim]
    let theta = &beta[0..cb_dim];
    let mut v_theta = vec![vec![0.0; cb_dim]; cb_dim];
    for a in 0..cb_dim {
        for b in 0..cb_dim {
            v_theta[a][b] = phi * xtwx_inv[a][b];
        }
    }

    let z95 = normal_quantile(0.975);
    let eval_point_and_se = |w_cb: &[f64; 9]| -> (f64, f64) {
        let log_rr: f64 = (0..cb_dim).map(|k| w_cb[k] * theta[k]).sum();
        let mut var = 0.0;
        for a in 0..cb_dim {
            for b in 0..cb_dim {
                var += w_cb[a] * v_theta[a][b] * w_cb[b];
            }
        }
        (log_rr, var.max(0.0).sqrt())
    };

    // Pre-sum lag basis across l = 0..=max_lag for cumulative exposure-response evaluation
    let mut sum_lag_basis = [0.0_f64; 3];
    for bl in &lag_mats {
        for m in 0..l_dim {
            sum_lag_basis[m] += bl[m];
        }
    }

    let grid_n = 33usize;
    let mut surface_exposures = Vec::with_capacity(grid_n);
    let mut cumulative_curve = Vec::with_capacity(grid_n);
    let mut surface_rr_grid = Vec::with_capacity(grid_n);
    let mut min_cum_rr = f64::INFINITY;
    let mut mme_x = x_ref;

    for g in 0..grid_n {
        let xg = x_min + (x_max - x_min) * (g as f64) / ((grid_n - 1) as f64);
        surface_exposures.push(xg);
        let bx = centered_exposure_basis(xg, x_min, x_max, x_ref, cfg.exposure_basis);

        let mut cum_w = [0.0_f64; 9];
        for j in 0..v_dim {
            for m in 0..l_dim {
                cum_w[j * l_dim + m] = bx[j] * sum_lag_basis[m];
            }
        }
        let (cum_lrr, cum_se) = eval_point_and_se(&cum_w);
        let cum_rr = cum_lrr.exp();
        if cum_rr < min_cum_rr {
            min_cum_rr = cum_rr;
            mme_x = xg;
        }
        cumulative_curve.push(CumulativeExposurePoint {
            exposure: xg,
            cumulative_rr: cum_rr,
            rr_ci_low_95: (cum_lrr - z95 * cum_se).exp(),
            rr_ci_high_95: (cum_lrr + z95 * cum_se).exp(),
        });

        let mut lag_row = Vec::with_capacity(max_lag + 1);
        for lag in 0..=max_lag {
            let bl = &lag_mats[lag];
            let mut lw = [0.0_f64; 9];
            for j in 0..v_dim {
                for m in 0..l_dim {
                    lw[j * l_dim + m] = bx[j] * bl[m];
                }
            }
            let (lrr, _) = eval_point_and_se(&lw);
            lag_row.push(lrr.exp());
        }
        surface_rr_grid.push(lag_row);
    }

    // Evaluate P95 cumulative RR and lag-by-lag slice
    let bx_p95 = centered_exposure_basis(x_p95, x_min, x_max, x_ref, cfg.exposure_basis);
    let mut p95_cum_w = [0.0_f64; 9];
    for j in 0..v_dim {
        for m in 0..l_dim {
            p95_cum_w[j * l_dim + m] = bx_p95[j] * sum_lag_basis[m];
        }
    }
    let (p95_lrr, p95_se) = eval_point_and_se(&p95_cum_w);

    let mut p95_lag_slice = Vec::with_capacity(max_lag + 1);
    for lag in 0..=max_lag {
        let bl = &lag_mats[lag];
        let mut lw = [0.0_f64; 9];
        for j in 0..v_dim {
            for m in 0..l_dim {
                lw[j * l_dim + m] = bx_p95[j] * bl[m];
            }
        }
        let (lrr, lse) = eval_point_and_se(&lw);
        p95_lag_slice.push(LagSlicePoint {
            lag,
            rr: lrr.exp(),
            rr_ci_low_95: (lrr - z95 * lse).exp(),
            rr_ci_high_95: (lrr + z95 * lse).exp(),
        });
    }

    // Backward/overall attributable events across fitted series: AN = sum_t y_t * (1 - exp(-eta_cb_t))
    let mut attr_events = 0.0;
    let mut total_events = 0.0;
    for i in 0..n_eff {
        let eta_cb: f64 = (0..cb_dim).map(|k| design_x[i][k] * theta[k]).sum();
        total_events += y_vec[i];
        if eta_cb > 0.0 {
            attr_events += mu_hat[i] * (1.0 - (-eta_cb).exp());
        }
    }
    let af_pct = if total_events > 0.0 {
        (attr_events / total_events) * 100.0
    } else {
        0.0
    };

    Ok(DlnmFitResult {
        title: dataset.title.clone(),
        exposure_unit: dataset.exposure_unit.clone(),
        effective_rows_fitted: n_eff,
        max_lag,
        reference_exposure: x_ref,
        minimum_risk_exposure: mme_x,
        p95_exposure: x_p95,
        p95_cumulative_rr: p95_lrr.exp(),
        p95_rr_ci_low_95: (p95_lrr - z95 * p95_se).exp(),
        p95_rr_ci_high_95: (p95_lrr + z95 * p95_se).exp(),
        dispersion_phi: phi,
        deviance,
        attributable_events: attr_events,
        attributable_fraction_pct: af_pct,
        cumulative_curve,
        p95_lag_slice,
        surface_exposures,
        surface_rr_grid,
    })
}

/// Parse daily/weekly environmental exposure + health outcome CSV/TSV (`date,exposure,outcome[,confounder]`).
pub fn parse_dlnm_exposure_csv(
    text: &str,
    title: &str,
    exposure_unit: &str,
) -> Result<DlnmDataset, EpiError> {
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let hdr = lines
        .next()
        .ok_or_else(|| EpiError::Schema("Empty DLNM CSV file".to_string()))?;
    let delim = if hdr.contains('\t') { '\t' } else { ',' };

    let mut observations = Vec::new();
    for line in lines {
        let cols: Vec<&str> = line.split(delim).map(|s| s.trim().trim_matches('"')).collect();
        if cols.len() < 3 {
            continue;
        }
        let Ok(exp) = cols[1].parse::<f64>() else { continue };
        let Ok(out_c) = cols[2].parse::<f64>() else { continue };
        let conf = cols.get(3).and_then(|s| s.parse::<f64>().ok());
        observations.push(DlnmObservation {
            date_label: cols[0].to_string(),
            exposure: exp,
            outcome_count: out_c.max(0.0),
            covariate_temp_or_confounder: conf,
        });
    }
    if observations.len() < 20 {
        return Err(EpiError::Schema(
            "DLNM CSV requires at least 20 valid observation rows".to_string(),
        ));
    }
    Ok(DlnmDataset {
        title: title.to_string(),
        exposure_unit: exposure_unit.to_string(),
        outcome_label: "Daily Cardiorespiratory ED Visits".to_string(),
        observations,
    })
}

/// Deterministic 365-day Pacific Northwest Wildfire Smoke PM2.5 (`ug/m3`) & Cardiorespiratory ED benchmark
/// exhibiting acute lag `0..4` day risk elevation and mild lag `8..11` day harvesting displacement.
pub fn builtin_wildfire_pm25_dlnm_dataset() -> DlnmDataset {
    let n = 365usize;
    let mut pm25 = vec![0.0_f64; n];
    // Deterministic LCG for realistic daily noise + late-summer wildfire smoke episodes
    let mut seed: u64 = 0x5EE1_2026_D1A4;
    let mut next_u01 = || -> f64 {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 11) as f64) / ((1u64 << 53) as f64)
    };

    for d in 0..n {
        let base = 7.5 + 2.5 * (2.0 * std::f64::consts::PI * (d as f64) / 365.25).cos();
        // Multi-day wildfire smoke plumes around days 215..224 and 248..258
        let plume1 = 92.0 * (-(((d as f64) - 219.0) / 3.8).powi(2)).exp();
        let plume2 = 128.0 * (-(((d as f64) - 252.0) / 4.2).powi(2)).exp();
        let plume3 = 48.0 * (-(((d as f64) - 190.0) / 2.8).powi(2)).exp();
        let jitter = (next_u01() - 0.5) * 3.6;
        pm25[d] = (base + plume1 + plume2 + plume3 + jitter).max(2.0);
    }

    // Distributed lag response kernel: strong acute positive effect at lags 0..4, mild harvesting dip at lags 8..11
    let true_lag_weights: [f64; 15] = [
        0.0028, 0.0024, 0.0018, 0.0011, 0.0006, 0.0002, 0.0, -0.0002, -0.0004, -0.0004, -0.0003,
        -0.0001, 0.0, 0.0, 0.0,
    ];

    let mut observations = Vec::with_capacity(n);
    for d in 0..n {
        let seasonal_base = 44.0
            + 8.5 * (2.0 * std::f64::consts::PI * ((d as f64) + 15.0) / 365.25).cos();
        let mut log_rr = 0.0;
        for (lag, &w) in true_lag_weights.iter().enumerate() {
            if d >= lag {
                let excess_pm = (pm25[d - lag] - 8.0).max(0.0);
                // Slightly supralinear/saturating non-linear smoke curve
                log_rr += w * excess_pm.powf(1.06);
            }
        }
        let mu = seasonal_base * log_rr.exp();
        // Approximate Poisson draw via normal perturbation rounded >= 0
        let z = (next_u01() + next_u01() + next_u01() - 1.5) * 2.0;
        let count = (mu + z * mu.sqrt()).round().max(1.0);
        let temp_c = 13.0 - 9.0 * (2.0 * std::f64::consts::PI * (d as f64) / 365.25).cos();
        observations.push(DlnmObservation {
            date_label: format!("2026-D{:03}", d + 1),
            exposure: (pm25[d] * 10.0).round() / 10.0,
            outcome_count: count,
            covariate_temp_or_confounder: Some((temp_c * 10.0).round() / 10.0),
        });
    }

    DlnmDataset {
        title: "Pacific NW Wildfire PM2.5 & Cardiorespiratory ED Visits (365d)".to_string(),
        exposure_unit: "ug/m3 PM2.5".to_string(),
        outcome_label: "Daily Cardiorespiratory ED Visits".to_string(),
        observations,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dlnm_cross_basis_acute_smoke_risk() {
        let ds = builtin_wildfire_pm25_dlnm_dataset();
        let cfg = DlnmConfig::default();
        let res = fit_dlnm_quasi_poisson(&ds, &cfg).expect("fit DLNM");
        assert_eq!(res.effective_rows_fitted, 365 - 14);
        assert!(res.p95_cumulative_rr > 1.05);
        assert!(res.p95_lag_slice[0].rr > res.p95_lag_slice[10].rr);
        assert!(res.attributable_events > 10.0);
    }
}
