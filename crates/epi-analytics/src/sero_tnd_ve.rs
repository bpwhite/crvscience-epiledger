use crate::math::normal_quantile;
use epi_core::EpiError;
use serde::{Deserialize, Serialize};

/// Sero-Catalytic Force-of-Infection (`FOI`) & Test-Negative Design (`TND`) Vaccine Effectiveness (`VE`):
/// 1. Griffiths/Muench age-structured catalytic and reverse-catalytic seroprevalence models estimating
///    annual force of infection (`lambda`), antibody seroreversion rate (`omega`), implied `R0`, and
///    population-weighted empirical susceptibility (`S_eff`) to calibrate `seir.rs`.
/// 2. Observational Test-Negative Design (`TND`) multivariable logistic regression (`VE = 1 - OR_adj`)
///    adjusting for calendar time, age stratum, and healthcare-seeking propensity, plus waning strata.
///
/// References:
/// - Griffiths (1974) / Muench age-structured catalytic FOI models: https://pmc.ncbi.nlm.nih.gov/articles/PMC4520443/
/// - Reverse-catalytic seroreversion / waning serosurvey models: https://pmc.ncbi.nlm.nih.gov/articles/PMC8517721/
/// - Coupling catalytic FOI with dynamic SEIR vaccination models: https://pmc.ncbi.nlm.nih.gov/articles/PMC4786291/
/// - Test-Negative Design (TND) VE estimation: https://pmc.ncbi.nlm.nih.gov/articles/PMC5013887/
/// - TND confounder control for calendar week & care-seeking: https://pmc.ncbi.nlm.nih.gov/articles/PMC7604521/
/// - Time-since-vaccination waning & recurrent TND extensions: https://pmc.ncbi.nlm.nih.gov/articles/PMC12706598/

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CatalyticModelKind {
    /// Lifelong immunity `pi(a) = 1 - exp(-lambda * a)`.
    IrreversibleCatalytic,
    /// Seroreversion / antibody waning `pi(a) = [lambda / (lambda + omega)] * (1 - exp(-(lambda + omega) * a))`.
    ReversibleWaningCatalytic,
}

impl CatalyticModelKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::IrreversibleCatalytic => "Irreversible Catalytic (Lifelong IgG)",
            Self::ReversibleWaningCatalytic => "Reverse-Catalytic (FOI + Antibody Waning)",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgeSeroprevalenceStratum {
    pub age_label: String,
    pub age_midpoint_years: f64,
    pub n_tested: u64,
    pub n_seropositive: u64,
    pub population_weight: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerosurveyDataset {
    pub title: String,
    pub strata: Vec<AgeSeroprevalenceStratum>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CatalyticCurvePoint {
    pub age_years: f64,
    pub fitted_seroprevalence: f64,
    pub ci_low_95: f64,
    pub ci_high_95: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CatalyticFitResult {
    pub title: String,
    pub model_kind: CatalyticModelKind,
    /// Annual Force of Infection `lambda` (`year^-1`).
    pub foi_lambda_per_year: f64,
    pub foi_ci_low_95: f64,
    pub foi_ci_high_95: f64,
    /// Annual seroreversion / waning rate `omega` (`year^-1`, `0.0` for irreversible).
    pub waning_omega_per_year: f64,
    pub mean_antibody_duration_years: Option<f64>,
    pub log_likelihood: f64,
    pub aic: f64,
    /// Population-weighted immune seropositive proportion across strata.
    pub weighted_seropositive_fraction: f64,
    /// Empirical susceptible fraction `1 - weighted_seropositive_fraction` for `seir.rs` calibration.
    pub effective_susceptible_fraction: f64,
    /// Implied endemic `R0 ≈ 1 + life_expectancy * lambda` (Anderson-May catalytic approximation, L = 77y).
    pub implied_endemic_r0: f64,
    pub age_curve: Vec<CatalyticCurvePoint>,
}

/// Fit an age-structured catalytic or reverse-catalytic seroprevalence model via Binomial MLE.
pub fn fit_catalytic_seroprevalence(
    dataset: &SerosurveyDataset,
    kind: CatalyticModelKind,
) -> Result<CatalyticFitResult, EpiError> {
    if dataset.strata.len() < 3 {
        return Err(EpiError::Schema(
            "Sero-catalytic model requires at least 3 age strata".to_string(),
        ));
    }

    let binom_ll = |lam: f64, omg: f64| -> f64 {
        if lam <= 1e-7 || omg < 0.0 {
            return f64::NEG_INFINITY;
        }
        let mut ll = 0.0;
        for s in &dataset.strata {
            let a = s.age_midpoint_years.max(0.25);
            let p = if omg <= 1e-12 {
                1.0 - (-lam * a).exp()
            } else {
                (lam / (lam + omg)) * (1.0 - (-(lam + omg) * a).exp())
            }
            .clamp(1e-9, 1.0 - 1e-9);
            let k = s.n_seropositive as f64;
            let n = s.n_tested as f64;
            ll += k * p.ln() + (n - k).max(0.0) * (1.0 - p).ln();
        }
        ll
    };

    let z95 = normal_quantile(0.975);

    let (best_lam, best_omg, se_lam) = match kind {
        CatalyticModelKind::IrreversibleCatalytic => {
            // Coarse log-grid bracket over [0.0005, 1.50] followed by 25 Newton-Raphson iterations
            let mut lam = 0.05_f64;
            let mut top_ll = f64::NEG_INFINITY;
            for g in 0..120 {
                let cand = 0.0005 * (1.50_f64 / 0.0005).powf((g as f64) / 119.0);
                let ll = binom_ll(cand, 0.0);
                if ll > top_ll {
                    top_ll = ll;
                    lam = cand;
                }
            }

            let mut fisher_info = 1.0_f64;
            for _ in 0..25 {
                let mut g1 = 0.0;
                let mut info = 0.0;
                for s in &dataset.strata {
                    let a = s.age_midpoint_years.max(0.25);
                    let e = (-lam * a).exp();
                    let p = (1.0 - e).clamp(1e-8, 1.0 - 1e-8);
                    let dp = a * e;
                    let k = s.n_seropositive as f64;
                    let n = s.n_tested as f64;
                    g1 += (k / p - (n - k) / (1.0 - p)) * dp;
                    // Expected Fisher information for Binomial: sum n * (dp/dlam)^2 / (p * (1 - p))
                    info += n * dp * dp / (p * (1.0 - p));
                }
                fisher_info = info.max(1e-9);
                let step = g1 / fisher_info;
                lam = (lam + step.clamp(-0.5 * lam, 0.8 * lam)).clamp(1e-5, 3.0);
                if step.abs() < 1e-8 {
                    break;
                }
            }
            (lam, 0.0, (1.0 / fisher_info).sqrt())
        }
        CatalyticModelKind::ReversibleWaningCatalytic => {
            // 2D log-grid search + coordinate / Nelder-style refinement over (lambda, omega)
            let mut best_l = 0.08_f64;
            let mut best_o = 0.02_f64;
            let mut top_ll = f64::NEG_INFINITY;
            for gl in 0..55 {
                let cand_l = 0.001 * (1.20_f64 / 0.001).powf((gl as f64) / 54.0);
                for go in 0..45 {
                    let cand_o = if go == 0 {
                        0.0
                    } else {
                        0.0005 * (0.50_f64 / 0.0005).powf(((go - 1) as f64) / 43.0)
                    };
                    let ll = binom_ll(cand_l, cand_o);
                    if ll > top_ll {
                        top_ll = ll;
                        best_l = cand_l;
                        best_o = cand_o;
                    }
                }
            }

            let mut step_l = best_l * 0.20;
            let mut step_o = (best_o.max(0.002)) * 0.20;
            for _ in 0..40 {
                let mut improved = false;
                for &(dl, dw) in &[
                    (1.0, 0.0),
                    (-1.0, 0.0),
                    (0.0, 1.0),
                    (0.0, -1.0),
                    (0.7, 0.7),
                    (0.7, -0.7),
                    (-0.7, 0.7),
                    (-0.7, -0.7),
                ] {
                    let cl = (best_l + dl * step_l).clamp(1e-5, 3.0);
                    let co = (best_o + dw * step_o).clamp(0.0, 1.5);
                    let ll = binom_ll(cl, co);
                    if ll > top_ll + 1e-10 {
                        top_ll = ll;
                        best_l = cl;
                        best_o = co;
                        improved = true;
                    }
                }
                if !improved {
                    step_l *= 0.5;
                    step_o *= 0.5;
                    if step_l < 1e-7 && step_o < 1e-7 {
                        break;
                    }
                }
            }

            // Compute 2x2 expected Fisher Information at (best_l, best_o) for marginal SE(lambda)
            let mut i_ll = 1e-6;
            let mut i_lo = 0.0;
            let mut i_oo = 1e-6;
            let eps = 1e-5;
            for s in &dataset.strata {
                let a = s.age_midpoint_years.max(0.25);
                let eval_p = |l: f64, o: f64| -> f64 {
                    if o <= 1e-12 {
                        1.0 - (-l * a).exp()
                    } else {
                        (l / (l + o)) * (1.0 - (-(l + o) * a).exp())
                    }
                };
                let p = eval_p(best_l, best_o).clamp(1e-7, 1.0 - 1e-7);
                let dp_dl = (eval_p(best_l + eps, best_o) - eval_p((best_l - eps).max(1e-7), best_o))
                    / (2.0 * eps);
                let dp_do = (eval_p(best_l, best_o + eps) - eval_p(best_l, (best_o - eps).max(0.0)))
                    / (best_o + eps - (best_o - eps).max(0.0)).max(1e-7);
                let w = (s.n_tested as f64) / (p * (1.0 - p));
                i_ll += w * dp_dl * dp_dl;
                i_lo += w * dp_dl * dp_do;
                i_oo += w * dp_do * dp_do;
            }
            let det = (i_ll * i_oo - i_lo * i_lo).max(1e-9);
            let var_l = (i_oo / det).max(1e-10);
            (best_l, best_o, var_l.sqrt())
        }
    };

    let final_ll = binom_ll(best_lam, best_omg);
    let k_params = match kind {
        CatalyticModelKind::IrreversibleCatalytic => 1.0,
        CatalyticModelKind::ReversibleWaningCatalytic => 2.0,
    };
    let aic = 2.0 * k_params - 2.0 * final_ll;

    let eval_prev = |lam: f64, omg: f64, age: f64| -> f64 {
        if omg <= 1e-12 {
            1.0 - (-lam * age).exp()
        } else {
            (lam / (lam + omg)) * (1.0 - (-(lam + omg) * age).exp())
        }
        .clamp(0.0, 1.0)
    };

    let lam_lo = (best_lam - z95 * se_lam).max(1e-5);
    let lam_hi = best_lam + z95 * se_lam;

    let mut weighted_sero = 0.0;
    let mut sum_w = 0.0;
    for s in &dataset.strata {
        let w = s.population_weight.max(1e-6);
        weighted_sero += w * eval_prev(best_lam, best_omg, s.age_midpoint_years);
        sum_w += w;
    }
    let sero_frac = (weighted_sero / sum_w.max(1e-9)).clamp(0.0, 1.0);
    let eff_susc = 1.0 - sero_frac;

    let max_age = dataset
        .strata
        .iter()
        .map(|s| s.age_midpoint_years)
        .fold(80.0_f64, f64::max);
    let mut age_curve = Vec::with_capacity(41);
    for g in 0..=40 {
        let age = (max_age * (g as f64) / 40.0).max(0.25);
        age_curve.push(CatalyticCurvePoint {
            age_years: age,
            fitted_seroprevalence: eval_prev(best_lam, best_omg, age),
            ci_low_95: eval_prev(lam_lo, best_omg, age),
            ci_high_95: eval_prev(lam_hi, best_omg, age),
        });
    }

    Ok(CatalyticFitResult {
        title: dataset.title.clone(),
        model_kind: kind,
        foi_lambda_per_year: best_lam,
        foi_ci_low_95: lam_lo,
        foi_ci_high_95: lam_hi,
        waning_omega_per_year: best_omg,
        mean_antibody_duration_years: if best_omg > 1e-5 {
            Some(1.0 / best_omg)
        } else {
            None
        },
        log_likelihood: final_ll,
        aic,
        weighted_seropositive_fraction: sero_frac,
        effective_susceptible_fraction: eff_susc,
        implied_endemic_r0: 1.0 + 77.0 * best_lam,
        age_curve,
    })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TndCellRecord {
    /// Time-since-vaccination group (`0` = unvaccinated, `1` = 14-89d, `2` = 90-179d, `3` = 180d+).
    pub vax_window_group: u8,
    pub age_senior: bool,
    /// Epidemic calendar epoch (`0..3`) controlling time-varying pathogen circulation & vax rollout.
    pub calendar_epoch: u8,
    /// High healthcare-seeking / chronic-condition indicator (`0` or `1`).
    pub high_care_seeking: bool,
    pub test_positive_cases: u64,
    pub test_negative_controls: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TndDataset {
    pub title: String,
    pub cells: Vec<TndCellRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TndWaningStratumVe {
    pub label: String,
    pub cases_vax: u64,
    pub controls_vax: u64,
    pub adjusted_ve_pct: f64,
    pub ve_ci_low_95: f64,
    pub ve_ci_high_95: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TndVeResult {
    pub title: String,
    pub total_cases: u64,
    pub total_controls: u64,
    pub crude_unadjusted_ve_pct: f64,
    pub crude_ve_ci_low_95: f64,
    pub crude_ve_ci_high_95: f64,
    /// Multivariable logistic regression VE adjusting for calendar epoch, age, and care-seeking propensity.
    pub adjusted_overall_ve_pct: f64,
    pub adjusted_ve_ci_low_95: f64,
    pub adjusted_ve_ci_high_95: f64,
    /// Bias removed by confounder adjustment (`adjusted_ve - crude_ve` in percentage points).
    pub confounding_bias_removed_pp: f64,
    pub waning_strata: Vec<TndWaningStratumVe>,
}

fn solve_linear_and_invert(a_in: &[Vec<f64>], b_in: &[f64]) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = b_in.len();
    let mut aug = vec![vec![0.0; 2 * n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a_in[i][j];
        }
        aug[i][n + i] = 1.0;
        aug[i][2 * n] = b_in[i];
    }
    for col in 0..n {
        let mut piv = col;
        let mut best = aug[col][col].abs();
        for r in (col + 1)..n {
            if aug[r][col].abs() > best {
                best = aug[r][col].abs();
                piv = r;
            }
        }
        if best < 1e-13 {
            return None;
        }
        if piv != col {
            aug.swap(col, piv);
        }
        let diag = aug[col][col];
        for j in col..=(2 * n) {
            aug[col][j] /= diag;
        }
        for r in 0..n {
            if r != col {
                let f = aug[r][col];
                if f.abs() > 0.0 {
                    for j in col..=(2 * n) {
                        aug[r][j] -= f * aug[col][j];
                    }
                }
            }
        }
    }
    let sol = (0..n).map(|i| aug[i][2 * n]).collect();
    let inv = (0..n)
        .map(|i| (0..n).map(|j| aug[i][n + j]).collect())
        .collect();
    Some((sol, inv))
}

fn fit_grouped_binomial_irls(
    rows: &[(Vec<f64>, f64, f64)], // (design_x, positives_y, total_trials_n)
) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let p = rows.first()?.0.len();
    let mut beta = vec![0.0_f64; p];
    let mut cov = vec![vec![0.0; p]; p];

    for _ in 0..40 {
        let mut grad = vec![0.0_f64; p];
        let mut info = vec![vec![0.0_f64; p]; p];

        for (x, y, n) in rows {
            if *n <= 0.0 {
                continue;
            }
            let eta: f64 = x
                .iter()
                .zip(&beta)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                .clamp(-14.0, 14.0);
            let prob = 1.0 / (1.0 + (-eta).exp());
            let resid = *y - *n * prob;
            let w = (*n) * prob * (1.0 - prob);
            for a in 0..p {
                grad[a] += resid * x[a];
                for b in a..p {
                    info[a][b] += w * x[a] * x[b];
                }
            }
        }
        for a in 0..p {
            for b in 0..a {
                info[a][b] = info[b][a];
            }
            let pen = if a == 0 { 1e-8 } else { 1e-5 };
            grad[a] -= pen * beta[a];
            info[a][a] += pen;
        }
        let (step, inv) = solve_linear_and_invert(&info, &grad)?;
        cov = inv;
        let max_s = step.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        for a in 0..p {
            beta[a] += step[a];
        }
        if max_s < 1e-7 {
            break;
        }
    }
    Some((beta, cov))
}

/// Estimate Crude and Confounder-Adjusted Test-Negative Design Vaccine Effectiveness (`VE = 1 - OR_adj`)
/// plus time-since-vaccination waning strata (`14-89d`, `90-179d`, `180d+`).
pub fn estimate_tnd_vaccine_effectiveness(ds: &TndDataset) -> Result<TndVeResult, EpiError> {
    if ds.cells.is_empty() {
        return Err(EpiError::Schema("Empty TND study dataset".to_string()));
    }

    let z95 = normal_quantile(0.975);
    let mut a_v_pos = 0u64;
    let mut b_v_neg = 0u64;
    let mut c_u_pos = 0u64;
    let mut d_u_neg = 0u64;

    for c in &ds.cells {
        if c.vax_window_group > 0 {
            a_v_pos += c.test_positive_cases;
            b_v_neg += c.test_negative_controls;
        } else {
            c_u_pos += c.test_positive_cases;
            d_u_neg += c.test_negative_controls;
        }
    }

    // Woolf-Haldane 0.5 continuity-corrected crude OR
    let af = a_v_pos as f64 + 0.5;
    let bf = b_v_neg as f64 + 0.5;
    let cf = c_u_pos as f64 + 0.5;
    let df = d_u_neg as f64 + 0.5;
    let ln_or_crude = ((af * df) / (bf * cf)).ln();
    let se_crude = (1.0 / af + 1.0 / bf + 1.0 / cf + 1.0 / df).sqrt();
    let crude_ve = (1.0 - ln_or_crude.exp()) * 100.0;
    let crude_ve_lo = (1.0 - (ln_or_crude + z95 * se_crude).exp()) * 100.0;
    let crude_ve_hi = (1.0 - (ln_or_crude - z95 * se_crude).exp()) * 100.0;

    // Model 1: Overall adjusted VE [intercept, any_vax, age_senior, cal_e1, cal_e2, cal_e3, high_seek]
    let overall_rows: Vec<(Vec<f64>, f64, f64)> = ds
        .cells
        .iter()
        .map(|c| {
            let x = vec![
                1.0,
                if c.vax_window_group > 0 { 1.0 } else { 0.0 },
                if c.age_senior { 1.0 } else { 0.0 },
                if c.calendar_epoch == 1 { 1.0 } else { 0.0 },
                if c.calendar_epoch == 2 { 1.0 } else { 0.0 },
                if c.calendar_epoch >= 3 { 1.0 } else { 0.0 },
                if c.high_care_seeking { 1.0 } else { 0.0 },
            ];
            let y = c.test_positive_cases as f64;
            let n = (c.test_positive_cases + c.test_negative_controls) as f64;
            (x, y, n)
        })
        .collect();

    let (b_ov, cov_ov) = fit_grouped_binomial_irls(&overall_rows)
        .ok_or_else(|| EpiError::Schema("Singular TND logistic regression matrix".to_string()))?;
    let beta_vax = b_ov[1];
    let se_vax = cov_ov[1][1].max(0.0).sqrt();
    let adj_ve = (1.0 - beta_vax.exp()) * 100.0;
    let adj_ve_lo = (1.0 - (beta_vax + z95 * se_vax).exp()) * 100.0;
    let adj_ve_hi = (1.0 - (beta_vax - z95 * se_vax).exp()) * 100.0;

    // Model 2: Time-since-vaccination waning strata [intercept, w1(14-89d), w2(90-179d), w3(180d+), age, cal1, cal2, cal3, seek]
    let waning_rows: Vec<(Vec<f64>, f64, f64)> = ds
        .cells
        .iter()
        .map(|c| {
            let x = vec![
                1.0,
                if c.vax_window_group == 1 { 1.0 } else { 0.0 },
                if c.vax_window_group == 2 { 1.0 } else { 0.0 },
                if c.vax_window_group >= 3 { 1.0 } else { 0.0 },
                if c.age_senior { 1.0 } else { 0.0 },
                if c.calendar_epoch == 1 { 1.0 } else { 0.0 },
                if c.calendar_epoch == 2 { 1.0 } else { 0.0 },
                if c.calendar_epoch >= 3 { 1.0 } else { 0.0 },
                if c.high_care_seeking { 1.0 } else { 0.0 },
            ];
            let y = c.test_positive_cases as f64;
            let n = (c.test_positive_cases + c.test_negative_controls) as f64;
            (x, y, n)
        })
        .collect();

    let mut waning_strata = Vec::new();
    if let Some((b_w, cov_w)) = fit_grouped_binomial_irls(&waning_rows) {
        let labels = [
            (1u8, "14-89 days post-dose"),
            (2u8, "90-179 days post-dose"),
            (3u8, "180+ days post-dose"),
        ];
        for (idx, (grp, lbl)) in labels.into_iter().enumerate() {
            let cases_g: u64 = ds
                .cells
                .iter()
                .filter(|c| c.vax_window_group == grp)
                .map(|c| c.test_positive_cases)
                .sum();
            let ctrls_g: u64 = ds
                .cells
                .iter()
                .filter(|c| c.vax_window_group == grp)
                .map(|c| c.test_negative_controls)
                .sum();
            if cases_g + ctrls_g > 0 {
                let bw = b_w[1 + idx];
                let sew = cov_w[1 + idx][1 + idx].max(0.0).sqrt();
                waning_strata.push(TndWaningStratumVe {
                    label: lbl.to_string(),
                    cases_vax: cases_g,
                    controls_vax: ctrls_g,
                    adjusted_ve_pct: (1.0 - bw.exp()) * 100.0,
                    ve_ci_low_95: (1.0 - (bw + z95 * sew).exp()) * 100.0,
                    ve_ci_high_95: (1.0 - (bw - z95 * sew).exp()) * 100.0,
                });
            }
        }
    }

    Ok(TndVeResult {
        title: ds.title.clone(),
        total_cases: a_v_pos + c_u_pos,
        total_controls: b_v_neg + d_u_neg,
        crude_unadjusted_ve_pct: crude_ve,
        crude_ve_ci_low_95: crude_ve_lo,
        crude_ve_ci_high_95: crude_ve_hi,
        adjusted_overall_ve_pct: adj_ve,
        adjusted_ve_ci_low_95: adj_ve_lo,
        adjusted_ve_ci_high_95: adj_ve_hi,
        confounding_bias_removed_pp: adj_ve - crude_ve,
        waning_strata,
    })
}

pub fn parse_serosurvey_csv(text: &str, title: &str) -> Result<SerosurveyDataset, EpiError> {
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let hdr = lines
        .next()
        .ok_or_else(|| EpiError::Schema("Empty serosurvey CSV".to_string()))?;
    let delim = if hdr.contains('\t') { '\t' } else { ',' };
    let mut strata = Vec::new();
    for line in lines {
        let cols: Vec<&str> = line.split(delim).map(|s| s.trim().trim_matches('"')).collect();
        if cols.len() < 4 {
            continue;
        }
        let Ok(age) = cols[1].parse::<f64>() else { continue };
        let Ok(n_test) = cols[2].parse::<u64>() else { continue };
        let Ok(n_pos) = cols[3].parse::<u64>() else { continue };
        let w = cols.get(4).and_then(|s| s.parse().ok()).unwrap_or(1.0);
        if n_test > 0 {
            strata.push(AgeSeroprevalenceStratum {
                age_label: cols[0].to_string(),
                age_midpoint_years: age,
                n_tested: n_test,
                n_seropositive: n_pos.min(n_test),
                population_weight: w,
            });
        }
    }
    if strata.len() < 3 {
        return Err(EpiError::Schema(
            "Serosurvey CSV requires at least 3 valid age strata".to_string(),
        ));
    }
    Ok(SerosurveyDataset {
        title: title.to_string(),
        strata,
    })
}

/// Built-in 9-stratum cross-sectional IgG seroprevalence survey (`lambda ≈ 0.085/yr`, `omega ≈ 0.015/yr`).
pub fn builtin_age_serosurvey_dataset() -> SerosurveyDataset {
    let raw: &[(&str, f64, u64, u64, f64)] = &[
        ("1-4y", 2.5, 420, 74, 0.055),
        ("5-9y", 7.5, 510, 219, 0.068),
        ("10-14y", 12.5, 490, 291, 0.072),
        ("15-24y", 20.0, 640, 468, 0.138),
        ("25-34y", 30.0, 710, 569, 0.136),
        ("35-44y", 40.0, 680, 568, 0.163),
        ("45-54y", 50.0, 620, 524, 0.135),
        ("55-64y", 60.0, 580, 491, 0.098),
        ("65+y", 74.0, 650, 549, 0.135),
    ];
    SerosurveyDataset {
        title: "Cross-Sectional Age-Stratified IgG Seroprevalence Survey (n=5,300)".to_string(),
        strata: raw
            .iter()
            .map(|(lbl, age, nt, np, w)| AgeSeroprevalenceStratum {
                age_label: (*lbl).to_string(),
                age_midpoint_years: *age,
                n_tested: *nt,
                n_seropositive: *np,
                population_weight: *w,
            })
            .collect(),
    }
}

/// Built-in Multi-Center Acute Symptomatic Test-Negative Design (`TND`) Vaccine Effectiveness Dataset
/// demonstrating calendar-epoch & care-seeking confounding adjustment and waning across `14-89d`, `90-179d`, `180d+`.
pub fn builtin_tnd_vaccine_study_dataset() -> TndDataset {
    let mut cells = Vec::new();
    // True conditional odds ratios: w1(14-89d)=0.28 (VE=72%), w2(90-179d)=0.46 (VE=54%), w3(180d+)=0.68 (VE=32%)
    // Confounding: vaccinated uptake is much higher in peak epidemic epoch (epoch 2-3) & seniors, biasing crude VE downward!
    let true_or = [1.0, 0.28, 0.46, 0.68];
    let epoch_attack_odds = [0.12, 0.28, 0.58, 0.42];
    let vax_dist_low_epoch = [0.64, 0.18, 0.11, 0.07];
    let vax_dist_high_epoch = [0.28, 0.26, 0.25, 0.21];

    for epoch in 0..4u8 {
        for &senior in &[false, true] {
            for &seek in &[false, true] {
                let base_n: f64 = if senior { 340.0 } else { 460.0 } * if seek { 0.7 } else { 1.0 };
                let vdist = if epoch >= 2 || senior {
                    &vax_dist_high_epoch
                } else {
                    &vax_dist_low_epoch
                };
                for grp in 0..4u8 {
                    let ctrls = (base_n * vdist[grp as usize]).round().max(4.0);
                    let odds = epoch_attack_odds[epoch as usize]
                        * true_or[grp as usize]
                        * if senior { 1.32 } else { 1.0 }
                        * if seek { 1.18 } else { 1.0 };
                    let cases = (ctrls * odds).round().max(1.0);
                    cells.push(TndCellRecord {
                        vax_window_group: grp,
                        age_senior: senior,
                        calendar_epoch: epoch,
                        high_care_seeking: seek,
                        test_positive_cases: cases as u64,
                        test_negative_controls: ctrls as u64,
                    });
                }
            }
        }
    }

    TndDataset {
        title: "Multi-Site Symptomatic Test-Negative Design (TND) VE Study".to_string(),
        cells,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_catalytic_and_reversible_foi() {
        let ds = builtin_age_serosurvey_dataset();
        let irrev = fit_catalytic_seroprevalence(&ds, CatalyticModelKind::IrreversibleCatalytic)
            .expect("irreversible fit");
        assert!(irrev.foi_lambda_per_year > 0.04 && irrev.foi_lambda_per_year < 0.15);
        assert!(irrev.effective_susceptible_fraction > 0.10 && irrev.effective_susceptible_fraction < 0.45);

        let rev = fit_catalytic_seroprevalence(&ds, CatalyticModelKind::ReversibleWaningCatalytic)
            .expect("reversible fit");
        assert!(rev.waning_omega_per_year > 0.005);
        assert!(rev.log_likelihood >= irrev.log_likelihood);
    }

    #[test]
    fn test_tnd_confounder_adjustment_and_waning() {
        let ds = builtin_tnd_vaccine_study_dataset();
        let res = estimate_tnd_vaccine_effectiveness(&ds).expect("TND VE");
        // Adjusted VE removes downward calendar-time confounding bias and recovers clear waning gradient
        assert!(res.adjusted_overall_ve_pct > res.crude_unadjusted_ve_pct + 5.0);
        assert_eq!(res.waning_strata.len(), 3);
        assert!(res.waning_strata[0].adjusted_ve_pct > res.waning_strata[1].adjusted_ve_pct);
        assert!(res.waning_strata[1].adjusted_ve_pct > res.waning_strata[2].adjusted_ve_pct);
    }
}
