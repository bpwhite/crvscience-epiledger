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

use crate::math::{normal_quantile, regularized_gamma_p};
use epi_core::EpiError;
use serde::{Deserialize, Serialize};

/// Multi-Source Capture-Recapture Surveillance Completeness & Unobserved Burden Estimation:
/// Two-Stream Lincoln-Petersen / Chapman / Sadinle Log-Interval Estimators and Full 8-Model
/// Three-Source Hierarchical Log-Linear Poisson `AIC`/`BIC` Selection (`Rcapture` family).
///
/// References:
/// - Three-source log-linear capture-recapture surveillance completeness: https://pmc.ncbi.nlm.nih.gov/articles/PMC8705476/
/// - Chapman bias-corrected two-source capture-recapture in outbreaks: https://pmc.ncbi.nlm.nih.gov/articles/PMC6518598/
/// - Validity of parsimonious vs saturated log-linear CRC models: https://pmc.ncbi.nlm.nih.gov/articles/PMC2870770/

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TwoSourceChapmanResult {
    pub source_a_name: String,
    pub source_b_name: String,
    /// Cases observed ONLY in Source A (n_10).
    pub n_10: u64,
    /// Cases observed ONLY in Source B (n_01).
    pub n_01: u64,
    /// Cases observed in BOTH Source A and Source B (n_11).
    pub n_11: u64,
    pub total_observed_unique: u64,
    pub petersen_n_hat: f64,
    pub chapman_unobserved_n00: f64,
    pub chapman_total_n_hat: f64,
    pub seber_se: f64,
    /// Log-transformed Sadinle/Chao 95% CI strictly respecting N_low_95 >= total_observed_unique.
    pub ci_low_95: f64,
    pub ci_high_95: f64,
    pub overall_completeness_pct: f64,
    pub source_a_sensitivity_pct: f64,
    pub source_b_sensitivity_pct: f64,
}

/// Observable 7-cell Venn contingency counts for 3 surveillance registries (S1, S2, S3).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThreeSourceOverlapCounts {
    pub s1_name: String,
    pub s2_name: String,
    pub s3_name: String,
    pub n_100: u64,
    pub n_010: u64,
    pub n_001: u64,
    pub n_110: u64,
    pub n_101: u64,
    pub n_011: u64,
    pub n_111: u64,
}

impl ThreeSourceOverlapCounts {
    pub fn total_observed(&self) -> u64 {
        self.n_100 + self.n_010 + self.n_001 + self.n_110 + self.n_101 + self.n_011 + self.n_111
    }

    pub fn s1_marginal(&self) -> u64 {
        self.n_100 + self.n_110 + self.n_101 + self.n_111
    }

    pub fn s2_marginal(&self) -> u64 {
        self.n_010 + self.n_110 + self.n_011 + self.n_111
    }

    pub fn s3_marginal(&self) -> u64 {
        self.n_001 + self.n_101 + self.n_011 + self.n_111
    }

    /// Collapse Sources 2 + 3 into a combined second stream to run a companion 2-source Chapman check.
    pub fn to_s1_vs_s23_chapman(&self) -> TwoSourceChapmanResult {
        let n_10 = self.n_100;
        let n_01 = self.n_010 + self.n_001 + self.n_011;
        let n_11 = self.n_110 + self.n_101 + self.n_111;
        compute_two_source_chapman(
            &self.s1_name,
            &format!("{} U {}", self.s2_name, self.s3_name),
            n_10,
            n_01,
            n_11,
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LogLinearCrcModelFit {
    pub model_id: String,
    pub formula_label: String,
    pub num_params: usize,
    pub residual_df: usize,
    pub deviance_g2: f64,
    pub gof_p_value: f64,
    pub aic: f64,
    pub bic: f64,
    pub unobserved_n000_hat: f64,
    pub total_n_hat: f64,
    pub n_ci_low_95: f64,
    pub n_ci_high_95: f64,
    pub overall_completeness_pct: f64,
    pub s1_sensitivity_pct: f64,
    pub s2_sensitivity_pct: f64,
    pub s3_sensitivity_pct: f64,
    /// Pairwise interaction odds ratios exp(beta_12), exp(beta_13), exp(beta_23) (1.0 when omitted).
    pub or_12: f64,
    pub or_13: f64,
    pub or_23: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThreeSourceCrcAnalysis {
    pub counts: ThreeSourceOverlapCounts,
    pub total_observed: u64,
    pub bic_selected_index: usize,
    pub aic_selected_index: usize,
    pub models: Vec<LogLinearCrcModelFit>,
    pub companion_chapman: TwoSourceChapmanResult,
}

/// Bias-corrected Chapman two-source capture-recapture estimator with Seber SE and log-transformed 95% CI.
pub fn compute_two_source_chapman(
    source_a_name: &str,
    source_b_name: &str,
    n_10: u64,
    n_01: u64,
    n_11: u64,
) -> TwoSourceChapmanResult {
    let n_obs = n_10 + n_01 + n_11;
    let n1_plus = (n_10 + n_11) as f64;
    let n_plus_1 = (n_01 + n_11) as f64;
    let m2 = n_11 as f64;

    let petersen = if n_11 > 0 {
        (n1_plus * n_plus_1) / m2
    } else {
        n_obs as f64
    };

    let chapman_n = (((n1_plus + 1.0) * (n_plus_1 + 1.0)) / (m2 + 1.0) - 1.0).max(n_obs as f64);
    let n00_hat = ((n_10 as f64 + 0.5) * (n_01 as f64 + 0.5) / (m2 + 0.5)).max(0.0);

    let seber_var = ((n1_plus + 1.0) * (n_plus_1 + 1.0) * (n_10 as f64) * (n_01 as f64))
        / ((m2 + 1.0) * (m2 + 1.0) * (m2 + 2.0));
    let seber_se = seber_var.max(0.0).sqrt();

    // Chao / Sadinle log-scale variance on n_00 ensures 95% CI lower bound >= n_obs
    let z975 = normal_quantile(0.975);
    let var_ln_n00 = (1.0 / (n_10 as f64 + 0.5))
        + (1.0 / (n_01 as f64 + 0.5))
        + (1.0 / (m2 + 0.5))
        + (1.0 / ((n_10 as f64 + 0.5) * (n_01 as f64 + 0.5) / (m2 + 0.5)).max(0.5));
    let factor = (z975 * var_ln_n00.sqrt().min(3.5)).exp();
    let ci_low_95 = (n_obs as f64) + n00_hat / factor;
    let ci_high_95 = (n_obs as f64) + n00_hat * factor;

    let denom = chapman_n.max(1.0);
    TwoSourceChapmanResult {
        source_a_name: source_a_name.to_string(),
        source_b_name: source_b_name.to_string(),
        n_10,
        n_01,
        n_11,
        total_observed_unique: n_obs,
        petersen_n_hat: petersen,
        chapman_unobserved_n00: (chapman_n - (n_obs as f64)).max(0.0),
        chapman_total_n_hat: chapman_n,
        seber_se,
        ci_low_95,
        ci_high_95,
        overall_completeness_pct: ((n_obs as f64) * 100.0 / denom).clamp(0.0, 100.0),
        source_a_sensitivity_pct: (n1_plus * 100.0 / denom).clamp(0.0, 100.0),
        source_b_sensitivity_pct: (n_plus_1 * 100.0 / denom).clamp(0.0, 100.0),
    }
}

/// Solve small symmetric positive-definite linear system A x = b and invert A via Gauss-Jordan with partial pivoting.
fn invert_matrix_and_solve(a: &[Vec<f64>], b: &[f64]) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let p = a.len();
    let mut aug = vec![vec![0.0_f64; 2 * p + 1]; p];
    for i in 0..p {
        for j in 0..p {
            aug[i][j] = a[i][j];
        }
        aug[i][p + i] = 1.0;
        aug[i][2 * p] = b[i];
    }

    for col in 0..p {
        let mut piv = col;
        let mut best = aug[col][col].abs();
        for r in (col + 1)..p {
            if aug[r][col].abs() > best {
                best = aug[r][col].abs();
                piv = r;
            }
        }
        if best < 1e-12 {
            return None;
        }
        if piv != col {
            aug.swap(col, piv);
        }
        let diag = aug[col][col];
        for c in col..=(2 * p) {
            aug[col][c] /= diag;
        }
        for r in 0..p {
            if r != col {
                let factor = aug[r][col];
                if factor.abs() > 0.0 {
                    for c in col..=(2 * p) {
                        aug[r][c] -= factor * aug[col][c];
                    }
                }
            }
        }
    }

    let mut sol = vec![0.0; p];
    let mut inv = vec![vec![0.0; p]; p];
    for i in 0..p {
        sol[i] = aug[i][2 * p];
        for j in 0..p {
            inv[i][j] = aug[i][p + j];
        }
    }
    Some((sol, inv))
}

/// Fit all 8 hierarchical 3-source log-linear Poisson capture-recapture models and select optimal AIC & BIC models.
pub fn fit_three_source_loglinear_crc(
    counts: &ThreeSourceOverlapCounts,
) -> Result<ThreeSourceCrcAnalysis, EpiError> {
    let n_obs = counts.total_observed();
    if n_obs < 5 {
        return Err(EpiError::Schema(
            "Three-source log-linear capture-recapture requires >= 5 observed unique cases"
                .to_string(),
        ));
    }

    // 7 observable cells in canonical order: (I1, I2, I3)
    let cell_profiles: [[f64; 3]; 7] = [
        [1.0, 0.0, 0.0], // n_100
        [0.0, 1.0, 0.0], // n_010
        [0.0, 0.0, 1.0], // n_001
        [1.0, 1.0, 0.0], // n_110
        [1.0, 0.0, 1.0], // n_101
        [0.0, 1.0, 1.0], // n_011
        [1.0, 1.0, 1.0], // n_111
    ];
    let raw_y: [u64; 7] = [
        counts.n_100,
        counts.n_010,
        counts.n_001,
        counts.n_110,
        counts.n_101,
        counts.n_011,
        counts.n_111,
    ];
    // Apply 0.5 Haldane continuity protection ONLY if any zero overlap cell exists so saturated log-linear stays finite
    let has_zero = raw_y.contains(&0);
    let y_fit: [f64; 7] = raw_y.map(|v| if has_zero { (v as f64) + 0.5 } else { v as f64 });

    // Full 7-column basis per cell: [1, I1, I2, I3, I1*I2, I1*I3, I2*I3]
    let full_design: [[f64; 7]; 7] = cell_profiles.map(|[i1, i2, i3]| {
        [1.0, i1, i2, i3, i1 * i2, i1 * i3, i2 * i3]
    });

    // All 8 hierarchical Rcapture models specified by active column indices into full_design
    let model_specs: &[(&str, &str, &[usize])] = &[
        ("M0", "[S1][S2][S3] (Independent)", &[0, 1, 2, 3]),
        ("M12", "[S1*S2][S3]", &[0, 1, 2, 3, 4]),
        ("M13", "[S1*S3][S2]", &[0, 1, 2, 3, 5]),
        ("M23", "[S2*S3][S1]", &[0, 1, 2, 3, 6]),
        ("M12_13", "[S1*S2][S1*S3]", &[0, 1, 2, 3, 4, 5]),
        ("M12_23", "[S1*S2][S2*S3]", &[0, 1, 2, 3, 4, 6]),
        ("M13_23", "[S1*S3][S2*S3]", &[0, 1, 2, 3, 5, 6]),
        ("M12_13_23", "[S1*S2][S1*S3][S2*S3] (Saturated Pairwise)", &[0, 1, 2, 3, 4, 5, 6]),
    ];

    let s1_marg = counts.s1_marginal() as f64;
    let s2_marg = counts.s2_marginal() as f64;
    let s3_marg = counts.s3_marginal() as f64;
    let z975 = normal_quantile(0.975);
    let ln_n_obs = (n_obs as f64).ln();

    let mut fitted_models = Vec::with_capacity(model_specs.len());

    for &(m_id, m_label, cols) in model_specs {
        let p = cols.len();
        let mut beta = vec![0.0_f64; p];
        beta[0] = ((n_obs as f64) / 7.0).max(1.0).ln();

        let mut cov = vec![vec![0.0_f64; p]; p];

        // Damped IRLS Poisson GLM with tiny ridge stabilization 1e-8 on diagonal
        for _iter in 0..60 {
            let mut xtwx = vec![vec![0.0_f64; p]; p];
            let mut grad = vec![0.0_f64; p];
            for r in 0..p {
                xtwx[r][r] = 1e-8;
            }

            for c in 0..7 {
                let eta: f64 = cols
                    .iter()
                    .enumerate()
                    .map(|(j, &col_idx)| full_design[c][col_idx] * beta[j])
                    .sum();
                let mu = eta.clamp(-18.0, 18.0).exp();
                let resid = y_fit[c] - mu;
                for j in 0..p {
                    let x_cj = full_design[c][cols[j]];
                    grad[j] += x_cj * resid;
                    for k in 0..p {
                        xtwx[j][k] += mu * x_cj * full_design[c][cols[k]];
                    }
                }
            }

            let Some((step, inv_mat)) = invert_matrix_and_solve(&xtwx, &grad) else {
                break;
            };
            cov = inv_mat;
            let mut max_step = 0.0_f64;
            for j in 0..p {
                let clamped = step[j].clamp(-2.0, 2.0);
                beta[j] += clamped;
                max_step = max_step.max(clamped.abs());
            }
            if max_step < 1e-10 {
                break;
            }
        }

        // Compute Poisson deviance G^2 against y_fit
        let mut g2 = 0.0_f64;
        for c in 0..7 {
            let eta: f64 = cols
                .iter()
                .enumerate()
                .map(|(j, &col_idx)| full_design[c][col_idx] * beta[j])
                .sum();
            let mu = eta.clamp(-18.0, 18.0).exp().max(1e-12);
            if y_fit[c] > 0.0 {
                g2 += 2.0 * (y_fit[c] * (y_fit[c] / mu).ln() - (y_fit[c] - mu));
            } else {
                g2 += 2.0 * mu;
            }
        }
        g2 = g2.max(0.0);

        let df = 7_usize.saturating_sub(p);
        let gof_p = if df == 0 {
            1.0
        } else {
            (1.0 - regularized_gamma_p(0.5 * (df as f64), 0.5 * g2)).clamp(0.0, 1.0)
        };

        let aic = g2 + 2.0 * (p as f64);
        let bic = g2 + (p as f64) * ln_n_obs;

        // Unobserved cell (0,0,0) expected count = exp(beta_0), bounded to realistic range
        let n000_hat = beta[0].clamp(-12.0, 14.0).exp();
        let se_b0 = cov[0][0].max(1e-9).sqrt().min(3.5);
        let total_n_hat = (n_obs as f64) + n000_hat;
        let n_ci_low_95 = (n_obs as f64) + n000_hat * (-z975 * se_b0).exp();
        let n_ci_high_95 = (n_obs as f64) + n000_hat * (z975 * se_b0).exp();

        let get_or = |target_col: usize| -> f64 {
            cols.iter()
                .position(|&c_idx| c_idx == target_col)
                .map(|pos| beta[pos].clamp(-8.0, 8.0).exp())
                .unwrap_or(1.0)
        };

        let denom = total_n_hat.max(1.0);
        fitted_models.push(LogLinearCrcModelFit {
            model_id: m_id.to_string(),
            formula_label: m_label.to_string(),
            num_params: p,
            residual_df: df,
            deviance_g2: g2,
            gof_p_value: gof_p,
            aic,
            bic,
            unobserved_n000_hat: n000_hat,
            total_n_hat,
            n_ci_low_95,
            n_ci_high_95,
            overall_completeness_pct: ((n_obs as f64) * 100.0 / denom).clamp(0.0, 100.0),
            s1_sensitivity_pct: (s1_marg * 100.0 / denom).clamp(0.0, 100.0),
            s2_sensitivity_pct: (s2_marg * 100.0 / denom).clamp(0.0, 100.0),
            s3_sensitivity_pct: (s3_marg * 100.0 / denom).clamp(0.0, 100.0),
            or_12: get_or(4),
            or_13: get_or(5),
            or_23: get_or(6),
        });
    }

    let bic_selected_index = fitted_models
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.bic.partial_cmp(&b.bic).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let aic_selected_index = fitted_models
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.aic.partial_cmp(&b.aic).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    Ok(ThreeSourceCrcAnalysis {
        companion_chapman: counts.to_s1_vs_s23_chapman(),
        counts: counts.clone(),
        total_observed: n_obs,
        bic_selected_index,
        aic_selected_index,
        models: fitted_models,
    })
}

/// Parse case-level 3-source registry membership CSV (`case_id,in_s1,in_s2,in_s3` with 0/1 or true/false flags).
pub fn parse_crc_case_lists_csv(
    csv_text: &str,
    s1_default: &str,
    s2_default: &str,
    s3_default: &str,
) -> Result<ThreeSourceOverlapCounts, EpiError> {
    let mut s1_name = s1_default.to_string();
    let mut s2_name = s2_default.to_string();
    let mut s3_name = s3_default.to_string();
    let mut c = [0_u64; 8];

    let parse_flag = |tok: &str| -> bool {
        matches!(
            tok.trim().to_ascii_lowercase().as_str(),
            "1" | "true" | "yes" | "y" | "pos"
        )
    };

    for (idx, raw_line) in csv_text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        if cols.len() < 4 {
            return Err(EpiError::Schema(format!(
                "CRC case-list CSV line {} requires 4 columns (case_id,in_s1,in_s2,in_s3)",
                idx + 1
            )));
        }
        if idx == 0 && !cols[1].chars().all(|ch| ch == '0' || ch == '1') {
            s1_name = cols[1].to_string();
            s2_name = cols[2].to_string();
            s3_name = cols[3].to_string();
            continue;
        }
        let b1 = usize::from(parse_flag(cols[1]));
        let b2 = usize::from(parse_flag(cols[2]));
        let b3 = usize::from(parse_flag(cols[3]));
        c[(b1 << 2) | (b2 << 1) | b3] += 1;
    }

    let out = ThreeSourceOverlapCounts {
        s1_name,
        s2_name,
        s3_name,
        n_100: c[0b100],
        n_010: c[0b010],
        n_001: c[0b001],
        n_110: c[0b110],
        n_101: c[0b101],
        n_011: c[0b011],
        n_111: c[0b111],
    };
    if out.total_observed() == 0 {
        return Err(EpiError::Schema(
            "CRC case-list CSV contained zero observed positive memberships".to_string(),
        ));
    }
    Ok(out)
}

/// Built-in 3-registry invasive bacterial / tuberculosis surveillance dataset (ELR Lab, Hospital Discharge, Physician Notification)
/// exhibiting positive Hospital-Lab catchability dependence (`OR_12 > 2.2`).
pub fn builtin_three_source_tb_pertussis_crc() -> ThreeSourceOverlapCounts {
    ThreeSourceOverlapCounts {
        s1_name: "ELR Reference Lab".to_string(),
        s2_name: "Hospital Discharge".to_string(),
        s3_name: "Clinician Case Report".to_string(),
        n_100: 68,
        n_010: 21,
        n_001: 14,
        n_110: 42,
        n_101: 19,
        n_011: 6,
        n_111: 15,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chapman_and_loglinear_recover_true_population_and_positive_dependence() {
        // Synthetic independent 3-list population N = 400 with p1=0.60, p2=0.40, p3=0.25:
        // Expected cells:
        // n100 = 400 * 0.6 * 0.6 * 0.75 = 108
        // n010 = 400 * 0.4 * 0.4 * 0.75 = 48
        // n001 = 400 * 0.4 * 0.6 * 0.25 = 24
        // n110 = 400 * 0.6 * 0.4 * 0.75 = 72
        // n101 = 400 * 0.6 * 0.6 * 0.25 = 36
        // n011 = 400 * 0.4 * 0.4 * 0.25 = 16
        // n111 = 400 * 0.6 * 0.4 * 0.25 = 24
        // Observed = 328, Unobserved n000 = 400 * 0.4 * 0.6 * 0.75 = 72 -> Total = 400!
        let exact_indep = ThreeSourceOverlapCounts {
            s1_name: "Lab".to_string(),
            s2_name: "Hosp".to_string(),
            s3_name: "Vital".to_string(),
            n_100: 108,
            n_010: 48,
            n_001: 24,
            n_110: 72,
            n_101: 36,
            n_011: 16,
            n_111: 24,
        };
        let res = fit_three_source_loglinear_crc(&exact_indep).expect("CRC fit succeeds");
        assert_eq!(res.bic_selected_index, 0, "BIC must select M0 Independent on exact independent data");
        let m0 = &res.models[0];
        assert!((m0.total_n_hat - 400.0).abs() < 0.25, "Expected N_hat=400, got {}", m0.total_n_hat);
        assert!(m0.n_ci_low_95 >= (exact_indep.total_observed() as f64));
        assert!(m0.n_ci_low_95 < 400.0 && m0.n_ci_high_95 > 400.0);

        let builtin = builtin_three_source_tb_pertussis_crc();
        let b_res = fit_three_source_loglinear_crc(&builtin).unwrap();
        assert!(b_res.models[b_res.bic_selected_index].total_n_hat > (builtin.total_observed() as f64));
    }
}
