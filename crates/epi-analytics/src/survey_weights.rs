use crate::math::normal_quantile;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Complex survey Taylor-series linearization (delta method) estimator for stratified
/// multi-stage cluster surveys such as CDC Continuous NHANES (`SDMVSTRA` / `SDMVPSU` / `WTMEC2YR`)
/// and CDC BRFSS (`_STSTR` / `_PSU` / `_LLCPWT`).
///
/// Follows NCHS Variance Estimation Guidelines (https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx):
/// - Retains full stratum/PSU design structure during subpopulation (`in_subpop`) domain estimation (`NOMCAR`).
/// - Supports conservative grand-mean centering for lonely/singleton PSUs (`n_h == 1`).

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LonelyPsuStrategy {
    /// Center singleton PSU residuals at the grand PSU mean across all strata (conservative NCHS/R `lonely.psu="adjust"`).
    GrandMeanCenter,
    /// Contribute zero variance from singleton-PSU strata (`lonely.psu="remove"`).
    ZeroContribution,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SurveyObservation {
    pub stratum: i64,
    pub psu: i64,
    pub weight: f64,
    pub value: f64,
    pub in_subpop: bool,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ComplexSurveyEstimate {
    pub weighted_mean_or_prop: f64,
    pub standard_error: f64,
    pub ci_low_95: f64,
    pub ci_high_95: f64,
    pub degrees_of_freedom: usize,
    pub n_strata: usize,
    pub n_psus: usize,
    pub lonely_strata_count: usize,
    pub subpop_sample_size: usize,
    pub sum_of_weights: f64,
    pub design_effect_deff: f64,
}

/// Compute survey-weighted mean (or prevalence proportion when `value` is `0.0` or `1.0`)
/// and Taylor-series linearized standard error across stratified PSUs.
pub fn compute_taylor_linearized_mean(
    obs: &[SurveyObservation],
    lonely_strategy: LonelyPsuStrategy,
) -> Option<ComplexSurveyEstimate> {
    // Filter strictly positive weights per NCHS rule while keeping out-of-subpopulation rows!
    let valid: Vec<&SurveyObservation> = obs
        .iter()
        .filter(|o| o.weight.is_finite() && o.weight > 0.0 && o.value.is_finite())
        .collect();
    if valid.is_empty() {
        return None;
    }

    let mut sum_w_sub = 0.0;
    let mut sum_wy_sub = 0.0;
    let mut n_sub = 0usize;

    for o in &valid {
        if o.in_subpop {
            sum_w_sub += o.weight;
            sum_wy_sub += o.weight * o.value;
            n_sub += 1;
        }
    }

    if n_sub < 2 || sum_w_sub <= 0.0 {
        return None;
    }

    let y_bar = sum_wy_sub / sum_w_sub;

    // Linearized Taylor score for ratio estimator y_bar = Y_hat / W_hat:
    // z_i = d_i * w_i * (y_i - y_bar) / W_hat
    // Summed within each (stratum h, PSU alpha) -> z_{h, alpha}
    let mut strata_psu_totals: BTreeMap<i64, BTreeMap<i64, f64>> = BTreeMap::new();
    for o in &valid {
        let z_i = if o.in_subpop {
            o.weight * (o.value - y_bar) / sum_w_sub
        } else {
            0.0
        };
        let psu_map = strata_psu_totals.entry(o.stratum).or_default();
        *psu_map.entry(o.psu).or_insert(0.0) += z_i;
    }

    let n_strata = strata_psu_totals.len();
    let total_psus: usize = strata_psu_totals.values().map(|m| m.len()).sum();

    let all_psu_scores: Vec<f64> = strata_psu_totals
        .values()
        .flat_map(|m| m.values().copied())
        .collect();
    let grand_psu_mean = if all_psu_scores.is_empty() {
        0.0
    } else {
        all_psu_scores.iter().sum::<f64>() / (all_psu_scores.len() as f64)
    };

    let mut var_tsl = 0.0;
    let mut lonely_count = 0usize;

    for psu_map in strata_psu_totals.values() {
        let m_h = psu_map.len();
        if m_h >= 2 {
            let m_f = m_h as f64;
            let z_h_bar = psu_map.values().sum::<f64>() / m_f;
            let sum_sq: f64 = psu_map
                .values()
                .map(|&z| (z - z_h_bar).powi(2))
                .sum();
            var_tsl += (m_f / (m_f - 1.0)) * sum_sq;
        } else if m_h == 1 {
            lonely_count += 1;
            if lonely_strategy == LonelyPsuStrategy::GrandMeanCenter {
                let z_only = *psu_map.values().next().unwrap_or(&0.0);
                var_tsl += (z_only - grand_psu_mean).powi(2);
            }
        }
    }

    // Simple Random Sampling (SRS) variance denominator for Design Effect (DEFF)
    let mut weighted_var_num = 0.0;
    for o in &valid {
        if o.in_subpop {
            weighted_var_num += o.weight * (o.value - y_bar).powi(2);
        }
    }
    let s2_sub = weighted_var_num / sum_w_sub;
    let var_srs = (s2_sub / (n_sub as f64)).max(1e-18);
    let deff = if s2_sub > 1e-15 {
        (var_tsl / var_srs).max(0.05)
    } else {
        1.0
    };

    let se = var_tsl.max(0.0).sqrt();
    let df = total_psus.saturating_sub(n_strata).max(1);
    // Cornish-Fisher Student-t critical value for 97.5% quantile with `df` degrees of freedom
    let z975 = normal_quantile(0.975);
    let df_f = df as f64;
    let t_crit = z975 + (z975.powi(3) + z975) / (4.0 * df_f)
        + (5.0 * z975.powi(5) + 16.0 * z975.powi(3) + 3.0 * z975) / (96.0 * df_f * df_f);

    Some(ComplexSurveyEstimate {
        weighted_mean_or_prop: y_bar,
        standard_error: se,
        ci_low_95: y_bar - t_crit * se,
        ci_high_95: y_bar + t_crit * se,
        degrees_of_freedom: df,
        n_strata,
        n_psus: total_psus,
        lonely_strata_count: lonely_count,
        subpop_sample_size: n_sub,
        sum_of_weights: sum_w_sub,
        design_effect_deff: deff,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_taylor_linearization_clustered_design_effect() {
        // Construct 4 strata x 2 PSUs = 8 PSUs with intra-PSU clustering -> DEFF > 1.0
        let mut rows = Vec::new();
        for s in 1..=4 {
            for p in 1..=2 {
                let cluster_val = if p == 1 { 0.80 } else { 0.20 };
                for k in 0..10 {
                    rows.push(SurveyObservation {
                        stratum: s,
                        psu: p,
                        weight: 1500.0 + (k as f64) * 50.0,
                        value: cluster_val,
                        in_subpop: true,
                    });
                }
            }
        }
        let est =
            compute_taylor_linearized_mean(&rows, LonelyPsuStrategy::GrandMeanCenter).unwrap();
        assert!((est.weighted_mean_or_prop - 0.50).abs() < 0.01);
        assert_eq!(est.degrees_of_freedom, 4);
        assert!(est.design_effect_deff > 2.0);
    }
}
