use crate::math::{gamma_quantile, spearman_rho};
use epi_core::{AgeStratumCount, JurisdictionRecord};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AgeStandardizedResult {
    pub crude_rate_per_100k: f64,
    pub age_adjusted_rate_per_100k: f64,
    pub ci_low_95: f64,
    pub ci_high_95: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LagCorrelationPoint {
    pub lag_steps: i32,
    pub spearman_rho: f64,
    pub paired_n: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatialClusterScore {
    pub fips: String,
    pub name: String,
    pub observed_cases: f64,
    pub expected_cases: f64,
    pub relative_risk: f64,
    pub kulldorff_llr: f64,
}

/// Direct age standardization against standard population weights (e.g. US 2000 Standard Million)
/// with Fay-Feuer Gamma exact 95% confidence bounds.
pub fn compute_direct_age_standardized_rate(
    strata: &[AgeStratumCount],
) -> Option<AgeStandardizedResult> {
    if strata.is_empty() {
        return None;
    }

    let mut total_cases = 0.0;
    let mut total_pop = 0.0;
    let mut adj_rate = 0.0;
    let mut var_rate = 0.0;
    let mut max_w_over_n: f64 = 0.0;
    let mut weight_sum = 0.0;

    for s in strata {
        if s.population <= 0.0 || s.standard_weight <= 0.0 {
            continue;
        }
        total_cases += s.cases.max(0.0);
        total_pop += s.population;
        weight_sum += s.standard_weight;

        let w_over_n = s.standard_weight / s.population;
        max_w_over_n = max_w_over_n.max(w_over_n);
        adj_rate += w_over_n * s.cases.max(0.0);
        var_rate += w_over_n * w_over_n * s.cases.max(0.0);
    }

    if total_pop <= 0.0 || weight_sum <= 0.0 {
        return None;
    }

    let norm_rate = adj_rate / weight_sum;
    let norm_var = var_rate / (weight_sum * weight_sum);
    let norm_max_w = max_w_over_n / weight_sum;

    // Fay & Feuer (1997) Gamma interval for directly standardized rates
    let (ci_low, ci_high) = if norm_rate > 0.0 && norm_var > 0.0 {
        let a_low = (norm_rate * norm_rate) / norm_var;
        let b_low = norm_rate / norm_var;
        let low = gamma_quantile(0.025, a_low, b_low);

        let r_plus = norm_rate + norm_max_w;
        let v_plus = norm_var + norm_max_w * norm_max_w;
        let a_high = (r_plus * r_plus) / v_plus;
        let b_high = r_plus / v_plus;
        let high = gamma_quantile(0.975, a_high, b_high);
        (low * 100_000.0, high * 100_000.0)
    } else {
        (0.0, 0.0)
    };

    Some(AgeStandardizedResult {
        crude_rate_per_100k: (total_cases / total_pop) * 100_000.0,
        age_adjusted_rate_per_100k: norm_rate * 100_000.0,
        ci_low_95: ci_low,
        ci_high_95: ci_high,
    })
}

/// Evaluate Spearman cross-correlation across `lag in -max_lag..=max_lag` between a leading signal
/// (such as NWSS wastewater percentile) and clinical notifications (`lag > 0` means wastewater leads clinical).
pub fn compute_wastewater_lead_lag(
    wastewater: &[f64],
    clinical: &[f64],
    max_lag: i32,
) -> Vec<LagCorrelationPoint> {
    let n = wastewater.len().min(clinical.len());
    let mut out = Vec::new();
    if n < 5 {
        return out;
    }

    for lag in -max_lag..=max_lag {
        let mut w_slice = Vec::new();
        let mut c_slice = Vec::new();
        for t in 0..n {
            let c_idx = (t as i32) + lag;
            if c_idx >= 0 && (c_idx as usize) < n {
                w_slice.push(wastewater[t]);
                c_slice.push(clinical[c_idx as usize]);
            }
        }
        if let Some(rho) = spearman_rho(&w_slice, &c_slice) {
            out.push(LagCorrelationPoint {
                lag_steps: lag,
                spearman_rho: rho,
                paired_n: w_slice.len(),
            });
        }
    }
    out
}

/// Rank jurisdictions by Kulldorff spatial Poisson excess Log-Likelihood Ratio (`LLR`):
/// `LLR = c * ln(c / E) + (C - c) * ln((C - c) / (C - E))` when `c > E`, else `0.0`.
pub fn compute_kulldorff_spatial_clusters(
    records: &[JurisdictionRecord],
    use_2025: bool,
) -> Vec<SpatialClusterScore> {
    let total_pop: f64 = records.iter().map(|r| r.census_pop_2025 as f64).sum();
    let total_cases: f64 = records
        .iter()
        .map(|r| {
            if let Some(custom) = r.active_custom_value {
                custom.max(0.0)
            } else if use_2025 {
                f64::from(r.cases_2025)
            } else {
                f64::from(r.cases_2024)
            }
        })
        .sum();

    if total_pop <= 0.0 || total_cases <= 0.0 {
        return Vec::new();
    }

    let mut scores: Vec<SpatialClusterScore> = records
        .iter()
        .map(|r| {
            let obs = if let Some(custom) = r.active_custom_value {
                custom.max(0.0)
            } else if use_2025 {
                f64::from(r.cases_2025)
            } else {
                f64::from(r.cases_2024)
            };
            let exp = total_cases * ((r.census_pop_2025 as f64) / total_pop);
            let rr = if exp > 0.0 { obs / exp } else { 0.0 };

            let llr = if obs > exp && exp > 0.0 && (total_cases - obs) > 0.0 && (total_cases - exp) > 0.0 {
                obs * (obs / exp).ln()
                    + (total_cases - obs) * ((total_cases - obs) / (total_cases - exp)).ln()
            } else {
                0.0
            };

            SpatialClusterScore {
                fips: r.fips.clone(),
                name: r.name.clone(),
                observed_cases: obs,
                expected_cases: exp,
                relative_risk: rr,
                kulldorff_llr: llr,
            }
        })
        .collect();

    scores.sort_by(|a, b| {
        b.kulldorff_llr
            .partial_cmp(&a.kulldorff_llr)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    scores
}
