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

//! Multi-jurisdiction differential-privacy `R_t` & competing-risk hospital surge
//! comparative risk digest (`JSON` / `CSV` export).

use crate::cfr_and_surge::{
    project_competing_risk_hospital_surge, HospitalSurgeConfig,
};
use crate::renewal_rt::{compute_cori_rt, CoriRtConfig};
use epi_core::{EpiError, SurveillanceCampaignPreset};
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;

/// Input surveillance time series and hospital capacity metadata for one jurisdiction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JurisdictionSeriesInput {
    pub jurisdiction_id: String,
    pub label: String,
    pub population: f64,
    pub staffed_icu_beds: f64,
    pub weekly_counts: Vec<f64>,
}

/// Configuration for multi-jurisdiction epsilon-DP `R_t` & surge comparative risk synthesis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DpRiskDigestConfig {
    /// Per-jurisdiction differential privacy budget (`> 0.005`).
    pub epsilon_per_jurisdiction: f64,
    /// L1 sensitivity of weekly count queries (`>= 0.01`).
    pub sensitivity_l1: f64,
    /// Optional deterministic SplitMix64 seed for reproducible DP noise; `None` draws OS time entropy.
    pub deterministic_seed: Option<u64>,
    /// Cori et al. (2013) renewal-equation `R_t` configuration.
    pub cori_cfg: CoriRtConfig,
    /// Competing-risk hospital/ICU surge configuration (`staffed_icu_beds` overridden per jurisdiction when positive).
    pub surge_cfg: HospitalSurgeConfig,
    /// Small-cell suppression threshold (`9.0` for NCHS `<= 9` flagging).
    pub suppression_threshold: f64,
}

impl Default for DpRiskDigestConfig {
    fn default() -> Self {
        Self {
            epsilon_per_jurisdiction: 0.75,
            sensitivity_l1: 1.0,
            deterministic_seed: Some(2026_1127),
            cori_cfg: CoriRtConfig::default(),
            surge_cfg: HospitalSurgeConfig::default(),
            suppression_threshold: 9.0,
        }
    }
}

/// Comparative multi-jurisdiction alert tier derived from DP-perturbed `R_t` and ICU strain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DpRiskTier {
    Low,
    Moderate,
    Elevated,
    Critical,
}

impl DpRiskTier {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Low => "LOW",
            Self::Moderate => "MODERATE",
            Self::Elevated => "ELEVATED",
            Self::Critical => "CRITICAL",
        }
    }
}

/// Differentially private `R_t` & competing-risk surge summary for a single jurisdiction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JurisdictionDpRiskSummary {
    pub rank: usize,
    pub jurisdiction_id: String,
    pub label: String,
    pub population: f64,
    pub staffed_icu_beds: f64,
    pub weeks_observed: usize,
    pub suppressed_weeks_raw: usize,
    pub epsilon_spent: f64,
    pub dp_empirical_mae: f64,
    pub dp_theoretical_rmse: f64,
    pub dp_total_cases: f64,
    pub dp_recent_4wk_cases: f64,
    pub dp_recent_4wk_per_100k: f64,
    pub dp_growth_ratio_4wk: f64,
    pub dp_rt_latest: Option<f64>,
    pub dp_rt_cri_low_95: Option<f64>,
    pub dp_rt_cri_high_95: Option<f64>,
    pub supercritical_rt_95: bool,
    pub dp_peak_ward_census: f64,
    pub dp_peak_icu_census: f64,
    pub dp_peak_icu_utilization_pct: f64,
    pub dp_icu_bed_days_deficit: f64,
    pub dp_projected_hospital_deaths: f64,
    pub comparative_risk_score: f64,
    pub risk_tier: DpRiskTier,
}

/// Cross-jurisdiction comparative DP risk digest and privacy-budget ledger summary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiJurisdictionDpRiskDigest {
    pub jurisdiction_count: usize,
    pub epsilon_per_jurisdiction: f64,
    pub total_epsilon_spent: f64,
    pub sensitivity_l1: f64,
    pub deterministic_seed: Option<u64>,
    pub highest_risk_jurisdiction: Option<String>,
    pub mean_dp_rt: Option<f64>,
    pub max_icu_utilization_pct: f64,
    pub elevated_or_critical_count: usize,
    pub jurisdictions: Vec<JurisdictionDpRiskSummary>,
}

#[inline]
fn splitmix64_next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn fnv1a_64(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

/// Exact Two-Sided Geometric (`Discrete Laplace`) draw with parameter `alpha = exp(-eps / sens)`.
fn sample_discrete_laplace(state: &mut u64, alpha: f64) -> i64 {
    if !alpha.is_finite() || alpha <= 1e-15 {
        return 0;
    }
    let u_raw = splitmix64_next(state);
    let sign_bit = splitmix64_next(state) & 1;
    let u = (((u_raw >> 11) as f64) + 0.5) / ((1u64 << 53) as f64);
    let p0 = (1.0 - alpha) / (1.0 + alpha);
    if u <= p0 {
        return 0;
    }
    let sign = if sign_bit == 0 { 1i64 } else { -1i64 };
    let u_tail = ((u - p0) / (1.0 - p0)).clamp(1e-15, 1.0 - 1e-15);
    let mag = 1 + (u_tail.ln() / alpha.ln()).floor().max(0.0) as i64;
    sign * mag
}

fn json_escape(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len() + 8);
    for ch in raw.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 || c.is_control() => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out
}

fn csv_escape(raw: &str) -> String {
    if raw.contains(',') || raw.contains('"') || raw.contains('\n') || raw.contains('\r') {
        format!("\"{}\"", raw.replace('"', "\"\""))
    } else {
        raw.to_string()
    }
}

fn fmt_opt_f64_json(v: Option<f64>, decimals: usize) -> String {
    match v {
        Some(x) if x.is_finite() => format!("{x:.decimals$}"),
        _ => "null".to_string(),
    }
}

fn fmt_opt_f64_csv(v: Option<f64>, decimals: usize) -> String {
    match v {
        Some(x) if x.is_finite() => format!("{x:.decimals$}"),
        _ => String::new(),
    }
}

/// Build an 8-jurisdiction benchmark cohort from a built-in `SurveillanceCampaignPreset`.
pub fn builtin_multi_jurisdiction_dp_inputs_for_campaign(
    preset: SurveillanceCampaignPreset,
) -> Vec<JurisdictionSeriesInput> {
    let base_ds = preset.build_weekly_dataset();
    let base_counts: Vec<f64> = base_ds.points.iter().map(|p| p.count.max(0.0)).collect();
    let base_sum: f64 = base_counts.iter().sum::<f64>().max(1.0);

    let mut jurisdictions = preset.build_jurisdictions();
    jurisdictions.sort_by(|a, b| {
        b.rate_2025
            .partial_cmp(&a.rate_2025)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.abbr.cmp(&b.abbr))
    });

    let selected = jurisdictions.into_iter().take(8).collect::<Vec<_>>();
    selected
        .into_iter()
        .enumerate()
        .map(|(idx, j)| {
            let pop = (j.census_pop_2025.max(100_000)) as f64;
            let annual_target = (j.cases_2025 as f64).max(40.0);
            let scale = (annual_target * 3.2) / base_sum;
            let phase_shift = idx % 3;
            let weekly_counts: Vec<f64> = base_counts
                .iter()
                .enumerate()
                .map(|(w, &c)| {
                    let mod_wave = 1.0 + 0.14 * ((w + phase_shift * 5) as f64 * 0.25).sin();
                    (c * scale * mod_wave).round().max(0.0)
                })
                .collect();
            let staffed_icu_beds = ((pop / 100_000.0) * 14.5).clamp(24.0, 2_800.0).round();

            JurisdictionSeriesInput {
                jurisdiction_id: format!("US-{}", j.abbr),
                label: j.name,
                population: pop,
                staffed_icu_beds,
                weekly_counts,
            }
        })
        .collect()
}

/// Default 8-jurisdiction benchmark cohort (US Measles 2022-2026 focal jurisdictions).
pub fn builtin_multi_jurisdiction_dp_inputs() -> Vec<JurisdictionSeriesInput> {
    builtin_multi_jurisdiction_dp_inputs_for_campaign(SurveillanceCampaignPreset::UsMeasles2022_2026)
}

/// Compute a multi-jurisdiction epsilon-DP `R_t` & competing-risk surge comparative risk digest.
pub fn compute_multi_jurisdiction_dp_risk_digest(
    inputs: &[JurisdictionSeriesInput],
    cfg: &DpRiskDigestConfig,
) -> Result<MultiJurisdictionDpRiskDigest, EpiError> {
    if inputs.is_empty() {
        return Err(EpiError::Schema(
            "Multi-jurisdiction DP risk digest requires at least 1 jurisdiction".to_string(),
        ));
    }
    let eps = cfg.epsilon_per_jurisdiction;
    if !eps.is_finite() || eps <= 0.005 {
        return Err(EpiError::Schema(
            "Differential privacy epsilon_per_jurisdiction must be finite and > 0.005".to_string(),
        ));
    }
    let sens = cfg.sensitivity_l1;
    if !sens.is_finite() || sens <= 0.0 {
        return Err(EpiError::Schema(
            "Differential privacy sensitivity_l1 must be finite and > 0.0".to_string(),
        ));
    }

    let alpha = (-eps / sens.max(0.01)).exp().clamp(0.0, 0.999_999_9);
    let theoretical_rmse = if alpha > 0.0 {
        (2.0 * alpha).sqrt() / (1.0 - alpha).max(1e-9)
    } else {
        0.0
    };

    let base_seed = cfg.deterministic_seed.unwrap_or_else(|| {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(0x517C_C1B7_2722_0A95)
    });

    let mut summaries = Vec::with_capacity(inputs.len());
    for (j_idx, inp) in inputs.iter().enumerate() {
        if inp.weekly_counts.is_empty() {
            return Err(EpiError::Schema(format!(
                "Jurisdiction '{}' has zero weekly observations",
                inp.jurisdiction_id
            )));
        }

        let mut rng_state = base_seed
            ^ fnv1a_64(inp.jurisdiction_id.as_bytes())
            ^ ((j_idx as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));

        let mut dp_weekly = Vec::with_capacity(inp.weekly_counts.len());
        let mut abs_err_sum = 0.0_f64;
        let mut suppressed_raw = 0usize;

        for &raw in &inp.weekly_counts {
            let clean = if raw.is_finite() { raw.max(0.0) } else { 0.0 };
            if clean > 0.0 && clean <= cfg.suppression_threshold {
                suppressed_raw += 1;
            }
            let noise = sample_discrete_laplace(&mut rng_state, alpha) as f64;
            let noisy = (clean + noise).max(0.0);
            abs_err_sum += (noisy - clean).abs();
            dp_weekly.push(noisy);
        }

        let n_wks = dp_weekly.len();
        let empirical_mae = abs_err_sum / (n_wks as f64);
        let dp_total_cases: f64 = dp_weekly.iter().sum();

        let recent_start = n_wks.saturating_sub(4);
        let prior_start = n_wks.saturating_sub(8);
        let dp_recent_4wk_cases: f64 = dp_weekly[recent_start..].iter().sum();
        let dp_prior_4wk_cases: f64 = if prior_start < recent_start {
            dp_weekly[prior_start..recent_start].iter().sum()
        } else {
            dp_recent_4wk_cases
        };

        let dp_growth_ratio_4wk = if dp_prior_4wk_cases > 1e-6 {
            (dp_recent_4wk_cases / dp_prior_4wk_cases).clamp(0.0, 20.0)
        } else if dp_recent_4wk_cases > 0.0 {
            2.0
        } else {
            1.0
        };

        let pop_safe = if inp.population.is_finite() && inp.population > 0.0 {
            inp.population
        } else {
            100_000.0
        };
        let dp_recent_4wk_per_100k = (dp_recent_4wk_cases / pop_safe) * 100_000.0;

        let rt_series = compute_cori_rt(&dp_weekly, &cfg.cori_cfg);
        let latest_rt = rt_series.estimates.iter().rev().flatten().next().copied();
        let (dp_rt_latest, dp_rt_cri_low_95, dp_rt_cri_high_95) = match latest_rt {
            Some(est) => (Some(est.mean), Some(est.cri_low_95), Some(est.cri_high_95)),
            None => (None, None, None),
        };
        let supercritical_rt_95 = dp_rt_cri_low_95.map(|lo| lo > 1.0).unwrap_or(false);

        // Disaggregate the most recent 12 weeks (or available tail) into daily cases for surge convolution
        let tail_wks = n_wks.min(12);
        let tail_slice = &dp_weekly[n_wks - tail_wks..];
        let mut daily_cases = Vec::with_capacity(tail_wks * 7);
        let mut daily_low = Vec::with_capacity(tail_wks * 7);
        let mut daily_high = Vec::with_capacity(tail_wks * 7);
        let noise_band = (1.96 * theoretical_rmse / 7.0).max(0.0);
        for &wk in tail_slice {
            let d = (wk / 7.0).max(0.0);
            for _ in 0..7 {
                daily_cases.push(d);
                daily_low.push((d - noise_band).max(0.0));
                daily_high.push(d + noise_band);
            }
        }

        let mut j_surge_cfg = cfg.surge_cfg.clone();
        if inp.staffed_icu_beds.is_finite() && inp.staffed_icu_beds > 0.0 {
            j_surge_cfg.staffed_icu_beds = inp.staffed_icu_beds;
            j_surge_cfg.staffed_ward_beds = (inp.staffed_icu_beds * 5.2).max(20.0);
            j_surge_cfg.available_ventilators = (inp.staffed_icu_beds * 0.58).max(4.0);
        }

        let surge = project_competing_risk_hospital_surge(
            &daily_cases,
            &daily_low,
            &daily_high,
            &j_surge_cfg,
        )?;

        let rt_component = dp_rt_latest.unwrap_or(1.0).clamp(0.0, 3.5) * 14.0;
        let supercrit_bonus = if supercritical_rt_95 { 12.0 } else { 0.0 };
        let growth_component = (dp_growth_ratio_4wk - 1.0).clamp(0.0, 3.0) * 8.0;
        let icu_component = surge.peak_icu_utilization_pct.clamp(0.0, 150.0) * 0.28;
        let inc_component = dp_recent_4wk_per_100k.clamp(0.0, 80.0) * 0.12;
        let score = (rt_component + supercrit_bonus + growth_component + icu_component + inc_component)
            .clamp(0.0, 100.0);

        let risk_tier = if score >= 70.0 || surge.peak_icu_utilization_pct >= 90.0 {
            DpRiskTier::Critical
        } else if score >= 45.0 || supercritical_rt_95 || surge.peak_icu_utilization_pct >= 65.0 {
            DpRiskTier::Elevated
        } else if score >= 25.0 || dp_rt_latest.unwrap_or(0.0) >= 1.05 {
            DpRiskTier::Moderate
        } else {
            DpRiskTier::Low
        };

        summaries.push(JurisdictionDpRiskSummary {
            rank: 0,
            jurisdiction_id: inp.jurisdiction_id.clone(),
            label: inp.label.clone(),
            population: pop_safe,
            staffed_icu_beds: j_surge_cfg.staffed_icu_beds,
            weeks_observed: n_wks,
            suppressed_weeks_raw: suppressed_raw,
            epsilon_spent: eps,
            dp_empirical_mae: empirical_mae,
            dp_theoretical_rmse: theoretical_rmse,
            dp_total_cases,
            dp_recent_4wk_cases,
            dp_recent_4wk_per_100k,
            dp_growth_ratio_4wk,
            dp_rt_latest,
            dp_rt_cri_low_95,
            dp_rt_cri_high_95,
            supercritical_rt_95,
            dp_peak_ward_census: surge.peak_ward_census,
            dp_peak_icu_census: surge.peak_icu_census,
            dp_peak_icu_utilization_pct: surge.peak_icu_utilization_pct,
            dp_icu_bed_days_deficit: surge.icu_bed_days_deficit,
            dp_projected_hospital_deaths: surge.total_projected_hospital_deaths,
            comparative_risk_score: score,
            risk_tier,
        });
    }

    summaries.sort_by(|a, b| {
        b.comparative_risk_score
            .partial_cmp(&a.comparative_risk_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.dp_peak_icu_utilization_pct
                    .partial_cmp(&a.dp_peak_icu_utilization_pct)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.jurisdiction_id.cmp(&b.jurisdiction_id))
    });

    for (idx, row) in summaries.iter_mut().enumerate() {
        row.rank = idx + 1;
    }

    let rt_vals: Vec<f64> = summaries.iter().filter_map(|s| s.dp_rt_latest).collect();
    let mean_dp_rt = if rt_vals.is_empty() {
        None
    } else {
        Some(rt_vals.iter().sum::<f64>() / (rt_vals.len() as f64))
    };
    let max_icu_utilization_pct = summaries
        .iter()
        .map(|s| s.dp_peak_icu_utilization_pct)
        .fold(0.0_f64, f64::max);
    let elevated_or_critical_count = summaries
        .iter()
        .filter(|s| matches!(s.risk_tier, DpRiskTier::Elevated | DpRiskTier::Critical))
        .count();
    let highest_risk_jurisdiction = summaries.first().map(|s| s.jurisdiction_id.clone());
    let total_epsilon_spent = eps * (summaries.len() as f64);

    Ok(MultiJurisdictionDpRiskDigest {
        jurisdiction_count: summaries.len(),
        epsilon_per_jurisdiction: eps,
        total_epsilon_spent,
        sensitivity_l1: sens,
        deterministic_seed: cfg.deterministic_seed,
        highest_risk_jurisdiction,
        mean_dp_rt,
        max_icu_utilization_pct,
        elevated_or_critical_count,
        jurisdictions: summaries,
    })
}

/// Export a `MultiJurisdictionDpRiskDigest` to RFC 8259 JSON.
pub fn export_dp_risk_digest_json(digest: &MultiJurisdictionDpRiskDigest) -> String {
    let mut out = String::with_capacity(2048);
    out.push_str("{\n");
    let _ = writeln!(out, "  \"jurisdiction_count\": {},", digest.jurisdiction_count);
    let _ = writeln!(
        out,
        "  \"epsilon_per_jurisdiction\": {:.4},",
        digest.epsilon_per_jurisdiction
    );
    let _ = writeln!(out, "  \"total_epsilon_spent\": {:.4},", digest.total_epsilon_spent);
    let _ = writeln!(out, "  \"sensitivity_l1\": {:.4},", digest.sensitivity_l1);
    let seed_str = digest
        .deterministic_seed
        .map(|s| s.to_string())
        .unwrap_or_else(|| "null".to_string());
    let _ = writeln!(out, "  \"deterministic_seed\": {seed_str},");
    let highest_str = digest
        .highest_risk_jurisdiction
        .as_deref()
        .map(|s| format!("\"{}\"", json_escape(s)))
        .unwrap_or_else(|| "null".to_string());
    let _ = writeln!(out, "  \"highest_risk_jurisdiction\": {highest_str},");
    let _ = writeln!(
        out,
        "  \"mean_dp_rt\": {},",
        fmt_opt_f64_json(digest.mean_dp_rt, 4)
    );
    let _ = writeln!(
        out,
        "  \"max_icu_utilization_pct\": {:.2},",
        digest.max_icu_utilization_pct
    );
    let _ = writeln!(
        out,
        "  \"elevated_or_critical_count\": {},",
        digest.elevated_or_critical_count
    );
    out.push_str("  \"jurisdictions\": [\n");

    for (idx, j) in digest.jurisdictions.iter().enumerate() {
        let comma = if idx + 1 < digest.jurisdictions.len() {
            ","
        } else {
            ""
        };
        out.push_str("    {\n");
        let _ = writeln!(out, "      \"rank\": {},", j.rank);
        let _ = writeln!(
            out,
            "      \"jurisdiction_id\": \"{}\",",
            json_escape(&j.jurisdiction_id)
        );
        let _ = writeln!(out, "      \"label\": \"{}\",", json_escape(&j.label));
        let _ = writeln!(out, "      \"population\": {:.0},", j.population);
        let _ = writeln!(out, "      \"staffed_icu_beds\": {:.1},", j.staffed_icu_beds);
        let _ = writeln!(out, "      \"weeks_observed\": {},", j.weeks_observed);
        let _ = writeln!(
            out,
            "      \"suppressed_weeks_raw\": {},",
            j.suppressed_weeks_raw
        );
        let _ = writeln!(out, "      \"epsilon_spent\": {:.4},", j.epsilon_spent);
        let _ = writeln!(out, "      \"dp_empirical_mae\": {:.4},", j.dp_empirical_mae);
        let _ = writeln!(
            out,
            "      \"dp_theoretical_rmse\": {:.4},",
            j.dp_theoretical_rmse
        );
        let _ = writeln!(out, "      \"dp_total_cases\": {:.1},", j.dp_total_cases);
        let _ = writeln!(
            out,
            "      \"dp_recent_4wk_cases\": {:.1},",
            j.dp_recent_4wk_cases
        );
        let _ = writeln!(
            out,
            "      \"dp_recent_4wk_per_100k\": {:.3},",
            j.dp_recent_4wk_per_100k
        );
        let _ = writeln!(
            out,
            "      \"dp_growth_ratio_4wk\": {:.4},",
            j.dp_growth_ratio_4wk
        );
        let _ = writeln!(
            out,
            "      \"dp_rt_latest\": {},",
            fmt_opt_f64_json(j.dp_rt_latest, 4)
        );
        let _ = writeln!(
            out,
            "      \"dp_rt_cri_low_95\": {},",
            fmt_opt_f64_json(j.dp_rt_cri_low_95, 4)
        );
        let _ = writeln!(
            out,
            "      \"dp_rt_cri_high_95\": {},",
            fmt_opt_f64_json(j.dp_rt_cri_high_95, 4)
        );
        let _ = writeln!(
            out,
            "      \"supercritical_rt_95\": {},",
            j.supercritical_rt_95
        );
        let _ = writeln!(
            out,
            "      \"dp_peak_ward_census\": {:.2},",
            j.dp_peak_ward_census
        );
        let _ = writeln!(
            out,
            "      \"dp_peak_icu_census\": {:.2},",
            j.dp_peak_icu_census
        );
        let _ = writeln!(
            out,
            "      \"dp_peak_icu_utilization_pct\": {:.2},",
            j.dp_peak_icu_utilization_pct
        );
        let _ = writeln!(
            out,
            "      \"dp_icu_bed_days_deficit\": {:.2},",
            j.dp_icu_bed_days_deficit
        );
        let _ = writeln!(
            out,
            "      \"dp_projected_hospital_deaths\": {:.2},",
            j.dp_projected_hospital_deaths
        );
        let _ = writeln!(
            out,
            "      \"comparative_risk_score\": {:.2},",
            j.comparative_risk_score
        );
        let _ = writeln!(out, "      \"risk_tier\": \"{}\"", j.risk_tier.as_str());
        let _ = writeln!(out, "    }}{comma}");
    }

    out.push_str("  ]\n}\n");
    out
}

/// Export a `MultiJurisdictionDpRiskDigest` to RFC 4180 CSV.
pub fn export_dp_risk_digest_csv(digest: &MultiJurisdictionDpRiskDigest) -> String {
    let mut out = String::with_capacity(1024);
    out.push_str("rank,jurisdiction_id,label,population,staffed_icu_beds,weeks_observed,suppressed_weeks_raw,epsilon_spent,dp_empirical_mae,dp_theoretical_rmse,dp_total_cases,dp_recent_4wk_cases,dp_recent_4wk_per_100k,dp_growth_ratio_4wk,dp_rt_latest,dp_rt_cri_low_95,dp_rt_cri_high_95,supercritical_rt_95,dp_peak_ward_census,dp_peak_icu_census,dp_peak_icu_utilization_pct,dp_icu_bed_days_deficit,dp_projected_hospital_deaths,comparative_risk_score,risk_tier\n");
    for j in &digest.jurisdictions {
        let _ = writeln!(
            out,
            "{},{},{},{:.0},{:.1},{},{},{:.4},{:.4},{:.4},{:.1},{:.1},{:.3},{:.4},{},{},{},{},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{}",
            j.rank,
            csv_escape(&j.jurisdiction_id),
            csv_escape(&j.label),
            j.population,
            j.staffed_icu_beds,
            j.weeks_observed,
            j.suppressed_weeks_raw,
            j.epsilon_spent,
            j.dp_empirical_mae,
            j.dp_theoretical_rmse,
            j.dp_total_cases,
            j.dp_recent_4wk_cases,
            j.dp_recent_4wk_per_100k,
            j.dp_growth_ratio_4wk,
            fmt_opt_f64_csv(j.dp_rt_latest, 4),
            fmt_opt_f64_csv(j.dp_rt_cri_low_95, 4),
            fmt_opt_f64_csv(j.dp_rt_cri_high_95, 4),
            j.supercritical_rt_95,
            j.dp_peak_ward_census,
            j.dp_peak_icu_census,
            j.dp_peak_icu_utilization_pct,
            j.dp_icu_bed_days_deficit,
            j.dp_projected_hospital_deaths,
            j.comparative_risk_score,
            j.risk_tier.as_str(),
        );
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_builtin_multi_jurisdiction_dp_risk_digest_ranks_and_exports() {
        let inputs = builtin_multi_jurisdiction_dp_inputs();
        assert_eq!(inputs.len(), 8);
        let cfg = DpRiskDigestConfig::default();
        let digest = compute_multi_jurisdiction_dp_risk_digest(&inputs, &cfg)
            .expect("digest computation succeeds");

        assert_eq!(digest.jurisdiction_count, 8);
        assert!((digest.total_epsilon_spent - 6.0).abs() < 1e-9);
        assert!(digest.highest_risk_jurisdiction.is_some());

        for (i, j) in digest.jurisdictions.iter().enumerate() {
            assert_eq!(j.rank, i + 1);
            assert!(j.comparative_risk_score.is_finite());
            if i > 0 {
                assert!(
                    digest.jurisdictions[i - 1].comparative_risk_score >= j.comparative_risk_score
                        - 1e-9
                );
            }
        }

        let json = export_dp_risk_digest_json(&digest);
        let parsed: serde_json::Value = serde_json::from_str(&json).expect("valid RFC 8259 JSON");
        assert_eq!(parsed["jurisdiction_count"], 8);
        assert_eq!(
            parsed["jurisdictions"].as_array().map(|a| a.len()),
            Some(8)
        );

        let csv = export_dp_risk_digest_csv(&digest);
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 9, "header + 8 jurisdiction rows");
    }

    #[test]
    fn test_dp_risk_digest_deterministic_seed_reproducibility() {
        let inputs = builtin_multi_jurisdiction_dp_inputs();
        let cfg_a = DpRiskDigestConfig {
            deterministic_seed: Some(42),
            ..Default::default()
        };
        let cfg_b = DpRiskDigestConfig {
            deterministic_seed: Some(99),
            ..Default::default()
        };
        let d1 = compute_multi_jurisdiction_dp_risk_digest(&inputs, &cfg_a).expect("d1");
        let d2 = compute_multi_jurisdiction_dp_risk_digest(&inputs, &cfg_a).expect("d2");
        let d3 = compute_multi_jurisdiction_dp_risk_digest(&inputs, &cfg_b).expect("d3");

        assert_eq!(export_dp_risk_digest_csv(&d1), export_dp_risk_digest_csv(&d2));
        assert_ne!(
            export_dp_risk_digest_csv(&d1),
            export_dp_risk_digest_csv(&d3),
            "distinct DP seeds should perturb counts differently"
        );
    }

    #[test]
    fn test_dp_risk_digest_edge_cases_and_validation() {
        let cfg = DpRiskDigestConfig::default();
        assert!(compute_multi_jurisdiction_dp_risk_digest(&[], &cfg).is_err());

        let bad_eps = DpRiskDigestConfig {
            epsilon_per_jurisdiction: 0.001,
            ..Default::default()
        };
        let single = vec![JurisdictionSeriesInput {
            jurisdiction_id: "TEST-ZERO \"A\"".to_string(),
            label: "Zero, Tiny \"County\"".to_string(),
            population: 0.0,
            staffed_icu_beds: 0.0,
            weekly_counts: vec![0.0, 2.0, f64::NAN, 5.0, 0.0],
        }];
        assert!(compute_multi_jurisdiction_dp_risk_digest(&single, &bad_eps).is_err());

        let digest = compute_multi_jurisdiction_dp_risk_digest(&single, &cfg)
            .expect("zero/low-count/NaN edge case handled safely");
        assert_eq!(digest.jurisdiction_count, 1);
        assert_eq!(digest.jurisdictions[0].suppressed_weeks_raw, 2);
        assert!(digest.jurisdictions[0].dp_rt_latest.is_none());

        let json = export_dp_risk_digest_json(&digest);
        let parsed: serde_json::Value = serde_json::from_str(&json).expect("escaped JSON valid");
        assert_eq!(parsed["jurisdictions"][0]["jurisdiction_id"], "TEST-ZERO \"A\"");

        let csv = export_dp_risk_digest_csv(&digest);
        assert!(csv.contains("\"Zero, Tiny \"\"County\"\"\""));
    }
}
