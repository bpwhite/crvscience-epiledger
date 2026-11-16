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
use serde::{Deserialize, Serialize};

/// Pharmacovigilance 2x2 contingency disproportionality analyzer for openFDA FAERS reports.
///
/// Computes three standard regulatory signal-detection metrics (PMC12443087, PMC3970098):
/// 1. **PRR (Proportional Reporting Ratio)** + Yates chi-square (`PRR >= 2.0`, `chi2 >= 4.0`, `a >= 3`).
/// 2. **ROR (Reporting Odds Ratio)** + Woolf exact 95% confidence interval (`ror_ci_low_95 > 1.0`).
/// 3. **BCPNN Information Component (`IC` & `IC025`)** with Bayesian Gamma-Poisson shrinkage (`ic_025 > 0.0`).

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DisproportionalitySignal {
    pub drug_label: String,
    pub reaction_term: String,
    /// Target drug + target reaction co-occurrences (`a`)
    pub a_target_both: u64,
    /// Target drug + all other reactions (`b`)
    pub b_target_other_rx: u64,
    /// Comparator drugs + target reaction (`c`)
    pub c_other_drug_rx: u64,
    /// Comparator drugs + all other reactions (`d`)
    pub d_background: u64,
    pub prr: f64,
    pub prr_chi2_yates: f64,
    pub ror: f64,
    pub ror_ci_low_95: f64,
    pub ror_ci_high_95: f64,
    pub bcpnn_ic: f64,
    pub bcpnn_ic_025: f64,
    pub flagged_evans_prr: bool,
    pub flagged_bcpnn_ic: bool,
}

pub fn compute_disproportionality_2x2(
    drug_label: &str,
    reaction_term: &str,
    a_u64: u64,
    b_u64: u64,
    c_u64: u64,
    d_u64: u64,
) -> DisproportionalitySignal {
    // Apply Haldane-Anscombe +0.5 continuity guard only when a zero cell would divide by zero
    let a = a_u64 as f64;
    let b = b_u64 as f64;
    let c = c_u64 as f64;
    let d = d_u64 as f64;

    let a_safe = if a == 0.0 { 0.5 } else { a };
    let b_safe = if b == 0.0 { 0.5 } else { b };
    let c_safe = if c == 0.0 { 0.5 } else { c };
    let d_safe = if d == 0.0 { 0.5 } else { d };

    // 1. PRR = [a / (a + b)] / [c / (c + d)]
    let p1 = a_safe / (a_safe + b_safe);
    let p2 = c_safe / (c_safe + d_safe);
    let prr = p1 / p2.max(1e-12);

    let total_n = (a + b + c + d).max(1.0);
    let denom_chi = (a + b) * (c + d) * (a + c) * (b + d);
    let prr_chi2_yates = if denom_chi > 0.0 {
        let diff = ((a * d - b * c).abs() - 0.5 * total_n).max(0.0);
        (total_n * diff * diff) / denom_chi
    } else {
        0.0
    };

    // 2. ROR = (a * d) / (b * c) with Woolf log-scale 95% CI
    let ror = (a_safe * d_safe) / (b_safe * c_safe);
    let se_ln_ror = (1.0 / a_safe + 1.0 / b_safe + 1.0 / c_safe + 1.0 / d_safe).sqrt();
    let z = normal_quantile(0.975);
    let ror_ci_low_95 = (ror.ln() - z * se_ln_ror).exp();
    let ror_ci_high_95 = (ror.ln() + z * se_ln_ror).exp();

    // 3. BCPNN Information Component (Norén et al. 2013 shrinkage: IC = log2((a + 0.5) / (E + 0.5)))
    let expected_a = ((a + b) * (a + c)) / total_n;
    let bcpnn_ic = ((a + 0.5) / (expected_a + 0.5)).log2();
    // Posterior variance on log2 scale approximated via digamma/trigamma moment match:
    // SD(IC) ~= 1 / (ln(2) * sqrt(a + 0.5))
    let sd_ic = 1.0 / (std::f64::consts::LN_2 * (a + 0.5).sqrt());
    let bcpnn_ic_025 = bcpnn_ic - z * sd_ic;

    let flagged_evans_prr = a_u64 >= 3 && prr >= 2.0 && prr_chi2_yates >= 4.0;
    let flagged_bcpnn_ic = a_u64 >= 3 && bcpnn_ic_025 > 0.0;

    DisproportionalitySignal {
        drug_label: drug_label.to_string(),
        reaction_term: reaction_term.to_string(),
        a_target_both: a_u64,
        b_target_other_rx: b_u64,
        c_other_drug_rx: c_u64,
        d_background: d_u64,
        prr,
        prr_chi2_yates,
        ror,
        ror_ci_low_95,
        ror_ci_high_95,
        bcpnn_ic,
        bcpnn_ic_025,
        flagged_evans_prr,
        flagged_bcpnn_ic,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_faers_disproportionality_signal_detection() {
        let sig = compute_disproportionality_2x2(
            "Oseltamivir",
            "HALLUCINATION",
            180,
            4_820,
            1_200,
            993_800,
        );
        assert!(sig.prr > 25.0);
        assert!(sig.ror_ci_low_95 > 20.0);
        assert!(sig.bcpnn_ic_025 > 3.5);
        assert!(sig.flagged_evans_prr && sig.flagged_bcpnn_ic);
    }
}
