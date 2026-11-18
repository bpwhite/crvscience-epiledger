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

use crate::bitemporal::EpiStore;
use crate::pprl_bloom::hmac_sha256;
use epi_core::{EpiError, SourceKind, TimeSeriesDataset, TimeSeriesPoint};
use rusqlite::{params, Connection, OptionalExtension};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::time::{SystemTime, UNIX_EPOCH};

/// Discrete Laplace (Two-Sided Geometric) `epsilon`-Differential Privacy Mechanism,
/// Per-Vintage SQLite Privacy-Budget Ledger, and Tamper-Evident `SHA-256` Hash-Chained Audit Log.
///
/// References:
/// - Differential privacy for low-count public health surveillance data: https://pmc.ncbi.nlm.nih.gov/articles/PMC8662814/
/// - Dynamically adjusting case reporting policy to balance privacy and utility: https://pmc.ncbi.nlm.nih.gov/articles/PMC9006705/

pub const GENESIS_PREV_HASH: &str =
    "0000000000000000000000000000000000000000000000000000000000000000";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AuditLogEntry {
    pub seq_no: i64,
    pub recorded_at_epoch: i64,
    pub actor_tag: String,
    pub action_kind: String,
    pub target_ref: String,
    pub detail_summary: String,
    pub prev_hash_hex: String,
    pub entry_hash_hex: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AuditChainVerification {
    pub is_valid: bool,
    pub total_entries: usize,
    pub tip_hash_hex: String,
    pub broken_at_seq: Option<i64>,
    pub failure_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DpBudgetLedgerEntry {
    pub entry_id: i64,
    pub vintage_id: String,
    pub released_at_epoch: i64,
    pub epsilon_spent: f64,
    pub sensitivity_l1: f64,
    pub cumulative_epsilon_after: f64,
    pub max_budget_cap: f64,
    pub empirical_mae: f64,
    pub output_checksum_hex: String,
    pub entry_kind: String,
    pub mechanism: String,
    pub zcdp_rho_spent: f64,
    pub cumulative_rho_after: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DpReleaseConfig {
    /// Privacy loss parameter `epsilon > 0` for this query release (`e.g., 0.25 .. 2.0`).
    pub epsilon: f64,
    /// Global L1 sensitivity (`1.0` for disjoint weekly incident case histograms).
    pub sensitivity_l1: f64,
    /// Maximum cumulative `epsilon` allowed on a single vintage before refusing release.
    pub max_budget_cap: f64,
    /// Optional deterministic seed (`None` mixes OS timestamp + vintage ID via HMAC-SHA256).
    pub deterministic_seed: Option<u64>,
}

impl Default for DpReleaseConfig {
    fn default() -> Self {
        Self {
            epsilon: 0.50,
            sensitivity_l1: 1.0,
            max_budget_cap: 4.0,
            deterministic_seed: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DpReleaseOutcome {
    pub noisy_dataset: TimeSeriesDataset,
    pub epsilon_spent: f64,
    pub cumulative_epsilon_for_vintage: f64,
    pub remaining_budget: f64,
    /// Two-sided geometric decay parameter `alpha = exp(-epsilon / sensitivity_l1)`.
    pub geometric_alpha: f64,
    /// Exact analytical pre-clamping RMSE `sqrt(2 * alpha) / (1 - alpha)`.
    pub theoretical_rmse: f64,
    pub empirical_mae: f64,
    pub empirical_rmse: f64,
    pub audit_entry_hash: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ZcdpReleaseConfig {
    /// Zero-concentrated differential privacy parameter `rho > 0`.
    pub rho: f64,
    /// Target privacy failure probability `delta in (0, 1)` for `rho -> (eps, delta)` conversion.
    pub delta: f64,
    /// Global L2 sensitivity (`1.0` for disjoint weekly incident case histograms).
    pub sensitivity_l2: f64,
    /// Maximum cumulative `rho` allowed on a single vintage before refusing release.
    pub max_rho_cap: f64,
    /// Maximum cumulative converted `epsilon` allowed on a single vintage.
    pub max_budget_cap: f64,
    /// Optional deterministic seed (`None` mixes OS timestamp + vintage ID via HMAC-SHA256).
    pub deterministic_seed: Option<u64>,
}

impl Default for ZcdpReleaseConfig {
    fn default() -> Self {
        Self {
            rho: 0.05,
            delta: 1e-5,
            sensitivity_l2: 1.0,
            max_rho_cap: 0.50,
            max_budget_cap: 6.0,
            deterministic_seed: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ZcdpReleaseOutcome {
    pub noisy_dataset: TimeSeriesDataset,
    pub rho_spent: f64,
    pub cumulative_rho_for_vintage: f64,
    pub remaining_rho_budget: f64,
    pub converted_epsilon_delta: f64,
    pub delta: f64,
    pub sigma: f64,
    pub theoretical_rmse: f64,
    pub empirical_mae: f64,
    pub empirical_rmse: f64,
    pub audit_entry_hash: String,
}

/// Convert `rho`-zCDP to `(epsilon, delta)`-DP via Bun & Steinke (2016) Proposition 1.3:
/// `epsilon(delta) = rho + 2 * sqrt(rho * ln(1 / delta))`.
pub fn zcdp_rho_to_epsilon(rho: f64, delta: f64) -> f64 {
    if rho <= 0.0 || !rho.is_finite() || delta <= 0.0 || delta >= 1.0 || !delta.is_finite() {
        return 0.0;
    }
    rho + 2.0 * (rho * (1.0 / delta).ln()).sqrt()
}

fn now_epoch_secs() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(1_763_020_800)
}

pub fn compute_audit_entry_hash(
    seq_no: i64,
    recorded_at_epoch: i64,
    actor_tag: &str,
    action_kind: &str,
    target_ref: &str,
    detail_summary: &str,
    prev_hash_hex: &str,
) -> String {
    let mut hasher = Sha256::new();
    hasher.update(
        format!(
            "{seq_no}|{recorded_at_epoch}|{actor_tag}|{action_kind}|{target_ref}|{detail_summary}|{prev_hash_hex}"
        )
        .as_bytes(),
    );
    let digest = hasher.finalize();
    let mut out = String::with_capacity(64);
    for b in digest {
        out.push_str(&format!("{b:02x}"));
    }
    out
}

pub(crate) fn append_audit_event_on_conn(
    conn: &Connection,
    actor_tag: &str,
    action_kind: &str,
    target_ref: &str,
    detail_summary: &str,
) -> Result<AuditLogEntry, EpiError> {
    let last_row: Option<(i64, String)> = conn
        .query_row(
            "SELECT seq_no, entry_hash_hex FROM audit_hash_chain ORDER BY seq_no DESC LIMIT 1",
            [],
            |r| Ok((r.get(0)?, r.get(1)?)),
        )
        .optional()
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let (next_seq, prev_hash) = match last_row {
        Some((s, h)) => (s + 1, h),
        None => (1, GENESIS_PREV_HASH.to_string()),
    };
    let ts = now_epoch_secs();
    let entry_hash = compute_audit_entry_hash(
        next_seq,
        ts,
        actor_tag,
        action_kind,
        target_ref,
        detail_summary,
        &prev_hash,
    );

    conn.execute(
        "INSERT INTO audit_hash_chain
         (seq_no, recorded_at_epoch, actor_tag, action_kind, target_ref, detail_summary, prev_hash_hex, entry_hash_hex)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
        params![
            next_seq,
            ts,
            actor_tag,
            action_kind,
            target_ref,
            detail_summary,
            prev_hash,
            entry_hash
        ],
    )
    .map_err(|e| EpiError::Storage(e.to_string()))?;

    Ok(AuditLogEntry {
        seq_no: next_seq,
        recorded_at_epoch: ts,
        actor_tag: actor_tag.to_string(),
        action_kind: action_kind.to_string(),
        target_ref: target_ref.to_string(),
        detail_summary: detail_summary.to_string(),
        prev_hash_hex: prev_hash,
        entry_hash_hex: entry_hash,
    })
}

pub fn append_audit_event(
    store: &EpiStore,
    actor_tag: &str,
    action_kind: &str,
    target_ref: &str,
    detail_summary: &str,
) -> Result<AuditLogEntry, EpiError> {
    append_audit_event_on_conn(
        store.conn(),
        actor_tag,
        action_kind,
        target_ref,
        detail_summary,
    )
}

pub fn list_audit_events(store: &EpiStore, limit: usize) -> Result<Vec<AuditLogEntry>, EpiError> {
    let mut stmt = store
        .conn()
        .prepare(
            "SELECT seq_no, recorded_at_epoch, actor_tag, action_kind, target_ref, detail_summary, prev_hash_hex, entry_hash_hex
             FROM audit_hash_chain
             ORDER BY seq_no DESC
             LIMIT ?1",
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let rows = stmt
        .query_map(params![limit as i64], |r| {
            Ok(AuditLogEntry {
                seq_no: r.get(0)?,
                recorded_at_epoch: r.get(1)?,
                actor_tag: r.get(2)?,
                action_kind: r.get(3)?,
                target_ref: r.get(4)?,
                detail_summary: r.get(5)?,
                prev_hash_hex: r.get(6)?,
                entry_hash_hex: r.get(7)?,
            })
        })
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let mut out = Vec::new();
    for r in rows {
        out.push(r.map_err(|e| EpiError::Storage(e.to_string()))?);
    }
    Ok(out)
}

/// Recompute every `SHA-256` hash link from genesis (`seq_no = 1`) to chain tip to verify
/// zero row tampering, insertion, or deletion.
pub fn verify_audit_chain(store: &EpiStore) -> Result<AuditChainVerification, EpiError> {
    let mut stmt = store
        .conn()
        .prepare(
            "SELECT seq_no, recorded_at_epoch, actor_tag, action_kind, target_ref, detail_summary, prev_hash_hex, entry_hash_hex
             FROM audit_hash_chain
             ORDER BY seq_no ASC",
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let rows = stmt
        .query_map([], |r| {
            Ok(AuditLogEntry {
                seq_no: r.get(0)?,
                recorded_at_epoch: r.get(1)?,
                actor_tag: r.get(2)?,
                action_kind: r.get(3)?,
                target_ref: r.get(4)?,
                detail_summary: r.get(5)?,
                prev_hash_hex: r.get(6)?,
                entry_hash_hex: r.get(7)?,
            })
        })
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let mut expected_prev = GENESIS_PREV_HASH.to_string();
    let mut expected_seq = 1i64;
    let mut count = 0usize;

    for row_res in rows {
        let e = row_res.map_err(|err| EpiError::Storage(err.to_string()))?;
        count += 1;
        if e.seq_no != expected_seq {
            return Ok(AuditChainVerification {
                is_valid: false,
                total_entries: count,
                tip_hash_hex: e.entry_hash_hex,
                broken_at_seq: Some(e.seq_no),
                failure_reason: Some(format!(
                    "Sequence gap detected: expected seq_no {expected_seq}, found {}",
                    e.seq_no
                )),
            });
        }
        if e.prev_hash_hex != expected_prev {
            return Ok(AuditChainVerification {
                is_valid: false,
                total_entries: count,
                tip_hash_hex: e.entry_hash_hex,
                broken_at_seq: Some(e.seq_no),
                failure_reason: Some(format!(
                    "Previous hash mismatch at seq_no {}",
                    e.seq_no
                )),
            });
        }
        let recomputed = compute_audit_entry_hash(
            e.seq_no,
            e.recorded_at_epoch,
            &e.actor_tag,
            &e.action_kind,
            &e.target_ref,
            &e.detail_summary,
            &e.prev_hash_hex,
        );
        if recomputed != e.entry_hash_hex {
            return Ok(AuditChainVerification {
                is_valid: false,
                total_entries: count,
                tip_hash_hex: e.entry_hash_hex,
                broken_at_seq: Some(e.seq_no),
                failure_reason: Some(format!(
                    "Entry SHA-256 payload digest mismatch at seq_no {}",
                    e.seq_no
                )),
            });
        }
        expected_prev = e.entry_hash_hex;
        expected_seq += 1;
    }

    Ok(AuditChainVerification {
        is_valid: true,
        total_entries: count,
        tip_hash_hex: expected_prev,
        broken_at_seq: None,
        failure_reason: None,
    })
}

pub fn get_vintage_spent_epsilon(store: &EpiStore, vintage_id: &str) -> Result<f64, EpiError> {
    let spent: Option<f64> = store
        .conn()
        .query_row(
            "SELECT SUM(epsilon_spent)
             FROM dp_budget_ledger
             WHERE vintage_id = ?1
               AND entry_id > COALESCE(
                   (SELECT MAX(entry_id) FROM dp_budget_ledger
                    WHERE vintage_id = ?1
                      AND (entry_kind = 'RESET' OR output_checksum_hex LIKE 'RESET%')),
                   0
               )",
            params![vintage_id],
            |r| r.get(0),
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;
    Ok(spent.unwrap_or(0.0))
}

pub fn get_vintage_spent_rho(store: &EpiStore, vintage_id: &str) -> Result<f64, EpiError> {
    let spent: Option<f64> = store
        .conn()
        .query_row(
            "SELECT SUM(zcdp_rho_spent)
             FROM dp_budget_ledger
             WHERE vintage_id = ?1
               AND entry_id > COALESCE(
                   (SELECT MAX(entry_id) FROM dp_budget_ledger
                    WHERE vintage_id = ?1
                      AND (entry_kind = 'RESET' OR output_checksum_hex LIKE 'RESET%')),
                   0
               )",
            params![vintage_id],
            |r| r.get(0),
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;
    Ok(spent.unwrap_or(0.0))
}

pub fn list_dp_budget_entries(
    store: &EpiStore,
    limit: usize,
) -> Result<Vec<DpBudgetLedgerEntry>, EpiError> {
    let mut stmt = store
        .conn()
        .prepare(
            "SELECT entry_id, vintage_id, released_at_epoch, epsilon_spent, sensitivity_l1,
                    cumulative_epsilon_after, max_budget_cap, empirical_mae, output_checksum_hex,
                    entry_kind, mechanism, zcdp_rho_spent, cumulative_rho_after
             FROM dp_budget_ledger
             ORDER BY entry_id DESC
             LIMIT ?1",
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let rows = stmt
        .query_map(params![limit as i64], |r| {
            Ok(DpBudgetLedgerEntry {
                entry_id: r.get(0)?,
                vintage_id: r.get(1)?,
                released_at_epoch: r.get(2)?,
                epsilon_spent: r.get(3)?,
                sensitivity_l1: r.get(4)?,
                cumulative_epsilon_after: r.get(5)?,
                max_budget_cap: r.get(6)?,
                empirical_mae: r.get(7)?,
                output_checksum_hex: r.get(8)?,
                entry_kind: r.get(9)?,
                mechanism: r.get(10)?,
                zcdp_rho_spent: r.get(11)?,
                cumulative_rho_after: r.get(12)?,
            })
        })
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let mut out = Vec::new();
    for r in rows {
        out.push(r.map_err(|e| EpiError::Storage(e.to_string()))?);
    }
    Ok(out)
}

/// Append a governance `RESET` epoch marker into the strictly append-only `dp_budget_ledger`
/// and record an immutable `RESET_DP_BUDGET` event in `audit_hash_chain`.
/// Prior release rows are never deleted or mutated.
pub fn reset_vintage_dp_budget(store: &EpiStore, vintage_id: &str) -> Result<AuditLogEntry, EpiError> {
    let prior_eps = get_vintage_spent_epsilon(store, vintage_id)?;
    let prior_rho = get_vintage_spent_rho(store, vintage_id)?;
    let audit = append_audit_event(
        store,
        "privacy-officer",
        "RESET_DP_BUDGET",
        vintage_id,
        &format!(
            "Appended governance reset epoch (prior spent_epsilon={prior_eps:.4}, prior_rho={prior_rho:.4}); historical ledger rows preserved"
        ),
    )?;
    let ts = now_epoch_secs();
    let tag = format!("RESET:{}", &audit.entry_hash_hex[..8]);
    store
        .conn()
        .execute(
            "INSERT INTO dp_budget_ledger
             (vintage_id, released_at_epoch, epsilon_spent, sensitivity_l1, cumulative_epsilon_after,
              max_budget_cap, empirical_mae, output_checksum_hex, entry_kind, mechanism, zcdp_rho_spent, cumulative_rho_after)
             VALUES (?1, ?2, 0.0, 0.0, 0.0, 0.0, 0.0, ?3, 'RESET', 'GOVERNANCE_RESET', 0.0, 0.0)",
            params![vintage_id, ts, tag],
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;
    Ok(audit)
}

/// Sample an exact integer draw `Z` from the Two-Sided Geometric (`Discrete Laplace`) distribution
/// `P(Z = k) = ((1 - alpha) / (1 + alpha)) * alpha^{|k|}` where `alpha = exp(-epsilon / sensitivity)`
/// using a keyed `HMAC-SHA256` counter stream.
fn sample_discrete_laplace_ctr(key: &[u8], counter: u64, alpha: f64) -> i64 {
    if alpha <= 1e-15 {
        return 0;
    }
    let mac = hmac_sha256(key, &counter.to_le_bytes());
    let u_raw = u64::from_le_bytes(mac[0..8].try_into().unwrap());
    // Uniform open interval (0, 1) with 53-bit mantissa precision
    let u = (((u_raw >> 11) as f64) + 0.5) / ((1u64 << 53) as f64);

    // P(Z = 0) = (1 - alpha) / (1 + alpha)
    let p0 = (1.0 - alpha) / (1.0 + alpha);
    if u <= p0 {
        return 0;
    }
    let sign = if (mac[8] & 1) == 0 { 1i64 } else { -1i64 };
    // Conditional on |Z| >= 1, |Z| - 1 follows Geometric(1 - alpha): P(|Z| >= k) = alpha^{k - 1}
    let u_tail = ((u - p0) / (1.0 - p0)).clamp(1e-15, 1.0 - 1e-15);
    let mag = 1 + (u_tail.ln() / alpha.ln()).floor().max(0.0) as i64;
    sign * mag
}

/// Sample an exact integer draw `Z` from the Discrete Gaussian distribution `N_Z(0, sigma^2)`
/// (`P(Z = k) proportional to exp(-k^2 / (2 * sigma^2))`) via Canonne-Kamath-Steinke (2020)
/// Algorithm 3 exact rejection sampling over a Discrete Laplace(`t = floor(sigma) + 1`) proposal.
pub fn sample_discrete_gaussian_ctr(key: &[u8], counter_base: u64, sigma: f64) -> i64 {
    if !sigma.is_finite() || sigma <= 1e-6 {
        return 0;
    }
    let t = (sigma.floor() + 1.0).max(1.0);
    let alpha = (-1.0 / t).exp().clamp(0.0, 0.999_999_9);
    let sigma2 = sigma * sigma;
    let shift = sigma2 / t;

    for attempt in 0u64..4096 {
        let ctr_prop = (counter_base << 16) | (attempt << 1);
        let y = sample_discrete_laplace_ctr(key, ctr_prop, alpha);
        let diff = (y.abs() as f64) - shift;
        let accept_prob = (-(diff * diff) / (2.0 * sigma2)).exp().clamp(0.0, 1.0);

        let ctr_acc = (counter_base << 16) | (attempt << 1) | 1;
        let mac = hmac_sha256(key, &ctr_acc.to_le_bytes());
        let u_raw = u64::from_le_bytes(mac[0..8].try_into().unwrap());
        let u = (((u_raw >> 11) as f64) + 0.5) / ((1u64 << 53) as f64);
        if u <= accept_prob {
            return y;
        }
    }
    0
}

/// Release an `epsilon`-Differentially Private surveillance time series via the Discrete Laplace
/// mechanism, enforcing `max_budget_cap` in SQLite and recording a `SHA-256` hash-chained audit event.
pub fn release_differentially_private_series(
    store: &EpiStore,
    dataset: &TimeSeriesDataset,
    cfg: &DpReleaseConfig,
) -> Result<DpReleaseOutcome, EpiError> {
    let eps = cfg.epsilon;
    let sens = cfg.sensitivity_l1.max(0.01);
    if !eps.is_finite() || eps <= 0.005 {
        return Err(EpiError::Schema(
            "Differential privacy epsilon must be > 0.005".to_string(),
        ));
    }

    let prior_spent = get_vintage_spent_epsilon(store, &dataset.vintage_id)?;
    let new_cum = prior_spent + eps;
    if new_cum > cfg.max_budget_cap + 1e-9 {
        return Err(EpiError::Schema(format!(
            "DP privacy budget exhausted for '{}': prior={prior_spent:.3} + requested={eps:.3} exceeds cap={:.3}",
            dataset.vintage_id, cfg.max_budget_cap
        )));
    }

    let alpha = (-eps / sens).exp().clamp(0.0, 0.999_999);
    let theoretical_rmse = (2.0 * alpha).sqrt() / (1.0 - alpha).max(1e-9);

    let ts = now_epoch_secs();
    let seed_val = cfg
        .deterministic_seed
        .unwrap_or_else(|| (ts as u64) ^ 0xA5A5_5A5A_D901_2026);
    let key_material = format!(
        "DP_CTR_KEY:{}:{}:{:.6}",
        dataset.vintage_id, seed_val, eps
    );

    let mut noisy_points = Vec::with_capacity(dataset.points.len());
    let mut abs_err_sum = 0.0;
    let mut sq_err_sum = 0.0;

    for (idx, pt) in dataset.points.iter().enumerate() {
        let z = sample_discrete_laplace_ctr(key_material.as_bytes(), idx as u64, alpha);
        // Post-processing projection to non-negative integer counts preserves exact epsilon-DP
        let true_c = pt.count.max(0.0).round();
        let noisy_c = ((true_c as i64) + z).max(0) as f64;
        let diff = noisy_c - pt.count.max(0.0);
        abs_err_sum += diff.abs();
        sq_err_sum += diff * diff;

        noisy_points.push(TimeSeriesPoint {
            date_label: pt.date_label.clone(),
            epiweek: pt.epiweek,
            count: noisy_c,
            suppressed: false,
        });
    }

    let n_pts = dataset.points.len().max(1) as f64;
    let emp_mae = abs_err_sum / n_pts;
    let emp_rmse = (sq_err_sum / n_pts).sqrt();

    let out_hash = EpiStore::compute_dataset_checksum(&noisy_points);
    let noisy_dataset = TimeSeriesDataset {
        vintage_id: format!("{}_dp_e{:03}_{}", dataset.vintage_id, (eps * 100.0).round() as u32, out_hash),
        title: format!("{} [DP eps={:.2}]", dataset.title, eps),
        source_kind: SourceKind::LocalCsvTsv,
        source_url: dataset.source_url.clone(),
        ingested_at_utc: format!("epoch:{ts}"),
        default_si_preset: dataset.default_si_preset.clone(),
        checksum_hex: out_hash.clone(),
        points: noisy_points,
    };

    let prior_rho = get_vintage_spent_rho(store, &dataset.vintage_id)?;
    store
        .conn()
        .execute(
            "INSERT INTO dp_budget_ledger
             (vintage_id, released_at_epoch, epsilon_spent, sensitivity_l1, cumulative_epsilon_after,
              max_budget_cap, empirical_mae, output_checksum_hex, entry_kind, mechanism, zcdp_rho_spent, cumulative_rho_after)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, 'RELEASE', 'DISCRETE_LAPLACE', 0.0, ?9)",
            params![
                dataset.vintage_id,
                ts,
                eps,
                sens,
                new_cum,
                cfg.max_budget_cap,
                emp_mae,
                out_hash,
                prior_rho
            ],
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let audit = append_audit_event(
        store,
        "dp-mechanism",
        "DP_RELEASE",
        &dataset.vintage_id,
        &format!(
            "DiscreteLaplace eps={eps:.3} sens={sens:.1} cum_eps={new_cum:.3}/{:.3} mae={emp_mae:.2} out_hash={out_hash}",
            cfg.max_budget_cap
        ),
    )?;

    Ok(DpReleaseOutcome {
        noisy_dataset,
        epsilon_spent: eps,
        cumulative_epsilon_for_vintage: new_cum,
        remaining_budget: (cfg.max_budget_cap - new_cum).max(0.0),
        geometric_alpha: alpha,
        theoretical_rmse,
        empirical_mae: emp_mae,
        empirical_rmse: emp_rmse,
        audit_entry_hash: audit.entry_hash_hex,
    })
}

/// Release a `rho`-zCDP Differentially Private surveillance time series via the exact
/// Canonne-Kamath-Steinke (2020) Discrete Gaussian mechanism (`sigma = sensitivity_l2 / sqrt(2 * rho)`),
/// tracking additive `rho`-composition in `dp_budget_ledger` and appending `DP_ZCDP_RELEASE`.
pub fn release_discrete_gaussian_zcdp_series(
    store: &EpiStore,
    dataset: &TimeSeriesDataset,
    cfg: &ZcdpReleaseConfig,
) -> Result<ZcdpReleaseOutcome, EpiError> {
    let rho = cfg.rho;
    let delta = cfg.delta;
    let sens = cfg.sensitivity_l2.max(0.01);
    if !rho.is_finite() || rho <= 0.001 {
        return Err(EpiError::Schema(
            "zCDP privacy parameter rho must be > 0.001".to_string(),
        ));
    }
    if !delta.is_finite() || delta <= 0.0 || delta >= 1.0 {
        return Err(EpiError::Schema(
            "zCDP target delta must lie strictly in (0, 1)".to_string(),
        ));
    }

    let prior_rho = get_vintage_spent_rho(store, &dataset.vintage_id)?;
    let new_cum_rho = prior_rho + rho;
    if new_cum_rho > cfg.max_rho_cap + 1e-9 {
        return Err(EpiError::Schema(format!(
            "zCDP rho budget exhausted for '{}': prior={prior_rho:.4} + requested={rho:.4} exceeds cap={:.4}",
            dataset.vintage_id, cfg.max_rho_cap
        )));
    }

    let prior_eps = get_vintage_spent_epsilon(store, &dataset.vintage_id)?;
    let cum_zcdp_eps = zcdp_rho_to_epsilon(new_cum_rho, delta);
    let prev_zcdp_eps = zcdp_rho_to_epsilon(prior_rho, delta);
    let marginal_eps = (cum_zcdp_eps - prev_zcdp_eps).max(0.0);
    let new_cum_eps = prior_eps + marginal_eps;
    if new_cum_eps > cfg.max_budget_cap + 1e-9 {
        return Err(EpiError::Schema(format!(
            "Converted epsilon budget exhausted for '{}': prior={prior_eps:.3} + marginal={marginal_eps:.3} exceeds cap={:.3}",
            dataset.vintage_id, cfg.max_budget_cap
        )));
    }

    let sigma = sens / (2.0 * rho).sqrt();
    let ts = now_epoch_secs();
    let seed_val = cfg
        .deterministic_seed
        .unwrap_or_else(|| (ts as u64) ^ 0x5A17_C0D3_9E37_79B9);
    let key_material = format!(
        "ZCDP_DG_KEY:{}:{}:{:.6}",
        dataset.vintage_id, seed_val, rho
    );

    let mut noisy_points = Vec::with_capacity(dataset.points.len());
    let mut abs_err_sum = 0.0;
    let mut sq_err_sum = 0.0;

    for (idx, pt) in dataset.points.iter().enumerate() {
        let z = sample_discrete_gaussian_ctr(key_material.as_bytes(), idx as u64, sigma);
        let true_c = pt.count.max(0.0).round();
        let noisy_c = ((true_c as i64) + z).max(0) as f64;
        let diff = noisy_c - pt.count.max(0.0);
        abs_err_sum += diff.abs();
        sq_err_sum += diff * diff;

        noisy_points.push(TimeSeriesPoint {
            date_label: pt.date_label.clone(),
            epiweek: pt.epiweek,
            count: noisy_c,
            suppressed: false,
        });
    }

    let n_pts = dataset.points.len().max(1) as f64;
    let emp_mae = abs_err_sum / n_pts;
    let emp_rmse = (sq_err_sum / n_pts).sqrt();

    let out_hash = EpiStore::compute_dataset_checksum(&noisy_points);
    let noisy_dataset = TimeSeriesDataset {
        vintage_id: format!(
            "{}_zcdp_r{:04}_{}",
            dataset.vintage_id,
            (rho * 1000.0).round() as u32,
            out_hash
        ),
        title: format!(
            "{} [zCDP rho={:.3}, sigma={:.2}]",
            dataset.title, rho, sigma
        ),
        source_kind: SourceKind::LocalCsvTsv,
        source_url: dataset.source_url.clone(),
        ingested_at_utc: format!("epoch:{ts}"),
        default_si_preset: dataset.default_si_preset.clone(),
        checksum_hex: out_hash.clone(),
        points: noisy_points,
    };

    store
        .conn()
        .execute(
            "INSERT INTO dp_budget_ledger
             (vintage_id, released_at_epoch, epsilon_spent, sensitivity_l1, cumulative_epsilon_after,
              max_budget_cap, empirical_mae, output_checksum_hex, entry_kind, mechanism, zcdp_rho_spent, cumulative_rho_after)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, 'RELEASE', 'DISCRETE_GAUSSIAN_ZCDP', ?9, ?10)",
            params![
                dataset.vintage_id,
                ts,
                marginal_eps,
                sens,
                new_cum_eps,
                cfg.max_budget_cap,
                emp_mae,
                out_hash,
                rho,
                new_cum_rho
            ],
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let audit = append_audit_event(
        store,
        "dp-mechanism",
        "DP_ZCDP_RELEASE",
        &dataset.vintage_id,
        &format!(
            "DiscreteGaussian zCDP rho={rho:.4} cum_rho={new_cum_rho:.4}/{:.4} sigma={sigma:.3} conv_eps={cum_zcdp_eps:.3} (delta={delta:.1e}) mae={emp_mae:.2} out_hash={out_hash}",
            cfg.max_rho_cap
        ),
    )?;

    Ok(ZcdpReleaseOutcome {
        noisy_dataset,
        rho_spent: rho,
        cumulative_rho_for_vintage: new_cum_rho,
        remaining_rho_budget: (cfg.max_rho_cap - new_cum_rho).max(0.0),
        converted_epsilon_delta: cum_zcdp_eps,
        delta,
        sigma,
        theoretical_rmse: sigma,
        empirical_mae: emp_mae,
        empirical_rmse: emp_rmse,
        audit_entry_hash: audit.entry_hash_hex,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use epi_core::builtin_measles_weekly_dataset;

    #[test]
    fn test_dp_budget_cap_and_tamper_evident_audit_chain() {
        let mut store = EpiStore::open_in_memory().unwrap();
        let ds = builtin_measles_weekly_dataset();
        store.commit_vintage(&ds).unwrap();

        let cfg = DpReleaseConfig {
            epsilon: 1.25,
            sensitivity_l1: 1.0,
            max_budget_cap: 2.0,
            deterministic_seed: Some(42),
        };
        let r1 = release_differentially_private_series(&store, &ds, &cfg).unwrap();
        assert!((r1.cumulative_epsilon_for_vintage - 1.25).abs() < 1e-9);

        // Second release of 1.25 exceeds cap 2.0 -> must be rejected!
        let r2 = release_differentially_private_series(&store, &ds, &cfg);
        assert!(r2.is_err());

        // Audit chain must verify cleanly across COMMIT_VINTAGE + DP_RELEASE
        let ver = verify_audit_chain(&store).unwrap();
        assert!(ver.is_valid);
        assert!(ver.total_entries >= 2);
    }

    #[test]
    fn test_dp_budget_ledger_strictly_append_only_on_reset() {
        let mut store = EpiStore::open_in_memory().unwrap();
        let ds = builtin_measles_weekly_dataset();
        store.commit_vintage(&ds).unwrap();

        let cfg = DpReleaseConfig {
            epsilon: 0.75,
            sensitivity_l1: 1.0,
            max_budget_cap: 2.0,
            deterministic_seed: Some(77),
        };
        release_differentially_private_series(&store, &ds, &cfg).unwrap();
        release_differentially_private_series(&store, &ds, &cfg).unwrap();
        assert!((get_vintage_spent_epsilon(&store, &ds.vintage_id).unwrap() - 1.5).abs() < 1e-9);

        let audit = reset_vintage_dp_budget(&store, &ds.vintage_id).unwrap();
        assert_eq!(audit.action_kind, "RESET_DP_BUDGET");
        assert!((get_vintage_spent_epsilon(&store, &ds.vintage_id).unwrap() - 0.0).abs() < 1e-9);

        // Prior rows MUST remain in dp_budget_ledger alongside the new RESET row
        let entries = list_dp_budget_entries(&store, 20).unwrap();
        assert_eq!(entries.len(), 3);
        assert_eq!(entries[0].entry_kind, "RESET");
        assert_eq!(entries[1].entry_kind, "RELEASE");
        assert_eq!(entries[2].entry_kind, "RELEASE");

        // Post-reset release uses fresh budget epoch
        release_differentially_private_series(&store, &ds, &cfg).unwrap();
        assert!((get_vintage_spent_epsilon(&store, &ds.vintage_id).unwrap() - 0.75).abs() < 1e-9);
        assert_eq!(list_dp_budget_entries(&store, 20).unwrap().len(), 4);

        // SQLite triggers must block direct DELETE and UPDATE on dp_budget_ledger
        let del_res = store
            .conn()
            .execute("DELETE FROM dp_budget_ledger WHERE entry_id = 1", []);
        assert!(del_res.is_err(), "DELETE trigger must abort");
        let upd_res = store
            .conn()
            .execute("UPDATE dp_budget_ledger SET epsilon_spent = 0.0 WHERE entry_id = 1", []);
        assert!(upd_res.is_err(), "UPDATE trigger must abort");

        assert!(verify_audit_chain(&store).unwrap().is_valid);
    }

    #[test]
    fn test_discrete_gaussian_zcdp_moments_and_determinism() {
        let key = b"zcdp_moments_test_key";
        let sigma = 3.2;
        let n = 20_000u64;
        let mut sum = 0.0f64;
        let mut sum_sq = 0.0f64;
        for i in 0..n {
            let z = sample_discrete_gaussian_ctr(key, i, sigma);
            let z2 = sample_discrete_gaussian_ctr(key, i, sigma);
            assert_eq!(z, z2, "HMAC-SHA256 sampler must be deterministic");
            let zf = z as f64;
            sum += zf;
            sum_sq += zf * zf;
        }
        let mean = sum / (n as f64);
        let var = sum_sq / (n as f64) - mean * mean;
        assert!(mean.abs() < 0.15, "mean {mean} far from 0");
        let expected_var = sigma * sigma;
        assert!(
            (var - expected_var).abs() < 0.65,
            "empirical variance {var:.3} far from sigma^2={expected_var:.3}"
        );

        // Small sigma check (sigma = 0.6 < 1.0)
        let mut s_sum = 0.0;
        let mut s_sq = 0.0;
        for i in 0..10_000u64 {
            let z = sample_discrete_gaussian_ctr(key, 100_000 + i, 0.6) as f64;
            s_sum += z;
            s_sq += z * z;
        }
        let s_mean = s_sum / 10_000.0;
        let s_var = s_sq / 10_000.0 - s_mean * s_mean;
        assert!(s_mean.abs() < 0.08);
        assert!(s_var > 0.10 && s_var < 0.50);
    }

    #[test]
    fn test_zcdp_release_additive_rho_composition_and_audit_chain() {
        let mut store = EpiStore::open_in_memory().unwrap();
        let ds = builtin_measles_weekly_dataset();
        store.commit_vintage(&ds).unwrap();

        let cfg = ZcdpReleaseConfig {
            rho: 0.05,
            delta: 1e-5,
            sensitivity_l2: 1.0,
            max_rho_cap: 0.12,
            max_budget_cap: 6.0,
            deterministic_seed: Some(909),
        };
        let o1 = release_discrete_gaussian_zcdp_series(&store, &ds, &cfg).unwrap();
        assert!((o1.cumulative_rho_for_vintage - 0.05).abs() < 1e-9);
        assert!((o1.converted_epsilon_delta - zcdp_rho_to_epsilon(0.05, 1e-5)).abs() < 1e-9);

        let o2 = release_discrete_gaussian_zcdp_series(&store, &ds, &cfg).unwrap();
        assert!((o2.cumulative_rho_for_vintage - 0.10).abs() < 1e-9);
        assert!((get_vintage_spent_rho(&store, &ds.vintage_id).unwrap() - 0.10).abs() < 1e-9);
        assert!(
            (get_vintage_spent_epsilon(&store, &ds.vintage_id).unwrap()
                - zcdp_rho_to_epsilon(0.10, 1e-5))
            .abs()
                < 1e-6
        );

        // 0.10 + 0.05 = 0.15 > 0.12 -> must be refused
        let o3 = release_discrete_gaussian_zcdp_series(&store, &ds, &cfg);
        assert!(o3.is_err());

        let ledger = list_dp_budget_entries(&store, 10).unwrap();
        assert_eq!(ledger.len(), 2);
        assert_eq!(ledger[0].mechanism, "DISCRETE_GAUSSIAN_ZCDP");
        assert!(verify_audit_chain(&store).unwrap().is_valid);
    }

    #[test]
    fn test_dp_budget_ledger_blocks_insert_or_replace_bypass() {
        let mut store = EpiStore::open_in_memory().unwrap();
        let ds = builtin_measles_weekly_dataset();
        store.commit_vintage(&ds).unwrap();

        let lap_cfg = DpReleaseConfig {
            epsilon: 0.75,
            sensitivity_l1: 1.0,
            max_budget_cap: 5.0,
            deterministic_seed: Some(123),
        };
        release_differentially_private_series(&store, &ds, &lap_cfg).unwrap();

        let zcdp_cfg = ZcdpReleaseConfig {
            rho: 0.04,
            delta: 1e-5,
            sensitivity_l2: 1.0,
            max_rho_cap: 0.5,
            max_budget_cap: 5.0,
            deterministic_seed: Some(456),
        };
        release_discrete_gaussian_zcdp_series(&store, &ds, &zcdp_cfg).unwrap();

        let eps_before = get_vintage_spent_epsilon(&store, &ds.vintage_id).unwrap();
        let rho_before = get_vintage_spent_rho(&store, &ds.vintage_id).unwrap();
        let entries_before = list_dp_budget_entries(&store, 10).unwrap();
        assert_eq!(entries_before.len(), 2);
        let existing_id = entries_before[0].entry_id;
        let orig_checksum = entries_before[0].output_checksum_hex.clone();

        let replace_res = store.conn().execute(
            "INSERT OR REPLACE INTO dp_budget_ledger
             (entry_id, vintage_id, released_at_epoch, epsilon_spent, sensitivity_l1,
              cumulative_epsilon_after, max_budget_cap, empirical_mae, output_checksum_hex)
             VALUES (?1, ?2, 0, 0.0, 1.0, 0.0, 10.0, 0.0, 'RESET:tamper')",
            params![existing_id, ds.vintage_id],
        );
        match replace_res {
            Err(rusqlite::Error::SqliteFailure(err, _)) => {
                assert_eq!(err.code, rusqlite::ErrorCode::ConstraintViolation);
            }
            other => panic!("INSERT OR REPLACE must fail with ConstraintViolation, got {other:?}"),
        }

        let eps_after = get_vintage_spent_epsilon(&store, &ds.vintage_id).unwrap();
        let rho_after = get_vintage_spent_rho(&store, &ds.vintage_id).unwrap();
        assert!((eps_after - eps_before).abs() < 1e-12);
        assert!((rho_after - rho_before).abs() < 1e-12);

        let entries_after = list_dp_budget_entries(&store, 10).unwrap();
        assert_eq!(entries_after.len(), 2);
        assert_eq!(entries_after[0].entry_id, existing_id);
        assert_eq!(entries_after[0].output_checksum_hex, orig_checksum);
    }
}
