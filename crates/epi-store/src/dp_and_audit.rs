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
            "SELECT SUM(epsilon_spent) FROM dp_budget_ledger WHERE vintage_id = ?1",
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
                    cumulative_epsilon_after, max_budget_cap, empirical_mae, output_checksum_hex
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
            })
        })
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    let mut out = Vec::new();
    for r in rows {
        out.push(r.map_err(|e| EpiError::Storage(e.to_string()))?);
    }
    Ok(out)
}

/// Reset the accumulated `epsilon` ledger for `vintage_id` while recording an immutable
/// `RESET_DP_BUDGET` event in `audit_hash_chain` so governance resets remain permanently auditable.
pub fn reset_vintage_dp_budget(store: &EpiStore, vintage_id: &str) -> Result<AuditLogEntry, EpiError> {
    let prior = get_vintage_spent_epsilon(store, vintage_id)?;
    store
        .conn()
        .execute(
            "DELETE FROM dp_budget_ledger WHERE vintage_id = ?1",
            params![vintage_id],
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;

    append_audit_event(
        store,
        "privacy-officer",
        "RESET_DP_BUDGET",
        vintage_id,
        &format!("Cleared prior spent_epsilon={prior:.4} under explicit governance reset"),
    )
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

    store
        .conn()
        .execute(
            "INSERT INTO dp_budget_ledger
             (vintage_id, released_at_epoch, epsilon_spent, sensitivity_l1, cumulative_epsilon_after, max_budget_cap, empirical_mae, output_checksum_hex)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            params![
                dataset.vintage_id,
                ts,
                eps,
                sens,
                new_cum,
                cfg.max_budget_cap,
                emp_mae,
                out_hash
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
}
