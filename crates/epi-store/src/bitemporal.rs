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

use crate::dp_and_audit::append_audit_event_on_conn;
use epi_core::{
    EpiError, MmwrWeek, SourceKind, TimeSeriesDataset, TimeSeriesPoint, VintageMetadata,
};
use rusqlite::{params, Connection, OptionalExtension};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackfillRevisionDelta {
    pub date_label: String,
    pub older_count: f64,
    pub newer_count: f64,
    pub delta: f64,
}

#[derive(Debug, Clone)]
pub struct CachedHttpPayload {
    pub cache_key: String,
    pub fetched_epoch_secs: i64,
    pub ttl_secs: i64,
    pub is_fresh: bool,
    pub body_json: String,
}

/// Embedded SQLite WAL repository storing bitemporal surveillance vintages, lazy-cached API responses,
/// Differential Privacy `epsilon` budgets, and a tamper-evident `SHA-256` hash-chained audit log.
pub struct EpiStore {
    conn: Connection,
}

impl EpiStore {
    pub fn open_file<P: AsRef<Path>>(path: P) -> Result<Self, EpiError> {
        let conn = Connection::open(path).map_err(|e| EpiError::Storage(e.to_string()))?;
        let store = Self { conn };
        store.init_schema()?;
        Ok(store)
    }

    pub fn open_in_memory() -> Result<Self, EpiError> {
        let conn = Connection::open_in_memory().map_err(|e| EpiError::Storage(e.to_string()))?;
        let store = Self { conn };
        store.init_schema()?;
        Ok(store)
    }

    pub(crate) fn conn(&self) -> &Connection {
        &self.conn
    }

    fn init_schema(&self) -> Result<(), EpiError> {
        self.conn
            .execute_batch(
                "
                PRAGMA journal_mode = WAL;
                PRAGMA synchronous = NORMAL;
                PRAGMA foreign_keys = ON;

                CREATE TABLE IF NOT EXISTS vintages (
                    vintage_id       TEXT PRIMARY KEY,
                    title            TEXT NOT NULL,
                    source_kind      TEXT NOT NULL,
                    source_url       TEXT NOT NULL,
                    ingested_at_utc  TEXT NOT NULL,
                    default_si       TEXT NOT NULL,
                    row_count        INTEGER NOT NULL,
                    checksum_hex     TEXT NOT NULL
                );

                CREATE TABLE IF NOT EXISTS vintage_points (
                    vintage_id       TEXT NOT NULL,
                    seq_idx          INTEGER NOT NULL,
                    date_label       TEXT NOT NULL,
                    mmwr_year        INTEGER,
                    mmwr_week        INTEGER,
                    count_val        REAL NOT NULL,
                    suppressed       INTEGER NOT NULL DEFAULT 0,
                    PRIMARY KEY (vintage_id, seq_idx),
                    FOREIGN KEY (vintage_id) REFERENCES vintages(vintage_id) ON DELETE CASCADE
                );

                CREATE TABLE IF NOT EXISTS http_lazy_cache (
                    cache_key          TEXT PRIMARY KEY,
                    url                TEXT NOT NULL,
                    fetched_epoch_secs INTEGER NOT NULL,
                    ttl_secs           INTEGER NOT NULL,
                    etag               TEXT,
                    body_json          TEXT NOT NULL
                );

                CREATE TABLE IF NOT EXISTS dp_budget_ledger (
                    entry_id                 INTEGER PRIMARY KEY AUTOINCREMENT,
                    vintage_id               TEXT NOT NULL,
                    released_at_epoch        INTEGER NOT NULL,
                    epsilon_spent            REAL NOT NULL,
                    sensitivity_l1           REAL NOT NULL,
                    cumulative_epsilon_after REAL NOT NULL,
                    max_budget_cap           REAL NOT NULL,
                    empirical_mae            REAL NOT NULL,
                    output_checksum_hex      TEXT NOT NULL,
                    entry_kind               TEXT NOT NULL DEFAULT 'RELEASE',
                    mechanism                TEXT NOT NULL DEFAULT 'DISCRETE_LAPLACE',
                    zcdp_rho_spent           REAL NOT NULL DEFAULT 0.0,
                    cumulative_rho_after     REAL NOT NULL DEFAULT 0.0
                );

                CREATE TABLE IF NOT EXISTS audit_hash_chain (
                    seq_no            INTEGER PRIMARY KEY,
                    recorded_at_epoch INTEGER NOT NULL,
                    actor_tag         TEXT NOT NULL,
                    action_kind       TEXT NOT NULL,
                    target_ref        TEXT NOT NULL,
                    detail_summary    TEXT NOT NULL,
                    prev_hash_hex     TEXT NOT NULL,
                    entry_hash_hex    TEXT NOT NULL
                );
                ",
            )
            .map_err(|e| EpiError::Storage(e.to_string()))?;

        // Migrate pre-v2.0 dp_budget_ledger tables before creating triggers
        let mut cols = std::collections::HashSet::new();
        {
            let mut stmt = self
                .conn
                .prepare("PRAGMA table_info(dp_budget_ledger)")
                .map_err(|e| EpiError::Storage(e.to_string()))?;
            let iter = stmt
                .query_map([], |r| r.get::<_, String>(1))
                .map_err(|e| EpiError::Storage(e.to_string()))?;
            for c in iter {
                cols.insert(c.map_err(|e| EpiError::Storage(e.to_string()))?);
            }
        }
        let migrations = [
            ("entry_kind", "ALTER TABLE dp_budget_ledger ADD COLUMN entry_kind TEXT NOT NULL DEFAULT 'RELEASE'"),
            ("mechanism", "ALTER TABLE dp_budget_ledger ADD COLUMN mechanism TEXT NOT NULL DEFAULT 'DISCRETE_LAPLACE'"),
            ("zcdp_rho_spent", "ALTER TABLE dp_budget_ledger ADD COLUMN zcdp_rho_spent REAL NOT NULL DEFAULT 0.0"),
            ("cumulative_rho_after", "ALTER TABLE dp_budget_ledger ADD COLUMN cumulative_rho_after REAL NOT NULL DEFAULT 0.0"),
        ];
        for (col_name, ddl) in migrations {
            if !cols.contains(col_name) {
                self.conn
                    .execute(ddl, [])
                    .map_err(|e| EpiError::Storage(e.to_string()))?;
            }
        }

        // Strictly enforce append-only immutability on dp_budget_ledger (audit_hash_chain remains
        // testable for external tamper detection in headless_validation).
        self.conn
            .execute_batch(
                "
                CREATE TRIGGER IF NOT EXISTS dp_budget_ledger_no_delete
                BEFORE DELETE ON dp_budget_ledger
                BEGIN
                    SELECT RAISE(ABORT, 'dp_budget_ledger is strictly append-only: DELETE is forbidden');
                END;

                CREATE TRIGGER IF NOT EXISTS dp_budget_ledger_no_update
                BEFORE UPDATE ON dp_budget_ledger
                BEGIN
                    SELECT RAISE(ABORT, 'dp_budget_ledger is strictly append-only: UPDATE is forbidden');
                END;

                CREATE TRIGGER IF NOT EXISTS dp_budget_ledger_no_replace
                BEFORE INSERT ON dp_budget_ledger
                WHEN NEW.entry_id IS NOT NULL AND EXISTS (
                    SELECT 1 FROM dp_budget_ledger WHERE entry_id = NEW.entry_id
                )
                BEGIN
                    SELECT RAISE(ABORT, 'dp_budget_ledger is strictly append-only: REPLACE is forbidden');
                END;
                ",
            )
            .map_err(|e| EpiError::Storage(e.to_string()))?;
        Ok(())
    }

    pub fn compute_dataset_checksum(points: &[TimeSeriesPoint]) -> String {
        let mut hasher = Sha256::new();
        for p in points {
            hasher.update(p.date_label.as_bytes());
            hasher.update(b":");
            hasher.update(format!("{:.4};", p.count).as_bytes());
        }
        let digest = hasher.finalize();
        digest[..4].iter().map(|b| format!("{b:02x}")).collect()
    }

    pub fn commit_vintage(&mut self, dataset: &TimeSeriesDataset) -> Result<(), EpiError> {
        let tx = self
            .conn
            .transaction()
            .map_err(|e| EpiError::Storage(e.to_string()))?;

        tx.execute(
            "INSERT OR REPLACE INTO vintages
             (vintage_id, title, source_kind, source_url, ingested_at_utc, default_si, row_count, checksum_hex)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            params![
                dataset.vintage_id,
                dataset.title,
                dataset.source_kind.as_str(),
                dataset.source_url,
                dataset.ingested_at_utc,
                dataset.default_si_preset,
                dataset.points.len() as i64,
                dataset.checksum_hex,
            ],
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;

        tx.execute(
            "DELETE FROM vintage_points WHERE vintage_id = ?1",
            params![dataset.vintage_id],
        )
        .map_err(|e| EpiError::Storage(e.to_string()))?;

        {
            let mut stmt = tx
                .prepare(
                    "INSERT INTO vintage_points
                     (vintage_id, seq_idx, date_label, mmwr_year, mmwr_week, count_val, suppressed)
                     VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
                )
                .map_err(|e| EpiError::Storage(e.to_string()))?;

            for (idx, pt) in dataset.points.iter().enumerate() {
                let (yr, wk) = match pt.epiweek {
                    Some(w) => (Some(w.year), Some(i32::from(w.week))),
                    None => (None, None),
                };
                stmt.execute(params![
                    dataset.vintage_id,
                    idx as i64,
                    pt.date_label,
                    yr,
                    wk,
                    pt.count,
                    i32::from(pt.suppressed),
                ])
                .map_err(|e| EpiError::Storage(e.to_string()))?;
            }
        }

        tx.commit().map_err(|e| EpiError::Storage(e.to_string()))?;

        append_audit_event_on_conn(
            &self.conn,
            "epi-store",
            "COMMIT_VINTAGE",
            &dataset.vintage_id,
            &format!(
                "source={} rows={} sha256_32={}",
                dataset.source_kind.as_str(),
                dataset.points.len(),
                dataset.checksum_hex
            ),
        )?;

        Ok(())
    }

    pub fn list_vintages(&self) -> Result<Vec<VintageMetadata>, EpiError> {
        let mut stmt = self
            .conn
            .prepare(
                "SELECT vintage_id, title, source_kind, source_url, ingested_at_utc, row_count, checksum_hex
                 FROM vintages
                 ORDER BY ingested_at_utc DESC, vintage_id DESC",
            )
            .map_err(|e| EpiError::Storage(e.to_string()))?;

        let rows = stmt
            .query_map([], |row| {
                Ok(VintageMetadata {
                    vintage_id: row.get(0)?,
                    title: row.get(1)?,
                    source_kind: row.get(2)?,
                    source_url: row.get(3)?,
                    ingested_at_utc: row.get(4)?,
                    row_count: row.get::<_, i64>(5)? as usize,
                    checksum_hex: row.get(6)?,
                })
            })
            .map_err(|e| EpiError::Storage(e.to_string()))?;

        let mut out = Vec::new();
        for r in rows {
            out.push(r.map_err(|e| EpiError::Storage(e.to_string()))?);
        }
        Ok(out)
    }

    pub fn load_vintage(&self, vintage_id: &str) -> Result<Option<TimeSeriesDataset>, EpiError> {
        let meta = self
            .conn
            .query_row(
                "SELECT title, source_kind, source_url, ingested_at_utc, default_si, checksum_hex
                 FROM vintages WHERE vintage_id = ?1",
                params![vintage_id],
                |r| {
                    Ok((
                        r.get::<_, String>(0)?,
                        r.get::<_, String>(1)?,
                        r.get::<_, String>(2)?,
                        r.get::<_, String>(3)?,
                        r.get::<_, String>(4)?,
                        r.get::<_, String>(5)?,
                    ))
                },
            )
            .optional()
            .map_err(|e| EpiError::Storage(e.to_string()))?;

        let Some((title, kind_str, source_url, ingested_at_utc, default_si_preset, checksum_hex)) =
            meta
        else {
            return Ok(None);
        };

        let source_kind = SourceKind::from_label(&kind_str);

        let mut stmt = self
            .conn
            .prepare(
                "SELECT date_label, mmwr_year, mmwr_week, count_val, suppressed
                 FROM vintage_points
                 WHERE vintage_id = ?1
                 ORDER BY seq_idx ASC",
            )
            .map_err(|e| EpiError::Storage(e.to_string()))?;

        let pt_iter = stmt
            .query_map(params![vintage_id], |row| {
                let yr: Option<i32> = row.get(1)?;
                let wk: Option<i32> = row.get(2)?;
                let epiweek = match (yr, wk) {
                    (Some(y), Some(w)) if (1..=53).contains(&w) => MmwrWeek::new(y, w as u8),
                    _ => None,
                };
                Ok(TimeSeriesPoint {
                    date_label: row.get(0)?,
                    epiweek,
                    count: row.get(3)?,
                    suppressed: row.get::<_, i32>(4)? != 0,
                })
            })
            .map_err(|e| EpiError::Storage(e.to_string()))?;

        let mut points = Vec::new();
        for p in pt_iter {
            points.push(p.map_err(|e| EpiError::Storage(e.to_string()))?);
        }

        Ok(Some(TimeSeriesDataset {
            vintage_id: vintage_id.to_string(),
            title,
            source_kind,
            source_url,
            ingested_at_utc,
            default_si_preset,
            checksum_hex,
            points,
        }))
    }

    pub fn diff_vintages(
        &self,
        older_vintage_id: &str,
        newer_vintage_id: &str,
    ) -> Result<Vec<BackfillRevisionDelta>, EpiError> {
        let Some(older) = self.load_vintage(older_vintage_id)? else {
            return Err(EpiError::Schema(format!(
                "Older vintage '{older_vintage_id}' not found"
            )));
        };
        let Some(newer) = self.load_vintage(newer_vintage_id)? else {
            return Err(EpiError::Schema(format!(
                "Newer vintage '{newer_vintage_id}' not found"
            )));
        };

        let mut old_map = BTreeMap::new();
        for p in &older.points {
            old_map.insert(p.date_label.clone(), p.count);
        }

        let mut diffs = Vec::new();
        for p in &newer.points {
            if let Some(&old_val) = old_map.get(&p.date_label) {
                let delta = p.count - old_val;
                if delta.abs() > 1e-6 {
                    diffs.push(BackfillRevisionDelta {
                        date_label: p.date_label.clone(),
                        older_count: old_val,
                        newer_count: p.count,
                        delta,
                    });
                }
            }
        }
        Ok(diffs)
    }

    pub fn get_cached_http(
        &self,
        cache_key: &str,
        now_epoch_secs: i64,
    ) -> Result<Option<CachedHttpPayload>, EpiError> {
        let row = self
            .conn
            .query_row(
                "SELECT fetched_epoch_secs, ttl_secs, body_json
                 FROM http_lazy_cache WHERE cache_key = ?1",
                params![cache_key],
                |r| Ok((r.get::<_, i64>(0)?, r.get::<_, i64>(1)?, r.get::<_, String>(2)?)),
            )
            .optional()
            .map_err(|e| EpiError::Storage(e.to_string()))?;

        Ok(row.map(|(fetched, ttl, body_json)| CachedHttpPayload {
            cache_key: cache_key.to_string(),
            fetched_epoch_secs: fetched,
            ttl_secs: ttl,
            is_fresh: (now_epoch_secs - fetched) < ttl,
            body_json,
        }))
    }

    pub fn put_cached_http(
        &self,
        cache_key: &str,
        url: &str,
        now_epoch_secs: i64,
        ttl_secs: i64,
        etag: Option<&str>,
        body_json: &str,
    ) -> Result<(), EpiError> {
        self.conn
            .execute(
                "INSERT OR REPLACE INTO http_lazy_cache
                 (cache_key, url, fetched_epoch_secs, ttl_secs, etag, body_json)
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
                params![cache_key, url, now_epoch_secs, ttl_secs, etag, body_json],
            )
            .map_err(|e| EpiError::Storage(e.to_string()))?;
        Ok(())
    }

    pub fn prune_expired_http_cache(&self, older_than_epoch_secs: i64) -> Result<usize, EpiError> {
        let deleted = self
            .conn
            .execute(
                "DELETE FROM http_lazy_cache WHERE fetched_epoch_secs < ?1",
                params![older_than_epoch_secs],
            )
            .map_err(|e| EpiError::Storage(e.to_string()))?;
        Ok(deleted)
    }
}
