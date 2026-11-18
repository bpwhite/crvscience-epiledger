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

use epi_core::{
    builtin_pertussis_spatial_records, EpiError, MmwrWeek, SourceKind, TimeSeriesDataset,
    TimeSeriesPoint,
};
use epi_store::EpiStore;
use serde_json::Value;
use std::collections::{BTreeMap, HashMap};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SodaDatasetKind {
    /// CDC NNDSS Weekly Data (`x9gk-5huc`)
    NndssWeekly,
    /// CDC NWSS Public SARS-CoV-2 Wastewater Metric Data (`2ew6-ywp6`)
    NwssWastewater,
}

impl SodaDatasetKind {
    pub fn resource_url(&self) -> &'static str {
        match self {
            Self::NndssWeekly => "https://data.cdc.gov/resource/x9gk-5huc.json",
            Self::NwssWastewater => "https://data.cdc.gov/resource/2ew6-ywp6.json",
        }
    }
}

#[derive(Debug, Clone)]
pub struct SodaFetchOutcome {
    pub dataset: TimeSeriesDataset,
    pub state_values_by_fips: HashMap<String, f64>,
    pub raw_row_count: usize,
    pub used_stale_cache: bool,
}

#[derive(Debug, Clone)]
pub struct SodaQuerySpec {
    pub condition_filter: String,
    pub jurisdiction_filter: String,
    pub start_date_or_year: Option<String>,
    pub end_date_or_year: Option<String>,
    pub limit: usize,
    pub offset: usize,
    pub max_pages: usize,
}

impl Default for SodaQuerySpec {
    fn default() -> Self {
        Self {
            condition_filter: String::new(),
            jurisdiction_filter: String::new(),
            start_date_or_year: None,
            end_date_or_year: None,
            limit: 1000,
            offset: 0,
            max_pages: 5,
        }
    }
}

/// Escape single quotes and strip control characters for safe SoQL string literals.
pub fn escape_soql_literal(s: &str) -> String {
    s.chars()
        .filter(|c| !c.is_control())
        .flat_map(|c| {
            if c == '\'' {
                vec!['\'', '\'']
            } else {
                vec![c]
            }
        })
        .collect()
}

fn build_where_clause(kind: SodaDatasetKind, spec: &SodaQuerySpec) -> Option<String> {
    let mut clauses: Vec<String> = Vec::new();
    let cond = spec.condition_filter.trim();
    let jur = spec.jurisdiction_filter.trim();

    match kind {
        SodaDatasetKind::NndssWeekly => {
            if !cond.is_empty() {
                let esc = escape_soql_literal(&cond.to_ascii_uppercase());
                clauses.push(format!("upper(label) like '%{esc}%'"));
            }
            if !jur.is_empty() {
                let records = builtin_pertussis_spatial_records();
                if let Some(rec) = records.iter().find(|r| {
                    r.name.eq_ignore_ascii_case(jur) || r.abbr.eq_ignore_ascii_case(jur)
                }) {
                    clauses.push(format!(
                        "(upper(states)='{}' OR upper(states)='{}')",
                        escape_soql_literal(&rec.name.to_ascii_uppercase()),
                        escape_soql_literal(&rec.abbr.to_ascii_uppercase()),
                    ));
                } else {
                    let esc = escape_soql_literal(&jur.to_ascii_uppercase());
                    clauses.push(format!("upper(states) like '%{esc}%'"));
                }
            }
            if let Some(ref start) = spec.start_date_or_year {
                let st = start.trim();
                if !st.is_empty() {
                    let yr = st
                        .get(0..4)
                        .filter(|p| p.chars().all(|c| c.is_ascii_digit()))
                        .unwrap_or(st);
                    clauses.push(format!("year >= '{}'", escape_soql_literal(yr)));
                }
            }
            if let Some(ref end) = spec.end_date_or_year {
                let ed = end.trim();
                if !ed.is_empty() {
                    let yr = ed
                        .get(0..4)
                        .filter(|p| p.chars().all(|c| c.is_ascii_digit()))
                        .unwrap_or(ed);
                    clauses.push(format!("year <= '{}'", escape_soql_literal(yr)));
                }
            }
        }
        SodaDatasetKind::NwssWastewater => {
            if !jur.is_empty() {
                let records = builtin_pertussis_spatial_records();
                if let Some(rec) = records.iter().find(|r| {
                    r.name.eq_ignore_ascii_case(jur) || r.abbr.eq_ignore_ascii_case(jur)
                }) {
                    clauses.push(format!(
                        "(upper(wwtp_jurisdiction)='{}' OR upper(wwtp_jurisdiction)='{}')",
                        escape_soql_literal(&rec.name.to_ascii_uppercase()),
                        escape_soql_literal(&rec.abbr.to_ascii_uppercase()),
                    ));
                } else {
                    clauses.push(format!(
                        "upper(wwtp_jurisdiction)='{}'",
                        escape_soql_literal(&jur.to_ascii_uppercase())
                    ));
                }
            }
            if let Some(ref start) = spec.start_date_or_year {
                let st = start.trim();
                if !st.is_empty() {
                    clauses.push(format!("date_end >= '{}'", escape_soql_literal(st)));
                }
            }
            if let Some(ref end) = spec.end_date_or_year {
                let ed = end.trim();
                if !ed.is_empty() {
                    clauses.push(format!("date_end <= '{}'", escape_soql_literal(ed)));
                }
            }
        }
    }

    if clauses.is_empty() {
        None
    } else {
        Some(clauses.join(" AND "))
    }
}

pub struct CdcSodaClient {
    app_token: Option<String>,
    timeout: Duration,
}

impl Default for CdcSodaClient {
    fn default() -> Self {
        Self {
            app_token: None,
            timeout: Duration::from_secs(15),
        }
    }
}

impl CdcSodaClient {
    pub fn new(app_token: Option<String>) -> Self {
        Self {
            app_token,
            timeout: Duration::from_secs(15),
        }
    }

    pub fn build_query_url_with_spec(&self, kind: SodaDatasetKind, spec: &SodaQuerySpec) -> String {
        let limit = spec.limit.clamp(10, 5000);
        let mut params = vec![
            format!("$limit={limit}"),
            format!("$offset={}", spec.offset),
        ];
        if let Some(where_sql) = build_where_clause(kind, spec) {
            params.push(format!("$where={}", percent_encode_query(&where_sql)));
        }
        let order = match kind {
            SodaDatasetKind::NndssWeekly => "year DESC,week DESC,states ASC,:id",
            SodaDatasetKind::NwssWastewater => "date_end DESC,:id",
        };
        params.push(format!("$order={}", percent_encode_query(order)));
        format!("{}?{}", kind.resource_url(), params.join("&"))
    }

    pub fn build_query_url(
        &self,
        kind: SodaDatasetKind,
        condition_filter: &str,
        jurisdiction_filter: &str,
        limit: usize,
    ) -> String {
        let mut params = vec![format!("$limit={}", limit.clamp(10, 5000))];
        let q_combined = [condition_filter.trim(), jurisdiction_filter.trim()]
            .into_iter()
            .filter(|s| !s.is_empty())
            .collect::<Vec<_>>()
            .join(" ");
        if !q_combined.is_empty() {
            params.push(format!("$q={}", percent_encode_query(&q_combined)));
        }
        if kind == SodaDatasetKind::NwssWastewater {
            params.push("$order=date_end%20DESC".to_string());
        }
        format!("{}?{}", kind.resource_url(), params.join("&"))
    }

    /// Fetch from live CDC SODA with cron-free SQLite lazy caching (`ttl_secs`) and stale fallback on HTTP error.
    pub fn fetch_with_lazy_cache(
        &self,
        store: &mut EpiStore,
        kind: SodaDatasetKind,
        condition_filter: &str,
        jurisdiction_filter: &str,
        limit: usize,
        ttl_secs: i64,
    ) -> Result<SodaFetchOutcome, EpiError> {
        let spec = SodaQuerySpec {
            condition_filter: condition_filter.to_string(),
            jurisdiction_filter: jurisdiction_filter.to_string(),
            start_date_or_year: None,
            end_date_or_year: None,
            limit,
            offset: 0,
            max_pages: 1,
        };
        self.fetch_paginated_with_lazy_cache(store, kind, &spec, ttl_secs)
    }

    /// Fetch multi-page CDC SODA records using deterministic `$order` + `$offset` pagination
    /// and per-page SQLite lazy caching (`ttl_secs`) with stale fallback on HTTP errors.
    /// If any page fails without a stale cache entry, returns `Err(EpiError::Network)` rather than
    /// silently truncating the dataset (already-fetched pages remain cached in `http_lazy_cache`).
    pub fn fetch_paginated_with_lazy_cache(
        &self,
        store: &mut EpiStore,
        kind: SodaDatasetKind,
        spec: &SodaQuerySpec,
        ttl_secs: i64,
    ) -> Result<SodaFetchOutcome, EpiError> {
        let agent = ureq::AgentBuilder::new()
            .timeout_connect(self.timeout)
            .timeout_read(self.timeout)
            .build();
        self.fetch_paginated_with_transport(store, kind, spec, ttl_secs, |url| {
            let mut req = agent.get(url);
            if let Some(tok) = &self.app_token {
                if !tok.trim().is_empty() {
                    req = req.set("X-App-Token", tok.trim());
                }
            }
            let resp = req.call().map_err(|e| format!("CDC SODA request failed ({url}): {e}"))?;
            let etag = resp.header("ETag").map(str::to_string);
            let body = resp.into_string().map_err(|e| e.to_string())?;
            Ok((body, etag))
        })
    }

    pub(crate) fn fetch_paginated_with_transport<F>(
        &self,
        store: &mut EpiStore,
        kind: SodaDatasetKind,
        spec: &SodaQuerySpec,
        ttl_secs: i64,
        mut transport: F,
    ) -> Result<SodaFetchOutcome, EpiError>
    where
        F: FnMut(&str) -> Result<(String, Option<String>), String>,
    {
        let page_size = spec.limit.clamp(10, 5000);
        let max_pages = spec.max_pages.clamp(1, 25);
        let now_secs = now_unix_secs();

        let mut all_rows: Vec<Value> = Vec::new();
        let mut first_url = String::new();
        let mut used_stale_any = false;

        for page in 0..max_pages {
            let mut page_spec = spec.clone();
            page_spec.limit = page_size;
            page_spec.offset = spec.offset + page * page_size;
            let url = self.build_query_url_with_spec(kind, &page_spec);
            if page == 0 {
                first_url = url.clone();
            }
            let cache_key = format!("soda:{url}");

            let page_rows: Vec<Value> = if let Some(cached) =
                store.get_cached_http(&cache_key, now_secs)?.filter(|c| c.is_fresh)
            {
                serde_json::from_str(&cached.body_json)?
            } else {
                match transport(&url) {
                    Ok((body, etag)) => {
                        let parsed: Vec<Value> = serde_json::from_str(&body)?;
                        store.put_cached_http(
                            &cache_key,
                            &url,
                            now_secs,
                            ttl_secs,
                            etag.as_deref(),
                            &body,
                        )?;
                        parsed
                    }
                    Err(net_err) => {
                        if let Some(stale) = store.get_cached_http(&cache_key, now_secs)? {
                            used_stale_any = true;
                            serde_json::from_str(&stale.body_json)?
                        } else if page > 0 {
                            return Err(EpiError::Network(format!(
                                "CDC SODA multi-page fetch interrupted on page {} without stale cache: {net_err}",
                                page + 1
                            )));
                        } else {
                            return Err(EpiError::Network(net_err));
                        }
                    }
                }
            };

            let fetched_len = page_rows.len();
            all_rows.extend(page_rows);
            if fetched_len < page_size || all_rows.len() >= 50_000 {
                break;
            }
        }

        let outcome = parse_soda_rows(
            kind,
            &spec.condition_filter,
            &first_url,
            &all_rows,
            used_stale_any,
        )?;
        store.commit_vintage(&outcome.dataset)?;
        Ok(outcome)
    }
}

fn parse_soda_rows(
    kind: SodaDatasetKind,
    condition_filter: &str,
    url: &str,
    rows: &[Value],
    used_stale_cache: bool,
) -> Result<SodaFetchOutcome, EpiError> {
    if rows.is_empty() {
        return Err(EpiError::Schema(
            "CDC SODA returned 0 matching rows for query".to_string(),
        ));
    }

    let state_lookup = build_state_name_lookup();
    let mut week_buckets: BTreeMap<String, (f64, usize, Option<MmwrWeek>)> = BTreeMap::new();
    let mut state_buckets: HashMap<String, (f64, usize)> = HashMap::new();

    for row in rows {
        let Some(obj) = row.as_object() else {
            continue;
        };

        match kind {
            SodaDatasetKind::NndssWeekly => {
                let yr_opt = extract_i32(obj, &["year", "mmwr_year", "current_mmwr_year"]);
                let wk_opt = extract_i32(obj, &["week", "mmwr_week"]);
                let area = extract_str(obj, &["states", "reporting_area"]).to_uppercase();
                let flag = extract_str(obj, &["m1_flag", "current_week_flag"]);
                let count_val = extract_f64(obj, &["m1", "current_week", "cases"]).unwrap_or(0.0);

                if let (Some(yr), Some(wk)) = (yr_opt, wk_opt) {
                    if (1..=53).contains(&wk)
                        && area != "TOTAL"
                        && area != "TOTALS"
                        && area != "US RESIDENTS"
                        && flag.trim().is_empty()
                    {
                        let mmwr = MmwrWeek::new(yr, wk as u8);
                        let key = format!("{yr:04}-W{wk:02}");
                        let entry = week_buckets.entry(key).or_insert((0.0, 0, mmwr));
                        entry.0 += count_val;
                        entry.1 += 1;
                    }
                }

                if let Some(fips) = state_lookup.get(&area) {
                    let st_entry = state_buckets.entry(fips.clone()).or_insert((0.0, 0));
                    st_entry.0 += count_val;
                    st_entry.1 += 1;
                }
            }
            SodaDatasetKind::NwssWastewater => {
                let dt = extract_str(obj, &["date_end", "date_start"]);
                let dt_short = dt.chars().take(10).collect::<String>();
                let pct = extract_f64(obj, &["percentile", "detect_prop_15d"]).unwrap_or(0.0);
                if !dt_short.is_empty() && (0.0..=100.0).contains(&pct) {
                    let mmwr = epi_core::parse_flexible_date_to_mmwr(&dt_short);
                    let entry = week_buckets.entry(dt_short).or_insert((0.0, 0, mmwr));
                    entry.0 += pct;
                    entry.1 += 1;
                }
                let jur = extract_str(obj, &["wwtp_jurisdiction"]).to_uppercase();
                if let Some(fips) = state_lookup.get(&jur) {
                    if (0.0..=100.0).contains(&pct) {
                        let st_entry = state_buckets.entry(fips.clone()).or_insert((0.0, 0));
                        st_entry.0 += pct;
                        st_entry.1 += 1;
                    }
                }
            }
        }
    }

    // Fallback for NNDSS queries filtered solely to TOTAL national summary rows
    if kind == SodaDatasetKind::NndssWeekly && week_buckets.is_empty() {
        for row in rows {
            if let Some(obj) = row.as_object() {
                let yr_opt = extract_i32(obj, &["year", "mmwr_year"]);
                let wk_opt = extract_i32(obj, &["week", "mmwr_week"]);
                let count_val = extract_f64(obj, &["m1", "current_week"]).unwrap_or(0.0);
                if let (Some(yr), Some(wk)) = (yr_opt, wk_opt) {
                    if (1..=53).contains(&wk) {
                        let mmwr = MmwrWeek::new(yr, wk as u8);
                        let key = format!("{yr:04}-W{wk:02}");
                        let entry = week_buckets.entry(key).or_insert((0.0, 0, mmwr));
                        entry.0 += count_val;
                        entry.1 += 1;
                    }
                }
            }
        }
    }

    let points: Vec<TimeSeriesPoint> = week_buckets
        .into_iter()
        .map(|(date_label, (sum, n, epiweek))| {
            let count = match kind {
                SodaDatasetKind::NndssWeekly => sum,
                SodaDatasetKind::NwssWastewater => sum / (n.max(1) as f64),
            };
            TimeSeriesPoint {
                date_label,
                epiweek,
                count,
                suppressed: false,
            }
        })
        .collect();

    let mut state_values_by_fips = HashMap::new();
    for (fips, (sum, n)) in state_buckets {
        let val = match kind {
            SodaDatasetKind::NndssWeekly => sum,
            SodaDatasetKind::NwssWastewater => sum / (n.max(1) as f64),
        };
        state_values_by_fips.insert(fips, val);
    }

    let checksum_hex = EpiStore::compute_dataset_checksum(&points);
    let label_tag = if condition_filter.trim().is_empty() {
        match kind {
            SodaDatasetKind::NndssWeekly => "NNDSS All",
            SodaDatasetKind::NwssWastewater => "NWSS SARS-CoV-2 Percentile",
        }
    } else {
        condition_filter.trim()
    };

    let default_si = if label_tag.to_ascii_lowercase().contains("measles") {
        "measles"
    } else if label_tag.to_ascii_lowercase().contains("pertussis") {
        "pertussis"
    } else {
        "sarscov2"
    };

    let dataset = TimeSeriesDataset {
        vintage_id: format!("soda_{}_{}", now_unix_secs(), checksum_hex),
        title: format!("CDC SODA Live: {label_tag}"),
        source_kind: SourceKind::LiveCdcSoda,
        source_url: url.to_string(),
        ingested_at_utc: format!("epoch:{}", now_unix_secs()),
        default_si_preset: default_si.to_string(),
        checksum_hex,
        points,
    };

    Ok(SodaFetchOutcome {
        dataset,
        state_values_by_fips,
        raw_row_count: rows.len(),
        used_stale_cache,
    })
}

fn build_state_name_lookup() -> HashMap<String, String> {
    let mut map = HashMap::new();
    for rec in builtin_pertussis_spatial_records() {
        map.insert(rec.name.to_uppercase(), rec.fips.clone());
        map.insert(rec.abbr.to_uppercase(), rec.fips.clone());
    }
    map.insert("D.C.".to_string(), "11".to_string());
    map.insert("NEW YORK CITY".to_string(), "36".to_string());
    map.insert("NEW YORK STATE".to_string(), "36".to_string());
    map
}

fn extract_str(obj: &serde_json::Map<String, Value>, keys: &[&str]) -> String {
    for &k in keys {
        if let Some(v) = obj.get(k) {
            if let Some(s) = v.as_str() {
                return s.trim().to_string();
            }
        }
    }
    String::new()
}

fn extract_i32(obj: &serde_json::Map<String, Value>, keys: &[&str]) -> Option<i32> {
    for &k in keys {
        if let Some(v) = obj.get(k) {
            if let Some(n) = v.as_i64() {
                return Some(n as i32);
            }
            if let Some(s) = v.as_str() {
                if let Ok(parsed) = s.trim().parse::<i32>() {
                    return Some(parsed);
                }
            }
        }
    }
    None
}

fn extract_f64(obj: &serde_json::Map<String, Value>, keys: &[&str]) -> Option<f64> {
    for &k in keys {
        if let Some(v) = obj.get(k) {
            if let Some(n) = v.as_f64() {
                return Some(n);
            }
            if let Some(s) = v.as_str() {
                if let Ok(parsed) = s.trim().replace(',', "").parse::<f64>() {
                    return Some(parsed);
                }
            }
        }
    }
    None
}

fn percent_encode_query(s: &str) -> String {
    let mut out = String::with_capacity(s.len() * 2);
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' => out.push(b as char),
            b' ' => out.push_str("%20"),
            _ => out.push_str(&format!("%{b:02X}")),
        }
    }
    out
}

pub(crate) fn now_unix_secs() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(1_794_550_000)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_soql_where_escaping_and_url_builder() {
        assert_eq!(escape_soql_literal("O'Brien\n"), "O''Brien");

        let client = CdcSodaClient::default();
        let spec_nndss = SodaQuerySpec {
            condition_filter: "Measles' OR '1'='1".to_string(),
            jurisdiction_filter: "Oregon".to_string(),
            start_date_or_year: Some("2024-01-01".to_string()),
            end_date_or_year: Some("2026".to_string()),
            limit: 500,
            offset: 1000,
            max_pages: 3,
        };
        let url_nndss = client.build_query_url_with_spec(SodaDatasetKind::NndssWeekly, &spec_nndss);
        assert!(url_nndss.contains("$limit=500"));
        assert!(url_nndss.contains("$offset=1000"));
        assert!(url_nndss.contains("$where="));
        assert!(url_nndss.contains("O%27%27R") || url_nndss.contains("MEASLES%27%27%20OR%20%27%271%27%27%3D%27%271"));
        assert!(url_nndss.contains("year%20%3E%3D%20%272024%27"));
        assert!(url_nndss.contains("year%20%3C%3D%20%272026%27"));
        assert!(url_nndss.contains("$order=year%20DESC%2Cweek%20DESC%2Cstates%20ASC%2C%3Aid"));

        let spec_nwss = SodaQuerySpec {
            condition_filter: String::new(),
            jurisdiction_filter: "OR".to_string(),
            start_date_or_year: Some("2025-01-01".to_string()),
            end_date_or_year: None,
            limit: 250,
            offset: 0,
            max_pages: 2,
        };
        let url_nwss = client.build_query_url_with_spec(SodaDatasetKind::NwssWastewater, &spec_nwss);
        assert!(url_nwss.contains("OREGON"));
        assert!(url_nwss.contains("%27OR%27"));
        assert!(url_nwss.contains("date_end%20%3E%3D%20%272025-01-01%27"));
        assert!(url_nwss.contains("$order=date_end%20DESC%2C%3Aid"));
    }

    #[test]
    fn test_soda_paginated_fetch_mock_transport_cache_and_cap() {
        let mut store = EpiStore::open_in_memory().unwrap();
        let client = CdcSodaClient::default();
        let spec = SodaQuerySpec {
            condition_filter: "Measles".to_string(),
            jurisdiction_filter: String::new(),
            start_date_or_year: Some("2025".to_string()),
            end_date_or_year: None,
            limit: 10,
            offset: 0,
            max_pages: 5,
        };

        let make_row = |yr: &str, wk: u32, st: &str, m1: u32| -> String {
            format!(r#"{{"year":"{yr}","week":"{wk}","states":"{st}","m1":"{m1}","m1_flag":""}}"#)
        };

        let mut call_urls = Vec::new();
        let outcome = client
            .fetch_paginated_with_transport(&mut store, SodaDatasetKind::NndssWeekly, &spec, 3600, |url| {
                call_urls.push(url.to_string());
                if url.contains("$offset=0") {
                    let rows: Vec<String> = (1..=10).map(|w| make_row("2025", w, "OREGON", w)).collect();
                    Ok((format!("[{}]", rows.join(",")), None))
                } else if url.contains("$offset=10") {
                    let rows: Vec<String> = (11..=13).map(|w| make_row("2025", w, "OREGON", w)).collect();
                    Ok((format!("[{}]", rows.join(",")), None))
                } else {
                    Err("unexpected third page".to_string())
                }
            })
            .unwrap();

        assert_eq!(call_urls.len(), 2, "must fetch 2 pages and stop on short page");
        assert_eq!(outcome.raw_row_count, 13);
        assert_eq!(outcome.dataset.points.len(), 13);
        assert!(!outcome.used_stale_cache);

        // Fresh cache hit on repeat call -> 0 transport calls
        let mut hit_calls = 0usize;
        let cached_outcome = client
            .fetch_paginated_with_transport(&mut store, SodaDatasetKind::NndssWeekly, &spec, 3600, |_| {
                hit_calls += 1;
                Err("should not hit network on fresh cache".to_string())
            })
            .unwrap();
        assert_eq!(hit_calls, 0);
        assert_eq!(cached_outcome.raw_row_count, 13);

        // max_pages cap test (stops after 2 full pages of 10)
        let cap_spec = SodaQuerySpec {
            limit: 10,
            max_pages: 2,
            condition_filter: "Pertussis".to_string(),
            ..Default::default()
        };
        let mut cap_calls = 0usize;
        let cap_outcome = client
            .fetch_paginated_with_transport(&mut store, SodaDatasetKind::NndssWeekly, &cap_spec, 3600, |_| {
                cap_calls += 1;
                let rows: Vec<String> = (1..=10).map(|w| make_row("2025", w, "TEXAS", 5)).collect();
                Ok((format!("[{}]", rows.join(",")), None))
            })
            .unwrap();
        assert_eq!(cap_calls, 2);
        assert_eq!(cap_outcome.raw_row_count, 20);
    }

    #[test]
    fn test_nndss_state_abbreviation_expansion_and_partial_page_error() {
        let client = CdcSodaClient::default();

        // 2-letter postal abbreviation "OR" must expand to exact Oregon OR-clause, never LIKE '%OR%'
        let or_spec = SodaQuerySpec {
            condition_filter: "Pertussis".to_string(),
            jurisdiction_filter: "OR".to_string(),
            ..Default::default()
        };
        let url_or = client.build_query_url_with_spec(SodaDatasetKind::NndssWeekly, &or_spec);
        assert!(url_or.contains("upper%28states%29%3D%27OREGON%27"));
        assert!(url_or.contains("upper%28states%29%3D%27OR%27"));
        assert!(!url_or.contains("%25OR%25"), "must not emit LIKE '%OR%'");

        let dc_spec = SodaQuerySpec {
            jurisdiction_filter: "DC".to_string(),
            ..Default::default()
        };
        let url_dc = client.build_query_url_with_spec(SodaDatasetKind::NndssWeekly, &dc_spec);
        assert!(url_dc.contains("DISTRICT%20OF%20COLUMBIA"));
        assert!(url_dc.contains("%27DC%27"));

        // Non-state multi-region fallback keeps LIKE '%...%'
        let region_spec = SodaQuerySpec {
            jurisdiction_filter: "NEW ENGLAND".to_string(),
            ..Default::default()
        };
        let url_reg = client.build_query_url_with_spec(SodaDatasetKind::NndssWeekly, &region_spec);
        assert!(url_reg.contains("upper%28states%29%20like%20%27%25NEW%20ENGLAND%25%27"));

        // Multi-page partial failure without stale cache must return Err instead of silent truncation
        let mut store = EpiStore::open_in_memory().unwrap();
        let page_spec = SodaQuerySpec {
            condition_filter: "Measles".to_string(),
            jurisdiction_filter: "OR".to_string(),
            limit: 10,
            offset: 0,
            max_pages: 3,
            ..Default::default()
        };
        let make_row = |wk: u32| -> String {
            format!(r#"{{"year":"2025","week":"{wk}","states":"OREGON","m1":"3","m1_flag":""}}"#)
        };
        let err_res = client.fetch_paginated_with_transport(
            &mut store,
            SodaDatasetKind::NndssWeekly,
            &page_spec,
            3600,
            |url| {
                if url.contains("$offset=0") {
                    let rows: Vec<String> = (1..=10).map(make_row).collect();
                    Ok((format!("[{}]", rows.join(",")), None))
                } else {
                    Err("simulated page 2 timeout".to_string())
                }
            },
        );
        assert!(matches!(err_res, Err(EpiError::Network(_))));
    }
}
