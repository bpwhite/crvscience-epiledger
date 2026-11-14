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
        let url = self.build_query_url(kind, condition_filter, jurisdiction_filter, limit);
        let cache_key = format!("soda:{url}");
        let now_secs = now_unix_secs();

        if let Some(cached) = store.get_cached_http(&cache_key, now_secs)? {
            if cached.is_fresh {
                return parse_soda_payload(
                    kind,
                    condition_filter,
                    &url,
                    &cached.body_json,
                    false,
                );
            }
        }

        let agent = ureq::AgentBuilder::new()
            .timeout_connect(self.timeout)
            .timeout_read(self.timeout)
            .build();

        let mut req = agent.get(&url);
        if let Some(tok) = &self.app_token {
            if !tok.trim().is_empty() {
                req = req.set("X-App-Token", tok.trim());
            }
        }

        match req.call() {
            Ok(resp) => {
                let etag = resp.header("ETag").map(str::to_string);
                let body = resp
                    .into_string()
                    .map_err(|e| EpiError::Network(e.to_string()))?;
                store.put_cached_http(
                    &cache_key,
                    &url,
                    now_secs,
                    ttl_secs,
                    etag.as_deref(),
                    &body,
                )?;
                let outcome = parse_soda_payload(kind, condition_filter, &url, &body, false)?;
                store.commit_vintage(&outcome.dataset)?;
                Ok(outcome)
            }
            Err(net_err) => {
                if let Some(stale) = store.get_cached_http(&cache_key, now_secs)? {
                    return parse_soda_payload(
                        kind,
                        condition_filter,
                        &url,
                        &stale.body_json,
                        true,
                    );
                }
                Err(EpiError::Network(format!(
                    "CDC SODA request failed ({url}): {net_err}"
                )))
            }
        }
    }
}

fn parse_soda_payload(
    kind: SodaDatasetKind,
    condition_filter: &str,
    url: &str,
    body_json: &str,
    used_stale_cache: bool,
) -> Result<SodaFetchOutcome, EpiError> {
    let rows: Vec<Value> = serde_json::from_str(body_json)?;
    if rows.is_empty() {
        return Err(EpiError::Schema(
            "CDC SODA returned 0 matching rows for query".to_string(),
        ));
    }

    let state_lookup = build_state_name_lookup();
    let mut week_buckets: BTreeMap<String, (f64, usize, Option<MmwrWeek>)> = BTreeMap::new();
    let mut state_buckets: HashMap<String, (f64, usize)> = HashMap::new();

    for row in &rows {
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
        for row in &rows {
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
