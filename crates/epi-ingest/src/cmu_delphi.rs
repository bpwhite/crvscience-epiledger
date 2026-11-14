use crate::cdc_soda::now_unix_secs;
use epi_core::{EpiError, MmwrWeek, SourceKind, TimeSeriesDataset, TimeSeriesPoint};
use epi_store::EpiStore;
use serde_json::Value;
use std::collections::BTreeMap;
use std::time::Duration;

/// Native client for the Carnegie Mellon University (CMU) Delphi Epidata API:
/// - FluView ILINet (`https://cmu-delphi.github.io/delphi-epidata/api/fluview.html`)
/// - COVIDcast V4 (`https://cmu-delphi.github.io/delphi-epidata/api/covidcast.html`)
pub struct CmuDelphiClient {
    api_key: Option<String>,
    timeout: Duration,
}

impl Default for CmuDelphiClient {
    fn default() -> Self {
        Self {
            api_key: None,
            timeout: Duration::from_secs(15),
        }
    }
}

impl CmuDelphiClient {
    pub fn new(api_key: Option<String>) -> Self {
        Self {
            api_key,
            timeout: Duration::from_secs(15),
        }
    }

    /// Fetch weekly outpatient Influenza-Like Illness (`num_ili` counts or `wili * 100` index)
    /// from CMU Delphi's `fluview` endpoint (`region` = `nat`, `hhs1`..`hhs10`, or state code like `or`, `ca`, `tx`).
    pub fn fetch_fluview_ilinet(
        &self,
        store: &mut EpiStore,
        region: &str,
        epiweeks_range: &str,
        use_ili_counts: bool,
        ttl_secs: i64,
    ) -> Result<TimeSeriesDataset, EpiError> {
        let reg = if region.trim().is_empty() {
            "nat".to_string()
        } else {
            region.trim().to_ascii_lowercase()
        };
        let wks = if epiweeks_range.trim().is_empty() {
            "202340-202645".to_string()
        } else {
            epiweeks_range.trim().to_string()
        };

        let mut url = format!(
            "https://api.delphi.cmu.edu/epidata/fluview/?regions={reg}&epiweeks={wks}"
        );
        if let Some(k) = &self.api_key {
            if !k.trim().is_empty() {
                url.push_str(&format!("&auth={}", k.trim()));
            }
        }

        let cache_key = format!("delphi_fluview:{reg}:{wks}:{use_ili_counts}");
        let now_secs = now_unix_secs();

        if let Some(cached) = store.get_cached_http(&cache_key, now_secs)? {
            if cached.is_fresh {
                return parse_fluview_json(&reg, &url, &cached.body_json, use_ili_counts);
            }
        }

        let agent = ureq::AgentBuilder::new()
            .timeout_connect(self.timeout)
            .timeout_read(self.timeout)
            .build();

        match agent.get(&url).call() {
            Ok(resp) => {
                let body = resp
                    .into_string()
                    .map_err(|e| EpiError::Network(e.to_string()))?;
                store.put_cached_http(&cache_key, &url, now_secs, ttl_secs, None, &body)?;
                let ds = parse_fluview_json(&reg, &url, &body, use_ili_counts)?;
                store.commit_vintage(&ds)?;
                Ok(ds)
            }
            Err(err) => {
                if let Some(stale) = store.get_cached_http(&cache_key, now_secs)? {
                    return parse_fluview_json(&reg, &url, &stale.body_json, use_ili_counts);
                }
                Err(EpiError::Network(format!(
                    "CMU Delphi FluView request failed: {err}"
                )))
            }
        }
    }

    /// Fetch curated indicator signals from CMU Delphi's `covidcast` V4 endpoint.
    pub fn fetch_covidcast_signal(
        &self,
        store: &mut EpiStore,
        data_source: &str,
        signal: &str,
        geo_type: &str,
        geo_value: &str,
        time_values: &str,
        ttl_secs: i64,
    ) -> Result<TimeSeriesDataset, EpiError> {
        let mut url = format!(
            "https://api.delphi.cmu.edu/epidata/covidcast/?data_source={}&signals={}&time_type=week&geo_type={}&time_values={}&geo_value={}",
            data_source.trim(),
            signal.trim(),
            geo_type.trim(),
            time_values.trim(),
            geo_value.trim()
        );
        if let Some(k) = &self.api_key {
            if !k.trim().is_empty() {
                url.push_str(&format!("&auth={}", k.trim()));
            }
        }

        let cache_key = format!("delphi_cc:{data_source}:{signal}:{geo_type}:{geo_value}:{time_values}");
        let now_secs = now_unix_secs();

        if let Some(cached) = store.get_cached_http(&cache_key, now_secs)? {
            if cached.is_fresh {
                return parse_covidcast_json(data_source, signal, geo_value, &url, &cached.body_json);
            }
        }

        let agent = ureq::AgentBuilder::new()
            .timeout_connect(self.timeout)
            .timeout_read(self.timeout)
            .build();

        let body = agent
            .get(&url)
            .call()
            .map_err(|e| EpiError::Network(format!("CMU Delphi COVIDcast error: {e}")))?
            .into_string()
            .map_err(|e| EpiError::Network(e.to_string()))?;

        store.put_cached_http(&cache_key, &url, now_secs, ttl_secs, None, &body)?;
        let ds = parse_covidcast_json(data_source, signal, geo_value, &url, &body)?;
        store.commit_vintage(&ds)?;
        Ok(ds)
    }
}

fn parse_fluview_json(
    region: &str,
    url: &str,
    body_json: &str,
    use_ili_counts: bool,
) -> Result<TimeSeriesDataset, EpiError> {
    let root: Value = serde_json::from_str(body_json)?;
    let arr = root
        .get("epidata")
        .and_then(Value::as_array)
        .ok_or_else(|| EpiError::Schema("Delphi FluView JSON missing `epidata` array".to_string()))?;

    let mut buckets: BTreeMap<MmwrWeek, f64> = BTreeMap::new();
    for item in arr {
        let ew_int = item.get("epiweek").and_then(Value::as_i64).unwrap_or(0);
        let yr = (ew_int / 100) as i32;
        let wk = (ew_int % 100) as u8;
        let Some(mmwr) = MmwrWeek::new(yr, wk) else {
            continue;
        };

        let val = if use_ili_counts {
            item.get("num_ili").and_then(Value::as_f64).unwrap_or(0.0)
        } else {
            // Scale weighted % ILI by 100 for integer-like visual resolution if requested
            item.get("wili").and_then(Value::as_f64).unwrap_or(0.0) * 100.0
        };
        buckets.insert(mmwr, val.max(0.0));
    }

    if buckets.is_empty() {
        return Err(EpiError::Schema(
            "CMU Delphi FluView returned 0 weekly observations".to_string(),
        ));
    }

    let points: Vec<TimeSeriesPoint> = buckets
        .into_iter()
        .map(|(wk, count)| TimeSeriesPoint {
            date_label: wk.compact_label(),
            epiweek: Some(wk),
            count,
            suppressed: false,
        })
        .collect();

    let checksum_hex = EpiStore::compute_dataset_checksum(&points);
    Ok(TimeSeriesDataset {
        vintage_id: format!("delphi_fluview_{region}_{checksum_hex}"),
        title: format!(
            "CMU Delphi FluView ILINet ({}, {})",
            region.to_uppercase(),
            if use_ili_counts { "ILI Cases" } else { "wILI x100" }
        ),
        source_kind: SourceKind::LiveCmuDelphi,
        source_url: url.to_string(),
        ingested_at_utc: format!("epoch:{}", now_unix_secs()),
        default_si_preset: "influenza".to_string(),
        checksum_hex,
        points,
    })
}

fn parse_covidcast_json(
    data_source: &str,
    signal: &str,
    geo_value: &str,
    url: &str,
    body_json: &str,
) -> Result<TimeSeriesDataset, EpiError> {
    let root: Value = serde_json::from_str(body_json)?;
    let arr = root
        .get("epidata")
        .and_then(Value::as_array)
        .ok_or_else(|| EpiError::Schema("Delphi COVIDcast JSON missing `epidata` array".to_string()))?;

    let mut buckets: BTreeMap<String, (f64, Option<MmwrWeek>)> = BTreeMap::new();
    for item in arr {
        let tv = item
            .get("time_value")
            .map(|v| match v {
                Value::Number(n) => n.to_string(),
                Value::String(s) => s.clone(),
                _ => String::new(),
            })
            .unwrap_or_default();
        let val = item.get("value").and_then(Value::as_f64).unwrap_or(0.0);
        if tv.len() == 6 {
            if let (Ok(yr), Ok(wk)) = (tv[0..4].parse::<i32>(), tv[4..6].parse::<u8>()) {
                if let Some(mmwr) = MmwrWeek::new(yr, wk) {
                    buckets.insert(mmwr.compact_label(), (val.max(0.0), Some(mmwr)));
                    continue;
                }
            }
        }
        if !tv.is_empty() {
            buckets.insert(tv, (val.max(0.0), None));
        }
    }

    if buckets.is_empty() {
        return Err(EpiError::Schema(
            "CMU Delphi COVIDcast returned 0 signal observations".to_string(),
        ));
    }

    let points: Vec<TimeSeriesPoint> = buckets
        .into_iter()
        .map(|(date_label, (count, epiweek))| TimeSeriesPoint {
            date_label,
            epiweek,
            count,
            suppressed: false,
        })
        .collect();

    let checksum_hex = EpiStore::compute_dataset_checksum(&points);
    Ok(TimeSeriesDataset {
        vintage_id: format!("delphi_cc_{signal}_{geo_value}_{checksum_hex}"),
        title: format!("CMU Delphi {data_source}:{signal} ({geo_value})"),
        source_kind: SourceKind::LiveCmuDelphi,
        source_url: url.to_string(),
        ingested_at_utc: format!("epoch:{}", now_unix_secs()),
        default_si_preset: "sarscov2".to_string(),
        checksum_hex,
        points,
    })
}
