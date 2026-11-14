use crate::cdc_soda::now_unix_secs;
use epi_core::{EpiError, SourceKind, TimeSeriesDataset, TimeSeriesPoint};
use epi_store::EpiStore;
use serde_json::Value;
use std::collections::BTreeMap;
use std::time::Duration;

/// Keyless client for the World Health Organization Global Health Observatory (GHO) OData API.
/// Documentation: https://www.who.int/data/gho/info/gho-odata-api
pub struct WhoGhoClient {
    timeout: Duration,
}

impl Default for WhoGhoClient {
    fn default() -> Self {
        Self {
            timeout: Duration::from_secs(15),
        }
    }
}

impl WhoGhoClient {
    /// Fetch an annual or periodic WHO GHO indicator (e.g. `WHS3_62` Measles reported cases) for an ISO3 country (`USA`, `DEU`, etc.).
    pub fn fetch_indicator_series(
        &self,
        store: &mut EpiStore,
        indicator_code: &str,
        iso3_country: &str,
        ttl_secs: i64,
    ) -> Result<TimeSeriesDataset, EpiError> {
        let ind: String = indicator_code
            .trim()
            .chars()
            .filter(|c| c.is_ascii_alphanumeric() || *c == '_')
            .collect();
        let iso: String = iso3_country
            .trim()
            .to_ascii_uppercase()
            .chars()
            .filter(|c| c.is_ascii_alphabetic())
            .collect();

        let url = if iso.is_empty() {
            format!("https://ghoapi.azureedge.net/api/{ind}")
        } else {
            format!("https://ghoapi.azureedge.net/api/{ind}?$filter=SpatialDim%20eq%20'{iso}'")
        };

        let cache_key = format!("who_gho:{ind}:{iso}");
        let now_secs = now_unix_secs();

        if let Some(cached) = store.get_cached_http(&cache_key, now_secs)? {
            if cached.is_fresh {
                return parse_who_odata(&ind, &iso, &url, &cached.body_json);
            }
        }

        let agent = ureq::AgentBuilder::new()
            .timeout_connect(self.timeout)
            .timeout_read(self.timeout)
            .build();

        let body = agent
            .get(&url)
            .call()
            .map_err(|e| EpiError::Network(format!("WHO GHO OData request failed: {e}")))?
            .into_string()
            .map_err(|e| EpiError::Network(e.to_string()))?;

        store.put_cached_http(&cache_key, &url, now_secs, ttl_secs, None, &body)?;
        let ds = parse_who_odata(&ind, &iso, &url, &body)?;
        store.commit_vintage(&ds)?;
        Ok(ds)
    }
}

fn parse_who_odata(
    indicator: &str,
    iso3: &str,
    url: &str,
    body_json: &str,
) -> Result<TimeSeriesDataset, EpiError> {
    let root: Value = serde_json::from_str(body_json)?;
    let mut buckets: BTreeMap<String, f64> = BTreeMap::new();

    if let Some(arr) = root.get("value").and_then(Value::as_array) {
        for row in arr {
            let time_dim = row
                .get("TimeDim")
                .map(|v| match v {
                    Value::Number(n) => n.to_string(),
                    Value::String(s) => s.trim().to_string(),
                    _ => String::new(),
                })
                .unwrap_or_default();
            let num_val = row.get("NumericValue").and_then(Value::as_f64).unwrap_or(0.0);
            if !time_dim.is_empty() && num_val.is_finite() {
                *buckets.entry(time_dim).or_insert(0.0) += num_val;
            }
        }
    }

    if buckets.is_empty() {
        return Err(EpiError::Schema(
            "WHO GHO OData response contained 0 numeric observations".to_string(),
        ));
    }

    let points: Vec<TimeSeriesPoint> = buckets
        .into_iter()
        .map(|(date_label, count)| TimeSeriesPoint {
            date_label,
            epiweek: None,
            count,
            suppressed: false,
        })
        .collect();

    let checksum_hex = EpiStore::compute_dataset_checksum(&points);
    Ok(TimeSeriesDataset {
        vintage_id: format!("who_{indicator}_{iso3}_{checksum_hex}"),
        title: format!("WHO GHO {indicator} ({})", if iso3.is_empty() { "Global" } else { iso3 }),
        source_kind: SourceKind::LiveWhoGho,
        source_url: url.to_string(),
        ingested_at_utc: format!("epoch:{}", now_unix_secs()),
        default_si_preset: "measles".to_string(),
        checksum_hex,
        points,
    })
}
