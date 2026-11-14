use crate::cdc_soda::now_unix_secs;
use epi_core::{
    parse_flexible_date_to_mmwr, EpiError, MmwrWeek, ReactionFacet, RecallRecord, SourceKind,
    TimeSeriesDataset, TimeSeriesPoint,
};
use epi_store::EpiStore;
use serde_json::Value;
use std::collections::BTreeMap;
use std::time::Duration;

#[derive(Debug, Clone)]
pub struct FdaFaersOutcome {
    pub weekly_series: TimeSeriesDataset,
    pub top_reactions: Vec<ReactionFacet>,
    pub used_stale_cache: bool,
}

pub struct OpenFdaClient {
    api_key: Option<String>,
    timeout: Duration,
}

impl Default for OpenFdaClient {
    fn default() -> Self {
        Self {
            api_key: None,
            timeout: Duration::from_secs(15),
        }
    }
}

impl OpenFdaClient {
    pub fn new(api_key: Option<String>) -> Self {
        Self {
            api_key,
            timeout: Duration::from_secs(15),
        }
    }

    pub fn build_faers_series_url(&self, generic_drug: &str) -> String {
        let search = build_faers_search_clause(generic_drug);
        let mut url = format!("https://api.fda.gov/drug/event.json?search={search}&count=receivedate");
        if let Some(k) = &self.api_key {
            if !k.trim().is_empty() {
                url.push_str(&format!("&api_key={}", k.trim()));
            }
        }
        url
    }

    pub fn build_faers_reactions_url(&self, generic_drug: &str, limit: usize) -> String {
        let search = build_faers_search_clause(generic_drug);
        let mut url = format!(
            "https://api.fda.gov/drug/event.json?search={search}&count=patient.reaction.reactionmeddrapt.exact&limit={}",
            limit.clamp(1, 50)
        );
        if let Some(k) = &self.api_key {
            if !k.trim().is_empty() {
                url.push_str(&format!("&api_key={}", k.trim()));
            }
        }
        url
    }

    pub fn fetch_faers_with_cache(
        &self,
        store: &mut EpiStore,
        generic_drug: &str,
        ttl_secs: i64,
    ) -> Result<FdaFaersOutcome, EpiError> {
        let series_url = self.build_faers_series_url(generic_drug);
        let react_url = self.build_faers_reactions_url(generic_drug, 12);
        let cache_key = format!("fda_faers:{generic_drug}");
        let react_cache_key = format!("fda_react:{generic_drug}");
        let now_secs = now_unix_secs();

        if let Some(cached_s) = store.get_cached_http(&cache_key, now_secs)? {
            if cached_s.is_fresh {
                let react_body = store
                    .get_cached_http(&react_cache_key, now_secs)?
                    .map(|c| c.body_json)
                    .unwrap_or_default();
                return parse_faers_payloads(
                    generic_drug,
                    &series_url,
                    &cached_s.body_json,
                    &react_body,
                    false,
                );
            }
        }

        let agent = ureq::AgentBuilder::new()
            .timeout_connect(self.timeout)
            .timeout_read(self.timeout)
            .build();

        match agent.get(&series_url).call() {
            Ok(resp_s) => {
                let body_s = resp_s
                    .into_string()
                    .map_err(|e| EpiError::Network(e.to_string()))?;
                store.put_cached_http(&cache_key, &series_url, now_secs, ttl_secs, None, &body_s)?;

                let body_r = match agent.get(&react_url).call() {
                    Ok(r) => {
                        let txt = r.into_string().unwrap_or_default();
                        let _ = store.put_cached_http(
                            &react_cache_key,
                            &react_url,
                            now_secs,
                            ttl_secs,
                            None,
                            &txt,
                        );
                        txt
                    }
                    Err(_) => String::new(),
                };

                let outcome = parse_faers_payloads(
                    generic_drug,
                    &series_url,
                    &body_s,
                    &body_r,
                    false,
                )?;
                store.commit_vintage(&outcome.weekly_series)?;
                Ok(outcome)
            }
            Err(err) => {
                if let Some(stale) = store.get_cached_http(&cache_key, now_secs)? {
                    let react_body = store
                        .get_cached_http(&react_cache_key, now_secs)?
                        .map(|c| c.body_json)
                        .unwrap_or_default();
                    return parse_faers_payloads(
                        generic_drug,
                        &series_url,
                        &stale.body_json,
                        &react_body,
                        true,
                    );
                }
                Err(EpiError::Network(format!("openFDA FAERS error: {err}")))
            }
        }
    }

    pub fn fetch_drug_recalls(&self, limit: usize) -> Result<Vec<RecallRecord>, EpiError> {
        let url = format!(
            "https://api.fda.gov/drug/enforcement.json?sort=report_date:desc&limit={}",
            limit.clamp(1, 100)
        );
        let agent = ureq::AgentBuilder::new()
            .timeout_connect(self.timeout)
            .timeout_read(self.timeout)
            .build();

        let body = agent
            .get(&url)
            .call()
            .map_err(|e| EpiError::Network(e.to_string()))?
            .into_string()
            .map_err(|e| EpiError::Network(e.to_string()))?;

        let root: Value = serde_json::from_str(&body)?;
        let mut out = Vec::new();
        if let Some(arr) = root.get("results").and_then(Value::as_array) {
            for item in arr {
                out.push(RecallRecord {
                    recall_number: str_field(item, "recall_number"),
                    classification: str_field(item, "classification"),
                    status: str_field(item, "status"),
                    recall_initiation_date: format_yyyymmdd(&str_field(
                        item,
                        "recall_initiation_date",
                    )),
                    recalling_firm: str_field(item, "recalling_firm"),
                    state: str_field(item, "state"),
                    product_description: str_field(item, "product_description"),
                    reason_for_recall: str_field(item, "reason_for_recall"),
                });
            }
        }
        Ok(out)
    }
}

fn build_faers_search_clause(generic_drug: &str) -> String {
    let trimmed = generic_drug.trim();
    if trimmed.is_empty() {
        "receivedate:[20240101+TO+20261112]".to_string()
    } else {
        let clean: String = trimmed
            .chars()
            .filter(|c| c.is_ascii_alphanumeric() || *c == ' ' || *c == '-')
            .collect();
        let encoded = clean.replace(' ', "+");
        format!("receivedate:[20240101+TO+20261112]+AND+patient.drug.openfda.generic_name:\"{encoded}\"")
    }
}

fn parse_faers_payloads(
    generic_drug: &str,
    source_url: &str,
    series_json: &str,
    react_json: &str,
    used_stale_cache: bool,
) -> Result<FdaFaersOutcome, EpiError> {
    let root: Value = serde_json::from_str(series_json)?;
    let mut week_buckets: BTreeMap<MmwrWeek, f64> = BTreeMap::new();

    if let Some(arr) = root.get("results").and_then(Value::as_array) {
        for item in arr {
            let t_str = str_field(item, "time");
            let cnt = item.get("count").and_then(Value::as_f64).unwrap_or(0.0);
            if let Some(mmwr) = parse_flexible_date_to_mmwr(&t_str) {
                *week_buckets.entry(mmwr).or_insert(0.0) += cnt;
            }
        }
    }

    if week_buckets.is_empty() {
        return Err(EpiError::Schema(
            "openFDA FAERS returned no weekly aggregatable records".to_string(),
        ));
    }

    let points: Vec<TimeSeriesPoint> = week_buckets
        .into_iter()
        .map(|(wk, count)| TimeSeriesPoint {
            date_label: wk.compact_label(),
            epiweek: Some(wk),
            count,
            suppressed: false,
        })
        .collect();

    let mut top_reactions = Vec::new();
    if !react_json.trim().is_empty() {
        if let Ok(r_val) = serde_json::from_str::<Value>(react_json) {
            if let Some(arr) = r_val.get("results").and_then(Value::as_array) {
                for item in arr {
                    let term = str_field(item, "term");
                    let count = item.get("count").and_then(Value::as_u64).unwrap_or(0);
                    if !term.is_empty() {
                        top_reactions.push(ReactionFacet { term, count });
                    }
                }
            }
        }
    }

    let checksum_hex = EpiStore::compute_dataset_checksum(&points);
    let label = if generic_drug.trim().is_empty() {
        "All Drugs"
    } else {
        generic_drug.trim()
    };

    let weekly_series = TimeSeriesDataset {
        vintage_id: format!("fda_faers_{}_{}", now_unix_secs(), checksum_hex),
        title: format!("openFDA FAERS Weekly Reports: {label}"),
        source_kind: SourceKind::LiveOpenFda,
        source_url: source_url.to_string(),
        ingested_at_utc: format!("epoch:{}", now_unix_secs()),
        default_si_preset: "influenza".to_string(),
        checksum_hex,
        points,
    };

    Ok(FdaFaersOutcome {
        weekly_series,
        top_reactions,
        used_stale_cache,
    })
}

fn str_field(v: &Value, key: &str) -> String {
    v.get(key)
        .and_then(Value::as_str)
        .unwrap_or("")
        .trim()
        .to_string()
}

fn format_yyyymmdd(s: &str) -> String {
    if s.len() == 8 && s.bytes().all(|b| b.is_ascii_digit()) {
        format!("{}-{}-{}", &s[0..4], &s[4..6], &s[6..8])
    } else {
        s.to_string()
    }
}
