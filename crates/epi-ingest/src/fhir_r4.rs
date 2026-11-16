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

use crate::cdc_soda::now_unix_secs;
use epi_core::{
    builtin_pertussis_spatial_records, parse_flexible_date_to_mmwr, EpiError, MmwrWeek,
    SourceKind, TimeSeriesDataset, TimeSeriesPoint,
};
use epi_store::EpiStore;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

/// Parsed summary from an HL7 FHIR R4 Electronic Case Reporting (eCR) or Electronic Lab Reporting (ELR) Bundle.
/// Specifications: https://hl7.org/fhir/R4/bundle.html and https://hl7.org/fhir/R4/observation.html
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FhirIngestSummary {
    pub dataset: TimeSeriesDataset,
    pub state_counts_by_fips: HashMap<String, f64>,
    pub coding_tallies: Vec<(String, usize)>,
    pub total_resources_scanned: usize,
}

/// Parse an HL7 FHIR R4 `Bundle` JSON payload (`Observation`, `Condition`, and `Patient` resources)
/// into an aggregated MMWR weekly `TimeSeriesDataset` and state spatial map.
pub fn parse_fhir_r4_bundle_str(
    json_str: &str,
    code_filter: Option<&str>,
    source_label: &str,
) -> Result<FhirIngestSummary, EpiError> {
    let root: Value = serde_json::from_str(json_str)?;
    let res_type = root
        .get("resourceType")
        .and_then(Value::as_str)
        .unwrap_or("");
    if res_type != "Bundle" {
        return Err(EpiError::Schema(format!(
            "Expected HL7 FHIR R4 `resourceType: \"Bundle\"`, found '{res_type}'"
        )));
    }

    let entries = root
        .get("entry")
        .and_then(Value::as_array)
        .ok_or_else(|| EpiError::Schema("FHIR Bundle contains no `entry` array".to_string()))?;

    let state_lookup = build_fhir_state_map();
    let mut patient_state_by_ref: HashMap<String, String> = HashMap::new();

    // Pass 1: Index Patient -> US state FIPS from `Patient.address[].state`
    for entry in entries {
        let Some(res) = entry.get("resource") else {
            continue;
        };
        if res.get("resourceType").and_then(Value::as_str) == Some("Patient") {
            let pid = res.get("id").and_then(Value::as_str).unwrap_or("");
            if let Some(addrs) = res.get("address").and_then(Value::as_array) {
                for a in addrs {
                    let st_str = a
                        .get("state")
                        .and_then(Value::as_str)
                        .unwrap_or("")
                        .trim()
                        .to_uppercase();
                    if let Some(fips) = state_lookup.get(&st_str) {
                        if !pid.is_empty() {
                            patient_state_by_ref.insert(format!("Patient/{pid}"), fips.clone());
                            patient_state_by_ref.insert(pid.to_string(), fips.clone());
                        }
                    }
                }
            }
        }
    }

    let filter_lc = code_filter
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_ascii_lowercase);

    let mut week_counts: BTreeMap<MmwrWeek, f64> = BTreeMap::new();
    let mut state_counts: HashMap<String, f64> = HashMap::new();
    let mut code_map: HashMap<String, usize> = HashMap::new();

    // Pass 2: Extract Observation and Condition case/count events
    for entry in entries {
        let Some(res) = entry.get("resource") else {
            continue;
        };
        let rtype = res.get("resourceType").and_then(Value::as_str).unwrap_or("");
        if rtype != "Observation" && rtype != "Condition" {
            continue;
        }

        let (code_display, matched_filter) = extract_coding_label(res, filter_lc.as_deref());
        if !matched_filter {
            continue;
        }

        let date_str = if rtype == "Observation" {
            res.get("effectiveDateTime")
                .or_else(|| res.get("issued"))
                .or_else(|| {
                    res.get("effectivePeriod")
                        .and_then(|p| p.get("start"))
                })
                .and_then(Value::as_str)
                .unwrap_or("")
        } else {
            res.get("onsetDateTime")
                .or_else(|| res.get("recordedDate"))
                .and_then(Value::as_str)
                .unwrap_or("")
        };

        let Some(mmwr) = parse_flexible_date_to_mmwr(date_str) else {
            continue;
        };

        // Support aggregate public-health Observation (`valueInteger` / `valueQuantity.value`) or line-list (1.0 per case)
        let weight = res
            .get("valueInteger")
            .and_then(Value::as_f64)
            .or_else(|| {
                res.get("valueQuantity")
                    .and_then(|vq| vq.get("value"))
                    .and_then(Value::as_f64)
            })
            .unwrap_or(1.0)
            .max(0.0);

        *week_counts.entry(mmwr).or_insert(0.0) += weight;
        *code_map.entry(code_display).or_insert(0) += 1;

        if let Some(subj_ref) = res
            .get("subject")
            .and_then(|s| s.get("reference"))
            .and_then(Value::as_str)
        {
            if let Some(fips) = patient_state_by_ref.get(subj_ref) {
                *state_counts.entry(fips.clone()).or_insert(0.0) += weight;
            }
        }
    }

    if week_counts.is_empty() {
        return Err(EpiError::Schema(
            "No dated Observation or Condition resources matched in HL7 FHIR R4 Bundle".to_string(),
        ));
    }

    let points: Vec<TimeSeriesPoint> = week_counts
        .into_iter()
        .map(|(wk, count)| TimeSeriesPoint {
            date_label: wk.compact_label(),
            epiweek: Some(wk),
            count,
            suppressed: false,
        })
        .collect();

    let mut coding_tallies: Vec<(String, usize)> = code_map.into_iter().collect();
    coding_tallies.sort_by(|a, b| b.1.cmp(&a.1));

    let checksum_hex = EpiStore::compute_dataset_checksum(&points);
    let dataset = TimeSeriesDataset {
        vintage_id: format!("fhir_r4_{}_{}", now_unix_secs(), checksum_hex),
        title: format!("HL7 FHIR R4 Bundle: {source_label}"),
        source_kind: SourceKind::LocalHl7Fhir,
        source_url: format!("fhir://{source_label}"),
        ingested_at_utc: format!("epoch:{}", now_unix_secs()),
        default_si_preset: "measles".to_string(),
        checksum_hex,
        points,
    };

    Ok(FhirIngestSummary {
        dataset,
        state_counts_by_fips: state_counts,
        coding_tallies,
        total_resources_scanned: entries.len(),
    })
}

pub fn import_fhir_r4_bundle_file<P: AsRef<Path>>(
    path: P,
    code_filter: Option<&str>,
) -> Result<FhirIngestSummary, EpiError> {
    let path_ref = path.as_ref();
    let raw = std::fs::read_to_string(path_ref)?;
    let label = path_ref
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("eCR_Bundle");
    parse_fhir_r4_bundle_str(&raw, code_filter, label)
}

/// Export a canonical `TimeSeriesDataset` as a standards-compliant HL7 FHIR R4 `collection` Bundle
/// of aggregate `Observation` resources coded with LOINC `85354-9` / public health surveillance measure tags.
pub fn export_timeseries_as_fhir_r4_bundle(
    dataset: &TimeSeriesDataset,
    loinc_or_snomed_code: &str,
) -> Result<String, EpiError> {
    let entries: Vec<Value> = dataset
        .points
        .iter()
        .enumerate()
        .map(|(idx, pt)| {
            let eff_date = pt
                .epiweek
                .map(|w| w.sunday_iso_date())
                .unwrap_or_else(|| pt.date_label.clone());
            json!({
                "fullUrl": format!("urn:uuid:{}-obs-{}", dataset.checksum_hex, idx),
                "resource": {
                    "resourceType": "Observation",
                    "id": format!("obs-{idx}"),
                    "status": "final",
                    "category": [{
                        "coding": [{
                            "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory",
                            "display": "Public Health Surveillance Aggregate"
                        }]
                    }],
                    "code": {
                        "coding": [{
                            "system": "http://snomed.info/sct",
                            "code": loinc_or_snomed_code,
                            "display": dataset.title
                        }]
                    },
                    "effectiveDateTime": eff_date,
                    "valueQuantity": {
                        "value": pt.count,
                        "unit": "cases",
                        "system": "http://unitsofmeasure.org",
                        "code": "{cases}"
                    }
                }
            })
        })
        .collect();

    let bundle = json!({
        "resourceType": "Bundle",
        "id": dataset.vintage_id,
        "type": "collection",
        "timestamp": dataset.ingested_at_utc,
        "entry": entries
    });

    serde_json::to_string_pretty(&bundle).map_err(EpiError::from)
}

fn extract_coding_label(res: &Value, filter_lc: Option<&str>) -> (String, bool) {
    let mut labels = Vec::new();
    if let Some(codings) = res
        .get("code")
        .and_then(|c| c.get("coding"))
        .and_then(Value::as_array)
    {
        for c in codings {
            let code = c.get("code").and_then(Value::as_str).unwrap_or("");
            let disp = c.get("display").and_then(Value::as_str).unwrap_or("");
            if !code.is_empty() || !disp.is_empty() {
                labels.push(format!("{code} {disp}").trim().to_string());
            }
        }
    }
    let primary = labels
        .first()
        .cloned()
        .unwrap_or_else(|| "Uncoded Resource".to_string());

    let matched = match filter_lc {
        None => true,
        Some(flt) => labels.iter().any(|l| l.to_ascii_lowercase().contains(flt)),
    };
    (primary, matched)
}

fn build_fhir_state_map() -> HashMap<String, String> {
    let mut map = HashMap::new();
    for rec in builtin_pertussis_spatial_records() {
        map.insert(rec.abbr.to_uppercase(), rec.fips.clone());
        map.insert(rec.name.to_uppercase(), rec.fips.clone());
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use epi_core::builtin_measles_weekly_dataset;

    #[test]
    fn test_fhir_r4_roundtrip_bundle() {
        let mut ds = builtin_measles_weekly_dataset();
        ds.points.truncate(6);
        let fhir_json = export_timeseries_as_fhir_r4_bundle(&ds, "14189004").unwrap();
        let parsed = parse_fhir_r4_bundle_str(&fhir_json, Some("14189004"), "MeaslesTest").unwrap();
        assert_eq!(parsed.dataset.points.len(), 6);
        assert_eq!(parsed.total_resources_scanned, 6);
    }
}
