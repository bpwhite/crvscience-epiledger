use crate::cdc_soda::now_unix_secs;
use epi_core::{
    parse_flexible_date_to_mmwr, CellValue, EpiError, SourceKind, SuppressionPolicy,
    TimeSeriesDataset, TimeSeriesPoint,
};
use epi_store::EpiStore;
use serde_json::Value;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct ColumnMapping {
    pub date_or_key_col: Option<String>,
    pub value_col: Option<String>,
    pub delimiter: Option<u8>,
}

impl Default for ColumnMapping {
    fn default() -> Self {
        Self {
            date_or_key_col: None,
            value_col: None,
            delimiter: None,
        }
    }
}

/// Import a local CSV or TSV public dataset (CDC WONDER, Project Tycho, OWID, state tables) into a `TimeSeriesDataset`.
pub fn import_delimited_timeseries<P: AsRef<Path>>(
    path: P,
    title: Option<&str>,
    mapping: &ColumnMapping,
) -> Result<TimeSeriesDataset, EpiError> {
    let path_ref = path.as_ref();
    let delim = mapping.delimiter.unwrap_or_else(|| {
        if path_ref
            .extension()
            .and_then(|s| s.to_str())
            .map(|e| e.eq_ignore_ascii_case("tsv") || e.eq_ignore_ascii_case("txt"))
            .unwrap_or(false)
        {
            b'\t'
        } else {
            b','
        }
    });

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim)
        .flexible(true)
        .trim(csv::Trim::All)
        .from_path(path_ref)
        .map_err(|e| EpiError::Schema(format!("CSV open failed: {e}")))?;

    let headers = rdr
        .headers()
        .map_err(|e| EpiError::Schema(format!("CSV header error: {e}")))?
        .clone();

    if headers.len() < 2 {
        return Err(EpiError::Schema(
            "Delimited dataset requires at least 2 columns (date/week + count)".to_string(),
        ));
    }

    let key_idx = resolve_column_index(&headers, mapping.date_or_key_col.as_deref(), &[
        "date",
        "week_start_date",
        "epiweek",
        "week",
        "year",
        "period",
    ])
    .unwrap_or(0);

    let val_idx = resolve_column_index(&headers, mapping.value_col.as_deref(), &[
        "confirmed_cases",
        "cases",
        "count",
        "new_cases",
        "deaths",
        "value",
        "m1",
    ])
    .unwrap_or(1.min(headers.len() - 1));

    let mut buckets: BTreeMap<String, (f64, bool)> = BTreeMap::new();
    for result in rdr.records() {
        let rec = result.map_err(|e| EpiError::Schema(format!("CSV row error: {e}")))?;
        let key_raw = rec.get(key_idx).unwrap_or("").trim();
        let val_raw = rec.get(val_idx).unwrap_or("").trim();
        if key_raw.is_empty() {
            continue;
        }
        if val_raw.eq_ignore_ascii_case("suppressed")
            || val_raw.starts_with("SUPP")
            || val_raw == "-"
        {
            buckets.entry(key_raw.to_string()).or_insert((0.0, true));
            continue;
        }
        if let Ok(v) = val_raw.replace(',', "").parse::<f64>() {
            if v.is_finite() {
                let entry = buckets.entry(key_raw.to_string()).or_insert((0.0, false));
                entry.0 += v.max(0.0);
            }
        }
    }

    if buckets.len() < 3 {
        return Err(EpiError::Schema(
            "Fewer than 3 valid time points parsed from delimited file".to_string(),
        ));
    }

    let points: Vec<TimeSeriesPoint> = buckets
        .into_iter()
        .map(|(date_label, (count, suppressed))| TimeSeriesPoint {
            epiweek: parse_flexible_date_to_mmwr(&date_label),
            date_label,
            count,
            suppressed,
        })
        .collect();

    let checksum_hex = EpiStore::compute_dataset_checksum(&points);
    let stem = title
        .map(str::to_string)
        .or_else(|| {
            path_ref
                .file_stem()
                .and_then(|s| s.to_str())
                .map(str::to_string)
        })
        .unwrap_or_else(|| "Imported CSV Dataset".to_string());

    Ok(TimeSeriesDataset {
        vintage_id: format!("csv_{}_{}", now_unix_secs(), checksum_hex),
        title: stem,
        source_kind: SourceKind::LocalCsvTsv,
        source_url: format!("file://{}", path_ref.display()),
        ingested_at_utc: format!("epoch:{}", now_unix_secs()),
        default_si_preset: "measles".to_string(),
        checksum_hex,
        points,
    })
}

/// Import a local JSON dataset (either an exported `TimeSeriesDataset` or a generic JSON array of objects).
pub fn import_json_timeseries<P: AsRef<Path>>(path: P) -> Result<TimeSeriesDataset, EpiError> {
    let path_ref = path.as_ref();
    let content = std::fs::read_to_string(path_ref)?;

    if let Ok(mut direct) = serde_json::from_str::<TimeSeriesDataset>(&content) {
        direct.source_kind = SourceKind::LocalJson;
        direct.checksum_hex = EpiStore::compute_dataset_checksum(&direct.points);
        return Ok(direct);
    }

    let root: Value = serde_json::from_str(&content)?;
    let arr = root
        .as_array()
        .or_else(|| root.get("points").and_then(Value::as_array))
        .or_else(|| root.get("results").and_then(Value::as_array))
        .ok_or_else(|| EpiError::Schema("JSON file has no top-level or `points` array".to_string()))?;

    let mut buckets: BTreeMap<String, f64> = BTreeMap::new();
    for item in arr {
        let Some(obj) = item.as_object() else { continue };
        let key = ["date", "date_label", "week", "period", "time"]
            .iter()
            .find_map(|k| obj.get(*k).and_then(Value::as_str))
            .unwrap_or("")
            .trim()
            .to_string();
        let val = ["count", "cases", "value", "NumericValue"]
            .iter()
            .find_map(|k| obj.get(*k).and_then(Value::as_f64))
            .unwrap_or(f64::NAN);
        if !key.is_empty() && val.is_finite() {
            *buckets.entry(key).or_insert(0.0) += val.max(0.0);
        }
    }

    if buckets.len() < 3 {
        return Err(EpiError::Schema(
            "Fewer than 3 valid time points found in JSON array".to_string(),
        ));
    }

    let points: Vec<TimeSeriesPoint> = buckets
        .into_iter()
        .map(|(date_label, count)| TimeSeriesPoint {
            epiweek: parse_flexible_date_to_mmwr(&date_label),
            date_label,
            count,
            suppressed: false,
        })
        .collect();

    let checksum_hex = EpiStore::compute_dataset_checksum(&points);
    let title = path_ref
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("Imported JSON Dataset")
        .to_string();

    Ok(TimeSeriesDataset {
        vintage_id: format!("json_{}_{}", now_unix_secs(), checksum_hex),
        title,
        source_kind: SourceKind::LocalJson,
        source_url: format!("file://{}", path_ref.display()),
        ingested_at_utc: format!("epoch:{}", now_unix_secs()),
        default_si_preset: "measles".to_string(),
        checksum_hex,
        points,
    })
}

pub fn export_dataset_csv<P: AsRef<Path>>(
    dataset: &TimeSeriesDataset,
    path: P,
    policy: SuppressionPolicy,
) -> Result<(), EpiError> {
    let mut wtr = csv::Writer::from_path(path)
        .map_err(|e| EpiError::Schema(format!("CSV writer error: {e}")))?;
    wtr.write_record(["date_label", "mmwr_week", "count"])
        .map_err(|e| EpiError::Schema(e.to_string()))?;

    for pt in &dataset.points {
        let wk_str = pt.epiweek.map(|w| w.compact_label()).unwrap_or_default();
        let cell = CellValue::from_count(Some(pt.count), pt.suppressed, policy);
        wtr.write_record([&pt.date_label, &wk_str, &cell.display_string(0)])
            .map_err(|e| EpiError::Schema(e.to_string()))?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn export_dataset_json<P: AsRef<Path>>(
    dataset: &TimeSeriesDataset,
    path: P,
) -> Result<(), EpiError> {
    let mut f = File::create(path)?;
    let txt = serde_json::to_string_pretty(dataset)?;
    f.write_all(txt.as_bytes())?;
    Ok(())
}

fn resolve_column_index(
    headers: &csv::StringRecord,
    explicit: Option<&str>,
    candidates: &[&str],
) -> Option<usize> {
    if let Some(exp) = explicit {
        for (i, h) in headers.iter().enumerate() {
            if h.eq_ignore_ascii_case(exp.trim()) {
                return Some(i);
            }
        }
    }
    for cand in candidates {
        for (i, h) in headers.iter().enumerate() {
            if h.to_ascii_lowercase().contains(cand) {
                return Some(i);
            }
        }
    }
    None
}
