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

use crate::epiweek::MmwrWeek;
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum EpiError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("Invalid dataset or schema: {0}")]
    Schema(String),
    #[error("Upstream live feed error: {0}")]
    Network(String),
    #[error("Storage error: {0}")]
    Storage(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SourceKind {
    VerifiedSnapshot,
    LiveCdcSoda,
    LiveCmuDelphi,
    LiveOpenFda,
    LiveWhoGho,
    LocalCsvTsv,
    LocalJson,
    LocalSasXport,
    LocalHl7Fhir,
}

impl SourceKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::VerifiedSnapshot => "Verified Snapshot",
            Self::LiveCdcSoda => "Live CDC SODA",
            Self::LiveCmuDelphi => "Live CMU Delphi",
            Self::LiveOpenFda => "Live openFDA",
            Self::LiveWhoGho => "Live WHO GHO",
            Self::LocalCsvTsv => "Local CSV/TSV",
            Self::LocalJson => "Local JSON",
            Self::LocalSasXport => "Local SAS XPT",
            Self::LocalHl7Fhir => "Local HL7 FHIR R4",
        }
    }

    pub fn from_label(s: &str) -> Self {
        match s {
            "Live CDC SODA" => Self::LiveCdcSoda,
            "Live CMU Delphi" => Self::LiveCmuDelphi,
            "Live openFDA" => Self::LiveOpenFda,
            "Live WHO GHO" => Self::LiveWhoGho,
            "Local CSV/TSV" => Self::LocalCsvTsv,
            "Local JSON" => Self::LocalJson,
            "Local SAS XPT" => Self::LocalSasXport,
            "Local HL7 FHIR R4" => Self::LocalHl7Fhir,
            _ => Self::VerifiedSnapshot,
        }
    }
}

/// Standardized taxonomic or clinical coding identifier for cross-dataset harmonization.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum TaxonomicCode {
    SnomedCt(String),
    Icd10(String),
    NcbiTaxon(u32),
    Loinc(String),
    NndssLabel(String),
}

/// Statistical disclosure policy for small-cell suppression.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SuppressionPolicy {
    /// Display exact counts without small-cell masking.
    Unmasked,
    /// NCHS / CDC WONDER rule: suppress positive counts of 1..=9 and rates derived from them.
    NchsWonderLe9,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum CellValue {
    Exact(f64),
    Suppressed,
    Missing,
}

impl CellValue {
    pub fn from_count(count: Option<f64>, flag_suppressed: bool, policy: SuppressionPolicy) -> Self {
        if flag_suppressed {
            return Self::Suppressed;
        }
        match count {
            None => Self::Missing,
            Some(v) if !v.is_finite() => Self::Missing,
            Some(v) => {
                if policy == SuppressionPolicy::NchsWonderLe9 && v > 0.0 && v <= 9.0 {
                    Self::Suppressed
                } else {
                    Self::Exact(v)
                }
            }
        }
    }

    pub fn display_string(&self, decimals: usize) -> String {
        match self {
            Self::Exact(v) => format!("{v:.decimals$}"),
            Self::Suppressed => "SUPP (<=9)".to_string(),
            Self::Missing => "--".to_string(),
        }
    }
}

/// Canonical bitemporal surveillance record normalized across SODA, Delphi, openFDA, WHO, FHIR, and local imports.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurveillanceRecord {
    pub source_id: String,
    pub fips_or_iso: String,
    pub condition_code: String,
    pub epiweek: MmwrWeek,
    pub count: Option<f64>,
    pub population: Option<u64>,
    pub suppressed_flag: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeSeriesPoint {
    pub date_label: String,
    pub epiweek: Option<MmwrWeek>,
    pub count: f64,
    pub suppressed: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeSeriesDataset {
    pub vintage_id: String,
    pub title: String,
    pub source_kind: SourceKind,
    pub source_url: String,
    pub ingested_at_utc: String,
    pub default_si_preset: String,
    pub checksum_hex: String,
    pub points: Vec<TimeSeriesPoint>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JurisdictionRecord {
    pub fips: String,
    pub abbr: String,
    pub name: String,
    pub census_pop_2025: u64,
    pub cases_2024: u32,
    pub rate_2024: f64,
    pub cases_2025: u32,
    pub rate_2025: f64,
    pub active_custom_value: Option<f64>,
    pub note: Option<String>,
}

impl JurisdictionRecord {
    pub fn yoy_pct_change(&self, policy: SuppressionPolicy) -> Option<f64> {
        if policy == SuppressionPolicy::NchsWonderLe9
            && ((self.cases_2025 > 0 && self.cases_2025 <= 9)
                || (self.cases_2024 > 0 && self.cases_2024 <= 9))
        {
            return None;
        }
        if self.cases_2024 == 0 {
            return Some(0.0);
        }
        Some(((self.cases_2025 as f64 - self.cases_2024 as f64) / self.cases_2024 as f64) * 100.0)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerialIntervalPreset {
    pub id: &'static str,
    pub label: &'static str,
    pub mean_days: f64,
    pub sd_days: f64,
    pub reference_url: &'static str,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgeStratumCount {
    pub label: String,
    pub cases: f64,
    pub population: f64,
    pub standard_weight: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecallRecord {
    pub recall_number: String,
    pub classification: String,
    pub status: String,
    pub recall_initiation_date: String,
    pub recalling_firm: String,
    pub state: String,
    pub product_description: String,
    pub reason_for_recall: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionFacet {
    pub term: String,
    pub count: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VintageMetadata {
    pub vintage_id: String,
    pub title: String,
    pub source_kind: String,
    pub source_url: String,
    pub ingested_at_utc: String,
    pub row_count: usize,
    pub checksum_hex: String,
}
