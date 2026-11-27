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

#![forbid(unsafe_code)]

pub mod epiweek;
pub mod schema;
pub mod seed_data;

pub use epiweek::{parse_flexible_date_to_mmwr, MmwrWeek};
pub use schema::{
    AgeStratumCount, CellValue, EpiError, JurisdictionRecord, ReactionFacet, RecallRecord,
    SerialIntervalPreset, SourceKind, SuppressionPolicy, SurveillanceRecord, TaxonomicCode,
    TimeSeriesDataset, TimeSeriesPoint, VintageMetadata,
};
pub use seed_data::{
    builtin_measles_weekly_dataset, builtin_pertussis_spatial_records, pathogen_si_presets,
    us_2000_standard_population_weights, SurveillanceCampaignPreset,
};
