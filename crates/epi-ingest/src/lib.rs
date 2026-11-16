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

pub mod cdc_soda;
pub mod cmu_delphi;
pub mod fhir_r4;
pub mod geojson_spatial;
pub mod hl7_v2_elr;
pub mod local_files;
pub mod openfda;
pub mod sas_xport;
pub mod syndromic_essence;
pub mod who_gho;

pub use cdc_soda::{CdcSodaClient, SodaDatasetKind, SodaFetchOutcome};
pub use cmu_delphi::CmuDelphiClient;
pub use fhir_r4::{
    export_timeseries_as_fhir_r4_bundle, import_fhir_r4_bundle_file, parse_fhir_r4_bundle_str,
    FhirIngestSummary,
};
pub use geojson_spatial::{
    parse_geojson_centroids_file, parse_geojson_centroids_str, GeoFeatureCentroid,
};
pub use hl7_v2_elr::{
    builtin_hl7_v2_elr_and_nssp_batch, import_hl7_v2_file, parse_hl7_v2_batch_str,
    Hl7v2IngestOutcome, NsspSyndromeCategory, ParsedHl7v2Event,
};
pub use local_files::{
    export_dataset_csv, export_dataset_json, import_delimited_timeseries, import_json_timeseries,
    ColumnMapping,
};
pub use openfda::{FdaFaersOutcome, OpenFdaClient};
pub use sas_xport::{parse_sas_xport_bytes, SasVariableMeta, SasXportDataset, SasXportValue};
pub use syndromic_essence::{
    build_syndrome_weekly_dataset, builtin_ed_chief_complaint_encounters,
    classify_ed_encounters_essence, evaluate_custom_essence_query, match_affirmative_or_negated,
    normalize_chief_complaint, parse_ed_encounters_csv, ClassifiedEdEncounter, EdTriageEncounter,
    EssenceBatchClassification, EssenceSyndromeId, SyndromeSummaryRow,
};
pub use who_gho::WhoGhoClient;
