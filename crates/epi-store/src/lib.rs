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

pub mod bitemporal;
pub mod dp_and_audit;
pub mod pprl_bloom;

pub use bitemporal::{BackfillRevisionDelta, CachedHttpPayload, EpiStore};
pub use dp_and_audit::{
    append_audit_event, compute_audit_entry_hash, get_vintage_spent_epsilon,
    get_vintage_spent_rho, list_audit_events, list_dp_budget_entries,
    release_differentially_private_series, release_discrete_gaussian_zcdp_series,
    reset_vintage_dp_budget, sample_discrete_gaussian_ctr, verify_audit_chain,
    zcdp_rho_to_epsilon, AuditChainVerification, AuditLogEntry, DpBudgetLedgerEntry,
    DpReleaseConfig, DpReleaseOutcome, ZcdpReleaseConfig, ZcdpReleaseOutcome, GENESIS_PREV_HASH,
};
pub use pprl_bloom::{
    builtin_pprl_hospital_lab_records, encode_pprl_record, generate_synthetic_pprl_scale_benchmark,
    hmac_sha256, link_and_deduplicate_pprl, BloomFilter512, EncodedPprlRecord, PprlCandidatePair,
    PprlConfig, PprlDemographicInput, PprlLinkageSummary,
};
