# Headless Validation (`crvscience-epiledger` v2.0.0)

Deterministic, dry-lab validation suite for the `epi-core`, `epi-store`, `epi-ingest`, and `epi-analytics` workspace crates executed via compiled Rust integration tests with zero GUI or external network dependencies.

## Command

```bash
cargo test --workspace
# Or run only the headless validation integration suite:
cargo test -p epi-analytics --test headless_validation -- --nocapture
```

## Test Coverage

1. **T1 (`t1_discrete_laplace_dp_cell_suppression_and_hash_chain_tamper_verification`)**: Evaluates NCHS / CDC WONDER small-cell suppression (`<= 9`) on the built-in US measles MMWR weekly series, executes `(epsilon, 0)` Discrete Laplace differential privacy releases (`epsilon = 0.75` and cumulative `1.75` against a `2.0` cap), asserts rejection of budget-exceeding queries (`2.25 > 2.0`), and proves out-of-band SQLite tampering of `audit_hash_chain` (`seq_no = 1`) is detected deterministically by SHA-256 chain verification (`verify_audit_chain`).
2. **T2 (`t2_clk_bloom_pprl_precision_recall_f1_and_essence_negex`)**: Benchmarks 512-bit HMAC-SHA256 Cryptographic Long-Term Key (CLK) Bloom-filter Privacy-Preserving Record Linkage (`link_and_deduplicate_pprl`) with birth-year + 4-band LSH hash blocking across 80 perturbed demographic ELR/ADT records (computing pair precision, recall, F1, and candidate reduction %) alongside NSSP ESSENCE `NegEx` clinical negation & exclusion filtering (`classify_ed_encounters_essence`) asserting zero false positives across all 8 negated control encounters.
3. **T3 (`t3_cori_rt_posterior_cri_and_freyja_wastewater_deconvolution`)**: Computes Cori renewal `Rt` posterior mean, median, and 95% credible intervals (`compute_cori_rt`) on both the built-in measles series (exporting weekly raw, NCHS-suppressed, DP `epsilon=0.75`, and `Rt` estimates to [`cori_rt_and_dp_audit.csv`](cori_rt_and_dp_audit.csv)) and the NSSP ILI respiratory series via EM daily disaggregation (`used_daily_em_disaggregation = true`), plus Freyja L1 constrained-least-squares wastewater SNV lineage deconvolution (`deconvolve_wastewater_timeseries`).

## Generated Validation Artifacts

- [`epiledger_validation_summary.json`](epiledger_validation_summary.json) — Structured JSON telemetry for T1–T3 (mirrored to `app-preprints/releases/epiledger-differential-privacy-epidemiology-validation-telemetry.json`).
- [`cori_rt_and_dp_audit.csv`](cori_rt_and_dp_audit.csv) — Weekly measles surveillance audit combining raw counts, NCHS WONDER suppression, DP `epsilon = 0.75` counts, and Cori `Rt` 95% credible intervals.
