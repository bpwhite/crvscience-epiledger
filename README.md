# EpiLedger (`epiledger-rs`) — Native Rust Public Health Informatics Workbench

100% pure-Rust, local-first public health surveillance, differential-privacy governance, forecast evaluation, privacy-preserving clinical interoperability, environmental biostatistics, sero-epidemiology, antimicrobial resistance (`CLSI M39`), wastewater `SNV` deconvolution, superspreading network analytics, and phylodynamic genomics workbench.

Ships with **three execution interfaces** in a single `epiledger` binary:
1. **16-Tab Native `wgpu`/`egui` Desktop Workstation** (`epiledger` or `epiledger gui`)
2. **8-Screen Headless [`ratatui`](https://ratatui.rs/) + [`crossterm`](https://docs.rs/crossterm/) SSH Terminal UI** (`epiledger tui`) for remote Linux servers without `X11`/`Wayland`
3. **39 Scriptable Batch CLI Subcommands** (`epiledger <subcommand> ...`) + standalone zero-JS HTML5/SVG Situation Report (`SitRep`) export

## Workspace Crate Architecture (5 Crates, 55 Rust Modules)

| Crate | Path | Modules & Capabilities |
| :--- | :--- | :--- |
| `epi-core` | `crates/epi-core` | `epiweek`, `schema`, `seed_data`: Pure-Rust CDC MMWR week math (`MmwrWeek`), [NCHS CDC WONDER small-cell suppression (`<=9`)](https://wonder.cdc.gov/mcd-icd10-provisional.html), SNOMED-CT/ICD-10/LOINC coding enums, and embedded [CDC US Weekly Measles (2022–2026)](https://www.cdc.gov/measles/data-research/index.html) & State Pertussis baselines |
| `epi-store` | `crates/epi-store` | `bitemporal`, `dp_and_audit`, `pprl_bloom`: Bundled SQLite WAL (`rusqlite`) bitemporal vintage ledger, **Discrete Laplace `epsilon`-Differential Privacy budget ledger** ([Public Health DP](https://pmc.ncbi.nlm.nih.gov/articles/PMC8662814/)), **tamper-evident `SHA-256` hash-chained audit log**, and **LSH + birth-year hash-blocked 512-bit `HMAC-SHA256` bigram Bloom PPRL** + Fellegi-Sunter EM ([WA Bloom PPRL](https://pmc.ncbi.nlm.nih.gov/articles/PMC8761329/)) |
| `epi-ingest` | `crates/epi-ingest` | `cdc_soda`, `cmu_delphi`, `openfda`, `who_gho`, `hl7_v2_elr`, `syndromic_essence`, `fhir_r4`, `geojson_spatial`, `local_files`, `sas_xport`: Keyless connectors for [CDC NNDSS](https://data.cdc.gov/NNDSS/NNDSS-Weekly-Data/x9gk-5huc), [CDC NWSS Wastewater](https://data.cdc.gov/Public-Health-Surveillance/NWSS-Public-SARS-CoV-2-Wastewater-Metric-Data/2ew6-ywp6), [CMU Delphi](https://cmu-delphi.github.io/delphi-epidata/api/fluview.html), [openFDA FAERS](https://open.fda.gov/apis/drug/event/), and [WHO GHO](https://www.who.int/data/gho/info/gho-odata-api); [HL7 v2.5.1 ELR](https://www.hl7.org/implement/standards/product_brief.cfm?product_id=98); deterministic **[NSSP ESSENCE](https://pmc.ncbi.nlm.nih.gov/articles/PMC11556616/) chief-complaint & `ICD-10-CM` `NegEx` classifier** ([validation](https://pmc.ncbi.nlm.nih.gov/articles/PMC11426074/)); [HL7 FHIR R4](https://hl7.org/fhir/R4/bundle.html); GeoJSON; and pure-Rust **SAS Transport v5 (`.XPT`)** `IBM370` binary decoding |
| `epi-analytics` | `crates/epi-analytics` | `renewal_rt`, `aberration`, `nowcast`, `seir`, `forecast_scoring`, `superspreading_network`, `wastewater_demix`, `amr_antibiogram`, `capture_recapture`, `cfr_and_surge`, `sero_tnd_ve`, `dlnm_env`, `euromomo_mortality`, `lineage_logistics`, `phylo_outbreak`, `spatial_scan_zones`, `survey_weights`, `pharmacovigilance`, `rates_and_lag`, `math`: [Cori `Rt`](https://pmc.ncbi.nlm.nih.gov/articles/PMC3816335/), [`K = 11` `WIS`](https://pmc.ncbi.nlm.nih.gov/articles/PMC7880475/) + [Inverse-WIS ensembles](https://pmc.ncbi.nlm.nih.gov/articles/PMC8963571/), [Lloyd-Smith Negative Binomial `k` MLE](https://pmc.ncbi.nlm.nih.gov/articles/PMC1791715/) & [bipartite contact networks](https://pmc.ncbi.nlm.nih.gov/articles/PMC13303329/), [Freyja `L1` wastewater demixing](https://pmc.ncbi.nlm.nih.gov/articles/PMC12330420/), [`CLSI M39` antibiograms](https://pmc.ncbi.nlm.nih.gov/articles/PMC9580356/), [3-source log-linear & Chapman capture-recapture](https://pmc.ncbi.nlm.nih.gov/articles/PMC8705476/), [Nishiura right-truncated `cCFR`](https://pmc.ncbi.nlm.nih.gov/articles/PMC7529598/) & [competing-risk `ICU` surge](https://pmc.ncbi.nlm.nih.gov/articles/PMC9200328/), [Gasparrini DLNM](https://pmc.ncbi.nlm.nih.gov/articles/PMC2998707/), [EuroMOMO excess mortality](https://pmc.ncbi.nlm.nih.gov/articles/PMC12944855/), [Kulldorff spatial scan](https://pmc.ncbi.nlm.nih.gov/articles/PMC5816564/), [Taylor survey variance](https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx), and [FAERS `PRR`/`ROR`/`IC`](https://pmc.ncbi.nlm.nih.gov/articles/PMC12443087/) |
| `epi-app` | `crates/epi-app` | `main`, `gui`, `tui`, `cli`, `cartogram`, `genomic_tab`, `env_mortality_tab`, `elr_pprl_sero_tab`, `scoring_sitrep_tab`, `dp_audit_tab`, `network_tab`, `wastewater_amr_tab`, `crc_syndromic_tab`, `severity_surge_tab`, `sitrep_export`: 16-tab `egui` desktop app, 8-screen `ratatui` SSH TUI, HTML5/SVG `SitRep` exporter, and 39 CLI subcommands |

## Build & Run (Linux, macOS, Windows, & Headless SSH)

```bash
# Run all unit tests across the 5-crate workspace
cargo test --workspace

# 1. Launch the 16-tab native wgpu/egui desktop workstation
cargo run --release -p epi-app

# 2. Launch the 8-screen headless SSH Terminal UI (no X11/Wayland needed)
cargo run --release -p epi-app -- tui

# 3. Representative CLI examples across Phases 1-13
cargo run --release -p epi-app -- offspring-k builtin all
cargo run --release -p epi-app -- contact-network builtin outbreak_bipartite.dot
cargo run --release -p epi-app -- wastewater-demix builtin log2 80
cargo run --release -p epi-app -- antibiogram builtin dedup
cargo run --release -p epi-app -- capture-recapture builtin
cargo run --release -p epi-app -- syndromic-classify builtin
cargo run --release -p epi-app -- cfr-adjust builtin 13.5 26
cargo run --release -p epi-app -- hospital-surge builtin 55 32
cargo run --release -p epi-app -- dp-release cdc_measles_2022_2026_v1 0.50 4.0 measles_dp.csv
cargo run --release -p epi-app -- audit-verify 20
```

## Headless Validation

Run the headless validation suite exercising `(epsilon, 0)` Discrete Laplace differential privacy & SHA-256 hash-chain tamper detection, 512-bit CLK Bloom-filter PPRL + NSSP ESSENCE `NegEx`, Cori renewal `Rt` estimation, and Freyja L1 wastewater deconvolution:

```bash
cargo test -p epi-analytics --test headless_validation -- --nocapture
```

Validation artifacts (`epiledger_validation_summary.json` and `cori_rt_and_dp_audit.csv`) are documented in [`validation/README.md`](validation/README.md).

## Citation & CRV Science Apps Hub

- **CRV Science Apps Hub**: <https://crvscience.com/apps>
- **Preprint (`APP-PREPRINT-06`)**: [EpiLedger: Differentially Private Public Health Surveillance, Bloom-Filter Record Linkage, and Renewal Rt Estimation](https://crvscience.com/apps/preprints/epiledger-differential-privacy-epidemiology)
- **Concept DOI**: [`10.5281/zenodo.24222382`](https://doi.org/10.5281/zenodo.24222382)

```bibtex
@article{white2026epiledger,
  author  = {Bryan P. White},
  title   = {EpiLedger: Differentially Private Public Health Surveillance, Bloom-Filter Record Linkage, and Renewal Rt Estimation},
  year    = {2026},
  journal = {CRV Science Apps Preprints},
  note    = {APP-PREPRINT-06, v2.0.0},
  doi     = {10.5281/zenodo.24222382},
  url     = {https://crvscience.com/apps/preprints/epiledger-differential-privacy-epidemiology}
}
```

See [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

