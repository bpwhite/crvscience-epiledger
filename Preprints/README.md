# Preprints for EpiLedger (`crvscience-epiledger`)

**Author:** Bryan P. White (ORCID: [0000-0003-4191-3511](https://orcid.org/0000-0003-4191-3511))  
**Hub:** https://crvscience.com/apps

This directory contains the self-contained LaTeX source (`.tex`) and single-column PDF (`.pdf`) editions for both:

1. **Companion Application Preprint (`APP-PREPRINT-06`)**
   - **Title:** EpiLedger: Privacy-Safe Disease Tracking, Clinical Text Screening, and Outbreak Monitoring
   - **Files:** `epiledger-differential-privacy-epidemiology.tex`, `epiledger-differential-privacy-epidemiology.pdf`
   - **Canonical Zenodo DOI:** https://doi.org/10.5281/zenodo.24222382
   - **Figshare Mirror:** https://figshare.com/articles/preprint/EpiLedger_Privacy-Safe_Disease_Tracking_Clinical_Text_Screening_and_Outbreak_Monitoring/35013757
   - **Web Preprint:** https://crvscience.com/apps/preprints/epiledger-differential-privacy-epidemiology

2. **CRV Science Public Apps Network Overview (`APP-PREPRINT-01`)**
   - **Title:** The CRV Science App Network: Twenty Open-Source Rust Programs for Science, Data, and Simulation
   - **Files:** `crv-science-public-apps-network.tex`, `crv-science-public-apps-network.pdf`
   - **Canonical Zenodo DOI:** https://doi.org/10.5281/zenodo.24222398
   - **Figshare Mirror:** https://figshare.com/articles/preprint/The_CRV_Science_App_Network_Seventeen_Open-Source_Rust_Programs_for_Science_Data_and_Simulation/35013736
   - **Web Preprint:** https://crvscience.com/apps/preprints/crv-science-public-apps-network

## Compiling Standalone LaTeX

Both `.tex` files inline their bibliography so each compiles standalone with standard `pdflatex`:

```bash
pdflatex epiledger-differential-privacy-epidemiology.tex
pdflatex crv-science-public-apps-network.tex
```
