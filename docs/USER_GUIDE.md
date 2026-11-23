# EpiLedger User Guide (`epiledger-rs` v2.0.0)

EpiLedger is a 5-crate pure-Rust public health informatics workbench combining local SQLite bitemporal vintage storage, statistical disclosure control, epidemiological analytics, live API connectors, a `wgpu`/`egui` desktop UI, an 8-screen `ratatui` SSH TUI, and 39 CLI subcommands.

---

## 1. Navigation & Keyboard Shortcuts

The desktop UI organizes **16 analytical modules** under **5 Domain Pillars**:

| Pillar (`Alt+1..5`) | Sub-Tabs (Cycle with `Ctrl+Tab` / `Ctrl+Shift+Tab`) |
| :--- | :--- |
| **`Alt+1` Outbreak & Rt** | 1. Rt & Aberration · 2. SEIR & Delay Nowcast · 3. Superspreading & Network |
| **`Alt+2` Spatial & Syndromic** | 4. Kulldorff Spatial Scan · 5. Capture-Recapture & ESSENCE · 6. Env DLNM & EuroMOMO |
| **`Alt+3` Genomic, Wastewater & Lab** | 7. Genomic & Phylodynamics · 8. Wastewater & CLSI AMR · 9. ELR, PPRL & Sero/VE |
| **`Alt+4` Forecasting, Audit & Ledger** | 10. WIS Scoring & SitRep · 11. cCFR & ICU Surge · 12. DP Privacy & Audit · 13. Vintage Ledger & Diff |
| **`Alt+5` Connectors & Surveys** | 14. Live CDC / Delphi / WHO APIs · 15. Local CSV / .XPT / FHIR · 16. Survey (NHANES) & PRR (FAERS) |

Additional global shortcuts:
- **`Ctrl+K` / `Ctrl+F`**: Focus the top-bar search input; type a module, state, pathogen, or gene keyword and press **`Enter`** to jump directly to the matching module.
- **`F1`**: Toggle the **About EpiLedger** window.
- **`Escape`**: Close the About window, or clear the active search query.

---

## 2. Privacy Controls (`NCHS <=9 Mask` & Differential Privacy)

### NCHS / CDC WONDER Small-Cell Confidentiality
Toggle **`NCHS <=9 Mask`** in the top navigation bar to activate `SuppressionPolicy::NchsWonderLe9` globally:
- Epidemic curve (`epi_curve_plot`) zeros out low-count cells `1..=9` and updates the y-axis label.
- Spatial cartogram (`Kulldorff Spatial Scan`) renders `1..=9` states in purple-slate with diagonal hatching and a `*` badge (`TX*`), visually distinct from `0`-case states.
- Hover tooltips, the selected-state inspector panel, and the `Single-Jurisdiction Excess LLR Table` all mask `1..=9` counts and derived incidence rates as `SUPP (<=9) (SUPP/100k)`.

### Differential Privacy & Audit Ledger (Tab 12)
- Opening Tab 12 lazily runs the 600-patient hash-blocked PPRL scale benchmark on first visit (eliminating launch latency).
- Choose **Discrete Laplace ($\varepsilon$-DP)** or **Discrete Gaussian ($\rho$-zCDP)**, configure budget caps, and release noisy time series backed by 256-bit OS CSPRNG entropy (`getrandom`).
- Every release writes an immutable row to `dp_budget_ledger` (tracking $\varepsilon$, $\rho$, cumulative budgets, and per-row `zcdp_delta`) and extends `audit_hash_chain`.

---

## 3. Core Analytical Workflows

1. **Rt & Aberration (`Rt & Aberration`)**: Inspect synchronized `epi_curve_plot` and `cori_rt_plot` charts with shaded 95% CrI polygon bands, single-entry multi-wave legends, and baseline/active outbreak cutoff sliders (`52w Active` / `26w Active`).
2. **SEIR & Delay Nowcast (`SEIR & Delay Nowcast`)**: Couple parametric or vintage reporting-triangle delay nowcasting into parallel RK4 + Binomial tau-leaping SEIR/SEIRS projections. Sliders recompute on `drag_stopped()` / `lost_focus()` for smooth interaction.
3. **Kulldorff Spatial Scan (`Kulldorff Spatial Scan`)**: Click or hover any tile on the responsive 12×8 US cartogram with perceptual ramp and colorbar legend, or import custom GeoJSON polygons for $M=999$ Monte Carlo circular scan clustering.
4. **Exporting Situation Reports (`WIS Scoring & SitRep`)**: Export self-contained HTML5/SVG Situation Reports with zero external JS dependencies.

---

## 4. Headless Terminal UI (`tui`) & CLI Quick Reference

```bash
# Headless SSH terminal UI (8 interactive views)
epiledger tui

# Inspect SQLite vintages & run Cori Rt + Farrington GLM
epiledger vintages
epiledger analyze cdc_measles_2022_2026_v1

# Run DP release and verify audit chain integrity
epiledger dp-release cdc_measles_2022_2026_v1 0.50 4.0 measles_dp.csv
epiledger audit-verify 20

# Export masked standalone HTML5/SVG Situation Report
epiledger export-sitrep cdc_measles_2022_2026_v1 sitrep.html mask
```
