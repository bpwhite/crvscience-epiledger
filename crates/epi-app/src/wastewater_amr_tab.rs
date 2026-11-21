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

use eframe::egui::{self, Color32, RichText};
use egui_plot::{Bar, BarChart, HLine, Legend, Line, Plot, PlotPoints};
use epi_analytics::{
    builtin_hospital_amr_isolates, builtin_wastewater_barcode_and_samples,
    compute_clsi_m39_antibiogram, deconvolve_wastewater_timeseries, parse_amr_isolates_csv,
    parse_wastewater_snv_csv, AmrIsolateRecord, AntibiogramConfig, BarcodeLibrary,
    CumulativeAntibiogramResult, DepthWeighting, IsolateDeduplicationPolicy, LineageCountTable,
    LongitudinalWastewaterResult, WastewaterDemixConfig, WastewaterSampleInput,
};
use epi_core::{SurveillanceCampaignPreset, TimeSeriesDataset};

const LINEAGE_COLORS: [Color32; 6] = [
    Color32::from_rgb(56, 189, 248),
    Color32::from_rgb(52, 211, 153),
    Color32::from_rgb(251, 191, 36),
    Color32::from_rgb(244, 63, 94),
    Color32::from_rgb(167, 139, 250),
    Color32::from_rgb(251, 146, 60),
];

#[derive(Default)]
pub struct Phase10GuiAction {
    pub push_wastewater_to_genomic: Option<LineageCountTable>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase10SubView {
    WastewaterFreyjaDemix,
    ClsiM39Antibiogram,
}

pub struct WastewaterAmrWorkbenchState {
    pub sub_view: Phase10SubView,

    // Wastewater SNV Deconvolution state
    pub ww_library: BarcodeLibrary,
    pub ww_samples: Vec<WastewaterSampleInput>,
    pub ww_cfg: WastewaterDemixConfig,
    pub ww_result: LongitudinalWastewaterResult,
    pub ww_csv_path: String,
    pub selected_sample_idx: usize,

    // CLSI M39 Antibiogram state
    pub amr_records: Vec<AmrIsolateRecord>,
    pub amr_cfg: AntibiogramConfig,
    pub amr_result: CumulativeAntibiogramResult,
    pub amr_csv_path: String,
    pub selected_organism_idx: usize,
}

impl Default for WastewaterAmrWorkbenchState {
    fn default() -> Self {
        Self::new()
    }
}

impl WastewaterAmrWorkbenchState {
    pub fn new() -> Self {
        let (ww_library, ww_samples) = builtin_wastewater_barcode_and_samples();
        let ww_cfg = WastewaterDemixConfig::default();
        let ww_result = deconvolve_wastewater_timeseries(&ww_library, &ww_samples, &ww_cfg)
            .expect("Built-in wastewater SNV series is valid");
        let selected_sample_idx = ww_result.samples.len().saturating_sub(1);

        let amr_records = builtin_hospital_amr_isolates();
        let amr_cfg = AntibiogramConfig::default();
        let amr_result = compute_clsi_m39_antibiogram(&amr_records, &amr_cfg)
            .expect("Built-in clinical AMR isolates are valid");

        Self {
            sub_view: Phase10SubView::WastewaterFreyjaDemix,
            ww_library,
            ww_samples,
            ww_cfg,
            ww_result,
            ww_csv_path: String::new(),
            selected_sample_idx,
            amr_records,
            amr_cfg,
            amr_result,
            amr_csv_path: String::new(),
            selected_organism_idx: 0,
        }
    }

    pub fn recompute_wastewater(&mut self) -> Result<(), String> {
        let res = deconvolve_wastewater_timeseries(&self.ww_library, &self.ww_samples, &self.ww_cfg)
            .map_err(|e| e.to_string())?;
        if self.selected_sample_idx >= res.samples.len() {
            self.selected_sample_idx = res.samples.len().saturating_sub(1);
        }
        self.ww_result = res;
        Ok(())
    }

    pub fn recompute_antibiogram(&mut self) -> Result<(), String> {
        let res = compute_clsi_m39_antibiogram(&self.amr_records, &self.amr_cfg)
            .map_err(|e| e.to_string())?;
        if self.selected_organism_idx >= res.organisms.len() {
            self.selected_organism_idx = 0;
        }
        self.amr_result = res;
        Ok(())
    }

    pub fn sync_to_campaign(
        &mut self,
        campaign: SurveillanceCampaignPreset,
        active_ds: &TimeSeriesDataset,
    ) {
        let (mut lib, mut smps) = builtin_wastewater_barcode_and_samples();
        let (lin_names, sewershed, load_scale) = match campaign {
            SurveillanceCampaignPreset::UsMeasles2022_2026 => (
                vec![
                    "MeV-B3.1".to_string(),
                    "MeV-D8-TX".to_string(),
                    "MeV-D8-OH".to_string(),
                    "MeV-D8-Import".to_string(),
                ],
                "South-Plains-Sentinel-WRF",
                0.18,
            ),
            SurveillanceCampaignPreset::UsPertussis2022_2026 => (
                vec![
                    "Bp-ptxP3".to_string(),
                    "Bp-PRN-del".to_string(),
                    "Bp-fhaB2".to_string(),
                    "Bp-MacRes".to_string(),
                ],
                "Great-Lakes-Municipal-WRF",
                0.35,
            ),
            SurveillanceCampaignPreset::InfluenzaIliWastewater2022_2026 => (
                vec![
                    "IAV-H1N1pdm09".to_string(),
                    "IAV-H3N2-2a3a1".to_string(),
                    "IAV-H3N2-J2".to_string(),
                    "IBV-Victoria".to_string(),
                ],
                "US-NWSS-ILINet-Sewershed",
                0.85,
            ),
            SurveillanceCampaignPreset::SarsCoV2VariantCohort2022_2026 => (
                lib.lineage_names.clone(),
                "Corvallis-Willamette-WRF",
                1.0,
            ),
        };
        lib.lineage_names = lin_names;
        let tail_start = active_ds.points.len().saturating_sub(smps.len());
        let tail = &active_ds.points[tail_start..];
        for (i, s) in smps.iter_mut().enumerate() {
            s.sewershed = sewershed.to_string();
            s.viral_copies_per_liter = (s.viral_copies_per_liter * load_scale).max(1.2e3);
            if let Some(pt) = tail.get(i) {
                s.date_label = pt
                    .epiweek
                    .map(|w| w.compact_label())
                    .unwrap_or_else(|| pt.date_label.clone());
            }
        }
        self.ww_library = lib;
        self.ww_samples = smps;
        let _ = self.recompute_wastewater();
    }

    pub fn render_with_filter(
        &mut self,
        ui: &mut egui::Ui,
        status_banner: &mut String,
        search_filter: &str,
    ) -> Phase10GuiAction {
        let mut action = Phase10GuiAction::default();

        ui.horizontal_wrapped(|ui| {
            ui.heading("Environmental Wastewater SNV Demixing & CLSI M39 Hospital Antibiograms");
            ui.separator();
            ui.selectable_value(
                &mut self.sub_view,
                Phase10SubView::WastewaterFreyjaDemix,
                "Freyja L1 Wastewater Deconvolution",
            );
            ui.selectable_value(
                &mut self.sub_view,
                Phase10SubView::ClsiM39Antibiogram,
                "CLSI M39 Antibiograms & MDR/XDR",
            );
        });
        ui.add_space(4.0);

        match self.sub_view {
            Phase10SubView::WastewaterFreyjaDemix => {
                self.render_wastewater_panel(ui, status_banner, &mut action, search_filter)
            }
            Phase10SubView::ClsiM39Antibiogram => {
                self.render_antibiogram_panel(ui, status_banner, search_filter)
            }
        }

        action
    }

    fn render_wastewater_panel(
        &mut self,
        ui: &mut egui::Ui,
        status_banner: &mut String,
        action: &mut Phase10GuiAction,
        search_filter: &str,
    ) {
        let latest = self
            .ww_result
            .samples
            .get(self.selected_sample_idx)
            .or_else(|| self.ww_result.samples.last());

        ui.columns(4, |cols| {
            if let Some(smp) = latest {
                let top_lin = smp.lineages.iter().max_by(|a, b| {
                    a.abundance
                        .partial_cmp(&b.abundance)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                cols[0].group(|ui| {
                    ui.label(RichText::new("DOMINANT SEWAGE LINEAGE").small().weak());
                    if let Some(tl) = top_lin {
                        ui.heading(
                            RichText::new(format!("{} ({:.1}%)", tl.lineage, tl.abundance * 100.0))
                                .monospace()
                                .color(Color32::from_rgb(56, 189, 248)),
                        );
                        ui.label(
                            RichText::new(format!(
                                "95% Boot [{:.1}%, {:.1}%]",
                                tl.ci_low_95 * 100.0,
                                tl.ci_high_95 * 100.0
                            ))
                            .small()
                            .monospace(),
                        );
                    }
                });

                cols[1].group(|ui| {
                    ui.label(RichText::new("WEIGHTED L1 / RMSE FIT").small().weak());
                    ui.heading(
                        RichText::new(format!("L1 = {:.4}", smp.weighted_l1_loss))
                            .monospace()
                            .color(Color32::from_rgb(52, 211, 153)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "RMSE {:.4} | {}/{} SNVs",
                            smp.weighted_rmse, smp.sites_used, smp.sites_total_library
                        ))
                        .small()
                        .monospace(),
                    );
                });

                cols[2].group(|ui| {
                    ui.label(RichText::new("SEWERSHED VIRAL CONCENTRATION").small().weak());
                    ui.heading(
                        RichText::new(format!("{:.2e} gc/L", smp.viral_copies_per_liter))
                            .monospace()
                            .color(Color32::from_rgb(251, 191, 36)),
                    );
                    ui.label(
                        RichText::new(format!("{} ({})", smp.sewershed, smp.date_label))
                            .small()
                            .monospace(),
                    );
                });

                cols[3].group(|ui| {
                    ui.label(RichText::new("CRYPTIC UNCATALOGUED SNVs").small().weak());
                    let n_cryp = smp.cryptic_sites.len();
                    let col = if n_cryp > 0 {
                        Color32::from_rgb(244, 63, 94)
                    } else {
                        Color32::from_rgb(52, 211, 153)
                    };
                    ui.heading(
                        RichText::new(format!("{n_cryp} flagged"))
                            .monospace()
                            .color(col),
                    );
                    let detail = smp
                        .cryptic_sites
                        .first()
                        .map(|c| format!("Top: {} (+{:.1}%)", c.mutation_site, c.residual_excess * 100.0))
                        .unwrap_or_else(|| "All SNVs explained by library".to_string());
                    ui.label(RichText::new(detail).small().monospace());
                });
            }
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Depth Weighting:").strong());
                let prev_dw = self.ww_cfg.depth_weighting;
                egui::ComboBox::from_id_source("ww_depth_weight_combo")
                    .selected_text(match self.ww_cfg.depth_weighting {
                        DepthWeighting::Log2DepthPlusOne => "Freyja log2(depth + 1)",
                        DepthWeighting::SqrtDepth => "sqrt(depth)",
                        DepthWeighting::Uniform => "Uniform (unweighted L1)",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.ww_cfg.depth_weighting,
                            DepthWeighting::Log2DepthPlusOne,
                            "Freyja log2(depth + 1)",
                        );
                        ui.selectable_value(
                            &mut self.ww_cfg.depth_weighting,
                            DepthWeighting::SqrtDepth,
                            "sqrt(depth)",
                        );
                        ui.selectable_value(
                            &mut self.ww_cfg.depth_weighting,
                            DepthWeighting::Uniform,
                            "Uniform (unweighted L1)",
                        );
                    });

                ui.label("Bootstraps:");
                let boot_changed = ui
                    .add(egui::DragValue::new(&mut self.ww_cfg.bootstrap_replicates).range(0..=300))
                    .changed();

                if self.ww_cfg.depth_weighting != prev_dw || boot_changed {
                    if let Err(e) = self.recompute_wastewater() {
                        *status_banner = format!("Wastewater demix error: {e}");
                    }
                }

                ui.separator();
                ui.label("Inspect Sample:");
                if !self.ww_result.samples.is_empty() {
                    let cur_lbl = self
                        .ww_result
                        .samples
                        .get(self.selected_sample_idx)
                        .map(|s| format!("{} ({})", s.sample_id, s.date_label))
                        .unwrap_or_default();
                    egui::ComboBox::from_id_source("ww_sample_selector")
                        .selected_text(cur_lbl)
                        .show_ui(ui, |ui| {
                            for (idx, s) in self.ww_result.samples.iter().enumerate() {
                                ui.selectable_value(
                                    &mut self.selected_sample_idx,
                                    idx,
                                    format!("{} ({})", s.sample_id, s.date_label),
                                );
                            }
                        });
                }

                if ui
                    .button("Send Deconvolved Shares to Genomic Logistic Fit")
                    .on_hover_text("Convert longitudinal sewage abundances into a LineageCountTable in Tab 7")
                    .clicked()
                {
                    action.push_wastewater_to_genomic = Some(
                        self.ww_result
                            .to_lineage_count_table(self.ww_cfg.pseudo_total_reads_per_sample),
                    );
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("SNV Barcode CSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.ww_csv_path)
                        .desired_width(240.0)
                        .hint_text("path/to/wastewater_snv_barcodes.csv"),
                );
                if ui.button("Browse...").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("CSV / TSV", &["csv", "tsv", "txt"])
                        .pick_file()
                    {
                        self.ww_csv_path = p.display().to_string();
                    }
                }
                if ui.button("Import Wastewater CSV").clicked() {
                    match std::fs::read_to_string(self.ww_csv_path.trim()) {
                        Ok(raw) => match parse_wastewater_snv_csv(&raw) {
                            Ok((lib, smps)) => {
                                self.ww_library = lib;
                                self.ww_samples = smps;
                                match self.recompute_wastewater() {
                                    Ok(()) => {
                                        *status_banner = format!(
                                            "Deconvolved {} wastewater samples across {} target lineages.",
                                            self.ww_result.samples.len(),
                                            self.ww_result.lineage_names.len()
                                        );
                                    }
                                    Err(e) => *status_banner = format!("Deconvolution error: {e}"),
                                }
                            }
                            Err(e) => *status_banner = format!("CSV parse error: {e}"),
                        },
                        Err(e) => *status_banner = format!("File read error: {e}"),
                    }
                }
                if ui.button("Reset Built-In Willamette WRF Series").clicked() {
                    let (lib, smps) = builtin_wastewater_barcode_and_samples();
                    self.ww_library = lib;
                    self.ww_samples = smps;
                    let _ = self.recompute_wastewater();
                    *status_banner =
                        "Reset to built-in 10-week Willamette Valley wastewater amplicon series."
                            .to_string();
                }
            });
        });

        ui.add_space(6.0);

        ui.columns(2, |cols| {
            let n_ww = self.ww_result.samples.len().max(1);
            let mut ab_y_max = 0.0_f64;
            let mut ab_left_max = 0.0_f64;
            let mut ab_right_max = 0.0_f64;
            let mut cp_y_max = 0.0_f64;
            let mut cp_left_max = 0.0_f64;
            let mut cp_right_max = 0.0_f64;
            let left_end = (n_ww / 3).max(1);
            let right_start = n_ww.saturating_sub((n_ww / 3).max(1));
            for (idx, s) in self.ww_result.samples.iter().enumerate() {
                for l in &s.lineages {
                    let pct = l.abundance * 100.0;
                    ab_y_max = ab_y_max.max(pct);
                    cp_y_max = cp_y_max.max(l.viral_load_weighted_copies);
                    if idx < left_end {
                        ab_left_max = ab_left_max.max(pct);
                        cp_left_max = cp_left_max.max(l.viral_load_weighted_copies);
                    }
                    if idx >= right_start {
                        ab_right_max = ab_right_max.max(pct);
                        cp_right_max = cp_right_max.max(l.viral_load_weighted_copies);
                    }
                }
            }
            let ab_corner = if ab_left_max <= ab_right_max {
                egui_plot::Corner::LeftTop
            } else {
                egui_plot::Corner::RightTop
            };
            let cp_corner = if cp_left_max <= cp_right_max {
                egui_plot::Corner::LeftTop
            } else {
                egui_plot::Corner::RightTop
            };

            cols[0].group(|ui| {
                ui.label(
                    RichText::new("Longitudinal Sewage Lineage Abundance Trajectories (%)")
                        .strong(),
                );
                let mut plot = Plot::new("ww_longitudinal_abundance_plot")
                    .height(235.0)
                    .legend(Legend::default().position(ab_corner).background_alpha(0.85))
                    .x_axis_label("Sampling Week Index")
                    .y_axis_label("Deconvolved Lineage Share (%)")
                    .include_y(0.0);
                if ab_y_max > 0.0 {
                    plot = plot.include_y((ab_y_max * 1.72).min(100.0).max(ab_y_max * 1.15));
                }
                plot.show(ui, |plot_ui| {
                        for (k_idx, lin_name) in self.ww_result.lineage_names.iter().enumerate() {
                            let pts: Vec<[f64; 2]> = self
                                .ww_result
                                .samples
                                .iter()
                                .map(|s| {
                                    let p = s
                                        .lineages
                                        .get(k_idx)
                                        .map(|l| l.abundance * 100.0)
                                        .unwrap_or(0.0);
                                    [s.time_weeks, p]
                                })
                                .collect();
                            plot_ui.line(
                                Line::new(PlotPoints::from(pts))
                                    .name(lin_name)
                                    .color(LINEAGE_COLORS[k_idx % LINEAGE_COLORS.len()])
                                    .width(2.5_f32),
                            );
                        }
                    });
            });

            cols[1].group(|ui| {
                ui.label(
                    RichText::new("Viral-Load-Weighted Lineage Signal (Gene Copies / L)").strong(),
                );
                let mut plot = Plot::new("ww_viral_copies_plot")
                    .height(235.0)
                    .legend(Legend::default().position(cp_corner).background_alpha(0.85))
                    .x_axis_label("Sampling Week Index")
                    .y_axis_label("Lineage Copies / L")
                    .y_axis_min_width(58.0)
                    .include_y(0.0)
                    .y_axis_formatter(|mark, _range| {
                        let v = mark.value;
                        if v.abs() >= 1_000_000.0 {
                            format!("{:.1}M", v / 1_000_000.0)
                        } else if v.abs() >= 1_000.0 {
                            let k = v / 1_000.0;
                            if (k - k.round()).abs() < 0.05 {
                                format!("{:.0}k", k)
                            } else {
                                format!("{:.1}k", k)
                            }
                        } else {
                            format!("{:.0}", v)
                        }
                    });
                if cp_y_max > 0.0 {
                    plot = plot.include_y(cp_y_max * 1.24);
                }
                plot.show(ui, |plot_ui| {
                        for (k_idx, lin_name) in self.ww_result.lineage_names.iter().enumerate() {
                            let pts: Vec<[f64; 2]> = self
                                .ww_result
                                .samples
                                .iter()
                                .map(|s| {
                                    let cp = s
                                        .lineages
                                        .get(k_idx)
                                        .map(|l| l.viral_load_weighted_copies)
                                        .unwrap_or(0.0);
                                    [s.time_weeks, cp]
                                })
                                .collect();
                            plot_ui.line(
                                Line::new(PlotPoints::from(pts))
                                    .name(lin_name)
                                    .color(LINEAGE_COLORS[k_idx % LINEAGE_COLORS.len()])
                                    .width(2.3_f32),
                            );
                        }
                    });
            });
        });

        ui.add_space(6.0);

        if let Some(smp) = self.ww_result.samples.get(self.selected_sample_idx) {
            ui.columns(2, |cols| {
                cols[0].group(|ui| {
                    ui.label(
                        RichText::new(format!(
                            "Sample {} ({}) Deconvolved Simplex Shares & 95% Bootstrap CIs",
                            smp.sample_id, smp.date_label
                        ))
                        .strong(),
                    );
                    egui::Grid::new("ww_sample_detail_grid")
                        .striped(true)
                        .min_col_width(80.0)
                        .show(ui, |ui| {
                            ui.strong("Lineage");
                            ui.strong("Share (%)");
                            ui.strong("95% Site-Boot CI");
                            ui.strong("Weighted gc/L");
                            ui.end_row();

                            for l in &smp.lineages {
                                ui.monospace(&l.lineage);
                                ui.monospace(format!("{:.2}%", l.abundance * 100.0));
                                ui.monospace(format!(
                                    "[{:.2}%, {:.2}%]",
                                    l.ci_low_95 * 100.0,
                                    l.ci_high_95 * 100.0
                                ));
                                ui.monospace(format!("{:.2e}", l.viral_load_weighted_copies));
                                ui.end_row();
                            }
                        });
                });

                cols[1].group(|ui| {
                    ui.label(
                        RichText::new("Longitudinal Sample Summary & Cryptic Excess Audit")
                            .strong(),
                    );
                    egui::ScrollArea::vertical()
                        .id_source("ww_series_table_scroll")
                        .max_height(185.0)
                        .show(ui, |ui| {
                            egui::Grid::new("ww_series_table_grid")
                                .striped(true)
                                .min_col_width(64.0)
                                .show(ui, |ui| {
                                    ui.strong("Sample");
                                    ui.strong("Week");
                                    ui.strong("L1 Loss");
                                    ui.strong("Top Lineage");
                                    ui.strong("Cryptic Flags");
                                    ui.end_row();

                                    let q = search_filter.trim().to_ascii_lowercase();
                                    for s in &self.ww_result.samples {
                                        let top = s
                                            .lineages
                                            .iter()
                                            .max_by(|a, b| {
                                                a.abundance
                                                    .partial_cmp(&b.abundance)
                                                    .unwrap_or(std::cmp::Ordering::Equal)
                                            })
                                            .map(|l| {
                                                format!("{} ({:.0}%)", l.lineage, l.abundance * 100.0)
                                            })
                                            .unwrap_or_default();
                                        if !q.is_empty()
                                            && !s.sample_id.to_ascii_lowercase().contains(&q)
                                            && !s.date_label.to_ascii_lowercase().contains(&q)
                                            && !top.to_ascii_lowercase().contains(&q)
                                        {
                                            continue;
                                        }
                                        ui.monospace(&s.sample_id);
                                        ui.monospace(&s.date_label);
                                        ui.monospace(format!("{:.4}", s.weighted_l1_loss));
                                        ui.monospace(top);
                                        ui.monospace(format!("{}", s.cryptic_sites.len()));
                                        ui.end_row();
                                    }
                                });
                        });
                });
            });
        }
    }

    fn render_antibiogram_panel(
        &mut self,
        ui: &mut egui::Ui,
        status_banner: &mut String,
        search_filter: &str,
    ) {
        let ab = &self.amr_result;
        let low_n_orgs = ab
            .organism_summaries
            .iter()
            .filter(|o| !o.meets_clsi_m39_n30)
            .count();
        let total_mdr: usize = ab.organism_summaries.iter().map(|o| o.mdr_isolates).sum();

        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("INCLUDED DIAGNOSTIC ISOLATES").small().weak());
                ui.heading(
                    RichText::new(format!("{} isolates", ab.total_included_isolates))
                        .monospace()
                        .color(Color32::from_rgb(56, 189, 248)),
                );
                ui.label(
                    RichText::new(format!(
                        "{} repeat cultures excluded",
                        ab.total_duplicates_excluded
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("MAGIORAKOS MDR BURDEN").small().weak());
                let mdr_pct = (total_mdr as f64) * 100.0 / (ab.total_included_isolates.max(1) as f64);
                ui.heading(
                    RichText::new(format!("{total_mdr} ({mdr_pct:.1}% MDR)"))
                        .monospace()
                        .color(Color32::from_rgb(244, 63, 94)),
                );
                ui.label(
                    RichText::new(">=1 non-susceptible agent in >=3 classes")
                        .small()
                        .monospace(),
                );
            });

            cols[2].group(|ui| {
                ui.label(RichText::new("CLSI M39 n >= 30 RELIABILITY").small().weak());
                let ok_orgs = ab.organisms.len().saturating_sub(low_n_orgs);
                let col = if low_n_orgs == 0 {
                    Color32::from_rgb(52, 211, 153)
                } else {
                    Color32::from_rgb(251, 191, 36)
                };
                ui.heading(
                    RichText::new(format!("{ok_orgs}/{} species >=30", ab.organisms.len()))
                        .monospace()
                        .color(col),
                );
                ui.label(
                    RichText::new(format!("{low_n_orgs} low-denominator caution"))
                        .small()
                        .monospace(),
                );
            });

            cols[3].group(|ui| {
                ui.label(RichText::new("STEWARDSHIP EMPIRIC CUTOFF").small().weak());
                let flagged = ab
                    .cells
                    .iter()
                    .filter(|c| c.empiric_caution_below_threshold)
                    .count();
                ui.heading(
                    RichText::new(format!(
                        "{flagged} pairs < {:.0}% S",
                        self.amr_cfg.empiric_stewardship_threshold_pct
                    ))
                    .monospace()
                    .color(Color32::from_rgb(167, 139, 250)),
                );
                ui.label(
                    RichText::new(format!("Unit scope: {}", ab.facility_unit_filter))
                        .small()
                        .monospace(),
                );
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Deduplication Policy:").strong());
                let prev_pol = self.amr_cfg.deduplication;
                egui::ComboBox::from_id_source("amr_dedup_policy_combo")
                    .selected_text(match self.amr_cfg.deduplication {
                        IsolateDeduplicationPolicy::ClsiM39FirstIsolatePerPatient => {
                            "CLSI M39 (1st Isolate / Patient / Species)"
                        }
                        IsolateDeduplicationPolicy::AllIsolatesUnfiltered => {
                            "Unfiltered All Cultures (Shows Repeat Bias)"
                        }
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.amr_cfg.deduplication,
                            IsolateDeduplicationPolicy::ClsiM39FirstIsolatePerPatient,
                            "CLSI M39 (1st Isolate / Patient / Species)",
                        );
                        ui.selectable_value(
                            &mut self.amr_cfg.deduplication,
                            IsolateDeduplicationPolicy::AllIsolatesUnfiltered,
                            "Unfiltered All Cultures (Shows Repeat Bias)",
                        );
                    });

                ui.label("Unit Filter:");
                let prev_unit = self.amr_cfg.facility_unit_filter.clone();
                egui::ComboBox::from_id_source("amr_unit_filter_combo")
                    .selected_text(if self.amr_cfg.facility_unit_filter.is_empty() {
                        "All Hospital Units"
                    } else {
                        self.amr_cfg.facility_unit_filter.as_str()
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.amr_cfg.facility_unit_filter,
                            String::new(),
                            "All Hospital Units",
                        );
                        for unit in [
                            "ICU",
                            "Inpatient Med-Surg",
                            "Emergency",
                            "Outpatient Clinic",
                        ] {
                            ui.selectable_value(
                                &mut self.amr_cfg.facility_unit_filter,
                                unit.to_string(),
                                unit,
                            );
                        }
                    });

                ui.label("Empiric %S Cutoff:");
                let thr_changed = ui
                    .add(
                        egui::Slider::new(
                            &mut self.amr_cfg.empiric_stewardship_threshold_pct,
                            60.0..=95.0,
                        )
                        .suffix("%"),
                    )
                    .changed();

                if self.amr_cfg.deduplication != prev_pol
                    || self.amr_cfg.facility_unit_filter != prev_unit
                    || thr_changed
                {
                    if let Err(e) = self.recompute_antibiogram() {
                        *status_banner = format!("Antibiogram filter warning: {e}");
                    }
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("Plot Pathogen:");
                if !self.amr_result.organisms.is_empty() {
                    let cur_org = self
                        .amr_result
                        .organisms
                        .get(self.selected_organism_idx)
                        .cloned()
                        .unwrap_or_default();
                    egui::ComboBox::from_id_source("amr_org_plot_combo")
                        .selected_text(cur_org)
                        .show_ui(ui, |ui| {
                            for (i, org) in self.amr_result.organisms.iter().enumerate() {
                                ui.selectable_value(&mut self.selected_organism_idx, i, org);
                            }
                        });
                }

                ui.separator();
                ui.label("Susceptibility CSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.amr_csv_path)
                        .desired_width(210.0)
                        .hint_text("path/to/amr_isolates.csv"),
                );
                if ui.button("Browse...").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("CSV / TSV", &["csv", "tsv", "txt"])
                        .pick_file()
                    {
                        self.amr_csv_path = p.display().to_string();
                    }
                }
                if ui.button("Import AMR CSV").clicked() {
                    match std::fs::read_to_string(self.amr_csv_path.trim()) {
                        Ok(raw) => match parse_amr_isolates_csv(&raw) {
                            Ok(recs) => {
                                self.amr_records = recs;
                                match self.recompute_antibiogram() {
                                    Ok(()) => {
                                        *status_banner = format!(
                                            "Computed CLSI M39 antibiogram on {} deduplicated isolates ({} excluded).",
                                            self.amr_result.total_included_isolates,
                                            self.amr_result.total_duplicates_excluded
                                        );
                                    }
                                    Err(e) => *status_banner = format!("Antibiogram error: {e}"),
                                }
                            }
                            Err(e) => *status_banner = format!("AMR CSV error: {e}"),
                        },
                        Err(e) => *status_banner = format!("File read error: {e}"),
                    }
                }
                if ui.button("Reset Built-In Hospital Panel").clicked() {
                    self.amr_records = builtin_hospital_amr_isolates();
                    self.amr_cfg.facility_unit_filter.clear();
                    let _ = self.recompute_antibiogram();
                }
            });
        });

        ui.add_space(6.0);

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                let sel_org = self
                    .amr_result
                    .organisms
                    .get(self.selected_organism_idx)
                    .cloned()
                    .unwrap_or_default();
                ui.label(
                    RichText::new(format!(
                        "{sel_org} — % Susceptible & Empiric Cutoff ({:.0}%)",
                        self.amr_cfg.empiric_stewardship_threshold_pct
                    ))
                    .strong(),
                );
                let org_cells: Vec<_> = self
                    .amr_result
                    .cells
                    .iter()
                    .filter(|c| c.organism == sel_org)
                    .collect();
                let bars: Vec<Bar> = org_cells
                    .iter()
                    .enumerate()
                    .map(|(idx, c)| {
                        let col = if c.empiric_caution_below_threshold {
                            Color32::from_rgb(244, 63, 94)
                        } else {
                            Color32::from_rgb(52, 211, 153)
                        };
                        Bar::new(idx as f64, c.pct_susceptible)
                            .name(&c.drug_code)
                            .width(0.58)
                            .fill(col)
                    })
                    .collect();

                Plot::new("amr_pct_susceptible_plot")
                    .height(215.0)
                    .legend(Legend::default())
                    .x_axis_label("Tested Antimicrobial Index")
                    .y_axis_label("% Susceptible (CLSI M39)")
                    .show(ui, |plot_ui| {
                        plot_ui.bar_chart(BarChart::new(bars).name("%S"));
                        plot_ui.hline(
                            HLine::new(self.amr_cfg.empiric_stewardship_threshold_pct)
                                .name("Empiric Stewardship Threshold")
                                .color(Color32::from_rgb(251, 191, 36)),
                        );
                    });
            });

            cols[1].group(|ui| {
                ui.label(
                    RichText::new("Organism Deduplication, MDR/XDR & Sentinel Resistance Summary")
                        .strong(),
                );
                egui::Grid::new("amr_org_summary_grid")
                    .striped(true)
                    .min_col_width(58.0)
                    .show(ui, |ui| {
                        ui.strong("Organism");
                        ui.strong("N (Excl)");
                        ui.strong("M39 >=30");
                        ui.strong("MDR %");
                        ui.strong("XDR %");
                        ui.strong("Sentinel Phenotype");
                        ui.end_row();

                        for s in &self.amr_result.organism_summaries {
                            ui.monospace(&s.organism);
                            ui.monospace(format!(
                                "{} (-{})",
                                s.deduplicated_isolates, s.duplicate_isolates_excluded
                            ));
                            let m39_lbl = if s.meets_clsi_m39_n30 {
                                RichText::new("PASS").color(Color32::from_rgb(52, 211, 153))
                            } else {
                                RichText::new("LOW-N <30").color(Color32::from_rgb(251, 191, 36))
                            };
                            ui.label(m39_lbl);
                            ui.monospace(format!("{:.1}%", s.pct_mdr));
                            ui.monospace(format!("{:.1}%", s.pct_xdr));
                            ui.monospace(format!(
                                "{}: {:.1}%",
                                s.sentinel_phenotype_label, s.pct_sentinel_resistant
                            ));
                            ui.end_row();
                        }
                    });
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.label(
                RichText::new(
                    "CLSI M39 Cumulative Susceptibility Matrix (%S, Wilson 95% CI, MIC50 / MIC90 ug/mL)",
                )
                .strong(),
            );
            egui::ScrollArea::vertical()
                .id_source("amr_cells_matrix_scroll")
                .max_height(205.0)
                .show(ui, |ui| {
                    egui::Grid::new("amr_cells_matrix_grid")
                        .striped(true)
                        .min_col_width(64.0)
                        .show(ui, |ui| {
                            ui.strong("Organism");
                            ui.strong("Drug");
                            ui.strong("Class");
                            ui.strong("N Tested");
                            ui.strong("%S");
                            ui.strong("Wilson 95% CI");
                            ui.strong("MIC50 / MIC90");
                            ui.strong("Guidance Flag");
                            ui.end_row();

                            let q = search_filter.trim().to_ascii_lowercase();
                            for c in &self.amr_result.cells {
                                if !q.is_empty()
                                    && !c.organism.to_ascii_lowercase().contains(&q)
                                    && !c.drug_code.to_ascii_lowercase().contains(&q)
                                    && !c.drug_name.to_ascii_lowercase().contains(&q)
                                    && !c.drug_class.to_ascii_lowercase().contains(&q)
                                {
                                    continue;
                                }
                                ui.monospace(&c.organism);
                                ui.monospace(format!("{} ({})", c.drug_code, c.drug_name));
                                ui.label(&c.drug_class);
                                ui.monospace(format!("{}", c.isolates_tested));
                                let s_col = if c.empiric_caution_below_threshold {
                                    Color32::from_rgb(244, 63, 94)
                                } else {
                                    Color32::from_rgb(52, 211, 153)
                                };
                                ui.label(
                                    RichText::new(format!("{:.1}%", c.pct_susceptible))
                                        .monospace()
                                        .color(s_col),
                                );
                                ui.monospace(format!(
                                    "[{:.1}%, {:.1}%]",
                                    c.wilson_ci_low_95, c.wilson_ci_high_95
                                ));
                                let mic_str = match (c.mic50_ug_ml, c.mic90_ug_ml) {
                                    (Some(m50), Some(m90)) => format!("{m50:.2} / {m90:.2}"),
                                    _ => "--".to_string(),
                                };
                                ui.monospace(mic_str);
                                let note = if !c.meets_clsi_m39_n30 {
                                    "CAUTION: N < 30"
                                } else if c.empiric_caution_below_threshold {
                                    "AVOID MONO-EMPIRIC"
                                } else {
                                    "EMPIRIC OK"
                                };
                                ui.monospace(note);
                                ui.end_row();
                            }
                        });
                });
        });
    }
}
