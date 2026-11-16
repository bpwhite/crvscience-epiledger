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

use eframe::egui::{self, Color32, Pos2, Rect, RichText, Sense, Stroke, Vec2};
use egui_plot::{Legend, Line, Plot, PlotPoints, Points};
use epi_analytics::{
    builtin_outbreak_fasta_and_metadata, builtin_respiratory_lineage_table,
    export_outbreak_graphviz_dot, fit_multinomial_lineage_logistics, parse_fasta_and_metadata,
    parse_lineage_counts_csv, reconstruct_outbreak_tree, LineageBasisMode, LineageCountTable,
    LineageLogisticsConfig, LineageLogisticsResult, OutbreakIsolateInput, OutbreakTreeConfig,
    OutbreakTreeResult,
};
use epi_core::{SurveillanceCampaignPreset, TimeSeriesDataset};
use std::collections::HashMap;

const CLADE_PALETTE: [Color32; 8] = [
    Color32::from_rgb(56, 189, 248),
    Color32::from_rgb(52, 211, 153),
    Color32::from_rgb(251, 191, 36),
    Color32::from_rgb(244, 114, 182),
    Color32::from_rgb(167, 139, 250),
    Color32::from_rgb(248, 113, 113),
    Color32::from_rgb(45, 212, 191),
    Color32::from_rgb(251, 146, 60),
];

pub struct GenomicWorkbenchState {
    // Multinomial lineage logistics state
    pub lineage_table: LineageCountTable,
    pub lineage_cfg: LineageLogisticsConfig,
    pub lineage_csv_path: String,
    pub cached_lineage: Option<LineageLogisticsResult>,

    // Outbreak phylodynamic tree state
    pub outbreak_isolates: Vec<OutbreakIsolateInput>,
    pub outbreak_cfg: OutbreakTreeConfig,
    pub fasta_path_input: String,
    pub metadata_path_input: String,
    pub export_prefix_input: String,
    pub cached_tree: Option<OutbreakTreeResult>,
}

impl Default for GenomicWorkbenchState {
    fn default() -> Self {
        Self::new()
    }
}

impl GenomicWorkbenchState {
    pub fn new() -> Self {
        let lineage_table = builtin_respiratory_lineage_table();
        let lineage_cfg = LineageLogisticsConfig::default();
        let cached_lineage = fit_multinomial_lineage_logistics(&lineage_table, &lineage_cfg).ok();

        let (fa, meta) = builtin_outbreak_fasta_and_metadata();
        let outbreak_isolates = parse_fasta_and_metadata(&fa, &meta).unwrap_or_default();
        let outbreak_cfg = OutbreakTreeConfig::default();
        let cached_tree = reconstruct_outbreak_tree(&outbreak_isolates, &outbreak_cfg).ok();

        Self {
            lineage_table,
            lineage_cfg,
            lineage_csv_path: String::new(),
            cached_lineage,
            outbreak_isolates,
            outbreak_cfg,
            fasta_path_input: String::new(),
            metadata_path_input: String::new(),
            export_prefix_input: "outbreak_export".to_string(),
            cached_tree,
        }
    }

    pub fn recompute_lineage(&mut self) -> Result<(), String> {
        if self.lineage_cfg.reference_clade_index >= self.lineage_table.clade_names.len() {
            self.lineage_cfg.reference_clade_index = 0;
        }
        match fit_multinomial_lineage_logistics(&self.lineage_table, &self.lineage_cfg) {
            Ok(res) => {
                self.cached_lineage = Some(res);
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    }

    pub fn recompute_tree(&mut self) -> Result<(), String> {
        match reconstruct_outbreak_tree(&self.outbreak_isolates, &self.outbreak_cfg) {
            Ok(res) => {
                self.cached_tree = Some(res);
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    }

    pub fn sync_to_campaign(
        &mut self,
        campaign: SurveillanceCampaignPreset,
        active_ds: &TimeSeriesDataset,
    ) {
        let mut table = builtin_respiratory_lineage_table();
        let tail_start = active_ds.points.len().saturating_sub(20);
        let tail = &active_ds.points[tail_start..];
        let mut new_clades = match campaign {
            SurveillanceCampaignPreset::UsMeasles2022_2026 => Some((
                "US Measles Genotype N-450 Clade Cohort (Recent 20w)".to_string(),
                vec![
                    "B3.1 (Ref)".to_string(),
                    "D8-MVs/TX.25".to_string(),
                    "D8-MVs/OH.26".to_string(),
                    "D8-Global".to_string(),
                    "B3-Import".to_string(),
                ],
            )),
            SurveillanceCampaignPreset::UsPertussis2022_2026 => Some((
                "US Bordetella pertussis Allelic & PRN-Deficient Clades (Recent 20w)".to_string(),
                vec![
                    "ptxP3/PRN+".to_string(),
                    "ptxP3/PRN-del".to_string(),
                    "fhaB-2400".to_string(),
                    "23S-A2047G".to_string(),
                    "ptxA1-Minor".to_string(),
                ],
            )),
            SurveillanceCampaignPreset::InfluenzaIliWastewater2022_2026 => Some((
                "US Influenza Type/Subclade Genomic Cohort (Recent 20w)".to_string(),
                vec![
                    "A(H1N1)pdm09 6B.1A".to_string(),
                    "A(H3N2) 2a.3a.1".to_string(),
                    "A(H3N2) J.2".to_string(),
                    "B/Victoria V1A.3a.2".to_string(),
                    "A(H1N1) 5a.2a.1".to_string(),
                ],
            )),
            SurveillanceCampaignPreset::SarsCoV2VariantCohort2022_2026 => None,
        };
        if let Some((title, clade_names)) = new_clades.take() {
            table.title = title;
            table.clade_names = clade_names;
            if tail.len() == table.rows.len() {
                for (r, pt) in table.rows.iter_mut().zip(tail.iter()) {
                    r.period_label = pt
                        .epiweek
                        .map(|w| w.compact_label())
                        .unwrap_or_else(|| pt.date_label.clone());
                    r.total_cases = pt.count.max(12.0);
                }
            }
        }
        self.lineage_table = table;
        self.lineage_cfg.reference_clade_index = 0;
        let _ = self.recompute_lineage();
    }

    pub fn render(&mut self, ui: &mut egui::Ui, status_banner: &mut String) {
        self.render_with_filter(ui, status_banner, "");
    }

    pub fn render_with_filter(
        &mut self,
        ui: &mut egui::Ui,
        status_banner: &mut String,
        search_filter: &str,
    ) {
        self.render_lineage_section(ui, status_banner, search_filter);
        ui.add_space(12.0);
        ui.separator();
        ui.add_space(8.0);
        self.render_outbreak_tree_section(ui, status_banner);
    }

    fn render_lineage_section(
        &mut self,
        ui: &mut egui::Ui,
        status_banner: &mut String,
        search_filter: &str,
    ) {
        ui.heading("Multinomial Lineage Competition, Spline Dynamics & Variant-Specific Rt");
        ui.label(
            RichText::new(
                "Penalized Newton-Raphson IRLS multinomial logistic regression & 3-df natural cubic splines with joint variant Cori Rt decomposition.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        if let Some(ref res) = self.cached_lineage {
            let dominant = res.clade_summaries.iter().max_by(|a, b| {
                a.forecast_end_share
                    .partial_cmp(&b.forecast_end_share)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let fastest = res
                .clade_summaries
                .iter()
                .filter(|c| !c.is_reference)
                .max_by(|a, b| {
                    a.growth_advantage_per_week
                        .partial_cmp(&b.growth_advantage_per_week)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            let latest_obs_pt = res.timeline.iter().rev().find(|p| !p.is_forecast);

            ui.columns(4, |cols| {
                cols[0].group(|ui| {
                    ui.label(RichText::new("PROJECTED DOMINANT CLADE").small().weak());
                    if let Some(dom) = dominant {
                        ui.heading(
                            RichText::new(format!(
                                "{} ({:.1}%)",
                                dom.clade_name,
                                dom.forecast_end_share * 100.0
                            ))
                            .monospace()
                            .color(Color32::from_rgb(52, 211, 153)),
                        );
                        ui.label(
                            RichText::new(format!(
                                "Current fitted: {:.1}%",
                                dom.latest_fitted_share * 100.0
                            ))
                            .small()
                            .monospace(),
                        );
                    }
                });

                cols[1].group(|ui| {
                    ui.label(RichText::new("TOP GROWTH ADVANTAGE").small().weak());
                    if let Some(fast) = fastest {
                        ui.heading(
                            RichText::new(format!(
                                "{} {:+.2}/wk",
                                fast.clade_name, fast.growth_advantage_per_week
                            ))
                            .monospace()
                            .color(Color32::from_rgb(56, 189, 248)),
                        );
                        ui.label(
                            RichText::new(format!(
                                "Per-gen {:+.1}% vs {}",
                                fast.multiplicative_advantage_pct, res.reference_clade
                            ))
                            .small()
                            .monospace(),
                        );
                    }
                });

                cols[2].group(|ui| {
                    ui.label(RichText::new("SHANNON CLADE DIVERSITY").small().weak());
                    if let Some(pt) = latest_obs_pt {
                        ui.heading(
                            RichText::new(format!("{:.2} bits", pt.shannon_bits))
                                .monospace()
                                .color(Color32::from_rgb(251, 191, 36)),
                        );
                        ui.label(
                            RichText::new(format!(
                                "Effective circulating clades: {:.2}",
                                pt.effective_clades
                            ))
                            .small()
                            .monospace(),
                        );
                    }
                });

                cols[3].group(|ui| {
                    ui.label(RichText::new("MULTINOMIAL FIT DIAGNOSTICS").small().weak());
                    ui.heading(RichText::new(format!("AIC {:.1}", res.aic)).monospace());
                    ui.label(
                        RichText::new(format!(
                            "LogLik {:.1} | {} ({} iters)",
                            res.log_likelihood,
                            if res.converged { "Converged" } else { "Capped" },
                            res.iterations_used
                        ))
                        .small()
                        .monospace(),
                    );
                });
            });
        }

        ui.add_space(6.0);
        let mut need_refit = false;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("Basis:");
                if ui
                    .selectable_label(
                        self.lineage_cfg.basis_mode == LineageBasisMode::LogLinearFitness,
                        "Log-Linear Fitness",
                    )
                    .clicked()
                {
                    self.lineage_cfg.basis_mode = LineageBasisMode::LogLinearFitness;
                    need_refit = true;
                }
                if ui
                    .selectable_label(
                        self.lineage_cfg.basis_mode == LineageBasisMode::NaturalCubicSpline,
                        "Natural Cubic Spline (df=3)",
                    )
                    .clicked()
                {
                    self.lineage_cfg.basis_mode = LineageBasisMode::NaturalCubicSpline;
                    need_refit = true;
                }

                ui.separator();
                ui.label("Ref Clade:");
                let cur_ref_name = self
                    .lineage_table
                    .clade_names
                    .get(self.lineage_cfg.reference_clade_index)
                    .cloned()
                    .unwrap_or_default();
                egui::ComboBox::from_id_source("ref_clade_combo")
                    .selected_text(cur_ref_name)
                    .show_ui(ui, |ui| {
                        for (idx, cname) in self.lineage_table.clade_names.iter().enumerate() {
                            if ui
                                .selectable_value(
                                    &mut self.lineage_cfg.reference_clade_index,
                                    idx,
                                    cname,
                                )
                                .changed()
                            {
                                need_refit = true;
                            }
                        }
                    });

                ui.separator();
                if ui
                    .add(
                        egui::Slider::new(&mut self.lineage_cfg.forecast_weeks, 1..=16)
                            .text("Forecast Wks"),
                    )
                    .changed()
                {
                    need_refit = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.lineage_cfg.generation_time_days, 2.0..=16.0)
                            .text("Gen Time (d)"),
                    )
                    .changed()
                {
                    need_refit = true;
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("Lineage Count CSV/TSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.lineage_csv_path)
                        .desired_width(260.0)
                        .hint_text("path/to/clade_counts.csv (wide or tidy long)"),
                );
                if ui.button("Load Clade CSV").clicked() {
                    let p = self.lineage_csv_path.trim().to_string();
                    match std::fs::read_to_string(&p)
                        .map_err(|e| e.to_string())
                        .and_then(|txt| {
                            parse_lineage_counts_csv(&txt, &p).map_err(|e| e.to_string())
                        }) {
                        Ok(tbl) => {
                            self.lineage_table = tbl;
                            self.lineage_cfg.reference_clade_index = 0;
                            need_refit = true;
                            *status_banner = format!(
                                "Loaded clade count table '{}' ({} clades, {} periods).",
                                self.lineage_table.title,
                                self.lineage_table.clade_names.len(),
                                self.lineage_table.rows.len()
                            );
                        }
                        Err(err) => {
                            *status_banner = format!("Clade CSV import error: {err}");
                        }
                    }
                }
                if ui.button("Reset Built-in Benchmark").clicked() {
                    self.lineage_table = builtin_respiratory_lineage_table();
                    self.lineage_cfg.reference_clade_index = 0;
                    need_refit = true;
                    *status_banner =
                        "Restored built-in 5-clade respiratory surveillance benchmark.".to_string();
                }
            });
        });

        if need_refit {
            if let Err(err) = self.recompute_lineage() {
                *status_banner = format!("Lineage fit warning: {err}");
            }
        }

        if let Some(ref res) = self.cached_lineage {
            ui.columns(2, |cols| {
                cols[0].group(|ui| {
                    ui.strong("Clade Share Trajectories & Forecast (%)");
                    Plot::new("lineage_shares_plot")
                        .height(230.0)
                        .legend(Legend::default())
                        .show(ui, |plot_ui| {
                            for (k, cname) in res.clade_names.iter().enumerate() {
                                let col = CLADE_PALETTE[k % CLADE_PALETTE.len()];
                                let fit_vec: Vec<[f64; 2]> = res
                                    .timeline
                                    .iter()
                                    .map(|pt| [pt.step_index as f64, pt.fitted_shares[k] * 100.0])
                                    .collect();
                                plot_ui.line(
                                    Line::new(PlotPoints::from(fit_vec))
                                        .name(cname)
                                        .color(col)
                                        .width(2.0),
                                );

                                let obs_pts: Vec<[f64; 2]> = res
                                    .timeline
                                    .iter()
                                    .filter_map(|pt| {
                                        pt.observed_shares[k]
                                            .map(|s| [pt.step_index as f64, s * 100.0])
                                    })
                                    .collect();
                                if !obs_pts.is_empty() {
                                    plot_ui.points(
                                        Points::new(PlotPoints::from(obs_pts))
                                            .name(format!("{cname} obs"))
                                            .color(col)
                                            .radius(3.0),
                                    );
                                }
                            }
                        });
                });

                cols[1].group(|ui| {
                    ui.strong("Joint Variant-Specific Effective Reproduction Number R(t, k)");
                    Plot::new("variant_rt_plot")
                        .height(230.0)
                        .legend(Legend::default())
                        .show(ui, |plot_ui| {
                            let max_x = res.timeline.len().saturating_sub(1) as f64;
                            plot_ui.line(
                                Line::new(PlotPoints::from(vec![[0.0, 1.0], [max_x, 1.0]]))
                                    .name("Rt = 1.0 threshold")
                                    .color(Color32::from_gray(130))
                                    .width(1.0),
                            );
                            for (k, cname) in res.clade_names.iter().enumerate() {
                                let col = CLADE_PALETTE[k % CLADE_PALETTE.len()];
                                let rt_pts: Vec<[f64; 2]> = res
                                    .timeline
                                    .iter()
                                    .filter_map(|pt| {
                                        pt.variant_rt[k].map(|est| [pt.step_index as f64, est.mean])
                                    })
                                    .collect();
                                if !rt_pts.is_empty() {
                                    plot_ui.line(
                                        Line::new(PlotPoints::from(rt_pts))
                                            .name(format!("{cname} Rt"))
                                            .color(col)
                                            .width(2.0),
                                    );
                                }
                            }
                        });
                });
            });

            ui.group(|ui| {
                ui.strong("Clade Growth Advantage & Crossover Summary");
                egui::Grid::new("clade_summary_grid")
                    .striped(true)
                    .min_col_width(82.0)
                    .show(ui, |ui| {
                        ui.strong("Clade");
                        ui.strong("Sequences");
                        ui.strong("Latest Obs");
                        ui.strong("Forecast End");
                        ui.strong("Growth s/wk [95% CI]");
                        ui.strong("Per-Gen Adv");
                        ui.strong("Crossover Wk");
                        ui.strong("Latest Variant Rt");
                        ui.end_row();

                        let q = search_filter.trim().to_ascii_lowercase();
                        for c in &res.clade_summaries {
                            if !q.is_empty() && !c.clade_name.to_ascii_lowercase().contains(&q) {
                                continue;
                            }
                            ui.monospace(&c.clade_name);
                            ui.monospace(format!("{:.0}", c.total_sequences));
                            ui.monospace(format!("{:.1}%", c.latest_observed_share * 100.0));
                            ui.monospace(format!("{:.1}%", c.forecast_end_share * 100.0));
                            if c.is_reference {
                                ui.monospace("0.000 (Reference)");
                            } else {
                                ui.monospace(format!(
                                    "{:+.3} [{:+.3}, {:+.3}]",
                                    c.growth_advantage_per_week,
                                    c.growth_ci_low_95,
                                    c.growth_ci_high_95
                                ));
                            }
                            ui.monospace(format!("{:+.1}%", c.multiplicative_advantage_pct));
                            ui.monospace(
                                c.crossover_week_vs_ref
                                    .map(|w| format!("Wk {w:.1}"))
                                    .unwrap_or_else(|| "--".to_string()),
                            );
                            ui.monospace(
                                c.latest_variant_rt
                                    .map(|r| {
                                        format!(
                                            "{:.2} [{:.2}, {:.2}]",
                                            r.mean, r.cri_low_95, r.cri_high_95
                                        )
                                    })
                                    .unwrap_or_else(|| "--".to_string()),
                            );
                            ui.end_row();
                        }
                    });
            });
        }
    }

    fn render_outbreak_tree_section(&mut self, ui: &mut egui::Ui, status_banner: &mut String) {
        ui.heading("Joint Genomic (FASTA SNP / K2P) + Symptom-Onset Outbreak Transmission Tree");
        ui.label(
            RichText::new(
                "Outbreaker2 / seqTrack maximum-a-posteriori transmission forest with >=90% unambiguous site QC, unsampled intermediate host inference, and multi-introduction separation.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        let mut need_tree = false;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                if ui
                    .checkbox(
                        &mut self.outbreak_cfg.exclude_quarantined,
                        "Exclude Quarantined",
                    )
                    .changed()
                {
                    need_tree = true;
                }
                if ui
                    .add(
                        egui::Slider::new(
                            &mut self.outbreak_cfg.min_unambiguous_coverage,
                            0.50..=1.00,
                        )
                        .text("Min ACGT Coverage"),
                    )
                    .changed()
                {
                    need_tree = true;
                }
                if ui
                    .add(
                        egui::Slider::new(
                            &mut self.outbreak_cfg.generation_mean_days,
                            2.0..=21.0,
                        )
                        .text("Gen Mean (d)"),
                    )
                    .changed()
                {
                    need_tree = true;
                }
                if ui
                    .add(
                        egui::Slider::new(
                            &mut self.outbreak_cfg.mutation_rate_snps_per_day,
                            0.005..=0.25,
                        )
                        .text("Clock SNPs/d"),
                    )
                    .changed()
                {
                    need_tree = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.outbreak_cfg.max_kappa_hops, 1..=5)
                            .text("Max Hops k"),
                    )
                    .changed()
                {
                    need_tree = true;
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("FASTA:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.fasta_path_input)
                        .desired_width(170.0)
                        .hint_text("aligned_isolates.fasta"),
                );
                ui.label("Metadata CSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.metadata_path_input)
                        .desired_width(170.0)
                        .hint_text("onset_line_list.csv"),
                );
                if ui.button("Load FASTA + Line List").clicked() {
                    let fa_p = self.fasta_path_input.trim().to_string();
                    let me_p = self.metadata_path_input.trim().to_string();
                    match std::fs::read_to_string(&fa_p)
                        .and_then(|fa| std::fs::read_to_string(&me_p).map(|me| (fa, me)))
                        .map_err(|e| e.to_string())
                        .and_then(|(fa, me)| {
                            parse_fasta_and_metadata(&fa, &me).map_err(|e| e.to_string())
                        }) {
                        Ok(isos) => {
                            let cnt = isos.len();
                            self.outbreak_isolates = isos;
                            need_tree = true;
                            *status_banner =
                                format!("Loaded {cnt} aligned FASTA outbreak isolates.");
                        }
                        Err(err) => {
                            *status_banner = format!("Outbreak import error: {err}");
                        }
                    }
                }
                if ui.button("Reset Measles D8 Outbreak").clicked() {
                    let (fa, me) = builtin_outbreak_fasta_and_metadata();
                    if let Ok(isos) = parse_fasta_and_metadata(&fa, &me) {
                        self.outbreak_isolates = isos;
                        need_tree = true;
                        *status_banner =
                            "Restored built-in 14-isolate Measles D8 multi-chain outbreak."
                                .to_string();
                    }
                }

                ui.separator();
                ui.add(
                    egui::TextEdit::singleline(&mut self.export_prefix_input)
                        .desired_width(120.0)
                        .hint_text("export_prefix"),
                );
                if ui.button("Export .nwk + .dot").clicked() {
                    if let Some(ref tree) = self.cached_tree {
                        let pref = self.export_prefix_input.trim();
                        let nwk_path = format!("{pref}.nwk");
                        let dot_path = format!("{pref}.dot");
                        let dot_txt = export_outbreak_graphviz_dot(tree);
                        match std::fs::write(&nwk_path, &tree.newick_tree)
                            .and_then(|_| std::fs::write(&dot_path, dot_txt))
                        {
                            Ok(()) => {
                                *status_banner =
                                    format!("Exported Newick ({nwk_path}) and DOT ({dot_path}).");
                            }
                            Err(e) => {
                                *status_banner = format!("Tree export error: {e}");
                            }
                        }
                    }
                }
            });
        });

        if need_tree {
            if let Err(err) = self.recompute_tree() {
                *status_banner = format!("Outbreak tree warning: {err}");
            }
        }

        if let Some(ref tree) = self.cached_tree {
            ui.columns(4, |cols| {
                cols[0].group(|ui| {
                    ui.label(RichText::new("RETAINED ISOLATES (QC)").small().weak());
                    ui.heading(
                        RichText::new(format!(
                            "{} ({} bp)",
                            tree.retained_cases, tree.alignment_length
                        ))
                        .monospace(),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Quarantined={} | LowCov={}",
                            tree.quarantined_excluded_ids.len(),
                            tree.low_coverage_excluded_ids.len()
                        ))
                        .small()
                        .monospace()
                        .color(Color32::from_rgb(251, 191, 36)),
                    );
                });
                cols[1].group(|ui| {
                    ui.label(RichText::new("INDEX INTRODUCTIONS").small().weak());
                    ui.heading(
                        RichText::new(format!("{} clusters", tree.total_importations))
                            .monospace()
                            .color(Color32::from_rgb(45, 212, 191)),
                    );
                    ui.label(RichText::new("History-host roots").small().monospace());
                });
                cols[2].group(|ui| {
                    ui.label(RichText::new("UNSAMPLED INTERMEDIATES").small().weak());
                    ui.heading(
                        RichText::new(format!("{} hosts", tree.total_unsampled_intermediates))
                            .monospace()
                            .color(Color32::from_rgb(251, 146, 60)),
                    );
                    ui.label(RichText::new("Inferred multi-hop links").small().monospace());
                });
                cols[3].group(|ui| {
                    ui.label(RichText::new("CASE REPRODUCTION Ri").small().weak());
                    ui.heading(
                        RichText::new(format!("Mean {:.2}", tree.mean_case_ri))
                            .monospace()
                            .color(Color32::from_rgb(56, 189, 248)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Max secondary: {} cases",
                            tree.max_superspreader_cases
                        ))
                        .small()
                        .monospace(),
                    );
                });
            });

            if !tree.quarantined_excluded_ids.is_empty()
                || !tree.low_coverage_excluded_ids.is_empty()
            {
                let low_str: Vec<String> = tree
                    .low_coverage_excluded_ids
                    .iter()
                    .map(|(id, c)| format!("{id} ({:.0}%)", c * 100.0))
                    .collect();
                ui.label(
                    RichText::new(format!(
                        "Sequence QC Gate Audit — Excluded Quarantined: [{}] | Excluded Below {:.0}% Unambiguous Sites: [{}]",
                        if tree.quarantined_excluded_ids.is_empty() {
                            "none".to_string()
                        } else {
                            tree.quarantined_excluded_ids.join(", ")
                        },
                        self.outbreak_cfg.min_unambiguous_coverage * 100.0,
                        if low_str.is_empty() {
                            "none".to_string()
                        } else {
                            low_str.join(", ")
                        }
                    ))
                    .small()
                    .monospace()
                    .color(Color32::from_rgb(251, 191, 36)),
                );
            }

            ui.add_space(4.0);
            self.draw_transmission_network_canvas(ui, tree);

            ui.add_space(6.0);
            ui.columns(2, |cols| {
                cols[0].group(|ui| {
                    ui.strong("MAP Transmission Links & Posterior Ancestor Support");
                    egui::ScrollArea::vertical().max_height(210.0).show(ui, |ui| {
                        egui::Grid::new("outbreak_edges_grid")
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Infector");
                                ui.strong("Infectee");
                                ui.strong("Location");
                                ui.strong("Lag (d)");
                                ui.strong("SNPs");
                                ui.strong("Unsamp");
                                ui.strong("Support");
                                ui.end_row();

                                for e in &tree.edges {
                                    let inf = e
                                        .infector_id
                                        .clone()
                                        .unwrap_or_else(|| format!("[Import #{}]", e.cluster_id));
                                    ui.monospace(inf);
                                    ui.monospace(&e.infectee_id);
                                    ui.monospace(&e.infectee_location);
                                    ui.monospace(
                                        e.onset_lag_days
                                            .map(|d| format!("{d:.0}"))
                                            .unwrap_or_else(|| "--".to_string()),
                                    );
                                    ui.monospace(
                                        e.snp_distance
                                            .map(|s| s.to_string())
                                            .unwrap_or_else(|| "--".to_string()),
                                    );
                                    ui.monospace(format!("{}", e.unsampled_intermediates));
                                    ui.monospace(format!("{:.1}%", e.posterior_support * 100.0));
                                    ui.end_row();
                                }
                            });
                    });
                });

                cols[1].group(|ui| {
                    ui.strong("Case-Level Expected Reproduction Number (Ri) & Newick Forest");
                    egui::ScrollArea::vertical().max_height(135.0).show(ui, |ui| {
                        egui::Grid::new("case_ri_grid").striped(true).show(ui, |ui| {
                            ui.strong("Case ID");
                            ui.strong("Onset");
                            ui.strong("Location");
                            ui.strong("MAP Kids");
                            ui.strong("Expected Ri");
                            ui.end_row();
                            for c in &tree.case_offspring {
                                ui.monospace(&c.case_id);
                                ui.monospace(&c.onset_label);
                                ui.monospace(&c.location);
                                ui.monospace(format!("{}", c.map_secondary_cases));
                                ui.monospace(format!("{:.2}", c.expected_ri));
                                ui.end_row();
                            }
                        });
                    });
                    ui.separator();
                    ui.label(RichText::new("Newick Outbreak Representation:").small().weak());
                    let mut nwk_preview = tree.newick_tree.clone();
                    ui.add(
                        egui::TextEdit::multiline(&mut nwk_preview)
                            .desired_rows(2)
                            .font(egui::TextStyle::Monospace)
                            .interactive(false),
                    );
                });
            });
        }
    }

    fn draw_transmission_network_canvas(&self, ui: &mut egui::Ui, tree: &OutbreakTreeResult) {
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.strong("Chronological Transmission Forest Canvas (X = Symptom Onset Day, Y = Introduction Cluster)");
                ui.label(
                    RichText::new("Teal = Index Importation | Solid Sky = Direct Link | Orange Link = Unsampled Intermediate Host (k >= 2)")
                        .small()
                        .weak(),
                );
            });

            let desired = Vec2::new(ui.available_width().max(480.0), 220.0);
            let (rect, _resp) = ui.allocate_exact_size(desired, Sense::hover());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 6.0, Color32::from_rgb(15, 23, 42));

            let n = tree.edges.len();
            if n == 0 {
                return;
            }

            let max_day = tree
                .edges
                .iter()
                .map(|e| e.infectee_onset_day)
                .fold(1.0_f64, f64::max);

            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by(|&a, &b| {
                tree.edges[a]
                    .cluster_id
                    .cmp(&tree.edges[b].cluster_id)
                    .then_with(|| {
                        tree.edges[a]
                            .infectee_onset_day
                            .partial_cmp(&tree.edges[b].infectee_onset_day)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
            });

            let pad_x = 54.0_f32;
            let pad_y = 24.0_f32;
            let plot_r = Rect::from_min_max(
                Pos2::new(rect.left() + pad_x, rect.top() + pad_y),
                Pos2::new(rect.right() - pad_x, rect.bottom() - pad_y),
            );

            let mut pos_by_id: HashMap<String, Pos2> = HashMap::new();
            for (rank, &idx) in order.iter().enumerate() {
                let e = &tree.edges[idx];
                let fx = (e.infectee_onset_day / max_day).clamp(0.0, 1.0) as f32;
                let fy = if n > 1 {
                    (rank as f32) / ((n - 1) as f32)
                } else {
                    0.5
                };
                let pt = Pos2::new(
                    plot_r.left() + fx * plot_r.width(),
                    plot_r.top() + fy * plot_r.height(),
                );
                pos_by_id.insert(e.infectee_id.clone(), pt);
            }

            for e in &tree.edges {
                if let Some(ref inf_id) = e.infector_id {
                    if let (Some(&p_src), Some(&p_dst)) =
                        (pos_by_id.get(inf_id), pos_by_id.get(&e.infectee_id))
                    {
                        let edge_col = if e.unsampled_intermediates > 0 {
                            Color32::from_rgb(251, 146, 60)
                        } else {
                            Color32::from_rgb(56, 189, 248)
                        };
                        let width = (1.2 + 1.8 * e.posterior_support as f32).clamp(1.2, 3.0);
                        let mid_x = 0.5 * (p_src.x + p_dst.x);
                        painter.line_segment(
                            [p_src, Pos2::new(mid_x, p_dst.y)],
                            Stroke::new(width, edge_col),
                        );
                        painter.line_segment(
                            [Pos2::new(mid_x, p_dst.y), p_dst],
                            Stroke::new(width, edge_col),
                        );

                        if e.unsampled_intermediates > 0 {
                            let mid_pt = Pos2::new(mid_x, p_dst.y);
                            painter.circle_stroke(mid_pt, 4.5, Stroke::new(1.5, edge_col));
                        }
                    }
                }
            }

            for e in &tree.edges {
                if let Some(&pt) = pos_by_id.get(&e.infectee_id) {
                    let fill = if e.is_importation {
                        Color32::from_rgb(45, 212, 191)
                    } else {
                        CLADE_PALETTE[(e.cluster_id.saturating_sub(1)) % CLADE_PALETTE.len()]
                    };
                    painter.circle_filled(pt, 5.5, fill);
                    painter.circle_stroke(pt, 5.5, Stroke::new(1.0, Color32::WHITE));
                    painter.text(
                        Pos2::new(pt.x + 8.0, pt.y),
                        egui::Align2::LEFT_CENTER,
                        &e.infectee_id,
                        egui::FontId::monospace(10.5),
                        Color32::from_gray(230),
                    );
                }
            }
        });
    }
}
