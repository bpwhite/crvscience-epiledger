use eframe::egui::{self, Color32, RichText, Stroke};
use egui_plot::{Bar, BarChart, HLine, Legend, Line, Plot, PlotPoints};
use epi_analytics::{
    analyze_contact_tracing_network, builtin_superspreading_outbreak_records,
    export_bipartite_network_dot, fit_negative_binomial_offspring, parse_contact_tracing_csv,
    BipartiteNetworkAnalysis, ContactTracingRecord, OffspringDispersionResult, OffspringFitConfig,
    ZeroTruncationMode,
};

#[derive(Default)]
pub struct Phase9GuiAction {
    pub apply_si_to_cori_rt: Option<(f64, f64)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RightPlotKind {
    OffspringPmfComparison,
    ProfileLogLikelihood,
}

pub struct SuperspreadingWorkbenchState {
    pub records: Vec<ContactTracingRecord>,
    pub fit_cfg: OffspringFitConfig,
    pub analysis: BipartiteNetworkAnalysis,
    pub active_offspring: OffspringDispersionResult,
    pub using_custom_counts: bool,
    pub custom_counts_input: String,
    pub csv_path_input: String,
    pub dot_export_path: String,
    pub right_plot: RightPlotKind,
}

impl Default for SuperspreadingWorkbenchState {
    fn default() -> Self {
        Self::new()
    }
}

impl SuperspreadingWorkbenchState {
    pub fn new() -> Self {
        let records = builtin_superspreading_outbreak_records();
        let fit_cfg = OffspringFitConfig::default();
        let analysis = analyze_contact_tracing_network(&records, &fit_cfg)
            .expect("Built-in superspreading outbreak records are valid");
        let active_offspring = analysis.offspring_fit.clone();
        Self {
            records,
            fit_cfg,
            analysis,
            active_offspring,
            using_custom_counts: false,
            custom_counts_input:
                "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,2,2,6,9,14".to_string(),
            csv_path_input: String::new(),
            dot_export_path: "superspreading_bipartite.dot".to_string(),
            right_plot: RightPlotKind::OffspringPmfComparison,
        }
    }

    pub fn recompute_network(&mut self) -> Result<(), String> {
        let net = analyze_contact_tracing_network(&self.records, &self.fit_cfg)
            .map_err(|e| e.to_string())?;
        self.active_offspring = net.offspring_fit.clone();
        self.analysis = net;
        self.using_custom_counts = false;
        Ok(())
    }

    pub fn fit_custom_offspring_vector(&mut self) -> Result<(), String> {
        let mut parsed = Vec::new();
        for tok in self
            .custom_counts_input
            .split(|c: char| c == ',' || c.is_whitespace() || c == ';')
        {
            let t = tok.trim();
            if t.is_empty() {
                continue;
            }
            let val: u32 = t
                .parse()
                .map_err(|_| format!("Invalid non-negative integer offspring count '{t}'"))?;
            parsed.push(val);
        }
        let fit =
            fit_negative_binomial_offspring(&parsed, &self.fit_cfg).map_err(|e| e.to_string())?;
        self.active_offspring = fit;
        self.using_custom_counts = true;
        Ok(())
    }

    pub fn render(&mut self, ui: &mut egui::Ui, status_banner: &mut String) -> Phase9GuiAction {
        let mut action = Phase9GuiAction::default();

        ui.heading(
            "Superspreading Offspring Dispersion (k MLE) & Bipartite Contact-Tracing Networks",
        );
        ui.label(
            RichText::new(
                "Lloyd-Smith Negative Binomial offspring R0 and dispersion k MLE with profile-likelihood 95% CIs, Pareto 80/20 share, Matula-Beck k-core decomposition, and Brandes betweenness centrality.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        let off = &self.active_offspring;
        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("DISPERSION k [95% PROFILE CI]").small().weak());
                let k_col = if off.k_mle < 0.5 {
                    Color32::from_rgb(244, 63, 94)
                } else if off.k_mle < 1.5 {
                    Color32::from_rgb(251, 191, 36)
                } else {
                    Color32::from_rgb(52, 211, 153)
                };
                ui.heading(
                    RichText::new(format!("k = {:.3}", off.k_mle))
                        .monospace()
                        .color(k_col),
                );
                let hi_str = if off.k_ci_high_95.is_finite() {
                    format!("{:.3}", off.k_ci_high_95)
                } else {
                    "inf".to_string()
                };
                ui.label(
                    RichText::new(format!(
                        "CI [{:.3}, {}] | p={:.2e}",
                        off.k_ci_low_95, hi_str, off.lrt_p_value
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("PARETO 80/20 TRANSMISSION SHARE").small().weak());
                ui.heading(
                    RichText::new(format!(
                        "{:.1}% from top 20%",
                        off.prop_transmission_top_20_theoretical * 100.0
                    ))
                    .monospace()
                    .color(Color32::from_rgb(56, 189, 248)),
                );
                ui.label(
                    RichText::new(format!(
                        "Empirical {:.1}% | 80% tx by {:.1}% cases",
                        off.prop_transmission_top_20_empirical * 100.0,
                        off.cases_driving_80_pct_transmission * 100.0
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[2].group(|ui| {
                ui.label(RichText::new("R0 MLE & EXTINCTION PROBABILITY").small().weak());
                ui.heading(
                    RichText::new(format!(
                        "R0 = {:.2} | q_ext = {:.1}%",
                        off.r0_mle,
                        off.single_intro_extinction_prob * 100.0
                    ))
                    .monospace()
                    .color(Color32::from_rgb(167, 139, 250)),
                );
                ui.label(
                    RichText::new(format!(
                        "99% Poisson SSE > {} ({} obs)",
                        off.sse_99_threshold, off.observed_sse_count
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[3].group(|ui| {
                ui.label(RichText::new("BIPARTITE GRAPH & SERIAL INTERVAL").small().weak());
                ui.heading(
                    RichText::new(format!(
                        "{} cases / {} venues",
                        self.analysis.total_cases, self.analysis.total_venues
                    ))
                    .monospace()
                    .color(Color32::from_rgb(52, 211, 153)),
                );
                ui.label(
                    RichText::new(format!(
                        "Max k-core={} | SI={:.2}d +/- {:.2}d",
                        self.analysis.max_bipartite_k_core,
                        self.analysis.serial_interval.mean_days,
                        self.analysis.serial_interval.sd_days
                    ))
                    .small()
                    .monospace(),
                );
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Observation Regime:").strong());
                let prev_mode = self.fit_cfg.mode;
                egui::ComboBox::from_id_source("offspring_zero_trunc_mode")
                    .selected_text(match self.fit_cfg.mode {
                        ZeroTruncationMode::StandardAllCases => "All Index Cases (z >= 0)",
                        ZeroTruncationMode::ZeroTruncatedClustersOnly => {
                            "Zero-Truncated Clusters (z >= 1, ZTNB)"
                        }
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.fit_cfg.mode,
                            ZeroTruncationMode::StandardAllCases,
                            "All Index Cases (z >= 0)",
                        );
                        ui.selectable_value(
                            &mut self.fit_cfg.mode,
                            ZeroTruncationMode::ZeroTruncatedClustersOnly,
                            "Zero-Truncated Clusters (z >= 1, ZTNB)",
                        );
                    });
                if self.fit_cfg.mode != prev_mode {
                    let res = if self.using_custom_counts {
                        self.fit_custom_offspring_vector()
                    } else {
                        self.recompute_network()
                    };
                    if let Err(e) = res {
                        *status_banner = format!("Offspring refit warning: {e}");
                    }
                }

                if ui.button("Load Built-In Multi-Venue Outbreak").clicked() {
                    self.records = builtin_superspreading_outbreak_records();
                    if let Err(e) = self.recompute_network() {
                        *status_banner = format!("Built-in network error: {e}");
                    } else {
                        *status_banner = format!(
                            "Loaded built-in respiratory superspreading network ({} cases, {} venues, k={:.3}).",
                            self.analysis.total_cases,
                            self.analysis.total_venues,
                            self.active_offspring.k_mle
                        );
                    }
                }

                if ui
                    .button("Apply Empirical SI to Cori Rt")
                    .on_hover_text("Push positive-shifted empirical serial interval mean & SD into Tab 1 Cori Rt")
                    .clicked()
                {
                    action.apply_si_to_cori_rt = Some((
                        self.analysis.serial_interval.positive_gamma_mean_days,
                        self.analysis.serial_interval.positive_gamma_sd_days,
                    ));
                }

                ui.separator();
                ui.label("Right Chart:");
                ui.radio_value(
                    &mut self.right_plot,
                    RightPlotKind::OffspringPmfComparison,
                    "Offspring PMF",
                );
                ui.radio_value(
                    &mut self.right_plot,
                    RightPlotKind::ProfileLogLikelihood,
                    "Profile LL(k)",
                );
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("Contact CSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.csv_path_input)
                        .desired_width(210.0)
                        .hint_text("path/to/contact_pairs.csv"),
                );
                if ui.button("Import Line-List CSV").clicked() {
                    match std::fs::read_to_string(self.csv_path_input.trim()) {
                        Ok(raw) => match parse_contact_tracing_csv(&raw) {
                            Ok(recs) => {
                                self.records = recs;
                                match self.recompute_network() {
                                    Ok(()) => {
                                        *status_banner = format!(
                                            "Imported {} contact links ({} cases, {} venues, k={:.3}).",
                                            self.records.len(),
                                            self.analysis.total_cases,
                                            self.analysis.total_venues,
                                            self.active_offspring.k_mle
                                        );
                                    }
                                    Err(e) => *status_banner = format!("Network fit error: {e}"),
                                }
                            }
                            Err(e) => *status_banner = format!("CSV parse error: {e}"),
                        },
                        Err(e) => *status_banner = format!("File read error: {e}"),
                    }
                }

                ui.separator();
                ui.label("Custom z_i:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.custom_counts_input)
                        .desired_width(220.0),
                );
                if ui.button("Fit Offspring z_i").clicked() {
                    match self.fit_custom_offspring_vector() {
                        Ok(()) => {
                            *status_banner = format!(
                                "Fitted custom offspring vector (n={}, R0={:.2}, k={:.3}, top-20% tx={:.1}%).",
                                self.active_offspring.n_cases_evaluated,
                                self.active_offspring.r0_mle,
                                self.active_offspring.k_mle,
                                self.active_offspring.prop_transmission_top_20_theoretical * 100.0
                            );
                        }
                        Err(e) => *status_banner = format!("Custom z_i error: {e}"),
                    }
                }

                ui.separator();
                ui.add(
                    egui::TextEdit::singleline(&mut self.dot_export_path).desired_width(165.0),
                );
                if ui.button("Export .DOT").clicked() {
                    let dot = export_bipartite_network_dot(&self.analysis, &self.records);
                    match std::fs::write(self.dot_export_path.trim(), dot) {
                        Ok(()) => {
                            *status_banner = format!(
                                "Exported bipartite Graphviz .dot to '{}'.",
                                self.dot_export_path.trim()
                            );
                        }
                        Err(e) => *status_banner = format!("DOT export failed: {e}"),
                    }
                }
            });
        });

        ui.add_space(6.0);

        let off = &self.active_offspring;
        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.label(
                    RichText::new(
                        "Cumulative Transmission Concentration (Lloyd-Smith Gamma vs Empirical Lorenz)",
                    )
                    .strong(),
                );
                let emp_pts: Vec<[f64; 2]> = off
                    .empirical_lorenz_curve
                    .iter()
                    .map(|p| [p[0] * 100.0, p[1] * 100.0])
                    .collect();
                let theo_pts: Vec<[f64; 2]> = off
                    .theoretical_lorenz_curve
                    .iter()
                    .map(|p| [p[0] * 100.0, p[1] * 100.0])
                    .collect();
                let diag_pts = vec![[0.0, 0.0], [100.0, 100.0]];

                Plot::new("superspreading_lorenz_plot")
                    .height(235.0)
                    .legend(Legend::default())
                    .x_axis_label("Top % Most Infectious Cases")
                    .y_axis_label("Cumulative % of Secondary Transmissions")
                    .show(ui, |plot_ui| {
                        plot_ui.line(
                            Line::new(PlotPoints::from(theo_pts))
                                .name(format!("Lloyd-Smith NB Theory (k={:.2})", off.k_mle))
                                .color(Color32::from_rgb(56, 189, 248))
                                .width(2.4),
                        );
                        plot_ui.line(
                            Line::new(PlotPoints::from(emp_pts))
                                .name("Observed Line-List Lorenz")
                                .color(Color32::from_rgb(244, 63, 94))
                                .width(2.0),
                        );
                        plot_ui.line(
                            Line::new(PlotPoints::from(diag_pts))
                                .name("Homogeneous Equality (k -> inf)")
                                .stroke(Stroke::new(1.2, Color32::from_rgb(148, 163, 184))),
                        );
                    });
            });

            cols[1].group(|ui| match self.right_plot {
                RightPlotKind::OffspringPmfComparison => {
                    ui.label(
                        RichText::new(
                            "Offspring Distribution: Observed vs Negative Binomial vs Poisson",
                        )
                        .strong(),
                    );
                    let bars: Vec<Bar> = off
                        .histogram
                        .iter()
                        .map(|b| {
                            Bar::new(f64::from(b.offspring_count), b.observed_prob * 100.0)
                                .width(0.62)
                                .fill(Color32::from_rgba_unmultiplied(96, 165, 250, 155))
                        })
                        .collect();
                    let nb_line: Vec<[f64; 2]> = off
                        .histogram
                        .iter()
                        .map(|b| [f64::from(b.offspring_count), b.fitted_nb_prob * 100.0])
                        .collect();
                    let pois_line: Vec<[f64; 2]> = off
                        .histogram
                        .iter()
                        .map(|b| [f64::from(b.offspring_count), b.fitted_poisson_prob * 100.0])
                        .collect();

                    Plot::new("offspring_pmf_plot")
                        .height(235.0)
                        .legend(Legend::default())
                        .x_axis_label("Secondary Offspring Count z")
                        .y_axis_label("Probability (%)")
                        .show(ui, |plot_ui| {
                            plot_ui.bar_chart(BarChart::new(bars).name("Observed %"));
                            plot_ui.line(
                                Line::new(PlotPoints::from(nb_line))
                                    .name("Fitted Negative Binomial")
                                    .color(Color32::from_rgb(52, 211, 153))
                                    .width(2.4),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(pois_line))
                                    .name("Homogeneous Poisson")
                                    .color(Color32::from_rgb(251, 191, 36))
                                    .width(1.8),
                            );
                        });
                }
                RightPlotKind::ProfileLogLikelihood => {
                    ui.label(
                        RichText::new(
                            "Profile Log-Likelihood Relative to MLE (95% CI Cutoff at -1.92)",
                        )
                        .strong(),
                    );
                    let ll_pts: Vec<[f64; 2]> = off
                        .profile_ll_curve
                        .iter()
                        .map(|p| [p[0].log10(), p[1]])
                        .collect();
                    Plot::new("offspring_profile_ll_plot")
                        .height(235.0)
                        .legend(Legend::default())
                        .x_axis_label("log10(k)")
                        .y_axis_label("ell(k) - ell(k_MLE)")
                        .show(ui, |plot_ui| {
                            plot_ui.line(
                                Line::new(PlotPoints::from(ll_pts))
                                    .name("Relative Profile Log-Likelihood")
                                    .color(Color32::from_rgb(167, 139, 250))
                                    .width(2.4),
                            );
                            plot_ui.hline(
                                HLine::new(-1.920729)
                                    .name("95% Chi-Square Cutoff (-1.92)")
                                    .color(Color32::from_rgb(244, 63, 94)),
                            );
                        });
                }
            });
        });

        ui.add_space(6.0);

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.label(
                    RichText::new("Exposure Venue Amplification & Bipartite Centrality Ranking")
                        .strong(),
                );
                egui::ScrollArea::vertical()
                    .id_source("venue_amplification_scroll")
                    .max_height(215.0)
                    .show(ui, |ui| {
                        egui::Grid::new("venue_amplification_grid")
                            .striped(true)
                            .min_col_width(54.0)
                            .show(ui, |ui| {
                                ui.strong("Venue ID");
                                ui.strong("Setting");
                                ui.strong("Sec Tx");
                                ui.strong("Share");
                                ui.strong("SAR");
                                ui.strong("k-Core");
                                ui.strong("Betweenness");
                                ui.end_row();

                                for v in &self.analysis.venue_metrics {
                                    ui.monospace(&v.venue_id);
                                    ui.label(&v.venue_category);
                                    ui.monospace(format!("{}", v.secondary_transmissions));
                                    ui.monospace(format!(
                                        "{:.1}%",
                                        v.share_of_outbreak_transmissions * 100.0
                                    ));
                                    ui.monospace(
                                        v.secondary_attack_rate
                                            .map(|s| format!("{:.1}%", s * 100.0))
                                            .unwrap_or_else(|| "--".to_string()),
                                    );
                                    ui.monospace(format!("{}", v.bipartite_k_core));
                                    ui.monospace(format!("{:.4}", v.bipartite_betweenness));
                                    ui.end_row();
                                }
                            });
                    });
            });

            cols[1].group(|ui| {
                ui.label(
                    RichText::new(
                        "Case-Level Offspring, Generation Depth & Bridge Betweenness (Top 25)",
                    )
                    .strong(),
                );
                egui::ScrollArea::vertical()
                    .id_source("case_network_scroll")
                    .max_height(215.0)
                    .show(ui, |ui| {
                        egui::Grid::new("case_network_grid")
                            .striped(true)
                            .min_col_width(52.0)
                            .show(ui, |ui| {
                                ui.strong("Case ID");
                                ui.strong("Gen");
                                ui.strong("z_i");
                                ui.strong("Primary Venue");
                                ui.strong("Deg");
                                ui.strong("k-Core");
                                ui.strong("Betweenness");
                                ui.strong("Flag");
                                ui.end_row();

                                for c in self.analysis.case_metrics.iter().take(25) {
                                    ui.monospace(&c.case_id);
                                    ui.monospace(format!("G{}", c.generation_depth));
                                    ui.monospace(format!("{}", c.secondary_offspring));
                                    ui.monospace(&c.primary_venue_id);
                                    ui.monospace(format!("{}", c.bipartite_degree));
                                    ui.monospace(format!("{}", c.bipartite_k_core));
                                    ui.monospace(format!("{:.4}", c.bipartite_betweenness));
                                    let role = if c.is_superspreader_99 {
                                        RichText::new("SSE >99%").color(Color32::from_rgb(244, 63, 94))
                                    } else if c.bipartite_degree >= 2 {
                                        RichText::new("BRIDGE").color(Color32::from_rgb(251, 191, 36))
                                    } else if c.is_primary_index {
                                        RichText::new("INDEX").color(Color32::from_rgb(56, 189, 248))
                                    } else {
                                        RichText::new("LEAF").weak()
                                    };
                                    ui.label(role);
                                    ui.end_row();
                                }
                            });
                    });
            });
        });

        action
    }
}
