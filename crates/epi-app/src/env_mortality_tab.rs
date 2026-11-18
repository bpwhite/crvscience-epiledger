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

use eframe::egui::{self, Color32, Pos2, Rect, RichText, Sense, Vec2};
use egui_plot::{Legend, Line, Plot, PlotPoints, Points};
use epi_analytics::{
    builtin_all_cause_weekly_mortality_dataset, builtin_wildfire_pm25_dlnm_dataset,
    compute_euromomo_excess_mortality, fit_dlnm_quasi_poisson, parse_dlnm_exposure_csv,
    parse_weekly_mortality_csv, DlnmConfig, DlnmDataset, DlnmFitResult, ExcessMortalityConfig,
    ExcessMortalityResult, ExposureBasisKind, LagBasisKind, MortalityBaselineMode,
    WeeklyMortalityDataset,
};

pub struct EnvMortalityWorkbenchState {
    // Gasparrini Cross-Basis DLNM state
    pub dlnm_dataset: DlnmDataset,
    pub dlnm_cfg: DlnmConfig,
    pub dlnm_csv_path: String,
    pub dlnm_unit_input: String,
    pub cached_dlnm: Option<DlnmFitResult>,

    // EuroMOMO / World Mortality Excess Mortality state
    pub mort_dataset: WeeklyMortalityDataset,
    pub mort_cfg: ExcessMortalityConfig,
    pub mort_csv_path: String,
    pub cached_mort: Option<ExcessMortalityResult>,
}

impl Default for EnvMortalityWorkbenchState {
    fn default() -> Self {
        Self::new()
    }
}

impl EnvMortalityWorkbenchState {
    pub fn new() -> Self {
        let dlnm_dataset = builtin_wildfire_pm25_dlnm_dataset();
        let dlnm_cfg = DlnmConfig::default();
        let cached_dlnm = fit_dlnm_quasi_poisson(&dlnm_dataset, &dlnm_cfg).ok();

        let mort_dataset = builtin_all_cause_weekly_mortality_dataset();
        let mort_cfg = ExcessMortalityConfig::default();
        let cached_mort = compute_euromomo_excess_mortality(&mort_dataset, &mort_cfg).ok();

        Self {
            dlnm_dataset,
            dlnm_cfg,
            dlnm_csv_path: String::new(),
            dlnm_unit_input: "ug/m3 PM2.5".to_string(),
            cached_dlnm,
            mort_dataset,
            mort_cfg,
            mort_csv_path: String::new(),
            cached_mort,
        }
    }

    pub fn recompute_dlnm(&mut self) -> Result<(), String> {
        match fit_dlnm_quasi_poisson(&self.dlnm_dataset, &self.dlnm_cfg) {
            Ok(res) => {
                self.cached_dlnm = Some(res);
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    }

    pub fn recompute_mortality(&mut self) -> Result<(), String> {
        match compute_euromomo_excess_mortality(&self.mort_dataset, &self.mort_cfg) {
            Ok(res) => {
                self.cached_mort = Some(res);
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    }

    pub fn render(&mut self, ui: &mut egui::Ui, status_banner: &mut String) {
        self.render_dlnm_section(ui, status_banner);
        ui.add_space(12.0);
        ui.separator();
        ui.add_space(8.0);
        self.render_mortality_section(ui, status_banner);
    }

    fn render_dlnm_section(&mut self, ui: &mut egui::Ui, status_banner: &mut String) {
        ui.heading("Gasparrini Cross-Basis Distributed Lag Non-Linear Model (DLNM)");
        ui.label(
            RichText::new(
                "Bi-dimensional exposure x lag tensor cross-basis estimated via Quasi-Poisson IRLS for wildfire PM2.5, extreme heat, and precipitation health lags.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        if let Some(ref res) = self.cached_dlnm {
            ui.columns(4, |cols| {
                cols[0].group(|ui| {
                    ui.label(RichText::new("P95 CUMULATIVE RELATIVE RISK").small().weak());
                    let col = if res.p95_cumulative_rr > 1.0 {
                        Color32::from_rgb(244, 63, 94)
                    } else {
                        Color32::from_rgb(52, 211, 153)
                    };
                    ui.heading(
                        RichText::new(format!("{:.3}x", res.p95_cumulative_rr))
                            .monospace()
                            .color(col),
                    );
                    ui.label(
                        RichText::new(format!(
                            "95% CI [{:.3}, {:.3}] @ {:.1}",
                            res.p95_rr_ci_low_95, res.p95_rr_ci_high_95, res.p95_exposure
                        ))
                        .small()
                        .monospace(),
                    );
                });

                cols[1].group(|ui| {
                    ui.label(RichText::new("MINIMUM RISK EXPOSURE (MME)").small().weak());
                    ui.heading(
                        RichText::new(format!(
                            "{:.1} {}",
                            res.minimum_risk_exposure, res.exposure_unit
                        ))
                        .monospace()
                        .color(Color32::from_rgb(56, 189, 248)),
                    );
                    ui.label(
                        RichText::new(format!("Centered ref: {:.1}", res.reference_exposure))
                            .small()
                            .monospace(),
                    );
                });

                cols[2].group(|ui| {
                    ui.label(RichText::new("ATTRIBUTABLE BURDEN (AN / AF)").small().weak());
                    ui.heading(
                        RichText::new(format!("{:.1} events", res.attributable_events))
                            .monospace()
                            .color(Color32::from_rgb(251, 191, 36)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Attributable fraction: {:.2}%",
                            res.attributable_fraction_pct
                        ))
                        .small()
                        .monospace(),
                    );
                });

                cols[3].group(|ui| {
                    ui.label(RichText::new("QUASI-POISSON FIT").small().weak());
                    ui.heading(
                        RichText::new(format!("phi = {:.2}", res.dispersion_phi)).monospace(),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Deviance {:.1} ({} steps)",
                            res.deviance, res.effective_rows_fitted
                        ))
                        .small()
                        .monospace(),
                    );
                });
            });
        }

        ui.add_space(6.0);
        let mut dlnm_dirty = false;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("Exposure Basis:");
                for kind in [
                    ExposureBasisKind::NaturalCubicSpline3Df,
                    ExposureBasisKind::QuadraticThreshold,
                ] {
                    if ui
                        .selectable_label(self.dlnm_cfg.exposure_basis == kind, kind.label())
                        .clicked()
                    {
                        self.dlnm_cfg.exposure_basis = kind;
                        dlnm_dirty = true;
                    }
                }
                ui.separator();
                ui.label("Lag Basis:");
                for lkind in [
                    LagBasisKind::NaturalCubicLogLag3Df,
                    LagBasisKind::QuadraticPolynomial,
                ] {
                    if ui
                        .selectable_label(self.dlnm_cfg.lag_basis == lkind, lkind.label())
                        .clicked()
                    {
                        self.dlnm_cfg.lag_basis = lkind;
                        dlnm_dirty = true;
                    }
                }
                ui.separator();
                if ui
                    .add(egui::Slider::new(&mut self.dlnm_cfg.max_lag, 3..=28).text("Max Lag"))
                    .changed()
                {
                    dlnm_dirty = true;
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("Exposure CSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.dlnm_csv_path)
                        .desired_width(220.0)
                        .hint_text("date,exposure,outcome[,confounder]"),
                );
                if ui.button("Browse...").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("CSV / TSV", &["csv", "tsv", "txt"])
                        .pick_file()
                    {
                        self.dlnm_csv_path = p.display().to_string();
                    }
                }
                ui.label("Unit:");
                ui.add(egui::TextEdit::singleline(&mut self.dlnm_unit_input).desired_width(100.0));
                if ui.button("Load Exposure CSV").clicked() {
                    let p = self.dlnm_csv_path.trim().to_string();
                    let u = self.dlnm_unit_input.trim().to_string();
                    match std::fs::read_to_string(&p)
                        .map_err(|e| e.to_string())
                        .and_then(|txt| parse_dlnm_exposure_csv(&txt, &p, &u).map_err(|e| e.to_string()))
                    {
                        Ok(ds) => {
                            let n = ds.observations.len();
                            self.dlnm_dataset = ds;
                            dlnm_dirty = true;
                            *status_banner = format!("Loaded DLNM dataset '{p}' ({n} rows).");
                        }
                        Err(e) => {
                            *status_banner = format!("DLNM CSV error: {e}");
                        }
                    }
                }
                if ui.button("Reset Wildfire PM2.5 Benchmark").clicked() {
                    self.dlnm_dataset = builtin_wildfire_pm25_dlnm_dataset();
                    dlnm_dirty = true;
                    *status_banner =
                        "Restored built-in 365-day Pacific NW Wildfire PM2.5 benchmark.".to_string();
                }
            });
        });

        if dlnm_dirty {
            if let Err(err) = self.recompute_dlnm() {
                *status_banner = format!("DLNM fit warning: {err}");
            }
        }

        if let Some(ref res) = self.cached_dlnm {
            ui.columns(2, |cols| {
                cols[0].group(|ui| {
                    ui.strong("Overall Cumulative Exposure-Response RR(x | x_ref)");
                    Plot::new("dlnm_cum_rr_plot")
                        .height(215.0)
                        .legend(Legend::default())
                        .show(ui, |plot_ui| {
                            let x_lo = res
                                .cumulative_curve
                                .first()
                                .map(|p| p.exposure)
                                .unwrap_or(0.0);
                            let x_hi = res
                                .cumulative_curve
                                .last()
                                .map(|p| p.exposure)
                                .unwrap_or(1.0);
                            plot_ui.line(
                                Line::new(PlotPoints::from(vec![[x_lo, 1.0], [x_hi, 1.0]]))
                                    .name("RR = 1.0")
                                    .color(Color32::from_gray(120)),
                            );
                            let mid: Vec<[f64; 2]> = res
                                .cumulative_curve
                                .iter()
                                .map(|p| [p.exposure, p.cumulative_rr])
                                .collect();
                            let lo: Vec<[f64; 2]> = res
                                .cumulative_curve
                                .iter()
                                .map(|p| [p.exposure, p.rr_ci_low_95])
                                .collect();
                            let hi: Vec<[f64; 2]> = res
                                .cumulative_curve
                                .iter()
                                .map(|p| [p.exposure, p.rr_ci_high_95])
                                .collect();
                            plot_ui.line(
                                Line::new(PlotPoints::from(mid))
                                    .name("Cumulative RR")
                                    .color(Color32::from_rgb(244, 63, 94))
                                    .width(2.2),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(lo))
                                    .name("95% CI Low")
                                    .color(Color32::from_rgb(148, 163, 184))
                                    .width(1.0),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(hi))
                                    .name("95% CI High")
                                    .color(Color32::from_rgb(148, 163, 184))
                                    .width(1.0),
                            );
                        });
                });

                cols[1].group(|ui| {
                    ui.strong(format!(
                        "P95 Lag-Response Slice RR(x={:.1}, lag=0..{})",
                        res.p95_exposure, res.max_lag
                    ));
                    Plot::new("dlnm_lag_slice_plot")
                        .height(215.0)
                        .legend(Legend::default())
                        .show(ui, |plot_ui| {
                            plot_ui.line(
                                Line::new(PlotPoints::from(vec![
                                    [0.0, 1.0],
                                    [res.max_lag as f64, 1.0],
                                ]))
                                .name("RR = 1.0 (Harvesting < 1.0)")
                                .color(Color32::from_gray(120)),
                            );
                            let mid: Vec<[f64; 2]> = res
                                .p95_lag_slice
                                .iter()
                                .map(|p| [p.lag as f64, p.rr])
                                .collect();
                            let lo: Vec<[f64; 2]> = res
                                .p95_lag_slice
                                .iter()
                                .map(|p| [p.lag as f64, p.rr_ci_low_95])
                                .collect();
                            let hi: Vec<[f64; 2]> = res
                                .p95_lag_slice
                                .iter()
                                .map(|p| [p.lag as f64, p.rr_ci_high_95])
                                .collect();
                            plot_ui.line(
                                Line::new(PlotPoints::from(mid))
                                    .name("Lag-Specific RR")
                                    .color(Color32::from_rgb(56, 189, 248))
                                    .width(2.2),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(lo))
                                    .name("95% CI Low")
                                    .color(Color32::from_rgb(148, 163, 184))
                                    .width(1.0),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(hi))
                                    .name("95% CI High")
                                    .color(Color32::from_rgb(148, 163, 184))
                                    .width(1.0),
                            );
                        });
                });
            });

            self.draw_dlnm_surface_heatmap(ui, res);
        }
    }

    fn draw_dlnm_surface_heatmap(&self, ui: &mut egui::Ui, res: &DlnmFitResult) {
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.strong("Bi-Dimensional Cross-Basis RR(exposure, lag) Surface Heatmap");
                ui.label(
                    RichText::new(
                        "X = Lag (0..max_lag) | Y = Exposure Low->High | Blue < 1.0 (Displacement) | Slate = 1.0 | Crimson > 1.0 (Acute Excess)",
                    )
                    .small()
                    .weak(),
                );
            });

            let n_exp = res.surface_rr_grid.len();
            let n_lag = res.max_lag + 1;
            if n_exp == 0 || n_lag == 0 {
                return;
            }

            let desired = Vec2::new(ui.available_width().max(480.0), 145.0);
            let (rect, _resp) = ui.allocate_exact_size(desired, Sense::hover());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 4.0, Color32::from_rgb(15, 23, 42));

            let pad_l = 52.0_f32;
            let pad_r = 16.0_f32;
            let pad_t = 10.0_f32;
            let pad_b = 20.0_f32;
            let grid_r = Rect::from_min_max(
                Pos2::new(rect.left() + pad_l, rect.top() + pad_t),
                Pos2::new(rect.right() - pad_r, rect.bottom() - pad_b),
            );

            let cw = grid_r.width() / (n_lag as f32);
            let ch = grid_r.height() / (n_exp as f32);

            let mut max_dev = 0.02_f64;
            for row in &res.surface_rr_grid {
                for &rr in row {
                    max_dev = max_dev.max(rr.ln().abs());
                }
            }

            for (ei, row) in res.surface_rr_grid.iter().enumerate() {
                // Draw high exposures at top, low exposures at bottom
                let inv_row = n_exp - 1 - ei;
                let y0 = grid_r.top() + (inv_row as f32) * ch;
                for (li, &rr) in row.iter().enumerate() {
                    let x0 = grid_r.left() + (li as f32) * cw;
                    let cell = Rect::from_min_max(Pos2::new(x0, y0), Pos2::new(x0 + cw, y0 + ch));
                    let norm = (rr.ln() / max_dev).clamp(-1.0, 1.0) as f32;
                    let col = if norm >= 0.0 {
                        let r = (30.0 + 214.0 * norm) as u8;
                        let g = (41.0 + 22.0 * (1.0 - norm)) as u8;
                        let b = (59.0 + 35.0 * (1.0 - norm)) as u8;
                        Color32::from_rgb(r, g, b)
                    } else {
                        let a = -norm;
                        let r = (30.0 * (1.0 - a) + 14.0 * a) as u8;
                        let g = (41.0 * (1.0 - a) + 165.0 * a) as u8;
                        let b = (59.0 * (1.0 - a) + 233.0 * a) as u8;
                        Color32::from_rgb(r, g, b)
                    };
                    painter.rect_filled(cell, 0.0, col);
                }
            }

            painter.text(
                Pos2::new(rect.left() + 6.0, grid_r.top() + 6.0),
                egui::Align2::LEFT_TOP,
                format!("{:.0}", res.surface_exposures.last().copied().unwrap_or(0.0)),
                egui::FontId::monospace(10.0),
                Color32::from_gray(200),
            );
            painter.text(
                Pos2::new(rect.left() + 6.0, grid_r.bottom() - 4.0),
                egui::Align2::LEFT_BOTTOM,
                format!("{:.0}", res.surface_exposures.first().copied().unwrap_or(0.0)),
                egui::FontId::monospace(10.0),
                Color32::from_gray(200),
            );
            painter.text(
                Pos2::new(grid_r.left(), rect.bottom() - 3.0),
                egui::Align2::LEFT_BOTTOM,
                "Lag 0",
                egui::FontId::monospace(10.0),
                Color32::from_gray(200),
            );
            painter.text(
                Pos2::new(grid_r.right(), rect.bottom() - 3.0),
                egui::Align2::RIGHT_BOTTOM,
                format!("Lag {}", res.max_lag),
                egui::FontId::monospace(10.0),
                Color32::from_gray(200),
            );
        });
    }

    fn render_mortality_section(&mut self, ui: &mut egui::Ui, status_banner: &mut String) {
        ui.heading(
            "EuroMOMO / Serfling Harmonic Excess Mortality, 2/3-Power Box-Cox Z-Scores & YPLL-75",
        );
        ui.label(
            RichText::new(
                "Harmonic Quasi-Poisson baseline with 2/3-power Box-Cox skewness stabilization (Z > 2 substantial excess, Z > 4 extreme excess) & World Mortality Dataset baseline toggle.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        if let Some(ref res) = self.cached_mort {
            ui.columns(4, |cols| {
                cols[0].group(|ui| {
                    ui.label(RichText::new("EVAL WINDOW CUMULATIVE EXCESS").small().weak());
                    ui.heading(
                        RichText::new(format!("{:+.0} deaths", res.evaluation_cumulative_excess))
                            .monospace()
                            .color(Color32::from_rgb(244, 63, 94)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Mean P-score: {:+.1}%",
                            res.evaluation_mean_p_score_pct
                        ))
                        .small()
                        .monospace(),
                    );
                });

                cols[1].group(|ui| {
                    ui.label(RichText::new("BOX-COX Z-SCORE ALERTS").small().weak());
                    ui.heading(
                        RichText::new(format!(
                            "{} wks (Z > {:.1})",
                            res.substantial_alert_weeks, self.mort_cfg.z_substantial_threshold
                        ))
                        .monospace()
                        .color(Color32::from_rgb(251, 191, 36)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Extreme (Z > {:.1}): {} wks",
                            self.mort_cfg.z_extreme_threshold, res.extreme_alert_weeks
                        ))
                        .small()
                        .monospace(),
                    );
                });

                cols[2].group(|ui| {
                    ui.label(RichText::new("EXCESS YPLL-75").small().weak());
                    ui.heading(
                        RichText::new(format!("{:.0} life-yrs", res.evaluation_ypll_75))
                            .monospace()
                            .color(Color32::from_rgb(56, 189, 248)),
                    );
                    ui.label(RichText::new("Premature years lost <75y").small().monospace());
                });

                cols[3].group(|ui| {
                    ui.label(RichText::new("PEAK SURVEILLANCE Z-SCORE").small().weak());
                    ui.heading(
                        RichText::new(format!("Z = {:+.2}", res.peak_z_score))
                            .monospace()
                            .color(Color32::from_rgb(248, 113, 113)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Peak {} | phi={:.2}",
                            res.peak_period_label, res.dispersion_phi
                        ))
                        .small()
                        .monospace(),
                    );
                });
            });
        }

        ui.add_space(6.0);
        let mut mort_dirty = false;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("Baseline Regime:");
                for mode in [
                    MortalityBaselineMode::EuroMomoSpringAutumn,
                    MortalityBaselineMode::WorldMortalityFullYear,
                ] {
                    if ui
                        .selectable_label(self.mort_cfg.baseline_mode == mode, mode.label())
                        .clicked()
                    {
                        self.mort_cfg.baseline_mode = mode;
                        mort_dirty = true;
                    }
                }
                ui.separator();
                if ui
                    .add(
                        egui::Slider::new(&mut self.mort_cfg.evaluation_window_weeks, 12..=104)
                            .text("Eval Wks"),
                    )
                    .changed()
                {
                    mort_dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.mort_cfg.seasonal_harmonics, 1..=3)
                            .text("Harmonics"),
                    )
                    .changed()
                {
                    mort_dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.mort_cfg.z_substantial_threshold, 1.5..=3.5)
                            .text("Z Threshold"),
                    )
                    .changed()
                {
                    mort_dirty = true;
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("Weekly Mortality CSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.mort_csv_path)
                        .desired_width(260.0)
                        .hint_text("week,iso_week,deaths[,population,mean_age]"),
                );
                if ui.button("Browse...").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("CSV / TSV", &["csv", "tsv", "txt"])
                        .pick_file()
                    {
                        self.mort_csv_path = p.display().to_string();
                    }
                }
                if ui.button("Load Mortality CSV").clicked() {
                    let p = self.mort_csv_path.trim().to_string();
                    match std::fs::read_to_string(&p)
                        .map_err(|e| e.to_string())
                        .and_then(|txt| parse_weekly_mortality_csv(&txt, &p).map_err(|e| e.to_string()))
                    {
                        Ok(ds) => {
                            let n = ds.rows.len();
                            self.mort_dataset = ds;
                            mort_dirty = true;
                            *status_banner =
                                format!("Loaded weekly mortality dataset '{p}' ({n} weeks).");
                        }
                        Err(e) => {
                            *status_banner = format!("Mortality CSV error: {e}");
                        }
                    }
                }
                if ui.button("Reset 208-Wk Mortality Benchmark").clicked() {
                    self.mort_dataset = builtin_all_cause_weekly_mortality_dataset();
                    mort_dirty = true;
                    *status_banner =
                        "Restored built-in 208-week regional all-cause mortality benchmark."
                            .to_string();
                }
            });
        });

        if mort_dirty {
            if let Err(err) = self.recompute_mortality() {
                *status_banner = format!("Excess mortality warning: {err}");
            }
        }

        if let Some(ref res) = self.cached_mort {
            ui.columns(2, |cols| {
                cols[0].group(|ui| {
                    ui.strong("Weekly All-Cause Deaths vs. Harmonic Baseline & +2Z / +4Z Bands");
                    Plot::new("euromomo_counts_plot")
                        .height(225.0)
                        .legend(Legend::default())
                        .show(ui, |plot_ui| {
                            let obs: Vec<[f64; 2]> = res
                                .timeline
                                .iter()
                                .map(|p| [p.step_index as f64, p.observed_deaths])
                                .collect();
                            let base: Vec<[f64; 2]> = res
                                .timeline
                                .iter()
                                .map(|p| [p.step_index as f64, p.expected_baseline])
                                .collect();
                            let u2: Vec<[f64; 2]> = res
                                .timeline
                                .iter()
                                .map(|p| [p.step_index as f64, p.substantial_upper_2z])
                                .collect();
                            let u4: Vec<[f64; 2]> = res
                                .timeline
                                .iter()
                                .map(|p| [p.step_index as f64, p.extreme_upper_4z])
                                .collect();
                            let alerts: Vec<[f64; 2]> = res
                                .timeline
                                .iter()
                                .filter(|p| p.is_substantial_excess)
                                .map(|p| [p.step_index as f64, p.observed_deaths])
                                .collect();

                            plot_ui.line(
                                Line::new(PlotPoints::from(base))
                                    .name("Expected Baseline")
                                    .color(Color32::from_rgb(52, 211, 153))
                                    .width(2.0),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(u2))
                                    .name("Substantial (+2Z)")
                                    .color(Color32::from_rgb(251, 191, 36))
                                    .width(1.2),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(u4))
                                    .name("Extreme (+4Z)")
                                    .color(Color32::from_rgb(244, 63, 94))
                                    .width(1.2),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(obs))
                                    .name("Observed Deaths")
                                    .color(Color32::from_rgb(148, 163, 184))
                                    .width(1.5),
                            );
                            if !alerts.is_empty() {
                                plot_ui.points(
                                    Points::new(PlotPoints::from(alerts))
                                        .name("Z > Threshold")
                                        .color(Color32::from_rgb(251, 113, 133))
                                        .radius(3.8),
                                );
                            }
                        });
                });

                cols[1].group(|ui| {
                    ui.strong("EuroMOMO 2/3-Power Box-Cox Standardized Z-Score Series");
                    Plot::new("euromomo_z_plot")
                        .height(225.0)
                        .legend(Legend::default())
                        .show(ui, |plot_ui| {
                            let max_x = res.timeline.len().saturating_sub(1) as f64;
                            plot_ui.line(
                                Line::new(PlotPoints::from(vec![
                                    [0.0, self.mort_cfg.z_substantial_threshold],
                                    [max_x, self.mort_cfg.z_substantial_threshold],
                                ]))
                                .name("Z = +2 Threshold")
                                .color(Color32::from_rgb(251, 191, 36)),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(vec![
                                    [0.0, self.mort_cfg.z_extreme_threshold],
                                    [max_x, self.mort_cfg.z_extreme_threshold],
                                ]))
                                .name("Z = +4 Extreme")
                                .color(Color32::from_rgb(244, 63, 94)),
                            );
                            let zs: Vec<[f64; 2]> = res
                                .timeline
                                .iter()
                                .map(|p| [p.step_index as f64, p.box_cox_z_score])
                                .collect();
                            plot_ui.line(
                                Line::new(PlotPoints::from(zs))
                                    .name("Box-Cox Z-Score")
                                    .color(Color32::from_rgb(56, 189, 248))
                                    .width(1.8),
                            );
                        });
                });
            });

            ui.group(|ui| {
                ui.strong("Recent Surveillance Window Weeks (Trailing 16 Weeks)");
                egui::ScrollArea::vertical().max_height(185.0).show(ui, |ui| {
                    egui::Grid::new("mort_tail_grid").striped(true).show(ui, |ui| {
                        ui.strong("Period");
                        ui.strong("Observed");
                        ui.strong("Baseline");
                        ui.strong("Excess");
                        ui.strong("P-Score");
                        ui.strong("Box-Cox Z");
                        ui.strong("YPLL-75");
                        ui.strong("Signal");
                        ui.end_row();

                        let tail_start = res.timeline.len().saturating_sub(16);
                        for p in &res.timeline[tail_start..] {
                            ui.monospace(&p.period_label);
                            ui.monospace(format!("{:.0}", p.observed_deaths));
                            ui.monospace(format!("{:.1}", p.expected_baseline));
                            ui.monospace(format!("{:+.1}", p.excess_deaths));
                            ui.monospace(format!("{:+.1}%", p.p_score_pct));
                            ui.monospace(format!("{:+.2}", p.box_cox_z_score));
                            ui.monospace(format!("{:.0}", p.excess_ypll_75));
                            let status = if p.is_extreme_excess {
                                RichText::new("EXTREME (>4Z)").color(Color32::from_rgb(244, 63, 94))
                            } else if p.is_substantial_excess {
                                RichText::new("EXCESS (>2Z)").color(Color32::from_rgb(251, 191, 36))
                            } else {
                                RichText::new("Baseline").color(Color32::from_rgb(52, 211, 153))
                            };
                            ui.label(status.monospace());
                            ui.end_row();
                        }
                    });
                });
            });
        }
    }
}
