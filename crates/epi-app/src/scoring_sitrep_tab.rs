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

use crate::sitrep_export::{render_standalone_sitrep_html, SitRepContext};
use eframe::egui::{self, Color32, RichText};
use egui_plot::{Bar, BarChart, Legend, Line, Plot, PlotPoints};
use epi_analytics::{
    run_rolling_wis_backtest_and_ensemble, AberrationSeries, BacktestConfig, CoriRtConfig,
    CoriRtSeries, EnsembleBacktestResult, MultiZoneCluster, NowcastSeries, SeirConfig,
};
use epi_core::{SuppressionPolicy, TimeSeriesDataset};

pub struct ScoringSitrepWorkbenchState {
    pub cfg: BacktestConfig,
    pub sitrep_out_path: String,
    pub last_vintage_id: String,
    pub cached_backtest: Option<EnsembleBacktestResult>,
}

impl ScoringSitrepWorkbenchState {
    pub fn new(
        dataset: &TimeSeriesDataset,
        rt_cfg: &CoriRtConfig,
        seir_cfg: &SeirConfig,
    ) -> Self {
        let cfg = BacktestConfig {
            horizon_weeks: 4,
            num_origins: 14,
            si_mean_days: rt_cfg.si_mean_days,
            si_sd_days: rt_cfg.si_sd_days,
        };
        let counts: Vec<f64> = dataset.points.iter().map(|p| p.count).collect();
        let cached_backtest = run_rolling_wis_backtest_and_ensemble(&counts, &cfg, seir_cfg).ok();

        Self {
            cfg,
            sitrep_out_path: "epiledger_sitrep.html".to_string(),
            last_vintage_id: dataset.vintage_id.clone(),
            cached_backtest,
        }
    }

    pub fn sync_and_recompute(
        &mut self,
        dataset: &TimeSeriesDataset,
        rt_cfg: &CoriRtConfig,
        seir_cfg: &SeirConfig,
    ) -> Result<(), String> {
        self.cfg.si_mean_days = rt_cfg.si_mean_days;
        self.cfg.si_sd_days = rt_cfg.si_sd_days;
        self.last_vintage_id = dataset.vintage_id.clone();
        let counts: Vec<f64> = dataset.points.iter().map(|p| p.count).collect();
        match run_rolling_wis_backtest_and_ensemble(&counts, &self.cfg, seir_cfg) {
            Ok(res) => {
                self.cached_backtest = Some(res);
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn render(
        &mut self,
        ui: &mut egui::Ui,
        dataset: &TimeSeriesDataset,
        suppression_policy: SuppressionPolicy,
        rt_cfg: &CoriRtConfig,
        rt_series: &CoriRtSeries,
        ab_series: &AberrationSeries,
        nowcast: &NowcastSeries,
        seir_cfg: &SeirConfig,
        spatial_clusters: &[MultiZoneCluster],
        status_banner: &mut String,
    ) {
        if self.last_vintage_id != dataset.vintage_id {
            let _ = self.sync_and_recompute(dataset, rt_cfg, seir_cfg);
        }

        ui.heading(
            "Out-of-Sample Weighted Interval Score (K = 11 WIS), Inverse-WIS Ensemble Stacking & SitRep Export",
        );
        ui.label(
            RichText::new(
                "Bracher et al. 23-quantile WIS decomposition (Sharpness + Overprediction + Underprediction), Bosse log(y+1) WIS, causal Inverse-WIS stacking, and zero-JS HTML/SVG Situation Reports.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        if let Some(ref bt) = self.cached_backtest {
            let inv_sc = bt.scorecards.iter().find(|s| s.model_id == "inv_wis_ens");
            ui.columns(4, |cols| {
                cols[0].group(|ui| {
                    ui.label(RichText::new("TOP OUT-OF-SAMPLE MODEL").small().weak());
                    ui.heading(
                        RichText::new(&bt.best_model_label)
                            .monospace()
                            .color(Color32::from_rgb(52, 211, 153)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "{} origins x {}w horizon",
                            bt.origins_evaluated, bt.horizon_weeks
                        ))
                        .small()
                        .monospace(),
                    );
                });

                cols[1].group(|ui| {
                    ui.label(RichText::new("INVERSE-WIS RELATIVE SKILL").small().weak());
                    if let Some(sc) = inv_sc {
                        let col = if sc.relative_wis_vs_naive < 1.0 {
                            Color32::from_rgb(56, 189, 248)
                        } else {
                            Color32::from_rgb(251, 191, 36)
                        };
                        ui.heading(
                            RichText::new(format!("{:.2}x vs Naive", sc.relative_wis_vs_naive))
                                .monospace()
                                .color(col),
                        );
                        ui.label(
                            RichText::new(format!(
                                "Mean WIS {:.2} | Log-WIS {:.3}",
                                sc.mean_wis, sc.mean_log_wis
                            ))
                            .small()
                            .monospace(),
                        );
                    }
                });

                cols[2].group(|ui| {
                    ui.label(RichText::new("ENSEMBLE PI CALIBRATION").small().weak());
                    if let Some(sc) = inv_sc {
                        ui.heading(
                            RichText::new(format!("95% Cov: {:.1}%", sc.coverage_95_pct))
                                .monospace()
                                .color(Color32::from_rgb(251, 191, 36)),
                        );
                        ui.label(
                            RichText::new(format!("50% PI Coverage: {:.1}%", sc.coverage_50_pct))
                                .small()
                                .monospace(),
                        );
                    }
                });

                cols[3].group(|ui| {
                    ui.label(RichText::new("WK +1 STACKED FORECAST").small().weak());
                    if let Some(st1) = bt.forward_steps.first() {
                        let ens = &st1.inverse_wis_ensemble;
                        ui.heading(
                            RichText::new(format!("{:.1}", ens.median()))
                                .monospace()
                                .color(Color32::from_rgb(52, 211, 153)),
                        );
                        ui.label(
                            RichText::new(format!(
                                "95% PI [{:.0}, {:.0}]",
                                ens.pi95_low(),
                                ens.pi95_high()
                            ))
                            .small()
                            .monospace(),
                        );
                    }
                });
            });
        }

        ui.add_space(6.0);
        let mut need_recompute = false;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                if ui
                    .add(egui::Slider::new(&mut self.cfg.horizon_weeks, 1..=10).text("Horizon Wks"))
                    .changed()
                {
                    need_recompute = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.cfg.num_origins, 4..=28)
                            .text("Backtest Origins"),
                    )
                    .changed()
                {
                    need_recompute = true;
                }
                if ui.button("Re-Run Rolling WIS Backtest").clicked() {
                    need_recompute = true;
                }

                ui.separator();
                ui.label("SitRep HTML:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.sitrep_out_path)
                        .desired_width(185.0)
                        .hint_text("epiledger_sitrep.html"),
                );
                if ui
                    .button("Export Standalone HTML/SVG SitRep")
                    .clicked()
                {
                    if let Some(ref bt) = self.cached_backtest {
                        let ctx = SitRepContext {
                            dataset,
                            suppression_policy,
                            rt_series,
                            ab_series,
                            nowcast,
                            backtest: bt,
                            spatial_clusters,
                        };
                        let html = render_standalone_sitrep_html(&ctx);
                        let p = self.sitrep_out_path.trim();
                        match std::fs::write(p, html) {
                            Ok(()) => {
                                *status_banner = format!(
                                    "Exported standalone HTML5 + inline vector SVG Situation Report to '{p}'."
                                );
                            }
                            Err(e) => {
                                *status_banner = format!("SitRep export error: {e}");
                            }
                        }
                    }
                }
            });
        });

        if need_recompute {
            if let Err(err) = self.sync_and_recompute(dataset, rt_cfg, seir_cfg) {
                *status_banner = format!("WIS backtest warning: {err}");
            }
        }

        if let Some(ref bt) = self.cached_backtest {
            ui.columns(2, |cols| {
                cols[0].group(|ui| {
                    ui.strong(
                        "WIS Decomposition by Model (Sharpness + Overpred + Underpred)",
                    );
                    Plot::new("wis_decomp_plot")
                        .height(230.0)
                        .legend(Legend::default())
                        .show(ui, |plot_ui| {
                            let mut sharp_bars = Vec::new();
                            let mut over_bars = Vec::new();
                            let mut under_bars = Vec::new();
                            for (i, sc) in bt.scorecards.iter().enumerate() {
                                let x = i as f64;
                                sharp_bars.push(
                                    Bar::new(x - 0.22, sc.mean_sharpness)
                                        .width(0.20)
                                        .fill(Color32::from_rgb(56, 189, 248)),
                                );
                                over_bars.push(
                                    Bar::new(x, sc.mean_overprediction)
                                        .width(0.20)
                                        .fill(Color32::from_rgb(251, 191, 36)),
                                );
                                under_bars.push(
                                    Bar::new(x + 0.22, sc.mean_underprediction)
                                        .width(0.20)
                                        .fill(Color32::from_rgb(244, 63, 94)),
                                );
                            }
                            plot_ui.bar_chart(
                                BarChart::new(sharp_bars).name("Sharpness (Interval Width)"),
                            );
                            plot_ui.bar_chart(
                                BarChart::new(over_bars).name("Overprediction Penalty"),
                            );
                            plot_ui.bar_chart(
                                BarChart::new(under_bars).name("Underprediction Penalty"),
                            );
                        });
                    ui.label(
                        RichText::new(
                            "[0] Naive  [1] SEIR  [2] Renewal  [3] Trend  [4] Median Ens  [5] Inverse-WIS Ens",
                        )
                        .small()
                        .monospace()
                        .weak(),
                    );
                });

                cols[1].group(|ui| {
                    ui.strong("Live Forward Stacked Forecast (Inverse-WIS Pool + Constituents)");
                    Plot::new("ensemble_forward_plot")
                        .height(230.0)
                        .legend(Legend::default())
                        .show(ui, |plot_ui| {
                            let n_pts = dataset.points.len();
                            let tail = n_pts.min(18);
                            let start = n_pts - tail;

                            let obs_line: Vec<[f64; 2]> = (0..tail)
                                .map(|k| [-(tail as f64 - 1.0) + (k as f64), dataset.points[start + k].count])
                                .collect();
                            let last_val = dataset.points.last().map(|p| p.count).unwrap_or(0.0);

                            let mut ens_med = vec![[0.0, last_val]];
                            let mut ens_lo95 = vec![[0.0, last_val]];
                            let mut ens_hi95 = vec![[0.0, last_val]];
                            let mut seir_m = vec![[0.0, last_val]];
                            let mut ren_m = vec![[0.0, last_val]];

                            for st in &bt.forward_steps {
                                let x = st.horizon_week as f64;
                                ens_med.push([x, st.inverse_wis_ensemble.median()]);
                                ens_lo95.push([x, st.inverse_wis_ensemble.pi95_low()]);
                                ens_hi95.push([x, st.inverse_wis_ensemble.pi95_high()]);
                                seir_m.push([x, st.seir_median]);
                                ren_m.push([x, st.renewal_median]);
                            }

                            plot_ui.line(
                                Line::new(PlotPoints::from(obs_line))
                                    .name("Observed Incidence")
                                    .color(Color32::from_rgb(148, 163, 184))
                                    .width(1.8),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(ens_med))
                                    .name("Inverse-WIS Median")
                                    .color(Color32::from_rgb(52, 211, 153))
                                    .width(2.5),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(ens_lo95))
                                    .name("Ensemble 95% Low")
                                    .color(Color32::from_rgb(56, 189, 248))
                                    .width(1.0),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(ens_hi95))
                                    .name("Ensemble 95% High")
                                    .color(Color32::from_rgb(56, 189, 248))
                                    .width(1.0),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(seir_m))
                                    .name("SEIR Constituent")
                                    .color(Color32::from_rgb(251, 191, 36))
                                    .width(1.3),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(ren_m))
                                    .name("Renewal Constituent")
                                    .color(Color32::from_rgb(244, 114, 182))
                                    .width(1.3),
                            );
                        });
                });
            });

            ui.group(|ui| {
                ui.strong(format!(
                    "Out-of-Sample Proper Scoring Rule Leaderboard (Trained Weights: SEIR {:.1}%, Renewal {:.1}%, Trend {:.1}%)",
                    bt.trained_constituent_weights[0] * 100.0,
                    bt.trained_constituent_weights[1] * 100.0,
                    bt.trained_constituent_weights[2] * 100.0
                ));
                egui::Grid::new("wis_leaderboard_grid")
                    .striped(true)
                    .min_col_width(72.0)
                    .show(ui, |ui| {
                        ui.strong("Model");
                        ui.strong("Mean WIS");
                        ui.strong("Sharpness");
                        ui.strong("Overpred");
                        ui.strong("Underpred");
                        ui.strong("Log-WIS");
                        ui.strong("MAE");
                        ui.strong("Rel. WIS");
                        ui.strong("50% Cov");
                        ui.strong("95% Cov");
                        ui.strong("Stack Weight");
                        ui.end_row();

                        for sc in &bt.scorecards {
                            ui.monospace(&sc.model_label);
                            ui.monospace(format!("{:.2}", sc.mean_wis));
                            ui.monospace(format!("{:.2}", sc.mean_sharpness));
                            ui.monospace(format!("{:.2}", sc.mean_overprediction));
                            ui.monospace(format!("{:.2}", sc.mean_underprediction));
                            ui.monospace(format!("{:.3}", sc.mean_log_wis));
                            ui.monospace(format!("{:.2}", sc.mean_mae));
                            ui.monospace(format!("{:.2}x", sc.relative_wis_vs_naive));
                            ui.monospace(format!("{:.1}%", sc.coverage_50_pct));
                            ui.monospace(format!("{:.1}%", sc.coverage_95_pct));
                            ui.monospace(if sc.final_ensemble_weight > 0.0 {
                                format!("{:.1}%", sc.final_ensemble_weight * 100.0)
                            } else {
                                "--".to_string()
                            });
                            ui.end_row();
                        }
                    });
            });
        }
    }
}
