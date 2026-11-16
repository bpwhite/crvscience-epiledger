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

use eframe::egui::{self, Color32, RichText, Stroke};
use egui_plot::{HLine, Legend, Line, Plot, PlotPoints};
use epi_analytics::{
    builtin_emerging_wave_cases_and_deaths, compute_delay_adjusted_cfr,
    parse_cfr_incidence_deaths_csv, project_competing_risk_hospital_surge,
    project_hospital_surge_from_seir, CaseDeathObservation, DelayAdjustedCfrResult,
    DelayCfrConfig, HospitalSurgeConfig, HospitalSurgeResult, SeirForecastResult,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase12SubView {
    RightTruncatedCfr,
    CompetingRiskSurge,
}

pub struct SeveritySurgeWorkbenchState {
    pub sub_view: Phase12SubView,

    // Right-Truncated Delay-Adjusted cCFR state
    pub cfr_series: Vec<CaseDeathObservation>,
    pub cfr_trunc_cutoff_days: usize,
    pub cfr_cfg: DelayCfrConfig,
    pub cfr_result: DelayAdjustedCfrResult,
    pub cfr_csv_path: String,

    // Competing-Risk Hospital / ICU / Ventilator Surge state
    pub surge_cfg: HospitalSurgeConfig,
    pub drive_from_active_seir: bool,
    pub built_in_scale_multiplier: f64,
    pub surge_result: HospitalSurgeResult,
}

impl SeveritySurgeWorkbenchState {
    pub fn new(initial_seir: &SeirForecastResult) -> Self {
        let cfr_series = builtin_emerging_wave_cases_and_deaths();
        let cfr_trunc_cutoff_days = 26;
        let cfr_cfg = DelayCfrConfig::default();
        let cfr_result = compute_delay_adjusted_cfr(&cfr_series[..cfr_trunc_cutoff_days], &cfr_cfg)
            .expect("Built-in cCFR series is valid");

        let surge_cfg = HospitalSurgeConfig::default();
        let built_in_scale_multiplier = 6.0;
        let cases_d: Vec<f64> = cfr_series
            .iter()
            .map(|r| r.incident_cases * built_in_scale_multiplier)
            .collect();
        let low_d: Vec<f64> = cases_d.iter().map(|&v| v * 0.78).collect();
        let high_d: Vec<f64> = cases_d.iter().map(|&v| v * 1.24).collect();
        let surge_result =
            project_competing_risk_hospital_surge(&cases_d, &low_d, &high_d, &surge_cfg)
                .or_else(|_| project_hospital_surge_from_seir(initial_seir, &surge_cfg))
                .expect("Hospital surge initialization valid");

        Self {
            sub_view: Phase12SubView::RightTruncatedCfr,
            cfr_series,
            cfr_trunc_cutoff_days,
            cfr_cfg,
            cfr_result,
            cfr_csv_path: String::new(),
            surge_cfg,
            drive_from_active_seir: false,
            built_in_scale_multiplier,
            surge_result,
        }
    }

    pub fn recompute_cfr(&mut self) -> Result<(), String> {
        let end = self
            .cfr_trunc_cutoff_days
            .clamp(3, self.cfr_series.len().max(3))
            .min(self.cfr_series.len());
        self.cfr_result =
            compute_delay_adjusted_cfr(&self.cfr_series[..end], &self.cfr_cfg)
                .map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn recompute_surge(&mut self, active_seir: &SeirForecastResult) -> Result<(), String> {
        if self.drive_from_active_seir {
            self.surge_result = project_hospital_surge_from_seir(active_seir, &self.surge_cfg)
                .map_err(|e| e.to_string())?;
        } else {
            let cases_d: Vec<f64> = self
                .cfr_series
                .iter()
                .map(|r| r.incident_cases * self.built_in_scale_multiplier)
                .collect();
            let low_d: Vec<f64> = cases_d.iter().map(|&v| v * 0.78).collect();
            let high_d: Vec<f64> = cases_d.iter().map(|&v| v * 1.24).collect();
            self.surge_result =
                project_competing_risk_hospital_surge(&cases_d, &low_d, &high_d, &self.surge_cfg)
                    .map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    pub fn render(
        &mut self,
        ui: &mut egui::Ui,
        active_seir: &SeirForecastResult,
        status_banner: &mut String,
    ) {
        ui.horizontal_wrapped(|ui| {
            ui.heading("Right-Truncated Delay-Adjusted cCFR & Competing-Risk Hospital/ICU Surge");
            ui.separator();
            ui.selectable_value(
                &mut self.sub_view,
                Phase12SubView::RightTruncatedCfr,
                "Nishiura Right-Truncated cCFR",
            );
            ui.selectable_value(
                &mut self.sub_view,
                Phase12SubView::CompetingRiskSurge,
                "Competing-Risk Ward / ICU / Ventilator Surge",
            );
        });
        ui.add_space(4.0);

        match self.sub_view {
            Phase12SubView::RightTruncatedCfr => self.render_cfr_panel(ui, status_banner),
            Phase12SubView::CompetingRiskSurge => {
                self.render_surge_panel(ui, active_seir, status_banner)
            }
        }
    }

    fn render_cfr_panel(&mut self, ui: &mut egui::Ui, status_banner: &mut String) {
        let r = &self.cfr_result;

        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("DELAY-ADJUSTED cCFR [95% CI]").small().weak());
                ui.heading(
                    RichText::new(format!("{:.2}%", r.final_adjusted_ccfr_pct))
                        .monospace()
                        .color(Color32::from_rgb(56, 189, 248)),
                );
                ui.label(
                    RichText::new(format!(
                        "Wilson 95% [{:.2}%, {:.2}%]",
                        r.final_adjusted_ci_low_95, r.final_adjusted_ci_high_95
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("NAIVE CRUDE CFR (UNADJUSTED)").small().weak());
                ui.heading(
                    RichText::new(format!("{:.2}%", r.final_crude_cfr_pct))
                        .monospace()
                        .color(Color32::from_rgb(244, 63, 94)),
                );
                ui.label(
                    RichText::new(format!(
                        "{:.0} deaths / {:.0} raw cases",
                        r.final_cumulative_deaths, r.final_cumulative_cases
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[2].group(|ui| {
                ui.label(RichText::new("NISHIURA EFFECTIVE DENOMINATOR u_t").small().weak());
                ui.heading(
                    RichText::new(format!("{:.1} resolved", r.final_resolved_denominator))
                        .monospace()
                        .color(Color32::from_rgb(52, 211, 153)),
                );
                let resolved_pct = (r.final_resolved_denominator * 100.0)
                    / r.final_cumulative_cases.max(1.0);
                ui.label(
                    RichText::new(format!("{resolved_pct:.1}% outcome-resolved"))
                        .small()
                        .monospace(),
                );
            });

            cols[3].group(|ui| {
                ui.label(RichText::new("PEAK RIGHT-TRUNCATION BIAS").small().weak());
                ui.heading(
                    RichText::new(format!("{:.2}x downward bias", r.peak_truncation_bias_factor))
                        .monospace()
                        .color(Color32::from_rgb(251, 191, 36)),
                );
                ui.label(
                    RichText::new("Ratio cCFR / Crude CFR during growth")
                        .small()
                        .monospace(),
                );
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                let max_len = self.cfr_series.len().max(3);
                ui.label(RichText::new("Active-Wave Evaluation Horizon:").strong());
                let mut ch = ui
                    .add(
                        egui::Slider::new(&mut self.cfr_trunc_cutoff_days, 8..=max_len)
                            .suffix(" steps"),
                    )
                    .changed();
                ui.label("Onset-to-Death Mean:");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.cfr_cfg.delay_mean)
                            .speed(0.2)
                            .range(1.0..=45.0),
                    )
                    .changed();
                ui.label("SD:");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.cfr_cfg.delay_sd)
                            .speed(0.2)
                            .range(0.5..=30.0),
                    )
                    .changed();

                if ch {
                    if let Err(e) = self.recompute_cfr() {
                        *status_banner = format!("cCFR update warning: {e}");
                    }
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("Case-Death CSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.cfr_csv_path)
                        .desired_width(225.0)
                        .hint_text("path/to/cases_and_deaths.csv"),
                );
                if ui.button("Import CFR CSV").clicked() {
                    match std::fs::read_to_string(self.cfr_csv_path.trim()) {
                        Ok(raw) => match parse_cfr_incidence_deaths_csv(&raw) {
                            Ok(rows) => {
                                self.cfr_trunc_cutoff_days = rows.len();
                                self.cfr_series = rows;
                                let _ = self.recompute_cfr();
                                *status_banner = format!(
                                    "Computed Nishiura delay-adjusted cCFR across {} periods.",
                                    self.cfr_series.len()
                                );
                            }
                            Err(e) => *status_banner = format!("CFR CSV error: {e}"),
                        },
                        Err(e) => *status_banner = format!("File read error: {e}"),
                    }
                }
                if ui.button("Reset Built-In 60-Day Emerging Wave (Early Day 26)").clicked() {
                    self.cfr_series = builtin_emerging_wave_cases_and_deaths();
                    self.cfr_trunc_cutoff_days = 26;
                    self.cfr_cfg = DelayCfrConfig::default();
                    let _ = self.recompute_cfr();
                }
                if ui.button("Show Full Resolved Wave (Day 60)").clicked() {
                    self.cfr_trunc_cutoff_days = self.cfr_series.len();
                    let _ = self.recompute_cfr();
                }
            });
        });

        ui.add_space(6.0);

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.label(
                    RichText::new("Real-Time cCFR vs Downward-Biased Naive Crude CFR (%)").strong(),
                );
                let adj_pts: Vec<[f64; 2]> = self
                    .cfr_result
                    .trajectory
                    .iter()
                    .skip(4)
                    .map(|p| [p.step_index as f64, p.delay_adjusted_ccfr_pct])
                    .collect();
                let crude_pts: Vec<[f64; 2]> = self
                    .cfr_result
                    .trajectory
                    .iter()
                    .skip(4)
                    .map(|p| [p.step_index as f64, p.crude_cfr_pct])
                    .collect();

                Plot::new("ccfr_vs_crude_plot")
                    .height(225.0)
                    .legend(Legend::default())
                    .x_axis_label("Outbreak Day / Period Index")
                    .y_axis_label("Estimated Case Fatality Risk (%)")
                    .show(ui, |plot_ui| {
                        plot_ui.line(
                            Line::new(PlotPoints::from(adj_pts))
                                .name("Nishiura Delay-Adjusted cCFR (%)")
                                .color(Color32::from_rgb(56, 189, 248))
                                .width(2.5),
                        );
                        plot_ui.line(
                            Line::new(PlotPoints::from(crude_pts))
                                .name("Unadjusted Naive Crude CFR (%)")
                                .color(Color32::from_rgb(244, 63, 94))
                                .width(2.0),
                        );
                    });
            });

            cols[1].group(|ui| {
                ui.label(
                    RichText::new("Cumulative Raw Cases vs Outcome-Resolved Denominator u_t")
                        .strong(),
                );
                let cum_c_pts: Vec<[f64; 2]> = self
                    .cfr_result
                    .trajectory
                    .iter()
                    .map(|p| [p.step_index as f64, p.cumulative_cases])
                    .collect();
                let ut_pts: Vec<[f64; 2]> = self
                    .cfr_result
                    .trajectory
                    .iter()
                    .map(|p| [p.step_index as f64, p.resolved_denominator_ut])
                    .collect();

                Plot::new("ccfr_denominators_plot")
                    .height(225.0)
                    .legend(Legend::default())
                    .x_axis_label("Outbreak Day / Period Index")
                    .y_axis_label("Cumulative Cases")
                    .show(ui, |plot_ui| {
                        plot_ui.line(
                            Line::new(PlotPoints::from(cum_c_pts))
                                .name("Cumulative Reported Cases C_t")
                                .color(Color32::from_rgb(251, 191, 36))
                                .width(2.2),
                        );
                        plot_ui.line(
                            Line::new(PlotPoints::from(ut_pts))
                                .name("Nishiura Outcome-Resolved u_t")
                                .color(Color32::from_rgb(52, 211, 153))
                                .width(2.4),
                        );
                    });
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.label(
                RichText::new("Period-by-Period Right-Truncation Adjustment Ledger").strong(),
            );
            egui::ScrollArea::vertical()
                .id_source("cfr_trajectory_scroll")
                .max_height(195.0)
                .show(ui, |ui| {
                    egui::Grid::new("cfr_trajectory_grid")
                        .striped(true)
                        .min_col_width(68.0)
                        .show(ui, |ui| {
                            ui.strong("Period");
                            ui.strong("Cum Cases");
                            ui.strong("Resolved u_t");
                            ui.strong("Cum Deaths");
                            ui.strong("Naive CFR");
                            ui.strong("Adjusted cCFR [95% CI]");
                            ui.strong("Bias Factor");
                            ui.end_row();

                            for p in self.cfr_result.trajectory.iter().rev().take(30) {
                                ui.monospace(&p.period_label);
                                ui.monospace(format!("{:.0}", p.cumulative_cases));
                                ui.monospace(format!("{:.1}", p.resolved_denominator_ut));
                                ui.monospace(format!("{:.1}", p.cumulative_deaths));
                                ui.monospace(format!("{:.2}%", p.crude_cfr_pct));
                                ui.monospace(format!(
                                    "{:.2}% [{:.2}%, {:.2}%]",
                                    p.delay_adjusted_ccfr_pct,
                                    p.adjusted_ci_low_95,
                                    p.adjusted_ci_high_95
                                ));
                                ui.monospace(format!("{:.2}x", p.truncation_bias_factor));
                                ui.end_row();
                            }
                        });
                });
        });
    }

    fn render_surge_panel(
        &mut self,
        ui: &mut egui::Ui,
        active_seir: &SeirForecastResult,
        status_banner: &mut String,
    ) {
        let s = &self.surge_result;

        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("PEAK ACUTE FLOOR CENSUS").small().weak());
                let col = if s.first_ward_breach_day.is_some() {
                    Color32::from_rgb(244, 63, 94)
                } else {
                    Color32::from_rgb(52, 211, 153)
                };
                ui.heading(
                    RichText::new(format!(
                        "{:.1} / {:.0} beds ({:.0}%)",
                        s.peak_ward_census,
                        self.surge_cfg.staffed_ward_beds,
                        s.peak_ward_utilization_pct
                    ))
                    .monospace()
                    .color(col),
                );
                let note = s
                    .first_ward_breach_day
                    .map(|d| format!("Breach on Day {d} (Peak Day {})", s.peak_ward_day))
                    .unwrap_or_else(|| format!("Within capacity (Peak Day {})", s.peak_ward_day));
                ui.label(RichText::new(note).small().monospace());
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("PEAK TOTAL ICU CENSUS").small().weak());
                let col = if s.first_icu_breach_day.is_some() {
                    Color32::from_rgb(244, 63, 94)
                } else {
                    Color32::from_rgb(56, 189, 248)
                };
                ui.heading(
                    RichText::new(format!(
                        "{:.1} / {:.0} ICU ({:.0}%)",
                        s.peak_icu_census,
                        self.surge_cfg.staffed_icu_beds,
                        s.peak_icu_utilization_pct
                    ))
                    .monospace()
                    .color(col),
                );
                ui.label(
                    RichText::new(format!(
                        "Deficit: {:.1} ICU bed-days",
                        s.icu_bed_days_deficit
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[2].group(|ui| {
                ui.label(RichText::new("PEAK INVASIVE VENTILATORS").small().weak());
                let col = if s.first_vent_breach_day.is_some() {
                    Color32::from_rgb(244, 63, 94)
                } else {
                    Color32::from_rgb(251, 191, 36)
                };
                ui.heading(
                    RichText::new(format!(
                        "{:.1} / {:.0} IMV ({:.0}%)",
                        s.peak_vent_census,
                        self.surge_cfg.available_ventilators,
                        s.peak_vent_utilization_pct
                    ))
                    .monospace()
                    .color(col),
                );
                let note = s
                    .first_vent_breach_day
                    .map(|d| format!("Ventilator breach Day {d}"))
                    .unwrap_or_else(|| "Within ventilator reserve".to_string());
                ui.label(RichText::new(note).small().monospace());
            });

            cols[3].group(|ui| {
                ui.label(RichText::new("COMPETING-RISK HFR / DEATHS").small().weak());
                ui.heading(
                    RichText::new(format!(
                        "HFR = {:.2}% ({:.1} deaths)",
                        s.implied_hospital_fatality_pct, s.total_projected_hospital_deaths
                    ))
                    .monospace()
                    .color(Color32::from_rgb(167, 139, 250)),
                );
                ui.label(
                    RichText::new(format!(
                        "Total admissions: {:.1} | IFR {:.3}%",
                        s.total_projected_admissions, s.implied_infection_fatality_pct
                    ))
                    .small()
                    .monospace(),
                );
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            let mut ch = false;
            ui.horizontal_wrapped(|ui| {
                ch |= ui
                    .checkbox(
                        &mut self.drive_from_active_seir,
                        "Couple Directly to Active Tab 2 SEIR Forecast",
                    )
                    .changed();
                if self.drive_from_active_seir {
                    ch |= ui
                        .checkbox(
                            &mut self.surge_cfg.use_true_seir_infections,
                            "Use True Latent SEIR Infections",
                        )
                        .changed();
                } else {
                    ui.label("Wave Scale:");
                    ch |= ui
                        .add(
                            egui::DragValue::new(&mut self.built_in_scale_multiplier)
                                .speed(0.25)
                                .range(0.5..=50.0),
                        )
                        .changed();
                }

                ui.separator();
                ui.label("P(Hosp):");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.surge_cfg.prob_hospital_admission)
                            .speed(0.002)
                            .range(0.001..=0.50),
                    )
                    .changed();
                ui.label("P(ICU|Hosp):");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.surge_cfg.prob_icu_given_hosp)
                            .speed(0.01)
                            .range(0.01..=0.80),
                    )
                    .changed();
                ui.label("P(Vent|ICU):");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.surge_cfg.prob_vent_given_icu)
                            .speed(0.01)
                            .range(0.01..=0.95),
                    )
                    .changed();
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("Mean LOS (d) — Acute:");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.surge_cfg.ward_only_los_days)
                            .speed(0.2)
                            .range(1.0..=30.0),
                    )
                    .changed();
                ui.label("ICU:");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.surge_cfg.icu_nonvent_los_days)
                            .speed(0.2)
                            .range(1.0..=35.0),
                    )
                    .changed();
                ui.label("Vent:");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.surge_cfg.vent_los_days)
                            .speed(0.2)
                            .range(1.0..=45.0),
                    )
                    .changed();

                ui.separator();
                ui.label("Staffed Capacity — Acute Beds:");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.surge_cfg.staffed_ward_beds)
                            .speed(5.0)
                            .range(10.0..=50_000.0),
                    )
                    .changed();
                ui.label("ICU Beds:");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.surge_cfg.staffed_icu_beds)
                            .speed(2.0)
                            .range(2.0..=10_000.0),
                    )
                    .changed();
                ui.label("Ventilators:");
                ch |= ui
                    .add(
                        egui::DragValue::new(&mut self.surge_cfg.available_ventilators)
                            .speed(1.0)
                            .range(1.0..=10_000.0),
                    )
                    .changed();
            });

            if ch {
                if let Err(e) = self.recompute_surge(active_seir) {
                    *status_banner = format!("Surge update warning: {e}");
                }
            }
        });

        ui.add_space(6.0);

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.label(
                    RichText::new("Projected Daily Acute, ICU & Ventilator Census vs Capacity")
                        .strong(),
                );
                let w_pts: Vec<[f64; 2]> = self
                    .surge_result
                    .daily_projections
                    .iter()
                    .map(|p| [p.day_index as f64, p.ward_census_mean])
                    .collect();
                let icu_pts: Vec<[f64; 2]> = self
                    .surge_result
                    .daily_projections
                    .iter()
                    .map(|p| [p.day_index as f64, p.icu_census_mean])
                    .collect();
                let icu_hi: Vec<[f64; 2]> = self
                    .surge_result
                    .daily_projections
                    .iter()
                    .map(|p| [p.day_index as f64, p.icu_census_high_95])
                    .collect();
                let vent_pts: Vec<[f64; 2]> = self
                    .surge_result
                    .daily_projections
                    .iter()
                    .map(|p| [p.day_index as f64, p.vent_census_mean])
                    .collect();

                Plot::new("surge_daily_census_plot")
                    .height(225.0)
                    .legend(Legend::default())
                    .x_axis_label("Projection Day")
                    .y_axis_label("Occupied Beds / Ventilators")
                    .show(ui, |plot_ui| {
                        plot_ui.line(
                            Line::new(PlotPoints::from(w_pts))
                                .name("Acute Floor Census")
                                .color(Color32::from_rgb(52, 211, 153))
                                .width(2.2),
                        );
                        plot_ui.line(
                            Line::new(PlotPoints::from(icu_pts))
                                .name("Total ICU Census")
                                .color(Color32::from_rgb(56, 189, 248))
                                .width(2.5),
                        );
                        plot_ui.line(
                            Line::new(PlotPoints::from(icu_hi))
                                .name("ICU 95% Upper Bound")
                                .stroke(Stroke::new(1.2, Color32::from_rgb(125, 211, 252))),
                        );
                        plot_ui.line(
                            Line::new(PlotPoints::from(vent_pts))
                                .name("Invasive Ventilators")
                                .color(Color32::from_rgb(251, 191, 36))
                                .width(2.2),
                        );
                        plot_ui.hline(
                            HLine::new(self.surge_cfg.staffed_icu_beds)
                                .name("Staffed ICU Capacity")
                                .color(Color32::from_rgb(244, 63, 94)),
                        );
                    });
            });

            cols[1].group(|ui| {
                ui.label(
                    RichText::new(
                        "Admitted Cohort Competing-Risk State & Cumulative Incidence (CIF)",
                    )
                    .strong(),
                );
                let disch_pts: Vec<[f64; 2]> = self
                    .surge_result
                    .cohort_competing_cif
                    .iter()
                    .map(|c| [c.day_since_admission as f64, c.cif_discharged_alive * 100.0])
                    .collect();
                let dead_pts: Vec<[f64; 2]> = self
                    .surge_result
                    .cohort_competing_cif
                    .iter()
                    .map(|c| [c.day_since_admission as f64, c.cif_in_hospital_death * 100.0])
                    .collect();
                let in_hosp_pts: Vec<[f64; 2]> = self
                    .surge_result
                    .cohort_competing_cif
                    .iter()
                    .map(|c| {
                        [
                            c.day_since_admission as f64,
                            (c.prop_in_ward + c.prop_in_icu_nonvent + c.prop_on_ventilator) * 100.0,
                        ]
                    })
                    .collect();

                Plot::new("competing_risk_cif_plot")
                    .height(225.0)
                    .legend(Legend::default())
                    .x_axis_label("Days Since Hospital Admission")
                    .y_axis_label("Cohort Probability (%)")
                    .show(ui, |plot_ui| {
                        plot_ui.line(
                            Line::new(PlotPoints::from(in_hosp_pts))
                                .name("Still Hospitalized (%)")
                                .color(Color32::from_rgb(56, 189, 248))
                                .width(2.2),
                        );
                        plot_ui.line(
                            Line::new(PlotPoints::from(disch_pts))
                                .name("CIF Alive Discharge (%)")
                                .color(Color32::from_rgb(52, 211, 153))
                                .width(2.3),
                        );
                        plot_ui.line(
                            Line::new(PlotPoints::from(dead_pts))
                                .name("CIF In-Hospital Death (%)")
                                .color(Color32::from_rgb(244, 63, 94))
                                .width(2.3),
                        );
                    });
            });
        });
    }
}
