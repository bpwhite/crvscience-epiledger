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
    builtin_three_source_tb_pertussis_crc, compute_two_source_chapman,
    fit_three_source_loglinear_crc, parse_crc_case_lists_csv, ThreeSourceCrcAnalysis,
    ThreeSourceOverlapCounts, TwoSourceChapmanResult,
};
use epi_core::{SurveillanceCampaignPreset, TimeSeriesDataset};
use epi_ingest::{
    build_syndrome_weekly_dataset, builtin_ed_chief_complaint_encounters,
    classify_ed_encounters_essence, parse_ed_encounters_csv, EdTriageEncounter,
    EssenceBatchClassification, EssenceSyndromeId,
};

#[derive(Default)]
pub struct Phase11GuiAction {
    pub commit_syndrome_vintage: Option<TimeSeriesDataset>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase11SubView {
    CaptureRecapture,
    SyndromicEssence,
}

pub struct CrcSyndromicWorkbenchState {
    pub sub_view: Phase11SubView,

    // 3-Source Log-Linear & 2-Source Chapman Capture-Recapture state
    pub crc_counts: ThreeSourceOverlapCounts,
    pub crc_analysis: ThreeSourceCrcAnalysis,
    pub crc_csv_path: String,
    pub chapman_a_name: String,
    pub chapman_b_name: String,
    pub chapman_n10: u64,
    pub chapman_n01: u64,
    pub chapman_n11: u64,
    pub standalone_chapman: TwoSourceChapmanResult,

    // NSSP ESSENCE Chief-Complaint & ICD-10-CM Classifier state
    pub ed_encounters: Vec<EdTriageEncounter>,
    pub custom_essence_query: String,
    pub neg_window_tokens: usize,
    pub selected_syndrome: EssenceSyndromeId,
    pub essence_batch: EssenceBatchClassification,
    pub ed_csv_path: String,
}

impl Default for CrcSyndromicWorkbenchState {
    fn default() -> Self {
        Self::new()
    }
}

impl CrcSyndromicWorkbenchState {
    pub fn new() -> Self {
        let crc_counts = builtin_three_source_tb_pertussis_crc();
        let crc_analysis = fit_three_source_loglinear_crc(&crc_counts)
            .expect("Built-in 3-source CRC counts are valid");

        let chapman_a_name = "ELR Laboratory".to_string();
        let chapman_b_name = "Hospital Discharge".to_string();
        let chapman_n10 = 87;
        let chapman_n01 = 27;
        let chapman_n11 = 57;
        let standalone_chapman = compute_two_source_chapman(
            &chapman_a_name,
            &chapman_b_name,
            chapman_n10,
            chapman_n01,
            chapman_n11,
        );

        let ed_encounters = builtin_ed_chief_complaint_encounters();
        let custom_essence_query =
            "(^fever^ OR ^cough^ OR ^j10^) ANDNOT (^immunization^ OR ^heat pack^)".to_string();
        let neg_window_tokens = 5;
        let essence_batch = classify_ed_encounters_essence(
            &ed_encounters,
            &custom_essence_query,
            neg_window_tokens,
        );

        Self {
            sub_view: Phase11SubView::CaptureRecapture,
            crc_counts,
            crc_analysis,
            crc_csv_path: String::new(),
            chapman_a_name,
            chapman_b_name,
            chapman_n10,
            chapman_n01,
            chapman_n11,
            standalone_chapman,
            ed_encounters,
            custom_essence_query,
            neg_window_tokens,
            selected_syndrome: EssenceSyndromeId::IliRespiratory,
            essence_batch,
            ed_csv_path: String::new(),
        }
    }

    pub fn recompute_crc(&mut self) -> Result<(), String> {
        self.crc_analysis =
            fit_three_source_loglinear_crc(&self.crc_counts).map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn recompute_chapman(&mut self) {
        self.standalone_chapman = compute_two_source_chapman(
            &self.chapman_a_name,
            &self.chapman_b_name,
            self.chapman_n10,
            self.chapman_n01,
            self.chapman_n11,
        );
    }

    pub fn recompute_essence(&mut self) {
        self.essence_batch = classify_ed_encounters_essence(
            &self.ed_encounters,
            &self.custom_essence_query,
            self.neg_window_tokens,
        );
    }

    pub fn sync_to_campaign(
        &mut self,
        campaign: SurveillanceCampaignPreset,
        active_ds: &TimeSeriesDataset,
    ) {
        let base = builtin_three_source_tb_pertussis_crc();
        let recent_mean = active_ds
            .points
            .iter()
            .rev()
            .take(8)
            .map(|p| p.count)
            .sum::<f64>()
            / 8.0;
        let scale = (recent_mean / 220.0).clamp(0.45, 4.0);
        self.crc_counts = match campaign {
            SurveillanceCampaignPreset::UsPertussis2022_2026 => base,
            SurveillanceCampaignPreset::UsMeasles2022_2026 => ThreeSourceOverlapCounts {
                s1_name: "Measles PCR/IgM Lab".to_string(),
                s2_name: "ED Rash/Fever Discharge".to_string(),
                s3_name: "NNDSS Case Report".to_string(),
                n_100: ((base.n_100 as f64) * scale).round() as u64,
                n_010: ((base.n_010 as f64) * scale).round() as u64,
                n_001: ((base.n_001 as f64) * scale).round() as u64,
                n_110: ((base.n_110 as f64) * scale).round() as u64,
                n_101: ((base.n_101 as f64) * scale).round() as u64,
                n_011: ((base.n_011 as f64) * scale).round() as u64,
                n_111: ((base.n_111 as f64) * scale).round() as u64,
            },
            SurveillanceCampaignPreset::InfluenzaIliWastewater2022_2026 => ThreeSourceOverlapCounts {
                s1_name: "State Virology Lab".to_string(),
                s2_name: "ILINet Inpatient".to_string(),
                s3_name: "ESSENCE Respiratory".to_string(),
                n_100: ((base.n_100 as f64) * 2.2).round() as u64,
                n_010: ((base.n_010 as f64) * 2.2).round() as u64,
                n_001: ((base.n_001 as f64) * 2.2).round() as u64,
                n_110: ((base.n_110 as f64) * 2.2).round() as u64,
                n_101: ((base.n_101 as f64) * 2.2).round() as u64,
                n_011: ((base.n_011 as f64) * 2.2).round() as u64,
                n_111: ((base.n_111 as f64) * 2.2).round() as u64,
            },
            SurveillanceCampaignPreset::SarsCoV2VariantCohort2022_2026 => ThreeSourceOverlapCounts {
                s1_name: "ELR Nucleic Acid Test".to_string(),
                s2_name: "NHSN Admission".to_string(),
                s3_name: "ED COVID-Like Visit".to_string(),
                n_100: ((base.n_100 as f64) * 1.9).round() as u64,
                n_010: ((base.n_010 as f64) * 1.9).round() as u64,
                n_001: ((base.n_001 as f64) * 1.9).round() as u64,
                n_110: ((base.n_110 as f64) * 1.9).round() as u64,
                n_101: ((base.n_101 as f64) * 1.9).round() as u64,
                n_011: ((base.n_011 as f64) * 1.9).round() as u64,
                n_111: ((base.n_111 as f64) * 1.9).round() as u64,
            },
        };
        self.custom_essence_query = match campaign {
            SurveillanceCampaignPreset::UsMeasles2022_2026 => {
                "(^fever^ OR ^rash^ OR ^measles^ OR ^b05^) ANDNOT (^immunization^ OR ^vaccine^)"
            }
            SurveillanceCampaignPreset::UsPertussis2022_2026 => {
                "(^cough^ OR ^pertussis^ OR ^whooping^ OR ^a37^) ANDNOT (^immunization^ OR ^vaccine^)"
            }
            SurveillanceCampaignPreset::InfluenzaIliWastewater2022_2026 => {
                "(^fever^ OR ^cough^ OR ^j10^ OR ^influenza^) ANDNOT (^immunization^ OR ^heat pack^)"
            }
            SurveillanceCampaignPreset::SarsCoV2VariantCohort2022_2026 => {
                "(^fever^ OR ^cough^ OR ^covid^ OR ^u07^) ANDNOT (^immunization^ OR ^negative^)"
            }
        }
        .to_string();
        let _ = self.recompute_crc();
        self.recompute_essence();
    }

    pub fn render(&mut self, ui: &mut egui::Ui, status_banner: &mut String) -> Phase11GuiAction {
        let mut action = Phase11GuiAction::default();

        ui.horizontal_wrapped(|ui| {
            ui.heading("Capture-Recapture Completeness & NSSP ESSENCE Syndromic Classifier");
            ui.separator();
            ui.selectable_value(
                &mut self.sub_view,
                Phase11SubView::CaptureRecapture,
                "3-Source Log-Linear & Chapman CRC",
            );
            ui.selectable_value(
                &mut self.sub_view,
                Phase11SubView::SyndromicEssence,
                "NSSP ESSENCE Chief-Complaint & NegEx",
            );
        });
        ui.add_space(4.0);

        match self.sub_view {
            Phase11SubView::CaptureRecapture => self.render_crc_panel(ui, status_banner),
            Phase11SubView::SyndromicEssence => {
                self.render_essence_panel(ui, status_banner, &mut action)
            }
        }

        action
    }

    fn render_crc_panel(&mut self, ui: &mut egui::Ui, status_banner: &mut String) {
        let best = &self.crc_analysis.models[self.crc_analysis.bic_selected_index];

        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("BIC-SELECTED LOG-LINEAR MODEL").small().weak());
                ui.heading(
                    RichText::new(format!("{} ({})", best.model_id, best.formula_label))
                        .monospace()
                        .color(Color32::from_rgb(56, 189, 248)),
                );
                ui.label(
                    RichText::new(format!(
                        "BIC {:.2} | AIC {:.2} | G2={:.2} (df={})",
                        best.bic, best.aic, best.deviance_g2, best.residual_df
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("ESTIMATED TOTAL BURDEN N_HAT").small().weak());
                ui.heading(
                    RichText::new(format!(
                        "{:.1} cases (+{:.1} missed)",
                        best.total_n_hat, best.unobserved_n000_hat
                    ))
                    .monospace()
                    .color(Color32::from_rgb(251, 191, 36)),
                );
                ui.label(
                    RichText::new(format!(
                        "95% CI [{:.1}, {:.1}] | Obs={}",
                        best.n_ci_low_95, best.n_ci_high_95, self.crc_analysis.total_observed
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[2].group(|ui| {
                ui.label(RichText::new("SYSTEM REPORTING COMPLETENESS").small().weak());
                let col = if best.overall_completeness_pct >= 80.0 {
                    Color32::from_rgb(52, 211, 153)
                } else {
                    Color32::from_rgb(244, 63, 94)
                };
                ui.heading(
                    RichText::new(format!("{:.1}% ascertained", best.overall_completeness_pct))
                        .monospace()
                        .color(col),
                );
                ui.label(
                    RichText::new(format!(
                        "S1 {:.1}% | S2 {:.1}% | S3 {:.1}%",
                        best.s1_sensitivity_pct, best.s2_sensitivity_pct, best.s3_sensitivity_pct
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[3].group(|ui| {
                ui.label(RichText::new("2-STREAM CHAPMAN ESTIMATOR").small().weak());
                ui.heading(
                    RichText::new(format!(
                        "N_chap = {:.1} ({:.1}%)",
                        self.standalone_chapman.chapman_total_n_hat,
                        self.standalone_chapman.overall_completeness_pct
                    ))
                    .monospace()
                    .color(Color32::from_rgb(167, 139, 250)),
                );
                ui.label(
                    RichText::new(format!(
                        "95% CI [{:.1}, {:.1}] | SE={:.1}",
                        self.standalone_chapman.ci_low_95,
                        self.standalone_chapman.ci_high_95,
                        self.standalone_chapman.seber_se
                    ))
                    .small()
                    .monospace(),
                );
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("3-Source Venn Cells:").strong());
                let mut changed = false;
                ui.label("n100:");
                changed |= ui
                    .add(egui::DragValue::new(&mut self.crc_counts.n_100).range(0..=50_000))
                    .changed();
                ui.label("n010:");
                changed |= ui
                    .add(egui::DragValue::new(&mut self.crc_counts.n_010).range(0..=50_000))
                    .changed();
                ui.label("n001:");
                changed |= ui
                    .add(egui::DragValue::new(&mut self.crc_counts.n_001).range(0..=50_000))
                    .changed();
                ui.label("n110:");
                changed |= ui
                    .add(egui::DragValue::new(&mut self.crc_counts.n_110).range(0..=50_000))
                    .changed();
                ui.label("n101:");
                changed |= ui
                    .add(egui::DragValue::new(&mut self.crc_counts.n_101).range(0..=50_000))
                    .changed();
                ui.label("n011:");
                changed |= ui
                    .add(egui::DragValue::new(&mut self.crc_counts.n_011).range(0..=50_000))
                    .changed();
                ui.label("n111:");
                changed |= ui
                    .add(egui::DragValue::new(&mut self.crc_counts.n_111).range(0..=50_000))
                    .changed();

                if changed {
                    if let Err(e) = self.recompute_crc() {
                        *status_banner = format!("CRC refit warning: {e}");
                    }
                }

                if ui.button("Reset Built-In 3-Registry Cohort").clicked() {
                    self.crc_counts = builtin_three_source_tb_pertussis_crc();
                    let _ = self.recompute_crc();
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("2-Source Chapman:").strong());
                let mut ch2 = false;
                ui.label("A only (n10):");
                ch2 |= ui
                    .add(egui::DragValue::new(&mut self.chapman_n10).range(0..=100_000))
                    .changed();
                ui.label("B only (n01):");
                ch2 |= ui
                    .add(egui::DragValue::new(&mut self.chapman_n01).range(0..=100_000))
                    .changed();
                ui.label("Both A&B (n11):");
                ch2 |= ui
                    .add(egui::DragValue::new(&mut self.chapman_n11).range(0..=100_000))
                    .changed();
                if ch2 {
                    self.recompute_chapman();
                }

                ui.separator();
                ui.label("Case-List CSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.crc_csv_path)
                        .desired_width(205.0)
                        .hint_text("path/to/crc_case_lists.csv"),
                );
                if ui.button("Browse...").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("CSV / TSV", &["csv", "tsv", "txt"])
                        .pick_file()
                    {
                        self.crc_csv_path = p.display().to_string();
                    }
                }
                if ui.button("Import CRC CSV").clicked() {
                    match std::fs::read_to_string(self.crc_csv_path.trim()) {
                        Ok(raw) => match parse_crc_case_lists_csv(&raw, "Lab", "Hospital", "Vital")
                        {
                            Ok(cnts) => {
                                self.crc_counts = cnts;
                                match self.recompute_crc() {
                                    Ok(()) => {
                                        *status_banner = format!(
                                            "Fitted 8 log-linear CRC models on {} observed cases.",
                                            self.crc_analysis.total_observed
                                        );
                                    }
                                    Err(e) => *status_banner = format!("CRC model error: {e}"),
                                }
                            }
                            Err(e) => *status_banner = format!("CRC CSV error: {e}"),
                        },
                        Err(e) => *status_banner = format!("File read error: {e}"),
                    }
                }
            });
        });

        ui.add_space(6.0);

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.label(
                    RichText::new("Estimated Total Population N_hat Across All 8 Log-Linear Models")
                        .strong(),
                );
                let bic_idx = self.crc_analysis.bic_selected_index;
                let bars: Vec<Bar> = self
                    .crc_analysis
                    .models
                    .iter()
                    .enumerate()
                    .map(|(i, m)| {
                        let col = if i == bic_idx {
                            Color32::from_rgb(52, 211, 153)
                        } else {
                            Color32::from_rgb(96, 165, 250)
                        };
                        Bar::new(i as f64, m.total_n_hat)
                            .name(&m.model_id)
                            .width(0.6)
                            .fill(col)
                    })
                    .collect();

                Plot::new("crc_models_nhat_plot")
                    .height(215.0)
                    .legend(Legend::default())
                    .x_axis_label("Hierarchical Log-Linear Model Index (0=M0 .. 7=Saturated)")
                    .y_axis_label("Estimated Total Cases N_hat")
                    .show(ui, |plot_ui| {
                        plot_ui.bar_chart(BarChart::new(bars).name("Model N_hat"));
                        plot_ui.hline(
                            HLine::new(self.crc_analysis.total_observed as f64)
                                .name("Observed Unique Cases")
                                .color(Color32::from_rgb(244, 63, 94)),
                        );
                    });
            });

            cols[1].group(|ui| {
                let best_m = &self.crc_analysis.models[self.crc_analysis.bic_selected_index];
                ui.label(
                    RichText::new("Registry Sensitivity (%) Under BIC-Optimal Model").strong(),
                );
                let sens_bars = vec![
                    Bar::new(0.0, best_m.overall_completeness_pct)
                        .name("Combined 3-Source")
                        .width(0.55)
                        .fill(Color32::from_rgb(52, 211, 153)),
                    Bar::new(1.0, best_m.s1_sensitivity_pct)
                        .name(&self.crc_counts.s1_name)
                        .width(0.55)
                        .fill(Color32::from_rgb(56, 189, 248)),
                    Bar::new(2.0, best_m.s2_sensitivity_pct)
                        .name(&self.crc_counts.s2_name)
                        .width(0.55)
                        .fill(Color32::from_rgb(251, 191, 36)),
                    Bar::new(3.0, best_m.s3_sensitivity_pct)
                        .name(&self.crc_counts.s3_name)
                        .width(0.55)
                        .fill(Color32::from_rgb(167, 139, 250)),
                ];
                Plot::new("crc_sensitivity_bars_plot")
                    .height(215.0)
                    .legend(Legend::default())
                    .x_axis_label("Surveillance Stream (0=Combined, 1=S1, 2=S2, 3=S3)")
                    .y_axis_label("Estimated Sensitivity (%)")
                    .show(ui, |plot_ui| {
                        plot_ui.bar_chart(BarChart::new(sens_bars).name("Sensitivity %"));
                    });
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.label(
                RichText::new(
                    "8-Model Hierarchical Log-Linear Capture-Recapture Comparison (Rcapture Family)",
                )
                .strong(),
            );
            egui::ScrollArea::vertical()
                .id_source("crc_models_table_scroll")
                .max_height(205.0)
                .show(ui, |ui| {
                    egui::Grid::new("crc_models_table_grid")
                        .striped(true)
                        .min_col_width(58.0)
                        .show(ui, |ui| {
                            ui.strong("Model");
                            ui.strong("Interaction Structure");
                            ui.strong("df");
                            ui.strong("Deviance G2");
                            ui.strong("AIC");
                            ui.strong("BIC");
                            ui.strong("Unobs n000");
                            ui.strong("Total N_hat [95% CI]");
                            ui.strong("Completeness");
                            ui.strong("OR12 / OR13 / OR23");
                            ui.end_row();

                            for (idx, m) in self.crc_analysis.models.iter().enumerate() {
                                let is_bic = idx == self.crc_analysis.bic_selected_index;
                                let tag = if is_bic {
                                    format!("*{}*", m.model_id)
                                } else {
                                    m.model_id.clone()
                                };
                                ui.monospace(tag);
                                ui.label(&m.formula_label);
                                ui.monospace(format!("{}", m.residual_df));
                                ui.monospace(format!("{:.2}", m.deviance_g2));
                                ui.monospace(format!("{:.2}", m.aic));
                                ui.monospace(format!("{:.2}", m.bic));
                                ui.monospace(format!("{:.1}", m.unobserved_n000_hat));
                                ui.monospace(format!(
                                    "{:.1} [{:.1}, {:.1}]",
                                    m.total_n_hat, m.n_ci_low_95, m.n_ci_high_95
                                ));
                                ui.monospace(format!("{:.1}%", m.overall_completeness_pct));
                                ui.monospace(format!(
                                    "{:.2} / {:.2} / {:.2}",
                                    m.or_12, m.or_13, m.or_23
                                ));
                                ui.end_row();
                            }
                        });
                });
        });
    }

    fn render_essence_panel(
        &mut self,
        ui: &mut egui::Ui,
        status_banner: &mut String,
        action: &mut Phase11GuiAction,
    ) {
        let b = &self.essence_batch;

        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("ED TRIAGE VISITS SCANNED").small().weak());
                ui.heading(
                    RichText::new(format!("{} visits", b.total_encounters))
                        .monospace()
                        .color(Color32::from_rgb(56, 189, 248)),
                );
                ui.label(
                    RichText::new(format!(
                        "{} matched >=1 syndrome",
                        b.encounters_matching_any_syndrome
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("NEGEX / EXCLUSIONS BLOCKED").small().weak());
                ui.heading(
                    RichText::new(format!("{} false-positives blocked", b.total_negations_blocked))
                        .monospace()
                        .color(Color32::from_rgb(52, 211, 153)),
                );
                ui.label(
                    RichText::new(format!(
                        "Negation token window = {}",
                        self.neg_window_tokens
                    ))
                    .small()
                    .monospace(),
                );
            });

            cols[2].group(|ui| {
                ui.label(RichText::new("SELECTED SYNDROME VOLUME").small().weak());
                let sel_row = b
                    .syndrome_summaries
                    .iter()
                    .find(|r| r.syndrome == self.selected_syndrome);
                if let Some(sr) = sel_row {
                    ui.heading(
                        RichText::new(format!(
                            "{} ({:.1}%)",
                            sr.matched_visits, sr.pct_of_ed_visits
                        ))
                        .monospace()
                        .color(Color32::from_rgb(251, 191, 36)),
                    );
                    ui.label(RichText::new(&sr.code).small().monospace());
                }
            });

            cols[3].group(|ui| {
                ui.label(RichText::new("CUSTOM ESSENCE BOOLEAN HITS").small().weak());
                ui.heading(
                    RichText::new(format!("{} visits", b.custom_query_matches))
                        .monospace()
                        .color(Color32::from_rgb(167, 139, 250)),
                );
                ui.label(RichText::new("OR / ANDNOT CCDD evaluator").small().monospace());
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Custom ESSENCE Query:").strong());
                let q_ch = ui
                    .add(
                        egui::TextEdit::singleline(&mut self.custom_essence_query)
                            .desired_width(370.0),
                    )
                    .changed();
                ui.label("Negation Window:");
                let w_ch = ui
                    .add(egui::DragValue::new(&mut self.neg_window_tokens).range(1..=12))
                    .changed();
                if q_ch || w_ch {
                    self.recompute_essence();
                }

                egui::ComboBox::from_id_source("essence_syndrome_selector")
                    .selected_text(self.selected_syndrome.label())
                    .show_ui(ui, |ui| {
                        for syn in EssenceSyndromeId::ALL {
                            ui.selectable_value(&mut self.selected_syndrome, syn, syn.label());
                        }
                    });

                if ui
                    .button("Commit Syndrome Series to SQLite & Run Rt")
                    .clicked()
                {
                    let ds = build_syndrome_weekly_dataset(&self.essence_batch, self.selected_syndrome);
                    action.commit_syndrome_vintage = Some(ds);
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("ED Triage CSV:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.ed_csv_path)
                        .desired_width(230.0)
                        .hint_text("path/to/ed_chief_complaints.csv"),
                );
                if ui.button("Browse...").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("CSV / TSV", &["csv", "tsv", "txt"])
                        .pick_file()
                    {
                        self.ed_csv_path = p.display().to_string();
                    }
                }
                if ui.button("Import ED Triage CSV").clicked() {
                    match std::fs::read_to_string(self.ed_csv_path.trim()) {
                        Ok(raw) => match parse_ed_encounters_csv(&raw) {
                            Ok(encs) => {
                                self.ed_encounters = encs;
                                self.recompute_essence();
                                *status_banner = format!(
                                    "Classified {} ED encounters ({} syndrome matches, {} negations blocked).",
                                    self.essence_batch.total_encounters,
                                    self.essence_batch.encounters_matching_any_syndrome,
                                    self.essence_batch.total_negations_blocked
                                );
                            }
                            Err(e) => *status_banner = format!("ED CSV error: {e}"),
                        },
                        Err(e) => *status_banner = format!("File read error: {e}"),
                    }
                }
                if ui.button("Reset Built-In 68-Visit ED Feed").clicked() {
                    self.ed_encounters = builtin_ed_chief_complaint_encounters();
                    self.recompute_essence();
                }
            });
        });

        ui.add_space(6.0);

        let weekly_ds = build_syndrome_weekly_dataset(&self.essence_batch, self.selected_syndrome);
        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.label(
                    RichText::new(format!(
                        "Weekly MMWR Trajectory — {}",
                        self.selected_syndrome.label()
                    ))
                    .strong(),
                );
                let pts: Vec<[f64; 2]> = weekly_ds
                    .points
                    .iter()
                    .enumerate()
                    .map(|(i, p)| [i as f64, p.count])
                    .collect();
                Plot::new("essence_weekly_syndrome_plot")
                    .height(205.0)
                    .legend(Legend::default())
                    .x_axis_label("MMWR Week Index")
                    .y_axis_label("Classified ED Visits")
                    .show(ui, |plot_ui| {
                        plot_ui.line(
                            Line::new(PlotPoints::from(pts))
                                .name(self.selected_syndrome.code())
                                .color(Color32::from_rgb(56, 189, 248))
                                .width(2.5_f32),
                        );
                    });
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("NSSP ESSENCE Syndrome Summary & Blocked Negations").strong());
                egui::Grid::new("essence_summary_grid")
                    .striped(true)
                    .min_col_width(64.0)
                    .show(ui, |ui| {
                        ui.strong("Code");
                        ui.strong("Syndrome Definition");
                        ui.strong("Hits");
                        ui.strong("% ED");
                        ui.strong("NegEx Blocked");
                        ui.end_row();

                        for r in &self.essence_batch.syndrome_summaries {
                            ui.monospace(&r.code);
                            ui.label(&r.label);
                            ui.monospace(format!("{}", r.matched_visits));
                            ui.monospace(format!("{:.1}%", r.pct_of_ed_visits));
                            ui.monospace(format!("{}", r.negated_false_positives_blocked));
                            ui.end_row();
                        }
                    });
            });
        });

        ui.add_space(6.0);

        ui.group(|ui| {
            ui.label(
                RichText::new(
                    "Encounter-Level CCDD Classification & Negation Audit Log (Showing Negations First)",
                )
                .strong(),
            );
            egui::ScrollArea::vertical()
                .id_source("essence_encounters_audit_scroll")
                .max_height(200.0)
                .show(ui, |ui| {
                    egui::Grid::new("essence_encounters_audit_grid")
                        .striped(true)
                        .min_col_width(64.0)
                        .show(ui, |ui| {
                            ui.strong("Visit ID");
                            ui.strong("Date");
                            ui.strong("Normalized Chief Complaint");
                            ui.strong("Matched Syndromes");
                            ui.strong("Blocked Negations / Exclusions");
                            ui.strong("Custom");
                            ui.end_row();

                            for enc in self.essence_batch.classified_encounters.iter().take(36) {
                                ui.monospace(&enc.visit_id);
                                ui.monospace(&enc.visit_iso_date);
                                ui.label(&enc.normalized_cc_text);
                                let syn_str = if enc.matched_syndromes.is_empty() {
                                    "--".to_string()
                                } else {
                                    enc.matched_syndromes
                                        .iter()
                                        .map(|s| s.code())
                                        .collect::<Vec<_>>()
                                        .join(", ")
                                };
                                ui.monospace(syn_str);
                                let neg_str = if enc.suppressed_negated_terms.is_empty() {
                                    "--".to_string()
                                } else {
                                    enc.suppressed_negated_terms.join(", ")
                                };
                                ui.monospace(neg_str);
                                ui.monospace(if enc.custom_query_hit { "HIT" } else { "--" });
                                ui.end_row();
                            }
                        });
                });
        });
    }
}
