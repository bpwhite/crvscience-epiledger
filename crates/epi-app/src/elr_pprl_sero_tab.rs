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
use egui_plot::{Bar, BarChart, Legend, Line, Plot, PlotPoints, Points};
use epi_analytics::{
    builtin_age_serosurvey_dataset, builtin_tnd_vaccine_study_dataset,
    estimate_tnd_vaccine_effectiveness, fit_catalytic_seroprevalence, parse_serosurvey_csv,
    CatalyticFitResult, CatalyticModelKind, SerosurveyDataset, TndDataset, TndVeResult,
};
use epi_core::TimeSeriesDataset;
use epi_ingest::{
    builtin_hl7_v2_elr_and_nssp_batch, import_hl7_v2_file, parse_hl7_v2_batch_str,
    Hl7v2IngestOutcome,
};
use epi_store::{
    builtin_pprl_hospital_lab_records, encode_pprl_record, link_and_deduplicate_pprl,
    EncodedPprlRecord, PprlConfig, PprlDemographicInput, PprlLinkageSummary,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase6SubView {
    Hl7v2ElrAndNssp,
    BloomPprlDeduplication,
    SeroCatalyticAndTndVe,
}

#[derive(Default)]
pub struct Phase6GuiAction {
    pub commit_hl7_dataset: Option<TimeSeriesDataset>,
    /// `(initial_immune_fraction, waning_days)` pushed from Sero-Catalytic fit into `seir.rs`.
    pub apply_sero_to_seir: Option<(f64, f64)>,
}

pub struct ElrPprlSeroWorkbenchState {
    pub sub_view: Phase6SubView,

    // HL7 v2.5.1 ELR & NSSP Syndromic stream state
    pub hl7_raw_buffer: String,
    pub hl7_file_path: String,
    pub hl7_dedup_enabled: bool,
    pub cached_hl7: Option<Hl7v2IngestOutcome>,

    // 512-bit HMAC-SHA256 Bloom PPRL state
    pub pprl_cfg: PprlConfig,
    pub pprl_records: Vec<PprlDemographicInput>,
    pub pprl_encoded: Vec<EncodedPprlRecord>,
    pub cached_pprl: PprlLinkageSummary,

    // Sero-Catalytic FOI & Test-Negative Design VE state
    pub sero_dataset: SerosurveyDataset,
    pub sero_kind: CatalyticModelKind,
    pub sero_csv_path: String,
    pub cached_sero: Option<CatalyticFitResult>,
    pub tnd_dataset: TndDataset,
    pub cached_tnd: Option<TndVeResult>,
}

impl Default for ElrPprlSeroWorkbenchState {
    fn default() -> Self {
        Self::new()
    }
}

impl ElrPprlSeroWorkbenchState {
    pub fn new() -> Self {
        let pprl_cfg = PprlConfig::default();
        let hl7_raw_buffer = builtin_hl7_v2_elr_and_nssp_batch();
        let cached_hl7 =
            parse_hl7_v2_batch_str(&hl7_raw_buffer, &pprl_cfg, true, "BuiltinELR_NSSP").ok();

        let pprl_records = builtin_pprl_hospital_lab_records();
        let pprl_encoded: Vec<EncodedPprlRecord> = pprl_records
            .iter()
            .map(|r| encode_pprl_record(r, &pprl_cfg))
            .collect();
        let cached_pprl = link_and_deduplicate_pprl(&pprl_encoded, &pprl_cfg);

        let sero_dataset = builtin_age_serosurvey_dataset();
        let sero_kind = CatalyticModelKind::ReversibleWaningCatalytic;
        let cached_sero = fit_catalytic_seroprevalence(&sero_dataset, sero_kind).ok();

        let tnd_dataset = builtin_tnd_vaccine_study_dataset();
        let cached_tnd = estimate_tnd_vaccine_effectiveness(&tnd_dataset).ok();

        Self {
            sub_view: Phase6SubView::Hl7v2ElrAndNssp,
            hl7_raw_buffer,
            hl7_file_path: String::new(),
            hl7_dedup_enabled: true,
            cached_hl7,
            pprl_cfg,
            pprl_records,
            pprl_encoded,
            cached_pprl,
            sero_dataset,
            sero_kind,
            sero_csv_path: String::new(),
            cached_sero,
            tnd_dataset,
            cached_tnd,
        }
    }

    pub fn recompute_pprl_and_hl7(&mut self) {
        self.pprl_encoded = self
            .pprl_records
            .iter()
            .map(|r| encode_pprl_record(r, &self.pprl_cfg))
            .collect();
        self.cached_pprl = link_and_deduplicate_pprl(&self.pprl_encoded, &self.pprl_cfg);
        self.cached_hl7 = parse_hl7_v2_batch_str(
            &self.hl7_raw_buffer,
            &self.pprl_cfg,
            self.hl7_dedup_enabled,
            "WorkbenchHL7v2",
        )
        .ok();
    }

    pub fn recompute_sero(&mut self) -> Result<(), String> {
        match fit_catalytic_seroprevalence(&self.sero_dataset, self.sero_kind) {
            Ok(res) => {
                self.cached_sero = Some(res);
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    }

    pub fn render(&mut self, ui: &mut egui::Ui, status_banner: &mut String) -> Phase6GuiAction {
        let mut action = Phase6GuiAction::default();

        ui.horizontal_wrapped(|ui| {
            ui.strong("Clinical Ingestion, Privacy Linkage & Vaccine / Sero Sub-Module:");
            ui.selectable_value(
                &mut self.sub_view,
                Phase6SubView::Hl7v2ElrAndNssp,
                "A. HL7 v2.5.1 ELR & NSSP Syndromic",
            );
            ui.selectable_value(
                &mut self.sub_view,
                Phase6SubView::BloomPprlDeduplication,
                "B. 512-Bit HMAC-SHA256 Bloom PPRL",
            );
            ui.selectable_value(
                &mut self.sub_view,
                Phase6SubView::SeroCatalyticAndTndVe,
                "C. Sero-Catalytic FOI & Test-Negative VE",
            );
        });
        ui.separator();

        match self.sub_view {
            Phase6SubView::Hl7v2ElrAndNssp => {
                self.render_hl7_view(ui, status_banner, &mut action);
            }
            Phase6SubView::BloomPprlDeduplication => {
                self.render_pprl_view(ui, status_banner);
            }
            Phase6SubView::SeroCatalyticAndTndVe => {
                self.render_sero_tnd_view(ui, status_banner, &mut action);
            }
        }

        action
    }

    fn render_hl7_view(
        &mut self,
        ui: &mut egui::Ui,
        status_banner: &mut String,
        action: &mut Phase6GuiAction,
    ) {
        ui.heading("HL7 Version 2.5.1 Pipe-Delimited ELR (ORU^R01) & NSSP Syndromic (ADT) Parser");
        ui.label(
            RichText::new(
                "Zero-copy MSH/PID/OBR/OBX/DG1 stream decoder mapping LOINC/SNOMED-CT lab positives and Emergency Department syndromes with inline 512-bit Bloom PPRL deduplication.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        if let Some(ref out) = self.cached_hl7 {
            ui.columns(4, |cols| {
                cols[0].group(|ui| {
                    ui.label(RichText::new("HL7 v2.5.1 MESSAGES").small().weak());
                    ui.heading(
                        RichText::new(format!("{} parsed", out.total_messages_parsed)).monospace(),
                    );
                    ui.label(RichText::new("ORU^R01 + ADT batch").small().monospace());
                });
                cols[1].group(|ui| {
                    ui.label(RichText::new("RAW POSITIVE / SYNDROMIC").small().weak());
                    ui.heading(
                        RichText::new(format!("{} hits", out.positive_or_syndromic_hits))
                            .monospace()
                            .color(Color32::from_rgb(251, 191, 36)),
                    );
                    ui.label(RichText::new("Before patient dedup").small().monospace());
                });
                cols[2].group(|ui| {
                    ui.label(RichText::new("DEDUPLICATED INCIDENT CASES").small().weak());
                    ui.heading(
                        RichText::new(format!("{} cases", out.deduplicated_incident_cases))
                            .monospace()
                            .color(Color32::from_rgb(52, 211, 153)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Collapsed {} repeat swabs",
                            out.positive_or_syndromic_hits
                                .saturating_sub(out.deduplicated_incident_cases)
                        ))
                        .small()
                        .monospace(),
                    );
                });
                cols[3].group(|ui| {
                    ui.label(RichText::new("LEADING NSSP SYNDROME").small().weak());
                    if let Some((syn, cnt)) = out.syndrome_counts.first() {
                        ui.heading(
                            RichText::new(format!("{syn} ({cnt})"))
                                .monospace()
                                .color(Color32::from_rgb(56, 189, 248)),
                        );
                    } else {
                        ui.heading("--");
                    }
                    ui.label(
                        RichText::new(format!("{} MMWR weeks", out.dataset.points.len()))
                            .small()
                            .monospace(),
                    );
                });
            });
        }

        ui.add_space(6.0);
        let mut hl7_dirty = false;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                if ui
                    .checkbox(
                        &mut self.hl7_dedup_enabled,
                        "Inline 512-Bit Bloom PPRL Patient Deduplication",
                    )
                    .changed()
                {
                    hl7_dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.pprl_cfg.dice_match_threshold, 0.60..=0.96)
                            .text("Dice Match Cutoff"),
                    )
                    .changed()
                {
                    hl7_dirty = true;
                }
                ui.separator();
                ui.add(
                    egui::TextEdit::singleline(&mut self.hl7_file_path)
                        .desired_width(190.0)
                        .hint_text("batch_messages.hl7"),
                );
                if ui.button("Browse...").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("HL7 v2 Messages", &["hl7", "txt"])
                        .pick_file()
                    {
                        self.hl7_file_path = p.display().to_string();
                    }
                }
                if ui.button("Import .hl7 File").clicked() {
                    let p = self.hl7_file_path.trim().to_string();
                    match import_hl7_v2_file(&p, &self.pprl_cfg, self.hl7_dedup_enabled) {
                        Ok(out) => {
                            if let Ok(raw) = std::fs::read_to_string(&p) {
                                self.hl7_raw_buffer = raw;
                            }
                            *status_banner = format!(
                                "Imported {} HL7 v2.5.1 messages ({} deduplicated cases).",
                                out.total_messages_parsed, out.deduplicated_incident_cases
                            );
                            self.cached_hl7 = Some(out);
                        }
                        Err(e) => *status_banner = format!("HL7 v2.5.1 import error: {e}"),
                    }
                }
                if ui.button("Reset Built-in Batch").clicked() {
                    self.hl7_raw_buffer = builtin_hl7_v2_elr_and_nssp_batch();
                    hl7_dirty = true;
                }
                if let Some(ref out) = self.cached_hl7 {
                    if !out.dataset.points.is_empty()
                        && ui
                            .button("Commit Deduplicated Series to Active Vintage")
                            .clicked()
                    {
                        action.commit_hl7_dataset = Some(out.dataset.clone());
                    }
                }
            });
        });

        if hl7_dirty {
            self.recompute_pprl_and_hl7();
        }

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.horizontal(|ui| {
                    ui.strong("Live Pipe-Delimited HL7 v2.5.1 Stream Buffer");
                    if ui.button("Re-Parse Buffer").clicked() {
                        self.recompute_pprl_and_hl7();
                    }
                });
                egui::ScrollArea::vertical().max_height(225.0).show(ui, |ui| {
                    ui.add(
                        egui::TextEdit::multiline(&mut self.hl7_raw_buffer)
                            .desired_rows(11)
                            .desired_width(f32::INFINITY)
                            .font(egui::TextStyle::Monospace),
                    );
                });
            });

            cols[1].group(|ui| {
                ui.strong("NSSP Syndrome Distribution (Deduplicated Incident Hits)");
                if let Some(ref out) = self.cached_hl7 {
                    let bars: Vec<Bar> = out
                        .syndrome_counts
                        .iter()
                        .enumerate()
                        .map(|(i, (_lbl, cnt))| {
                            Bar::new(i as f64, *cnt as f64)
                                .width(0.6)
                                .fill(Color32::from_rgb(56, 189, 248))
                        })
                        .collect();
                    Plot::new("hl7_syndrome_plot")
                        .height(175.0)
                        .show(ui, |plot_ui| {
                            plot_ui.bar_chart(BarChart::new(bars).name("Incident Cases"));
                        });
                    ui.horizontal_wrapped(|ui| {
                        for (i, (lbl, cnt)) in out.syndrome_counts.iter().enumerate() {
                            ui.monospace(format!("[{i}] {lbl}: {cnt}  "));
                        }
                    });
                }
            });
        });

        if let Some(ref out) = self.cached_hl7 {
            ui.group(|ui| {
                ui.strong("Parsed HL7 v2.5.1 Clinical / Lab Line List (Zero-PII Bloom Fingerprints)");
                egui::ScrollArea::vertical().max_height(195.0).show(ui, |ui| {
                    egui::Grid::new("hl7_events_grid").striped(true).show(ui, |ui| {
                        ui.strong("Ctrl ID");
                        ui.strong("Type");
                        ui.strong("Date");
                        ui.strong("Facility");
                        ui.strong("Syndrome");
                        ui.strong("CLK-512 Prefix");
                        ui.strong("Entity");
                        ui.strong("Disposition");
                        ui.end_row();

                        for ev in &out.events {
                            ui.monospace(&ev.message_control_id);
                            ui.monospace(&ev.message_type);
                            ui.monospace(&ev.event_iso_date);
                            ui.monospace(&ev.sending_facility);
                            ui.monospace(ev.syndrome.label());
                            ui.monospace(&ev.bloom_fingerprint_prefix);
                            ui.monospace(if ev.linked_entity_id > 0 {
                                format!("P#{}", ev.linked_entity_id)
                            } else {
                                "--".to_string()
                            });
                            let disp = if !ev.is_positive_or_case {
                                RichText::new("Rule-Out / Neg").color(Color32::from_gray(150))
                            } else if ev.retained_after_pprl {
                                RichText::new("RETAINED CASE").color(Color32::from_rgb(52, 211, 153))
                            } else {
                                RichText::new("PPRL Duplicate").color(Color32::from_rgb(251, 191, 36))
                            };
                            ui.label(disp.monospace());
                            ui.end_row();
                        }
                    });
                });
            });
        }
    }

    fn render_pprl_view(&mut self, ui: &mut egui::Ui, status_banner: &mut String) {
        ui.heading(
            "512-Bit Keyed HMAC-SHA256 Unpadded Bigram Bloom Filters & Fellegi-Sunter EM Linkage",
        );
        ui.label(
            RichText::new(
                "Cross-institution patient deduplication without storing plaintext PII: unpadded character bigrams + POPCNT Dice similarity + 5-bin Fellegi-Sunter EM posterior weights.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        let sum = &self.cached_pprl;
        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("INPUT ENCOUNTERS").small().weak());
                ui.heading(RichText::new(format!("{}", sum.total_input_records)).monospace());
                ui.label(RichText::new("ADT + ELR line list").small().monospace());
            });
            cols[1].group(|ui| {
                ui.label(RichText::new("UNIQUE LINKED PATIENTS").small().weak());
                ui.heading(
                    RichText::new(format!("{}", sum.unique_linked_entities))
                        .monospace()
                        .color(Color32::from_rgb(56, 189, 248)),
                );
                ui.label(
                    RichText::new(format!(
                        "EM prior match rate: {:.1}%",
                        sum.em_prior_match_rate * 100.0
                    ))
                    .small()
                    .monospace(),
                );
            });
            cols[2].group(|ui| {
                ui.label(RichText::new("INCIDENT EPISODES").small().weak());
                ui.heading(
                    RichText::new(format!("{}", sum.incident_episodes_retained))
                        .monospace()
                        .color(Color32::from_rgb(52, 211, 153)),
                );
                ui.label(
                    RichText::new(format!(
                        "Window: {:.0} days",
                        self.pprl_cfg.episode_window_days
                    ))
                    .small()
                    .monospace(),
                );
            });
            cols[3].group(|ui| {
                ui.label(RichText::new("COLLAPSED DUPLICATES").small().weak());
                ui.heading(
                    RichText::new(format!("{}", sum.duplicate_encounters_collapsed))
                        .monospace()
                        .color(Color32::from_rgb(251, 191, 36)),
                );
                ui.label(RichText::new("Fuzzy typo-tolerant").small().monospace());
            });
        });

        ui.add_space(6.0);
        let mut pprl_dirty = false;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("HMAC Key:");
                if ui
                    .add(
                        egui::TextEdit::singleline(&mut self.pprl_cfg.secret_key)
                            .desired_width(210.0),
                    )
                    .changed()
                {
                    pprl_dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.pprl_cfg.dice_match_threshold, 0.55..=0.98)
                            .text("Dice Cutoff"),
                    )
                    .changed()
                {
                    pprl_dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.pprl_cfg.episode_window_days, 0.0..=120.0)
                            .text("Episode Window (d)"),
                    )
                    .changed()
                {
                    pprl_dirty = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut self.pprl_cfg.hardening_salted_bits, 0..=8)
                            .text("Hardening Bits"),
                    )
                    .changed()
                {
                    pprl_dirty = true;
                }
            });
        });

        if pprl_dirty {
            self.recompute_pprl_and_hl7();
            *status_banner = format!(
                "Re-encoded 512-bit HMAC-SHA256 Bloom filters: {} unique patients, {} incident episodes.",
                self.cached_pprl.unique_linked_entities,
                self.cached_pprl.incident_episodes_retained
            );
        }

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.strong("Pairwise 512-Bit Bloom Sørensen-Dice Similarity Matrix");
                let n = self.pprl_encoded.len();
                let desired = Vec2::new(ui.available_width().max(260.0), 220.0);
                let (rect, _) = ui.allocate_exact_size(desired, Sense::hover());
                let painter = ui.painter_at(rect);
                painter.rect_filled(rect, 4.0, Color32::from_rgb(15, 23, 42));

                if n > 0 {
                    let side = (rect.height() - 16.0).min(rect.width() - 16.0);
                    let x0 = rect.left() + 8.0;
                    let y0 = rect.top() + 8.0;
                    let cs = side / (n as f32);
                    for r in 0..n {
                        for c in 0..n {
                            let d = self.pprl_encoded[r]
                                .clk_bloom
                                .dice_similarity(&self.pprl_encoded[c].clk_bloom)
                                as f32;
                            let col = if r == c {
                                Color32::from_rgb(56, 189, 248)
                            } else if f64::from(d) >= self.pprl_cfg.dice_match_threshold {
                                Color32::from_rgb(52, 211, 153)
                            } else {
                                let v = (20.0 + 160.0 * d) as u8;
                                Color32::from_rgb(v / 2, v / 2, v)
                            };
                            let cell = Rect::from_min_max(
                                Pos2::new(x0 + (c as f32) * cs, y0 + (r as f32) * cs),
                                Pos2::new(
                                    x0 + ((c + 1) as f32) * cs - 1.0,
                                    y0 + ((r + 1) as f32) * cs - 1.0,
                                ),
                            );
                            painter.rect_filled(cell, 0.0, col);
                        }
                    }
                }
            });

            cols[1].group(|ui| {
                ui.strong("Top Candidate Pairs (Dice >= 0.45) & Fellegi-Sunter Posterior");
                egui::ScrollArea::vertical().max_height(220.0).show(ui, |ui| {
                    egui::Grid::new("pprl_pairs_grid").striped(true).show(ui, |ui| {
                        ui.strong("Rec A");
                        ui.strong("Rec B");
                        ui.strong("Dice");
                        ui.strong("Jaccard");
                        ui.strong("FS Bits");
                        ui.strong("P(Match)");
                        ui.strong("Link");
                        ui.end_row();

                        for cp in &self.cached_pprl.candidate_pairs {
                            ui.monospace(&cp.record_a);
                            ui.monospace(&cp.record_b);
                            ui.monospace(format!("{:.3}", cp.dice_similarity));
                            ui.monospace(format!("{:.3}", cp.jaccard_similarity));
                            ui.monospace(format!("{:+.2}", cp.fellegi_sunter_weight_bits));
                            ui.monospace(format!("{:.1}%", cp.posterior_match_prob * 100.0));
                            let lbl = if cp.is_linked_match {
                                RichText::new("MATCH").color(Color32::from_rgb(52, 211, 153))
                            } else {
                                RichText::new("Distinct").color(Color32::from_gray(150))
                            };
                            ui.label(lbl.monospace());
                            ui.end_row();
                        }
                    });
                });
            });
        });
    }

    fn render_sero_tnd_view(
        &mut self,
        ui: &mut egui::Ui,
        status_banner: &mut String,
        action: &mut Phase6GuiAction,
    ) {
        ui.columns(2, |cols| {
            // Left Column: Age-Structured Catalytic & Reverse-Catalytic FOI
            cols[0].group(|ui| {
                ui.heading("Age-Structured Sero-Catalytic Force of Infection (FOI)");
                ui.label(
                    RichText::new(
                        "Griffiths/Muench catalytic & reverse-catalytic Binomial MLE estimating annual FOI (lambda), seroreversion (omega), and SEIR initial immunity.",
                    )
                    .small()
                    .weak(),
                );

                let mut sero_dirty = false;
                ui.horizontal_wrapped(|ui| {
                    for k in [
                        CatalyticModelKind::IrreversibleCatalytic,
                        CatalyticModelKind::ReversibleWaningCatalytic,
                    ] {
                        if ui
                            .selectable_label(self.sero_kind == k, k.label())
                            .clicked()
                        {
                            self.sero_kind = k;
                            sero_dirty = true;
                        }
                    }
                });

                ui.horizontal_wrapped(|ui| {
                    ui.add(
                        egui::TextEdit::singleline(&mut self.sero_csv_path)
                            .desired_width(155.0)
                            .hint_text("serosurvey_strata.csv"),
                    );
                    if ui.button("Browse...").clicked() {
                        if let Some(p) = rfd::FileDialog::new()
                            .add_filter("CSV / TSV", &["csv", "tsv", "txt"])
                            .pick_file()
                        {
                            self.sero_csv_path = p.display().to_string();
                        }
                    }
                    if ui.button("Load Sero CSV").clicked() {
                        let p = self.sero_csv_path.trim().to_string();
                        match std::fs::read_to_string(&p)
                            .map_err(|e| e.to_string())
                            .and_then(|t| parse_serosurvey_csv(&t, &p).map_err(|e| e.to_string()))
                        {
                            Ok(ds) => {
                                self.sero_dataset = ds;
                                sero_dirty = true;
                            }
                            Err(e) => *status_banner = format!("Serosurvey CSV error: {e}"),
                        }
                    }
                    if ui.button("Reset Built-in").clicked() {
                        self.sero_dataset = builtin_age_serosurvey_dataset();
                        sero_dirty = true;
                    }
                });

                if sero_dirty {
                    if let Err(e) = self.recompute_sero() {
                        *status_banner = format!("Catalytic fit warning: {e}");
                    }
                }

                if let Some(ref fit) = self.cached_sero {
                    ui.separator();
                    ui.monospace(format!(
                        "Annual FOI lambda = {:.4}/yr [95% CI: {:.4}, {:.4}]  |  Implied R0 = {:.2}",
                        fit.foi_lambda_per_year,
                        fit.foi_ci_low_95,
                        fit.foi_ci_high_95,
                        fit.implied_endemic_r0
                    ));
                    ui.monospace(format!(
                        "Waning omega      = {:.4}/yr ({})  |  AIC = {:.1}",
                        fit.waning_omega_per_year,
                        fit.mean_antibody_duration_years
                            .map(|y| format!("mean {:.1} yrs", y))
                            .unwrap_or_else(|| "lifelong".to_string()),
                        fit.aic
                    ));
                    ui.monospace(format!(
                        "Weighted Immune   = {:.1}%  |  Empirical Susceptible S_eff = {:.1}%",
                        fit.weighted_seropositive_fraction * 100.0,
                        fit.effective_susceptible_fraction * 100.0
                    ));

                    if ui
                        .button("Push Empirical Immunity & Waning into SEIR Tab")
                        .clicked()
                    {
                        let wan_d = fit
                            .mean_antibody_duration_years
                            .map(|y| y * 365.25)
                            .unwrap_or(0.0);
                        action.apply_sero_to_seir =
                            Some((fit.weighted_seropositive_fraction, wan_d));
                    }

                    Plot::new("sero_catalytic_plot")
                        .height(215.0)
                        .legend(Legend::default())
                        .show(ui, |plot_ui| {
                            let mid: Vec<[f64; 2]> = fit
                                .age_curve
                                .iter()
                                .map(|p| [p.age_years, p.fitted_seroprevalence * 100.0])
                                .collect();
                            let lo: Vec<[f64; 2]> = fit
                                .age_curve
                                .iter()
                                .map(|p| [p.age_years, p.ci_low_95 * 100.0])
                                .collect();
                            let hi: Vec<[f64; 2]> = fit
                                .age_curve
                                .iter()
                                .map(|p| [p.age_years, p.ci_high_95 * 100.0])
                                .collect();
                            let obs: Vec<[f64; 2]> = self
                                .sero_dataset
                                .strata
                                .iter()
                                .map(|s| {
                                    [
                                        s.age_midpoint_years,
                                        (s.n_seropositive as f64 / s.n_tested.max(1) as f64) * 100.0,
                                    ]
                                })
                                .collect();

                            plot_ui.line(
                                Line::new(PlotPoints::from(mid))
                                    .name("Catalytic IgG %")
                                    .color(Color32::from_rgb(52, 211, 153))
                                    .width(2.2_f32),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(lo))
                                    .name("95% CI Low")
                                    .color(Color32::from_rgb(148, 163, 184))
                                    .width(1.0_f32),
                            );
                            plot_ui.line(
                                Line::new(PlotPoints::from(hi))
                                    .name("95% CI High")
                                    .color(Color32::from_rgb(148, 163, 184))
                                    .width(1.0_f32),
                            );
                            plot_ui.points(
                                Points::new(PlotPoints::from(obs))
                                    .name("Observed Serosurvey")
                                    .color(Color32::from_rgb(251, 191, 36))
                                    .radius(4.0_f32),
                            );
                        });
                }
            });

            // Right Column: Test-Negative Design (TND) Confounder-Adjusted VE
            cols[1].group(|ui| {
                ui.heading("Test-Negative Design (TND) Confounder-Adjusted VE");
                ui.label(
                    RichText::new(
                        "Multivariable logistic IRLS (VE = 1 - OR_adj) controlling for calendar epoch, senior age stratum, and healthcare-seeking propensity.",
                    )
                    .small()
                    .weak(),
                );

                if let Some(ref tnd) = self.cached_tnd {
                    ui.separator();
                    ui.monospace(format!(
                        "TND Study Cohort      : {} ({} stratified cells)",
                        self.tnd_dataset.title,
                        self.tnd_dataset.cells.len()
                    ));
                    ui.monospace(format!(
                        "Symptomatic Specimens : {} Test-Pos Cases  |  {} Test-Neg Controls",
                        tnd.total_cases, tnd.total_controls
                    ));
                    ui.monospace(format!(
                        "Crude Unadjusted VE   : {:>5.1}%  [95% CI: {:.1}%, {:.1}%]",
                        tnd.crude_unadjusted_ve_pct, tnd.crude_ve_ci_low_95, tnd.crude_ve_ci_high_95
                    ));
                    ui.monospace(format!(
                        "Adjusted Overall VE   : {:>5.1}%  [95% CI: {:.1}%, {:.1}%]  (Bias removed: {:+.1} pp)",
                        tnd.adjusted_overall_ve_pct,
                        tnd.adjusted_ve_ci_low_95,
                        tnd.adjusted_ve_ci_high_95,
                        tnd.confounding_bias_removed_pp
                    ));

                    let mut bars = Vec::new();
                    bars.push(
                        Bar::new(0.0, tnd.crude_unadjusted_ve_pct)
                            .width(0.55)
                            .fill(Color32::from_rgb(148, 163, 184)),
                    );
                    bars.push(
                        Bar::new(1.0, tnd.adjusted_overall_ve_pct)
                            .width(0.55)
                            .fill(Color32::from_rgb(52, 211, 153)),
                    );
                    for (idx, w) in tnd.waning_strata.iter().enumerate() {
                        bars.push(
                            Bar::new((2 + idx) as f64, w.adjusted_ve_pct)
                                .width(0.55)
                                .fill(Color32::from_rgb(56, 189, 248)),
                        );
                    }

                    Plot::new("tnd_ve_bar_plot")
                        .height(175.0)
                        .show(ui, |plot_ui| {
                            plot_ui.bar_chart(BarChart::new(bars).name("VE (%)"));
                        });

                    egui::Grid::new("tnd_waning_grid").striped(true).show(ui, |ui| {
                        ui.strong("Vaccination Stratum");
                        ui.strong("Vax Cases");
                        ui.strong("Vax Controls");
                        ui.strong("Adjusted VE [95% CI]");
                        ui.end_row();

                        for (idx, w) in tnd.waning_strata.iter().enumerate() {
                            ui.monospace(format!("[{}] {}", 2 + idx, w.label));
                            ui.monospace(format!("{}", w.cases_vax));
                            ui.monospace(format!("{}", w.controls_vax));
                            ui.monospace(format!(
                                "{:.1}% [{:.1}%, {:.1}%]",
                                w.adjusted_ve_pct, w.ve_ci_low_95, w.ve_ci_high_95
                            ));
                            ui.end_row();
                        }
                    });
                }
            });
        });
    }
}
