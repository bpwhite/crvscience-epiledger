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
use egui_plot::{Legend, Line, Plot, PlotPoints};
use epi_core::TimeSeriesDataset;
use epi_store::{
    encode_pprl_record, generate_synthetic_pprl_scale_benchmark, get_vintage_spent_epsilon,
    link_and_deduplicate_pprl, list_audit_events, list_dp_budget_entries,
    release_differentially_private_series, reset_vintage_dp_budget, verify_audit_chain,
    AuditChainVerification, AuditLogEntry, DpBudgetLedgerEntry, DpReleaseConfig, DpReleaseOutcome,
    EncodedPprlRecord, EpiStore, PprlConfig, PprlLinkageSummary,
};
use std::time::Instant;

#[derive(Default)]
pub struct Phase8GuiAction {
    pub commit_dp_vintage: Option<TimeSeriesDataset>,
}

pub struct PprlScaleBenchReport {
    pub num_unique_patients: usize,
    pub total_encounters: usize,
    pub used_blocking: bool,
    pub elapsed_ms: f64,
    pub summary: PprlLinkageSummary,
}

pub struct DpAuditWorkbenchState {
    pub dp_cfg: DpReleaseConfig,
    pub commit_on_release: bool,
    pub last_dp_outcome: Option<DpReleaseOutcome>,
    pub spent_epsilon: f64,
    pub dp_history: Vec<DpBudgetLedgerEntry>,

    pub chain_verification: AuditChainVerification,
    pub audit_entries: Vec<AuditLogEntry>,

    pub bench_unique_patients: usize,
    pub bench_use_blocking: bool,
    pub os_threads: usize,
    pub last_bench: Option<PprlScaleBenchReport>,
}

impl DpAuditWorkbenchState {
    pub fn new(store: &EpiStore, active_vintage_id: &str) -> Self {
        let spent_epsilon = get_vintage_spent_epsilon(store, active_vintage_id).unwrap_or(0.0);
        let dp_history = list_dp_budget_entries(store, 25).unwrap_or_default();
        let chain_verification = verify_audit_chain(store).unwrap_or(AuditChainVerification {
            is_valid: true,
            total_entries: 0,
            tip_hash_hex: String::new(),
            broken_at_seq: None,
            failure_reason: None,
        });
        let audit_entries = list_audit_events(store, 40).unwrap_or_default();
        let os_threads = std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(4);

        let mut st = Self {
            dp_cfg: DpReleaseConfig::default(),
            commit_on_release: false,
            last_dp_outcome: None,
            spent_epsilon,
            dp_history,
            chain_verification,
            audit_entries,
            bench_unique_patients: 600,
            bench_use_blocking: true,
            os_threads,
            last_bench: None,
        };
        st.run_pprl_benchmark();
        st
    }

    pub fn refresh_from_store(&mut self, store: &EpiStore, active_vintage_id: &str) {
        self.spent_epsilon = get_vintage_spent_epsilon(store, active_vintage_id).unwrap_or(0.0);
        self.dp_history = list_dp_budget_entries(store, 25).unwrap_or_default();
        if let Ok(ver) = verify_audit_chain(store) {
            self.chain_verification = ver;
        }
        self.audit_entries = list_audit_events(store, 40).unwrap_or_default();
    }

    pub fn trigger_dp_release_demo(
        &mut self,
        store: &EpiStore,
        active_dataset: &TimeSeriesDataset,
        status_banner: &mut String,
    ) {
        if let Ok(outcome) =
            release_differentially_private_series(store, active_dataset, &self.dp_cfg)
        {
            *status_banner = format!(
                "Released DP series (eps={:.2}, cum={:.2}/{:.2}, empirical MAE={:.2}, audit={}).",
                outcome.epsilon_spent,
                outcome.cumulative_epsilon_for_vintage,
                self.dp_cfg.max_budget_cap,
                outcome.empirical_mae,
                &outcome.audit_entry_hash[0..12]
            );
            self.last_dp_outcome = Some(outcome);
            self.refresh_from_store(store, &active_dataset.vintage_id);
        }
    }

    pub fn run_pprl_benchmark(&mut self) {
        let t0 = Instant::now();
        let raw = generate_synthetic_pprl_scale_benchmark(self.bench_unique_patients);
        let mut cfg = PprlConfig::default();
        cfg.use_hash_blocking = self.bench_use_blocking;
        let enc: Vec<EncodedPprlRecord> = raw.iter().map(|r| encode_pprl_record(r, &cfg)).collect();
        let summary = link_and_deduplicate_pprl(&enc, &cfg);
        let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;

        self.last_bench = Some(PprlScaleBenchReport {
            num_unique_patients: self.bench_unique_patients,
            total_encounters: raw.len(),
            used_blocking: self.bench_use_blocking,
            elapsed_ms,
            summary,
        });
    }

    pub fn render(
        &mut self,
        ui: &mut egui::Ui,
        store: &EpiStore,
        active_dataset: &TimeSeriesDataset,
        status_banner: &mut String,
    ) -> Phase8GuiAction {
        let mut action = Phase8GuiAction::default();
        self.spent_epsilon =
            get_vintage_spent_epsilon(store, &active_dataset.vintage_id).unwrap_or(0.0);

        ui.heading(
            "Discrete Laplace Differential Privacy Ledger, SHA-256 Audit Chain & Multi-Core Blocking",
        );
        ui.label(
            RichText::new(
                "Two-sided geometric epsilon-DP count protection with SQLite privacy-budget enforcement, genesis-to-tip SHA-256 audit verification, and LSH hash-blocked PPRL.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        let rem_budget = (self.dp_cfg.max_budget_cap - self.spent_epsilon).max(0.0);
        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("ACTIVE VINTAGE EPSILON SPENT").small().weak());
                let col = if rem_budget < self.dp_cfg.epsilon {
                    Color32::from_rgb(244, 63, 94)
                } else {
                    Color32::from_rgb(56, 189, 248)
                };
                ui.heading(
                    RichText::new(format!(
                        "{:.2} / {:.2}",
                        self.spent_epsilon, self.dp_cfg.max_budget_cap
                    ))
                    .monospace()
                    .color(col),
                );
                ui.label(
                    RichText::new(format!("Remaining budget: {:.2}", rem_budget))
                        .small()
                        .monospace(),
                );
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("LATEST DP RELEASE RMSE").small().weak());
                if let Some(ref out) = self.last_dp_outcome {
                    ui.heading(
                        RichText::new(format!("{:.2} cases", out.empirical_rmse))
                            .monospace()
                            .color(Color32::from_rgb(251, 191, 36)),
                    );
                    ui.label(
                        RichText::new(format!(
                            "Theory {:.2} | MAE {:.2}",
                            out.theoretical_rmse, out.empirical_mae
                        ))
                        .small()
                        .monospace(),
                    );
                } else {
                    let alpha = (-self.dp_cfg.epsilon / self.dp_cfg.sensitivity_l1.max(0.01)).exp();
                    let th = (2.0 * alpha).sqrt() / (1.0 - alpha).max(1e-9);
                    ui.heading(RichText::new(format!("~{th:.2} est")).monospace());
                    ui.label(RichText::new("Analytical pre-release RMSE").small().monospace());
                }
            });

            cols[2].group(|ui| {
                ui.label(RichText::new("SHA-256 AUDIT HASH CHAIN").small().weak());
                let (txt, col) = if self.chain_verification.is_valid {
                    (
                        format!("VALID ({} links)", self.chain_verification.total_entries),
                        Color32::from_rgb(52, 211, 153),
                    )
                } else {
                    ("TAMPER DETECTED".to_string(), Color32::from_rgb(244, 63, 94))
                };
                ui.heading(RichText::new(txt).monospace().color(col));
                let tip_short = if self.chain_verification.tip_hash_hex.len() >= 14 {
                    &self.chain_verification.tip_hash_hex[0..14]
                } else {
                    "--"
                };
                ui.label(RichText::new(format!("Tip: {tip_short}...")).small().monospace());
            });

            cols[3].group(|ui| {
                ui.label(RichText::new("SCOPED THREADS & BLOCKING").small().weak());
                ui.heading(
                    RichText::new(format!("{} OS threads", self.os_threads))
                        .monospace()
                        .color(Color32::from_rgb(56, 189, 248)),
                );
                if let Some(ref b) = self.last_bench {
                    ui.label(
                        RichText::new(format!(
                            "PPRL pruned {:.1}% ({:.1} ms)",
                            b.summary.blocking_pair_reduction_pct, b.elapsed_ms
                        ))
                        .small()
                        .monospace(),
                    );
                }
            });
        });

        ui.add_space(6.0);
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.add(
                    egui::Slider::new(&mut self.dp_cfg.epsilon, 0.05..=2.50)
                        .text("Query Epsilon"),
                );
                ui.add(
                    egui::Slider::new(&mut self.dp_cfg.sensitivity_l1, 1.0..=5.0)
                        .text("L1 Sensitivity"),
                );
                ui.add(
                    egui::Slider::new(&mut self.dp_cfg.max_budget_cap, 0.50..=10.0)
                        .text("Max Budget Cap"),
                );
                ui.checkbox(&mut self.commit_on_release, "Save as New SQLite Vintage");

                if ui.button("Execute Discrete Laplace DP Release").clicked() {
                    match release_differentially_private_series(store, active_dataset, &self.dp_cfg)
                    {
                        Ok(outcome) => {
                            *status_banner = format!(
                                "Released DP series (eps={:.2}, cum={:.2}/{:.2}, empirical MAE={:.2}, audit={}).",
                                outcome.epsilon_spent,
                                outcome.cumulative_epsilon_for_vintage,
                                self.dp_cfg.max_budget_cap,
                                outcome.empirical_mae,
                                &outcome.audit_entry_hash[0..12]
                            );
                            if self.commit_on_release {
                                action.commit_dp_vintage = Some(outcome.noisy_dataset.clone());
                            }
                            self.last_dp_outcome = Some(outcome);
                            self.refresh_from_store(store, &active_dataset.vintage_id);
                        }
                        Err(err) => {
                            *status_banner = format!("DP Budget Guard Blocked Release: {err}");
                        }
                    }
                }

                if ui.button("Governance Reset Vintage Budget").clicked() {
                    match reset_vintage_dp_budget(store, &active_dataset.vintage_id) {
                        Ok(ev) => {
                            self.refresh_from_store(store, &active_dataset.vintage_id);
                            *status_banner = format!(
                                "Cleared DP budget for '{}' (logged immutable audit seq #{}).",
                                active_dataset.vintage_id, ev.seq_no
                            );
                        }
                        Err(e) => *status_banner = format!("Reset error: {e}"),
                    }
                }
            });
        });

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.strong("True Surveillance Counts vs. Discrete Laplace DP Release");
                Plot::new("dp_overlay_plot")
                    .height(205.0)
                    .legend(Legend::default())
                    .show(ui, |plot_ui| {
                        let true_pts: Vec<[f64; 2]> = active_dataset
                            .points
                            .iter()
                            .enumerate()
                            .map(|(i, p)| [i as f64, p.count])
                            .collect();
                        plot_ui.line(
                            Line::new(PlotPoints::from(true_pts))
                                .name("True Vintage Counts")
                                .color(Color32::from_rgb(148, 163, 184))
                                .width(1.6),
                        );
                        if let Some(ref out) = self.last_dp_outcome {
                            let dp_pts: Vec<[f64; 2]> = out
                                .noisy_dataset
                                .points
                                .iter()
                                .enumerate()
                                .map(|(i, p)| [i as f64, p.count])
                                .collect();
                            plot_ui.line(
                                Line::new(PlotPoints::from(dp_pts))
                                    .name(format!("Discrete Laplace (eps={:.2})", out.epsilon_spent))
                                    .color(Color32::from_rgb(52, 211, 153))
                                    .width(1.8),
                            );
                        }
                    });
            });

            cols[1].group(|ui| {
                ui.strong("Sub-Quadratic Hash-Blocked PPRL Scale Benchmark");
                ui.horizontal_wrapped(|ui| {
                    ui.add(
                        egui::Slider::new(&mut self.bench_unique_patients, 100..=2_500)
                            .text("Unique Patients"),
                    );
                    ui.checkbox(&mut self.bench_use_blocking, "LSH + Birth-Year Hash Blocking");
                    if ui.button("Run Scale Benchmark").clicked() {
                        self.run_pprl_benchmark();
                    }
                });
                if let Some(ref b) = self.last_bench {
                    ui.separator();
                    ui.monospace(format!(
                        "Input Encounters (ADT + ELR) : {} records",
                        b.total_encounters
                    ));
                    ui.monospace(format!(
                        "Cartesian All-Pairs O(N^2)   : {} pairs",
                        b.summary.total_cartesian_pairs
                    ));
                    ui.monospace(format!(
                        "Blocked Candidate Pairs      : {} pairs ({:.1}% pruned)",
                        b.summary.evaluated_candidate_pairs, b.summary.blocking_pair_reduction_pct
                    ));
                    ui.monospace(format!(
                        "Recovered Unique Entities    : {} (expected {})",
                        b.summary.unique_linked_entities, b.num_unique_patients
                    ));
                    ui.monospace(format!(
                        "Wall-Clock Execution Latency : {:.2} ms",
                        b.elapsed_ms
                    ));
                }
            });
        });

        ui.add_space(6.0);
        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.horizontal(|ui| {
                    ui.strong("Tamper-Evident SHA-256 Audit Chain");
                    if ui.button("Re-Verify Genesis->Tip").clicked() {
                        self.refresh_from_store(store, &active_dataset.vintage_id);
                        *status_banner = format!(
                            "Audit hash chain verified: valid={}, {} entries.",
                            self.chain_verification.is_valid,
                            self.chain_verification.total_entries
                        );
                    }
                });
                egui::ScrollArea::vertical().max_height(190.0).show(ui, |ui| {
                    egui::Grid::new("audit_chain_grid").striped(true).show(ui, |ui| {
                        ui.strong("Seq");
                        ui.strong("Action");
                        ui.strong("Target");
                        ui.strong("Prev SHA256");
                        ui.strong("Entry SHA256");
                        ui.end_row();

                        for e in &self.audit_entries {
                            ui.monospace(format!("#{}", e.seq_no));
                            ui.monospace(&e.action_kind);
                            ui.monospace(&e.target_ref);
                            ui.monospace(&e.prev_hash_hex[0..10]);
                            ui.monospace(&e.entry_hash_hex[0..10]);
                            ui.end_row();
                        }
                    });
                });
            });

            cols[1].group(|ui| {
                ui.strong("SQLite Differential Privacy Budget Ledger");
                egui::ScrollArea::vertical().max_height(190.0).show(ui, |ui| {
                    egui::Grid::new("dp_ledger_grid").striped(true).show(ui, |ui| {
                        ui.strong("ID");
                        ui.strong("Vintage");
                        ui.strong("Eps");
                        ui.strong("Cumul");
                        ui.strong("Cap");
                        ui.strong("MAE");
                        ui.strong("Hash");
                        ui.end_row();

                        for d in &self.dp_history {
                            ui.monospace(format!("#{}", d.entry_id));
                            ui.monospace(&d.vintage_id);
                            ui.monospace(format!("{:.2}", d.epsilon_spent));
                            ui.monospace(format!("{:.2}", d.cumulative_epsilon_after));
                            ui.monospace(format!("{:.1}", d.max_budget_cap));
                            ui.monospace(format!("{:.2}", d.empirical_mae));
                            ui.monospace(&d.output_checksum_hex);
                            ui.end_row();
                        }
                    });
                });
            });
        });

        action
    }
}
