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
use epi_core::{parse_flexible_date_to_mmwr, CellValue, SuppressionPolicy, TimeSeriesDataset};
use epi_store::{
    encode_pprl_record, generate_synthetic_pprl_scale_benchmark, get_vintage_spent_epsilon,
    get_vintage_spent_rho, link_and_deduplicate_pprl, list_audit_events, list_dp_budget_entries,
    release_differentially_private_series, release_discrete_gaussian_zcdp_series,
    reset_vintage_dp_budget, verify_audit_chain, zcdp_rho_to_epsilon, AuditChainVerification,
    AuditLogEntry, DpBudgetLedgerEntry, DpReleaseConfig, DpReleaseOutcome, EncodedPprlRecord,
    EpiStore, PprlConfig, PprlLinkageSummary, ZcdpReleaseConfig, ZcdpReleaseOutcome,
};
use std::sync::Arc;
use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DpMechanismKind {
    DiscreteLaplace,
    DiscreteGaussianZcdp,
}

#[derive(Default)]
pub struct Phase8GuiAction {
    pub commit_dp_vintage: Option<TimeSeriesDataset>,
}

/// Build plot coordinates for true vintage counts while honoring `SuppressionPolicy`.
/// When `NchsWonderLe9` is active, counts in `1..=9` map to `[x, f64::NAN]` so suppressed
/// low-count cells are never leaked in plot curves or hover tooltips.
pub fn dp_overlay_true_points(
    dataset: &TimeSeriesDataset,
    policy: SuppressionPolicy,
) -> Vec<[f64; 2]> {
    dataset
        .points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            match CellValue::from_count(Some(p.count), p.suppressed, policy) {
                CellValue::Exact(v) => [i as f64, v],
                CellValue::Suppressed | CellValue::Missing => [i as f64, f64::NAN],
            }
        })
        .collect()
}

fn split_finite_segments(pts: &[[f64; 2]]) -> Vec<Vec<[f64; 2]>> {
    let mut segs: Vec<Vec<[f64; 2]>> = Vec::new();
    let mut cur: Vec<[f64; 2]> = Vec::new();
    for &p in pts {
        if p[1].is_finite() {
            cur.push(p);
        } else if !cur.is_empty() {
            segs.push(std::mem::take(&mut cur));
        }
    }
    if !cur.is_empty() {
        segs.push(cur);
    }
    segs
}

pub struct PprlScaleBenchReport {
    pub num_unique_patients: usize,
    pub total_encounters: usize,
    pub used_blocking: bool,
    pub elapsed_ms: f64,
    pub summary: PprlLinkageSummary,
}

pub struct DpAuditWorkbenchState {
    pub mechanism_kind: DpMechanismKind,
    pub dp_cfg: DpReleaseConfig,
    pub zcdp_cfg: ZcdpReleaseConfig,
    pub commit_on_release: bool,
    pub last_dp_outcome: Option<DpReleaseOutcome>,
    pub last_zcdp_outcome: Option<ZcdpReleaseOutcome>,
    pub spent_epsilon: f64,
    pub spent_rho: f64,
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
        let spent_rho = get_vintage_spent_rho(store, active_vintage_id).unwrap_or(0.0);
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

        let st = Self {
            mechanism_kind: DpMechanismKind::DiscreteLaplace,
            dp_cfg: DpReleaseConfig::default(),
            zcdp_cfg: ZcdpReleaseConfig::default(),
            commit_on_release: false,
            last_dp_outcome: None,
            last_zcdp_outcome: None,
            spent_epsilon,
            spent_rho,
            dp_history,
            chain_verification,
            audit_entries,
            bench_unique_patients: 600,
            bench_use_blocking: true,
            os_threads,
            last_bench: None,
        };
        st
    }

    pub fn refresh_from_store(&mut self, store: &EpiStore, active_vintage_id: &str) {
        self.spent_epsilon = get_vintage_spent_epsilon(store, active_vintage_id).unwrap_or(0.0);
        self.spent_rho = get_vintage_spent_rho(store, active_vintage_id).unwrap_or(0.0);
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
            self.last_zcdp_outcome = None;
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
        suppression_policy: SuppressionPolicy,
        status_banner: &mut String,
    ) -> Phase8GuiAction {
        let mut action = Phase8GuiAction::default();
        if self.last_bench.is_none() {
            self.run_pprl_benchmark();
        }
        self.spent_epsilon =
            get_vintage_spent_epsilon(store, &active_dataset.vintage_id).unwrap_or(0.0);
        self.spent_rho = get_vintage_spent_rho(store, &active_dataset.vintage_id).unwrap_or(0.0);

        ui.heading(
            "Discrete Laplace / Discrete Gaussian (zCDP) Differential Privacy Ledger & Audit Chain",
        );
        ui.label(
            RichText::new(
                "Two-sided geometric epsilon-DP & exact Discrete Gaussian rho-zCDP count protection with append-only SQLite privacy-budget ledger, SHA-256 audit verification, and LSH PPRL.",
            )
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        let rem_budget = (self.dp_cfg.max_budget_cap - self.spent_epsilon).max(0.0);
        let rem_rho = (self.zcdp_cfg.max_rho_cap - self.spent_rho).max(0.0);
        ui.columns(4, |cols| {
            cols[0].group(|ui| {
                match self.mechanism_kind {
                    DpMechanismKind::DiscreteLaplace => {
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
                            RichText::new(format!(
                                "Remaining: {:.2} | rho={:.3}",
                                rem_budget, self.spent_rho
                            ))
                            .small()
                            .monospace(),
                        );
                    }
                    DpMechanismKind::DiscreteGaussianZcdp => {
                        ui.label(RichText::new("ACTIVE VINTAGE zCDP RHO SPENT").small().weak());
                        let col = if rem_rho < self.zcdp_cfg.rho {
                            Color32::from_rgb(244, 63, 94)
                        } else {
                            Color32::from_rgb(56, 189, 248)
                        };
                        ui.heading(
                            RichText::new(format!(
                                "{:.3} / {:.3}",
                                self.spent_rho, self.zcdp_cfg.max_rho_cap
                            ))
                            .monospace()
                            .color(col),
                        );
                        ui.label(
                            RichText::new(format!(
                                "Rem rho={:.3} | eps={:.2}",
                                rem_rho, self.spent_epsilon
                            ))
                            .small()
                            .monospace(),
                        );
                    }
                }
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("LATEST DP RELEASE RMSE").small().weak());
                match self.mechanism_kind {
                    DpMechanismKind::DiscreteGaussianZcdp => {
                        if let Some(ref out) = self.last_zcdp_outcome {
                            ui.heading(
                                RichText::new(format!("{:.2} cases", out.empirical_rmse))
                                    .monospace()
                                    .color(Color32::from_rgb(251, 191, 36)),
                            );
                            ui.label(
                                RichText::new(format!(
                                    "sigma={:.2} | MAE {:.2}",
                                    out.sigma, out.empirical_mae
                                ))
                                .small()
                                .monospace(),
                            );
                        } else {
                            let sigma = self.zcdp_cfg.sensitivity_l2.max(0.01)
                                / (2.0 * self.zcdp_cfg.rho.max(1e-4)).sqrt();
                            let conv_eps =
                                zcdp_rho_to_epsilon(self.zcdp_cfg.rho, self.zcdp_cfg.delta);
                            ui.heading(RichText::new(format!("~{sigma:.2} est")).monospace());
                            ui.label(
                                RichText::new(format!("conv eps={conv_eps:.2} @ d=1e-5"))
                                    .small()
                                    .monospace(),
                            );
                        }
                    }
                    DpMechanismKind::DiscreteLaplace => {
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
                            let alpha =
                                (-self.dp_cfg.epsilon / self.dp_cfg.sensitivity_l1.max(0.01)).exp();
                            let th = (2.0 * alpha).sqrt() / (1.0 - alpha).max(1e-9);
                            ui.heading(RichText::new(format!("~{th:.2} est")).monospace());
                            ui.label(RichText::new("Analytical pre-release RMSE").small().monospace());
                        }
                    }
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
                ui.label("Mechanism:");
                ui.selectable_value(
                    &mut self.mechanism_kind,
                    DpMechanismKind::DiscreteLaplace,
                    "Discrete Laplace (eps-DP)",
                );
                ui.selectable_value(
                    &mut self.mechanism_kind,
                    DpMechanismKind::DiscreteGaussianZcdp,
                    "Discrete Gaussian (rho-zCDP)",
                );
                ui.separator();

                match self.mechanism_kind {
                    DpMechanismKind::DiscreteLaplace => {
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
                    }
                    DpMechanismKind::DiscreteGaussianZcdp => {
                        ui.add(
                            egui::Slider::new(&mut self.zcdp_cfg.rho, 0.005..=0.25)
                                .text("Query Rho (zCDP)"),
                        );
                        ui.add(
                            egui::Slider::new(&mut self.zcdp_cfg.sensitivity_l2, 1.0..=5.0)
                                .text("L2 Sensitivity"),
                        );
                        ui.add(
                            egui::Slider::new(&mut self.zcdp_cfg.max_rho_cap, 0.05..=2.0)
                                .text("Max Rho Cap"),
                        );
                    }
                }
                ui.checkbox(&mut self.commit_on_release, "Save as New SQLite Vintage");

                match self.mechanism_kind {
                    DpMechanismKind::DiscreteLaplace => {
                        if ui.button("Execute Discrete Laplace DP Release").clicked() {
                            match release_differentially_private_series(
                                store,
                                active_dataset,
                                &self.dp_cfg,
                            ) {
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
                                        action.commit_dp_vintage =
                                            Some(outcome.noisy_dataset.clone());
                                    }
                                    self.last_dp_outcome = Some(outcome);
                                    self.last_zcdp_outcome = None;
                                    self.refresh_from_store(store, &active_dataset.vintage_id);
                                }
                                Err(err) => {
                                    *status_banner =
                                        format!("DP Budget Guard Blocked Release: {err}");
                                }
                            }
                        }
                    }
                    DpMechanismKind::DiscreteGaussianZcdp => {
                        if ui.button("Execute Discrete Gaussian (zCDP) Release").clicked() {
                            match release_discrete_gaussian_zcdp_series(
                                store,
                                active_dataset,
                                &self.zcdp_cfg,
                            ) {
                                Ok(outcome) => {
                                    *status_banner = format!(
                                        "Released Discrete Gaussian zCDP series (rho={:.3}, sigma={:.2}, conv_eps={:.2}, MAE={:.2}, audit={}).",
                                        outcome.rho_spent,
                                        outcome.sigma,
                                        outcome.converted_epsilon_delta,
                                        outcome.empirical_mae,
                                        &outcome.audit_entry_hash[0..12]
                                    );
                                    if self.commit_on_release {
                                        action.commit_dp_vintage =
                                            Some(outcome.noisy_dataset.clone());
                                    }
                                    self.last_zcdp_outcome = Some(outcome);
                                    self.last_dp_outcome = None;
                                    self.refresh_from_store(store, &active_dataset.vintage_id);
                                }
                                Err(err) => {
                                    *status_banner =
                                        format!("zCDP Budget Guard Blocked Release: {err}");
                                }
                            }
                        }
                    }
                }

                if ui.button("Governance Reset Vintage Budget").clicked() {
                    match reset_vintage_dp_budget(store, &active_dataset.vintage_id) {
                        Ok(ev) => {
                            self.refresh_from_store(store, &active_dataset.vintage_id);
                            *status_banner = format!(
                                "Appended governance reset epoch for '{}' (immutable audit seq #{}).",
                                active_dataset.vintage_id, ev.seq_no
                            );
                        }
                        Err(e) => *status_banner = format!("Reset error: {e}"),
                    }
                }
            });
        });

        let x_labels: Arc<Vec<String>> = Arc::new(
            active_dataset
                .points
                .iter()
                .map(|p| {
                    p.epiweek
                        .map(|w| w.compact_label())
                        .or_else(|| {
                            parse_flexible_date_to_mmwr(&p.date_label).map(|w| w.compact_label())
                        })
                        .unwrap_or_else(|| p.date_label.clone())
                })
                .collect(),
        );
        let n_pts = x_labels.len();
        let x_labels_fmt = Arc::clone(&x_labels);
        let true_pts = dp_overlay_true_points(active_dataset, suppression_policy);
        let true_segs = split_finite_segments(&true_pts);
        let true_label = if suppression_policy == SuppressionPolicy::NchsWonderLe9 {
            "True Vintage Counts (NCHS <=9 Masked)"
        } else {
            "True Vintage Counts"
        };

        let mut y_max = true_pts
            .iter()
            .filter(|p| p[1].is_finite())
            .map(|p| p[1])
            .fold(0.0_f64, f64::max);
        let noisy_info: Option<(String, Vec<[f64; 2]>)> = if let Some(ref z) = self.last_zcdp_outcome {
            let pts: Vec<[f64; 2]> = z
                .noisy_dataset
                .points
                .iter()
                .enumerate()
                .map(|(i, p)| [i as f64, p.count])
                .collect();
            y_max = pts.iter().map(|p| p[1]).fold(y_max, f64::max);
            Some((
                format!("Discrete Gaussian (rho={:.3}, sigma={:.2})", z.rho_spent, z.sigma),
                pts,
            ))
        } else if let Some(ref out) = self.last_dp_outcome {
            let pts: Vec<[f64; 2]> = out
                .noisy_dataset
                .points
                .iter()
                .enumerate()
                .map(|(i, p)| [i as f64, p.count])
                .collect();
            y_max = pts.iter().map(|p| p[1]).fold(y_max, f64::max);
            Some((format!("Discrete Laplace (eps={:.2})", out.epsilon_spent), pts))
        } else {
            None
        };

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                let plot_title = match self.mechanism_kind {
                    DpMechanismKind::DiscreteLaplace => {
                        "True Surveillance Counts vs. Discrete Laplace DP Release"
                    }
                    DpMechanismKind::DiscreteGaussianZcdp => {
                        "True Surveillance Counts vs. Discrete Gaussian (zCDP) Release"
                    }
                };
                ui.strong(plot_title);
                let mut plot = Plot::new("dp_overlay_plot")
                    .height(225.0)
                    .legend(Legend::default().background_alpha(0.85))
                    .x_grid_spacer(move |input| {
                        let (min_x, max_x) = input.bounds;
                        let span = (max_x - min_x).abs().max(1.0);
                        let major_step: usize = if span > 180.0 {
                            26
                        } else if span > 80.0 {
                            13
                        } else if span > 32.0 {
                            4
                        } else if span > 14.0 {
                            2
                        } else {
                            1
                        };
                        let minor_step: usize = match major_step {
                            26 => 13,
                            13 => 4,
                            4 => 1,
                            _ => 1,
                        };
                        let start = (min_x.floor() as isize).max(0) as usize;
                        let end = ((max_x.ceil() as isize).max(0) as usize)
                            .min(n_pts.saturating_sub(1));
                        let mut marks = Vec::new();
                        for idx in start..=end {
                            if idx % major_step == 0 {
                                marks.push(egui_plot::GridMark {
                                    value: idx as f64,
                                    step_size: major_step as f64,
                                });
                            } else if minor_step < major_step && idx % minor_step == 0 {
                                marks.push(egui_plot::GridMark {
                                    value: idx as f64,
                                    step_size: minor_step as f64,
                                });
                            }
                        }
                        marks
                    })
                    .x_axis_formatter(move |mark, range| {
                        let span = (range.end() - range.start()).abs().max(1.0);
                        let major_step = if span > 180.0 {
                            26.0
                        } else if span > 80.0 {
                            13.0
                        } else if span > 32.0 {
                            4.0
                        } else if span > 14.0 {
                            2.0
                        } else {
                            1.0
                        };
                        if mark.step_size + 1e-6 < major_step {
                            return String::new();
                        }
                        let rounded = mark.value.round();
                        if (mark.value - rounded).abs() > 0.05 || rounded < 0.0 {
                            return String::new();
                        }
                        x_labels_fmt
                            .get(rounded as usize)
                            .cloned()
                            .unwrap_or_default()
                    })
                    .x_axis_label("MMWR Week (YYYY-Www)");
                if y_max > 0.0 {
                    plot = plot.include_y(y_max * 1.15);
                }
                plot.show(ui, |plot_ui| {
                    for seg in &true_segs {
                        plot_ui.line(
                            Line::new(PlotPoints::from(seg.clone()))
                                .name(true_label)
                                .color(Color32::from_rgb(148, 163, 184))
                                .width(1.6_f32),
                        );
                    }
                    if let Some((ref label, ref pts)) = noisy_info {
                        plot_ui.line(
                            Line::new(PlotPoints::from(pts.clone()))
                                .name(label)
                                .color(Color32::from_rgb(52, 211, 153))
                                .width(1.8_f32),
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
                        "{} Candidate Pairs      : {} pairs ({:.1}% pruned)",
                        if b.used_blocking { "Blocked" } else { "All" },
                        b.summary.evaluated_candidate_pairs,
                        b.summary.blocking_pair_reduction_pct
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
                ui.strong("SQLite Differential Privacy Budget Ledger (Append-Only)");
                egui::ScrollArea::vertical().max_height(190.0).show(ui, |ui| {
                    egui::Grid::new("dp_ledger_grid").striped(true).show(ui, |ui| {
                        ui.strong("ID");
                        ui.strong("Kind");
                        ui.strong("Mech");
                        ui.strong("Eps");
                        ui.strong("Rho");
                        ui.strong("Cumul");
                        ui.strong("MAE");
                        ui.strong("Hash");
                        ui.end_row();

                        for d in &self.dp_history {
                            let mech_short = match d.mechanism.as_str() {
                                "DISCRETE_GAUSSIAN_ZCDP" => "zCDP-DG",
                                "DISCRETE_LAPLACE" => "Laplace",
                                "GOVERNANCE_RESET" => "Reset",
                                other => other,
                            };
                            ui.monospace(format!("#{}", d.entry_id));
                            ui.monospace(&d.entry_kind);
                            ui.monospace(mech_short);
                            ui.monospace(format!("{:.2}", d.epsilon_spent));
                            ui.monospace(format!("{:.3}", d.zcdp_rho_spent));
                            ui.monospace(format!("{:.2}", d.cumulative_epsilon_after));
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

#[cfg(test)]
mod tests {
    use super::*;
    use epi_core::{SourceKind, TimeSeriesPoint};

    #[test]
    fn test_dp_overlay_true_points_nchs_le9_masking() {
        let ds = TimeSeriesDataset {
            vintage_id: "mask_test".to_string(),
            title: "Mask Test".to_string(),
            source_kind: SourceKind::LocalCsvTsv,
            source_url: "local://test".to_string(),
            ingested_at_utc: "2026-11-18T00:00:00Z".to_string(),
            default_si_preset: "measles".to_string(),
            checksum_hex: "deadbeef".to_string(),
            points: vec![
                TimeSeriesPoint {
                    date_label: "2026-W01".to_string(),
                    epiweek: None,
                    count: 0.0,
                    suppressed: false,
                },
                TimeSeriesPoint {
                    date_label: "2026-W02".to_string(),
                    epiweek: None,
                    count: 1.0,
                    suppressed: false,
                },
                TimeSeriesPoint {
                    date_label: "2026-W03".to_string(),
                    epiweek: None,
                    count: 9.0,
                    suppressed: false,
                },
                TimeSeriesPoint {
                    date_label: "2026-W04".to_string(),
                    epiweek: None,
                    count: 10.0,
                    suppressed: false,
                },
            ],
        };

        let masked = dp_overlay_true_points(&ds, SuppressionPolicy::NchsWonderLe9);
        assert_eq!(masked[0][1], 0.0, "0 must remain exact");
        assert!(masked[1][1].is_nan(), "1 must be masked (NaN)");
        assert!(masked[2][1].is_nan(), "9 must be masked (NaN)");
        assert_eq!(masked[3][1], 10.0, "10 must remain unmasked");

        let segs = split_finite_segments(&masked);
        assert_eq!(segs.len(), 2, "masking 1..=9 splits curve into 2 disjoint segments");
        assert_eq!(segs[0], vec![[0.0, 0.0]]);
        assert_eq!(segs[1], vec![[3.0, 10.0]]);

        let unmasked = dp_overlay_true_points(&ds, SuppressionPolicy::Unmasked);
        assert_eq!(unmasked.iter().map(|p| p[1]).collect::<Vec<_>>(), vec![0.0, 1.0, 9.0, 10.0]);
    }
}
