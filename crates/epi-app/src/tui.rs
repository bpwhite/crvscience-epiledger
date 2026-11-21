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
use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use crossterm::execute;
use crossterm::terminal::{
    disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen,
};
use epi_analytics::{
    analyze_contact_tracing_network, build_us_state_scan_zones,
    builtin_emerging_wave_cases_and_deaths, builtin_hospital_amr_isolates,
    builtin_superspreading_outbreak_records, builtin_three_source_tb_pertussis_crc,
    builtin_wastewater_barcode_and_samples, calibrate_and_project_seir, compute_aberrations,
    compute_clsi_m39_antibiogram, compute_cori_rt, compute_delay_adjusted_cfr,
    deconvolve_wastewater_timeseries, fit_three_source_loglinear_crc, nowcast_parametric_delay,
    project_hospital_surge_from_seir, run_kulldorff_circular_scan,
    run_rolling_wis_backtest_and_ensemble, AberrationConfig, AberrationSeries, AntibiogramConfig,
    BacktestConfig, BipartiteNetworkAnalysis, CoriRtConfig, CoriRtSeries,
    CumulativeAntibiogramResult, DelayAdjustedCfrResult, DelayCfrConfig, EnsembleBacktestResult,
    HospitalSurgeConfig, HospitalSurgeResult, LongitudinalWastewaterResult, MultiZoneCluster,
    NowcastSeries, OffspringFitConfig, ParametricNowcastConfig, SeirConfig, SeirForecastResult,
    ThreeSourceCrcAnalysis, WastewaterDemixConfig,
};
use epi_core::{
    builtin_measles_weekly_dataset, builtin_pertussis_spatial_records, CellValue, EpiError,
    SuppressionPolicy, TimeSeriesDataset, VintageMetadata,
};
use epi_ingest::{
    builtin_ed_chief_complaint_encounters, classify_ed_encounters_essence,
    EssenceBatchClassification,
};
use epi_store::{
    get_vintage_spent_epsilon, list_audit_events, release_differentially_private_series,
    verify_audit_chain, AuditChainVerification, AuditLogEntry, DpReleaseConfig, EpiStore,
};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::symbols;
use ratatui::text::{Line, Span};
use ratatui::widgets::{
    Axis, Block, Borders, Cell, Chart, Dataset, GraphType, Paragraph, Row, Table, Tabs,
};
use ratatui::{Frame, Terminal};
use std::io::{self, Stdout};
use std::path::Path;
use std::time::Duration;

const TUI_TAB_TITLES: [&str; 8] = [
    "1:Rt & EARS",
    "2:SEIR & ICU",
    "3:Spatial Scan",
    "4:Superspread",
    "5:Sewage & AMR",
    "6:CRC & ESSENCE",
    "7:WIS Stack",
    "8:DP & Audit",
];

struct TerminalRestoreGuard;

impl Drop for TerminalRestoreGuard {
    fn drop(&mut self) {
        let _ = disable_raw_mode();
        let _ = execute!(io::stdout(), LeaveAlternateScreen);
    }
}

pub struct EpiTuiState {
    store: EpiStore,
    vintages: Vec<VintageMetadata>,
    selected_vintage_idx: usize,
    active_dataset: TimeSeriesDataset,
    suppression_policy: SuppressionPolicy,
    active_tab: usize,
    status_line: String,

    rt_cfg: CoriRtConfig,
    ab_cfg: AberrationConfig,
    nowcast_cfg: ParametricNowcastConfig,
    seir_cfg: SeirConfig,
    surge_cfg: HospitalSurgeConfig,

    cached_rt: CoriRtSeries,
    cached_ab: AberrationSeries,
    cached_nowcast: NowcastSeries,
    cached_seir: SeirForecastResult,
    cached_cfr: DelayAdjustedCfrResult,
    cached_surge: HospitalSurgeResult,
    cached_clusters: Vec<MultiZoneCluster>,
    cached_net: BipartiteNetworkAnalysis,
    cached_ww: LongitudinalWastewaterResult,
    cached_amr: CumulativeAntibiogramResult,
    cached_crc: ThreeSourceCrcAnalysis,
    cached_essence: EssenceBatchClassification,
    cached_wis: EnsembleBacktestResult,

    spent_epsilon: f64,
    audit_ver: AuditChainVerification,
    audit_logs: Vec<AuditLogEntry>,
}

impl EpiTuiState {
    pub fn new(db_path: &Path) -> Result<Self, EpiError> {
        let mut store = EpiStore::open_file(db_path).or_else(|_| EpiStore::open_in_memory())?;
        let seed_ds = builtin_measles_weekly_dataset();
        if store.list_vintages()?.is_empty() {
            store.commit_vintage(&seed_ds)?;
        }
        let vintages = store.list_vintages()?;
        let active_dataset = store
            .load_vintage(&seed_ds.vintage_id)?
            .unwrap_or(seed_ds);

        let rt_cfg = CoriRtConfig::default();
        let ab_cfg = AberrationConfig::default();
        let nowcast_cfg = ParametricNowcastConfig::default();
        let seir_cfg = SeirConfig::default();
        let surge_cfg = HospitalSurgeConfig::default();

        let counts: Vec<f64> = active_dataset.points.iter().map(|p| p.count).collect();
        let cached_rt = compute_cori_rt(&counts, &rt_cfg);
        let cached_ab = compute_aberrations(&counts, &ab_cfg);
        let cached_nowcast = nowcast_parametric_delay(&counts, &nowcast_cfg);
        let cached_seir = calibrate_and_project_seir(&counts, true, Some(&nowcast_cfg), &seir_cfg);

        let wave = builtin_emerging_wave_cases_and_deaths();
        let cached_cfr = compute_delay_adjusted_cfr(&wave[..26], &DelayCfrConfig::default())?;
        let cached_surge = project_hospital_surge_from_seir(&cached_seir, &surge_cfg)?;

        let jurisdictions = builtin_pertussis_spatial_records();
        let zones = build_us_state_scan_zones(&jurisdictions, true);
        let cached_clusters = run_kulldorff_circular_scan(&zones, 0.25, 499, 8);

        let net_recs = builtin_superspreading_outbreak_records();
        let cached_net = analyze_contact_tracing_network(&net_recs, &OffspringFitConfig::default())?;

        let (ww_lib, ww_smps) = builtin_wastewater_barcode_and_samples();
        let cached_ww =
            deconvolve_wastewater_timeseries(&ww_lib, &ww_smps, &WastewaterDemixConfig::default())?;

        let amr_recs = builtin_hospital_amr_isolates();
        let cached_amr = compute_clsi_m39_antibiogram(&amr_recs, &AntibiogramConfig::default())?;

        let crc_cnts = builtin_three_source_tb_pertussis_crc();
        let cached_crc = fit_three_source_loglinear_crc(&crc_cnts)?;

        let ed_encs = builtin_ed_chief_complaint_encounters();
        let cached_essence = classify_ed_encounters_essence(
            &ed_encs,
            "(^fever^ OR ^cough^ OR ^j10^) ANDNOT (^immunization^ OR ^heat pack^)",
            5,
        );

        let cached_wis = run_rolling_wis_backtest_and_ensemble(
            &counts,
            &BacktestConfig::default(),
            &seir_cfg,
        )?;

        let spent_epsilon = get_vintage_spent_epsilon(&store, &active_dataset.vintage_id)?;
        let audit_ver = verify_audit_chain(&store)?;
        let audit_logs = list_audit_events(&store, 14)?;

        Ok(Self {
            store,
            vintages,
            selected_vintage_idx: 0,
            active_dataset,
            suppression_policy: SuppressionPolicy::Unmasked,
            active_tab: 0,
            status_line: "Keys: [Tab/1-8] Switch View | [v] Cycle Vintage | [m] NCHS <=9 Mask | [d] DP Release | [s] Export SitRep | [q] Quit".to_string(),
            rt_cfg,
            ab_cfg,
            nowcast_cfg,
            seir_cfg,
            surge_cfg,
            cached_rt,
            cached_ab,
            cached_nowcast,
            cached_seir,
            cached_cfr,
            cached_surge,
            cached_clusters,
            cached_net,
            cached_ww,
            cached_amr,
            cached_crc,
            cached_essence,
            cached_wis,
            spent_epsilon,
            audit_ver,
            audit_logs,
        })
    }

    fn reload_active_vintage(&mut self) {
        if let Ok(list) = self.store.list_vintages() {
            self.vintages = list;
        }
        if self.vintages.is_empty() {
            return;
        }
        self.selected_vintage_idx %= self.vintages.len();
        let vid = self.vintages[self.selected_vintage_idx].vintage_id.clone();
        if let Ok(Some(ds)) = self.store.load_vintage(&vid) {
            self.active_dataset = ds;
            let counts: Vec<f64> = self.active_dataset.points.iter().map(|p| p.count).collect();
            self.cached_rt = compute_cori_rt(&counts, &self.rt_cfg);
            self.cached_ab = compute_aberrations(&counts, &self.ab_cfg);
            self.cached_nowcast = nowcast_parametric_delay(&counts, &self.nowcast_cfg);
            self.cached_seir =
                calibrate_and_project_seir(&counts, true, Some(&self.nowcast_cfg), &self.seir_cfg);
            if let Ok(sg) = project_hospital_surge_from_seir(&self.cached_seir, &self.surge_cfg) {
                self.cached_surge = sg;
            }
            if let Ok(wis) = run_rolling_wis_backtest_and_ensemble(
                &counts,
                &BacktestConfig::default(),
                &self.seir_cfg,
            ) {
                self.cached_wis = wis;
            }
        }
        self.refresh_governance();
    }

    fn refresh_governance(&mut self) {
        self.spent_epsilon =
            get_vintage_spent_epsilon(&self.store, &self.active_dataset.vintage_id).unwrap_or(0.0);
        if let Ok(ver) = verify_audit_chain(&self.store) {
            self.audit_ver = ver;
        }
        if let Ok(logs) = list_audit_events(&self.store, 14) {
            self.audit_logs = logs;
        }
    }

    fn trigger_dp_release(&mut self) {
        let mut cfg = DpReleaseConfig::default();
        cfg.epsilon = 0.50;
        match release_differentially_private_series(&self.store, &self.active_dataset, &cfg) {
            Ok(out) => {
                if let Err(e) = self.store.commit_vintage(&out.noisy_dataset) {
                    self.status_line = format!("DP vintage commit failed: {e}");
                } else {
                    self.refresh_governance();
                    if let Ok(list) = self.store.list_vintages() {
                        self.vintages = list;
                    }
                    self.status_line = format!(
                        "Committed Discrete Laplace DP vintage '{}' (eps=0.50, empirical RMSE={:.2}, spent={:.2}).",
                        out.noisy_dataset.vintage_id, out.empirical_rmse, out.cumulative_epsilon_for_vintage
                    );
                }
            }
            Err(e) => {
                self.status_line = format!("DP budget guardrail blocked release: {e}");
            }
        }
    }

    fn export_sitrep(&mut self) {
        let out_path = "epiledger_tui_sitrep.html";
        let ctx = SitRepContext {
            dataset: &self.active_dataset,
            suppression_policy: self.suppression_policy,
            rt_series: &self.cached_rt,
            ab_series: &self.cached_ab,
            nowcast: &self.cached_nowcast,
            backtest: &self.cached_wis,
            spatial_clusters: &self.cached_clusters,
        };
        let html = render_standalone_sitrep_html(&ctx);
        match std::fs::write(out_path, html) {
            Ok(()) => {
                self.status_line =
                    format!("Exported standalone zero-JS HTML5/SVG SitRep to '{out_path}'.");
            }
            Err(e) => {
                self.status_line = format!("SitRep export error: {e}");
            }
        }
    }
}

/// Launch the full-screen headless SSH Terminal User Interface (`epiledger tui`).
pub fn run_tui(db_path: &Path) -> Result<(), EpiError> {
    let mut app = EpiTuiState::new(db_path)?;

    let orig_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |panic_info| {
        let _ = disable_raw_mode();
        let _ = execute!(io::stdout(), LeaveAlternateScreen);
        orig_hook(panic_info);
    }));

    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)?;
    let _guard = TerminalRestoreGuard;

    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;
    terminal.clear()?;

    let loop_res = run_event_loop(&mut terminal, &mut app);
    let _ = terminal.show_cursor();
    loop_res
}

fn run_event_loop(
    terminal: &mut Terminal<CrosstermBackend<Stdout>>,
    app: &mut EpiTuiState,
) -> Result<(), EpiError> {
    loop {
        terminal.draw(|frame| draw_ui(frame, app))?;

        if event::poll(Duration::from_millis(200))? {
            if let Event::Key(key) = event::read()? {
                if key.kind != KeyEventKind::Press {
                    continue;
                }
                match key.code {
                    KeyCode::Char('q') | KeyCode::Esc => break,
                    KeyCode::Tab | KeyCode::Right => {
                        app.active_tab = (app.active_tab + 1) % TUI_TAB_TITLES.len();
                    }
                    KeyCode::BackTab | KeyCode::Left => {
                        app.active_tab =
                            (app.active_tab + TUI_TAB_TITLES.len() - 1) % TUI_TAB_TITLES.len();
                    }
                    KeyCode::Char(c @ '1'..='8') => {
                        app.active_tab = (c as usize) - ('1' as usize);
                    }
                    KeyCode::Char('v') => {
                        if !app.vintages.is_empty() {
                            app.selected_vintage_idx =
                                (app.selected_vintage_idx + 1) % app.vintages.len();
                            app.reload_active_vintage();
                            app.status_line = format!(
                                "Switched active SQLite vintage to '{}'.",
                                app.active_dataset.vintage_id
                            );
                        }
                    }
                    KeyCode::Char('m') => {
                        app.suppression_policy = match app.suppression_policy {
                            SuppressionPolicy::Unmasked => SuppressionPolicy::NchsWonderLe9,
                            SuppressionPolicy::NchsWonderLe9 => SuppressionPolicy::Unmasked,
                        };
                        app.status_line = format!(
                            "Toggled small-cell disclosure policy to {:?}.",
                            app.suppression_policy
                        );
                    }
                    KeyCode::Char('d') => app.trigger_dp_release(),
                    KeyCode::Char('s') => app.export_sitrep(),
                    KeyCode::Char('r') => {
                        app.reload_active_vintage();
                        app.status_line =
                            "Recomputed analytics and re-verified SHA-256 audit chain.".to_string();
                    }
                    _ => {}
                }
            }
        }
    }
    Ok(())
}

fn draw_ui(f: &mut Frame<'_>, app: &EpiTuiState) {
    let root_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Length(3),
            Constraint::Min(14),
            Constraint::Length(2),
        ])
        .split(f.area());

    // 1. Header Banner
    let mask_lbl = if app.suppression_policy == SuppressionPolicy::NchsWonderLe9 {
        "NCHS <=9 MASKED"
    } else {
        "UNMASKED"
    };
    let chain_lbl = if app.audit_ver.is_valid {
        format!("AUDIT CHAIN VALID ({} links)", app.audit_ver.total_entries)
    } else {
        "AUDIT CHAIN TAMPERED".to_string()
    };
    let hdr = Paragraph::new(Line::from(vec![
        Span::styled(
            " EpiLedger SSH TUI ",
            Style::default()
                .fg(Color::Black)
                .bg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::raw("  Vintage: "),
        Span::styled(
            &app.active_dataset.vintage_id,
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        ),
        Span::raw(format!(" ({} pts)  |  [{}]  |  ", app.active_dataset.points.len(), mask_lbl)),
        Span::styled(
            chain_lbl,
            Style::default().fg(if app.audit_ver.is_valid {
                Color::Green
            } else {
                Color::Red
            }),
        ),
    ]))
    .block(Block::default().borders(Borders::ALL).title("Headless Public Health Workstation"));
    f.render_widget(hdr, root_chunks[0]);

    // 2. Navigation Tabs
    let tab_lines: Vec<Line<'_>> = TUI_TAB_TITLES.iter().map(|&t| Line::from(t)).collect();
    let tabs = Tabs::new(tab_lines)
        .select(app.active_tab)
        .block(Block::default().borders(Borders::ALL).title("Surveillance Modules"))
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD | Modifier::UNDERLINED),
        );
    f.render_widget(tabs, root_chunks[1]);

    // 3. Active Module Body
    match app.active_tab {
        0 => draw_rt_farrington_tab(f, app, root_chunks[2]),
        1 => draw_seir_surge_tab(f, app, root_chunks[2]),
        2 => draw_spatial_scan_tab(f, app, root_chunks[2]),
        3 => draw_superspreading_tab(f, app, root_chunks[2]),
        4 => draw_wastewater_amr_tab(f, app, root_chunks[2]),
        5 => draw_crc_essence_tab(f, app, root_chunks[2]),
        6 => draw_wis_ensemble_tab(f, app, root_chunks[2]),
        _ => draw_dp_audit_tab(f, app, root_chunks[2]),
    }

    // 4. Footer Keybindings & Status Bar
    let footer = Paragraph::new(Line::from(vec![
        Span::styled(" STATUS: ", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        Span::raw(&app.status_line),
    ]));
    f.render_widget(footer, root_chunks[3]);
}

fn draw_rt_farrington_tab(f: &mut Frame<'_>, app: &EpiTuiState, area: Rect) {
    let v_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(56), Constraint::Percentage(44)])
        .split(area);
    let top_cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
        .split(v_chunks[0]);

    let obs_pts: Vec<(f64, f64)> = app
        .active_dataset
        .points
        .iter()
        .enumerate()
        .map(|(i, p)| (i as f64, p.count))
        .collect();
    let thr_pts: Vec<(f64, f64)> = app
        .cached_ab
        .points
        .iter()
        .enumerate()
        .filter_map(|(i, a)| a.upper_threshold.map(|u| (i as f64, u)))
        .collect();
    let max_y = obs_pts
        .iter()
        .chain(thr_pts.iter())
        .map(|p| p.1)
        .fold(10.0_f64, f64::max)
        * 1.12;
    let x_hi = (obs_pts.len().max(1)) as f64;

    let inc_datasets = vec![
        Dataset::default()
            .name("Observed")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Cyan))
            .data(&obs_pts),
        Dataset::default()
            .name("Farrington Upper")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Red))
            .data(&thr_pts),
    ];
    let inc_chart = Chart::new(inc_datasets)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(format!("Weekly Incidence & Farrington Alerts ({} flagged)", app.cached_ab.total_alerts)),
        )
        .x_axis(Axis::default().bounds([0.0, x_hi]))
        .y_axis(Axis::default().bounds([0.0, max_y]));
    f.render_widget(inc_chart, top_cols[0]);

    let rt_pts: Vec<(f64, f64)> = app
        .cached_rt
        .estimates
        .iter()
        .enumerate()
        .filter_map(|(i, e)| e.map(|est| (i as f64, est.mean)))
        .collect();
    let one_line = vec![(0.0, 1.0), (x_hi, 1.0)];
    let rt_max = rt_pts.iter().map(|p| p.1).fold(2.5_f64, f64::max).min(8.0);
    let rt_datasets = vec![
        Dataset::default()
            .name("Cori Rt")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Green))
            .data(&rt_pts),
        Dataset::default()
            .name("Rt = 1.0")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Yellow))
            .data(&one_line),
    ];
    let rt_chart = Chart::new(rt_datasets)
        .block(Block::default().borders(Borders::ALL).title("Instantaneous Cori Rt"))
        .x_axis(Axis::default().bounds([0.0, x_hi]))
        .y_axis(Axis::default().bounds([0.0, rt_max]));
    f.render_widget(rt_chart, top_cols[1]);

    let n = app.active_dataset.points.len();
    let start = n.saturating_sub(8);
    let rows: Vec<Row<'_>> = (start..n)
        .map(|i| {
            let p = &app.active_dataset.points[i];
            let a = &app.cached_ab.points[i];
            let cv = CellValue::from_count(Some(p.count), p.suppressed, app.suppression_policy);
            let rt_s = app.cached_rt.estimates[i]
                .map(|e| format!("{:.2} [{:.2}, {:.2}]", e.mean, e.cri_low_95, e.cri_high_95))
                .unwrap_or_else(|| "--".to_string());
            Row::new(vec![
                Cell::from(p.date_label.clone()),
                Cell::from(cv.display_string(0)),
                Cell::from(
                    a.expected_baseline
                        .map(|v| format!("{v:.1}"))
                        .unwrap_or_else(|| "--".to_string()),
                ),
                Cell::from(
                    a.upper_threshold
                        .map(|v| format!("{v:.1}"))
                        .unwrap_or_else(|| "--".to_string()),
                ),
                Cell::from(rt_s),
                Cell::from(if a.is_alert { "ALERT" } else { "ok" }),
            ])
        })
        .collect();

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(14),
            Constraint::Length(12),
            Constraint::Length(12),
            Constraint::Length(12),
            Constraint::Length(22),
            Constraint::Length(10),
        ],
    )
    .header(
        Row::new(vec!["Period", "Observed", "Expected", "Upper Lim", "Cori Rt [95% CrI]", "Status"])
            .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(Block::default().borders(Borders::ALL).title("Recent Surveillance Ledger"));
    f.render_widget(tbl, v_chunks[1]);
}

fn draw_seir_surge_tab(f: &mut Frame<'_>, app: &EpiTuiState, area: Rect) {
    let v_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(60), Constraint::Percentage(40)])
        .split(area);
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(50), Constraint::Percentage(50)])
        .split(v_chunks[0]);

    let seir_m: Vec<(f64, f64)> = app
        .cached_seir
        .steps
        .iter()
        .map(|s| (s.week_offset as f64, s.reported_weekly_mean))
        .collect();
    let seir_hi: Vec<(f64, f64)> = app
        .cached_seir
        .steps
        .iter()
        .map(|s| (s.week_offset as f64, s.reported_pi_high_95))
        .collect();
    let y_max = seir_hi.iter().map(|p| p.1).fold(20.0_f64, f64::max) * 1.1;
    let x_max = (seir_m.len().max(1)) as f64;

    let seir_chart = Chart::new(vec![
        Dataset::default()
            .name("SEIR Mean")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Cyan))
            .data(&seir_m),
        Dataset::default()
            .name("95% Stochastic PI")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Magenta))
            .data(&seir_hi),
    ])
    .block(Block::default().borders(Borders::ALL).title(format!(
        "Scoped-Parallel SEIR Forecast (R0={:.2}, R_eff={:.2})",
        app.cached_seir.calibrated_r0, app.cached_seir.initial_r_eff
    )))
    .x_axis(Axis::default().bounds([0.0, x_max]))
    .y_axis(Axis::default().bounds([0.0, y_max]));
    f.render_widget(seir_chart, cols[0]);

    let ward_pts: Vec<(f64, f64)> = app
        .cached_surge
        .daily_projections
        .iter()
        .map(|d| (d.day_index as f64, d.ward_census_mean))
        .collect();
    let icu_pts: Vec<(f64, f64)> = app
        .cached_surge
        .daily_projections
        .iter()
        .map(|d| (d.day_index as f64, d.icu_census_mean))
        .collect();
    let vent_pts: Vec<(f64, f64)> = app
        .cached_surge
        .daily_projections
        .iter()
        .map(|d| (d.day_index as f64, d.vent_census_mean))
        .collect();
    let sg_ymax = ward_pts.iter().map(|p| p.1).fold(25.0_f64, f64::max) * 1.15;
    let sg_xmax = (ward_pts.len().max(1)) as f64;

    let sg_chart = Chart::new(vec![
        Dataset::default()
            .name("Acute Beds")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Green))
            .data(&ward_pts),
        Dataset::default()
            .name("ICU Census")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Cyan))
            .data(&icu_pts),
        Dataset::default()
            .name("Ventilators")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Yellow))
            .data(&vent_pts),
    ])
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title("Competing-Risk Hospital / ICU / Ventilator Census"),
    )
    .x_axis(Axis::default().bounds([0.0, sg_xmax]))
    .y_axis(Axis::default().bounds([0.0, sg_ymax]));
    f.render_widget(sg_chart, cols[1]);

    let rows = vec![
        Row::new(vec![
            Cell::from("Nishiura Delay-Adjusted cCFR (Active Wave)"),
            Cell::from(format!(
                "{:.2}% [95% CI: {:.2}%, {:.2}%]",
                app.cached_cfr.final_adjusted_ccfr_pct,
                app.cached_cfr.final_adjusted_ci_low_95,
                app.cached_cfr.final_adjusted_ci_high_95
            )),
            Cell::from(format!(
                "Naive CFR {:.2}% ({:.2}x right-truncation bias)",
                app.cached_cfr.final_crude_cfr_pct, app.cached_cfr.peak_truncation_bias_factor
            )),
        ]),
        Row::new(vec![
            Cell::from("Nowcast Imputed Unreported Cases"),
            Cell::from(format!("{:.1} cases", app.cached_nowcast.total_imputed_unreported)),
            Cell::from("Parametric Gamma reporting-delay adjustment"),
        ]),
        Row::new(vec![
            Cell::from("Competing-Risk Peak ICU / Ventilator Demand"),
            Cell::from(format!(
                "ICU {:.1} ({:.0}%) | IMV {:.1} ({:.0}%)",
                app.cached_surge.peak_icu_census,
                app.cached_surge.peak_icu_utilization_pct,
                app.cached_surge.peak_vent_census,
                app.cached_surge.peak_vent_utilization_pct
            )),
            Cell::from(format!(
                "Implied HFR {:.2}% | ICU deficit {:.1} bed-days",
                app.cached_surge.implied_hospital_fatality_pct,
                app.cached_surge.icu_bed_days_deficit
            )),
        ]),
    ];
    let tbl = Table::new(
        rows,
        [
            Constraint::Length(42),
            Constraint::Length(34),
            Constraint::Min(28),
        ],
    )
    .header(
        Row::new(vec!["Severity / Forecast Metric", "Estimate", "Clinical Interpretation"])
            .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(Block::default().borders(Borders::ALL).title("Nowcast, SEIR, cCFR & Hospital Capacity Summary"));
    f.render_widget(tbl, v_chunks[1]);
}

fn draw_spatial_scan_tab(f: &mut Frame<'_>, app: &EpiTuiState, area: Rect) {
    let rows: Vec<Row<'_>> = app
        .cached_clusters
        .iter()
        .enumerate()
        .map(|(idx, c)| {
            Row::new(vec![
                Cell::from(format!("#{}", idx + 1)),
                Cell::from(c.center_id.clone()),
                Cell::from(format!("{}", c.member_ids.len())),
                Cell::from(format!("{:.0} km", c.radius_km)),
                Cell::from(format!("{:.0}", c.observed_cases)),
                Cell::from(format!("{:.1}", c.expected_cases)),
                Cell::from(format!("{:.2}x", c.relative_risk)),
                Cell::from(format!("{:.2}", c.kulldorff_llr)),
                Cell::from(format!("{:.4}", c.monte_carlo_p_value)),
            ])
        })
        .collect();

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(6),
            Constraint::Length(10),
            Constraint::Length(8),
            Constraint::Length(12),
            Constraint::Length(10),
            Constraint::Length(10),
            Constraint::Length(9),
            Constraint::Length(10),
            Constraint::Length(10),
        ],
    )
    .header(
        Row::new(vec![
            "Rank", "Center", "Zones", "Radius", "Observed", "Expected", "RelRisk", "LLR", "MC p-val",
        ])
        .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title("Scoped-Parallel Multi-Zone Circular Kulldorff Spatial Scan (US Pertussis)"),
    );
    f.render_widget(tbl, area);
}

fn draw_superspreading_tab(f: &mut Frame<'_>, app: &EpiTuiState, area: Rect) {
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(46), Constraint::Percentage(54)])
        .split(area);

    let off = &app.cached_net.offspring_fit;
    let theo_pts: Vec<(f64, f64)> = off
        .theoretical_lorenz_curve
        .iter()
        .map(|p| (p[0] * 100.0, p[1] * 100.0))
        .collect();
    let diag = vec![(0.0, 0.0), (100.0, 100.0)];

    let chart = Chart::new(vec![
        Dataset::default()
            .name("Lloyd-Smith NB Lorenz")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Cyan))
            .data(&theo_pts),
        Dataset::default()
            .name("Homogeneous Equality")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::DarkGray))
            .data(&diag),
    ])
    .block(Block::default().borders(Borders::ALL).title(format!(
        "Offspring k={:.3} (Top 20% -> {:.1}% tx | q_ext={:.1}%)",
        off.k_mle,
        off.prop_transmission_top_20_theoretical * 100.0,
        off.single_intro_extinction_prob * 100.0
    )))
    .x_axis(Axis::default().bounds([0.0, 100.0]))
    .y_axis(Axis::default().bounds([0.0, 100.0]));
    f.render_widget(chart, cols[0]);

    let rows: Vec<Row<'_>> = app
        .cached_net
        .venue_metrics
        .iter()
        .map(|v| {
            Row::new(vec![
                Cell::from(v.venue_id.clone()),
                Cell::from(v.venue_category.clone()),
                Cell::from(format!("{}", v.secondary_transmissions)),
                Cell::from(format!("{:.1}%", v.share_of_outbreak_transmissions * 100.0)),
                Cell::from(format!("{}", v.bipartite_k_core)),
                Cell::from(format!("{:.3}", v.bipartite_betweenness)),
            ])
        })
        .collect();

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(18),
            Constraint::Length(16),
            Constraint::Length(7),
            Constraint::Length(9),
            Constraint::Length(7),
            Constraint::Length(11),
        ],
    )
    .header(
        Row::new(vec!["Venue ID", "Setting", "SecTx", "Share", "k-Core", "Betweenness"])
            .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(Block::default().borders(Borders::ALL).title("Bipartite Exposure Venue Centrality"));
    f.render_widget(tbl, cols[1]);
}

fn draw_wastewater_amr_tab(f: &mut Frame<'_>, app: &EpiTuiState, area: Rect) {
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(45), Constraint::Percentage(55)])
        .split(area);

    let jn1: Vec<(f64, f64)> = app
        .cached_ww
        .samples
        .iter()
        .map(|s| (s.time_weeks, s.lineages[0].abundance * 100.0))
        .collect();
    let xec: Vec<(f64, f64)> = app
        .cached_ww
        .samples
        .iter()
        .map(|s| (s.time_weeks, s.lineages.last().map(|l| l.abundance * 100.0).unwrap_or(0.0)))
        .collect();

    let chart = Chart::new(vec![
        Dataset::default()
            .name("JN.1 (%)")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Cyan))
            .data(&jn1),
        Dataset::default()
            .name("XEC (%)")
            .marker(symbols::Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Style::default().fg(Color::Red))
            .data(&xec),
    ])
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title("Freyja L1 Wastewater SNV Lineage Deconvolution"),
    )
    .x_axis(Axis::default().bounds([0.0, 9.5]))
    .y_axis(Axis::default().bounds([0.0, 100.0]));
    f.render_widget(chart, cols[0]);

    let rows: Vec<Row<'_>> = app
        .cached_amr
        .cells
        .iter()
        .take(14)
        .map(|c| {
            Row::new(vec![
                Cell::from(c.organism.clone()),
                Cell::from(c.drug_code.clone()),
                Cell::from(format!("{}", c.isolates_tested)),
                Cell::from(format!("{:.1}%", c.pct_susceptible)),
                Cell::from(if !c.meets_clsi_m39_n30 {
                    "LOW_N<30"
                } else if c.empiric_caution_below_threshold {
                    "AVOID_EMP"
                } else {
                    "OK"
                }),
            ])
        })
        .collect();

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(23),
            Constraint::Length(6),
            Constraint::Length(6),
            Constraint::Length(8),
            Constraint::Length(12),
        ],
    )
    .header(
        Row::new(vec!["Organism", "Drug", "N", "%S", "CLSI M39"])
            .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(Block::default().borders(Borders::ALL).title(format!(
        "CLSI M39 Antibiogram ({} repeat cultures excluded)",
        app.cached_amr.total_duplicates_excluded
    )));
    f.render_widget(tbl, cols[1]);
}

fn draw_crc_essence_tab(f: &mut Frame<'_>, app: &EpiTuiState, area: Rect) {
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(52), Constraint::Percentage(48)])
        .split(area);

    let crc_rows: Vec<Row<'_>> = app
        .cached_crc
        .models
        .iter()
        .enumerate()
        .map(|(i, m)| {
            let star = if i == app.cached_crc.bic_selected_index {
                "*"
            } else {
                " "
            };
            Row::new(vec![
                Cell::from(format!("{star}{}", m.model_id)),
                Cell::from(format!("{:.1}", m.bic)),
                Cell::from(format!("{:.1}", m.unobserved_n000_hat)),
                Cell::from(format!("{:.1}", m.total_n_hat)),
                Cell::from(format!("{:.1}%", m.overall_completeness_pct)),
            ])
        })
        .collect();

    let crc_tbl = Table::new(
        crc_rows,
        [
            Constraint::Length(12),
            Constraint::Length(9),
            Constraint::Length(10),
            Constraint::Length(11),
            Constraint::Length(10),
        ],
    )
    .header(
        Row::new(vec!["Model", "BIC", "Unobs 000", "Total N_hat", "Complete"])
            .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title("3-Source Log-Linear Capture-Recapture (* = BIC Optimal)"),
    );
    f.render_widget(crc_tbl, cols[0]);

    let ess_rows: Vec<Row<'_>> = app
        .cached_essence
        .syndrome_summaries
        .iter()
        .map(|r| {
            Row::new(vec![
                Cell::from(r.code.clone()),
                Cell::from(format!("{}", r.matched_visits)),
                Cell::from(format!("{:.1}%", r.pct_of_ed_visits)),
                Cell::from(format!("{}", r.negated_false_positives_blocked)),
            ])
        })
        .collect();

    let ess_tbl = Table::new(
        ess_rows,
        [
            Constraint::Length(18),
            Constraint::Length(8),
            Constraint::Length(9),
            Constraint::Length(12),
        ],
    )
    .header(
        Row::new(vec!["ESSENCE Code", "Hits", "ED %", "NegEx Blocked"])
            .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title("NSSP ESSENCE Chief-Complaint & Negation Classifier"),
    );
    f.render_widget(ess_tbl, cols[1]);
}

fn draw_wis_ensemble_tab(f: &mut Frame<'_>, app: &EpiTuiState, area: Rect) {
    let rows: Vec<Row<'_>> = app
        .cached_wis
        .scorecards
        .iter()
        .map(|sc| {
            Row::new(vec![
                Cell::from(sc.model_label.clone()),
                Cell::from(format!("{:.2}", sc.mean_wis)),
                Cell::from(format!("{:.3}", sc.mean_log_wis)),
                Cell::from(format!("{:.2}", sc.mean_sharpness)),
                Cell::from(format!("{:.1}%", sc.coverage_95_pct)),
                Cell::from(format!("{:.1}%", sc.final_ensemble_weight * 100.0)),
            ])
        })
        .collect();

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(30),
            Constraint::Length(12),
            Constraint::Length(12),
            Constraint::Length(12),
            Constraint::Length(12),
            Constraint::Length(14),
        ],
    )
    .header(
        Row::new(vec![
            "Forecast Model",
            "Mean WIS",
            "Log-WIS",
            "Sharpness",
            "95% Cov",
            "Inv-WIS Weight",
        ])
        .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title("K=11 Weighted Interval Score (WIS) Backtest & Causal Ensemble Stacking"),
    );
    f.render_widget(tbl, area);
}

fn draw_dp_audit_tab(f: &mut Frame<'_>, app: &EpiTuiState, area: Rect) {
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(44), Constraint::Percentage(56)])
        .split(area);

    let v_rows: Vec<Row<'_>> = app
        .vintages
        .iter()
        .enumerate()
        .map(|(i, v)| {
            let mark = if i == app.selected_vintage_idx { ">" } else { " " };
            Row::new(vec![
                Cell::from(format!("{mark}{}", v.vintage_id)),
                Cell::from(format!("{}", v.row_count)),
                Cell::from(v.checksum_hex.clone()),
            ])
        })
        .collect();

    let v_tbl = Table::new(
        v_rows,
        [
            Constraint::Length(28),
            Constraint::Length(8),
            Constraint::Length(12),
        ],
    )
    .header(
        Row::new(vec!["SQLite Vintage ID ([v] cycle)", "Pts", "Checksum"])
            .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(Block::default().borders(Borders::ALL).title(format!(
        "Bitemporal Ledger (Active DP eps spent = {:.2} | [d] Release)",
        app.spent_epsilon
    )));
    f.render_widget(v_tbl, cols[0]);

    let a_rows: Vec<Row<'_>> = app
        .audit_logs
        .iter()
        .map(|e| {
            let short_hash: String = e.entry_hash_hex.chars().take(12).collect();
            Row::new(vec![
                Cell::from(format!("#{}", e.seq_no)),
                Cell::from(e.action_kind.clone()),
                Cell::from(e.target_ref.clone()),
                Cell::from(short_hash),
            ])
        })
        .collect();

    let a_tbl = Table::new(
        a_rows,
        [
            Constraint::Length(6),
            Constraint::Length(18),
            Constraint::Length(24),
            Constraint::Length(14),
        ],
    )
    .header(
        Row::new(vec!["Seq", "Action", "Target", "SHA-256 Prefix"])
            .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
    )
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title("Genesis-to-Tip SHA-256 Tamper-Evident Audit Chain"),
    );
    f.render_widget(a_tbl, cols[1]);
}
