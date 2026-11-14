mod cartogram;
mod cli;
mod crc_syndromic_tab;
mod dp_audit_tab;
mod elr_pprl_sero_tab;
mod env_mortality_tab;
mod genomic_tab;
mod gui;
mod network_tab;
mod scoring_sitrep_tab;
mod severity_surge_tab;
mod sitrep_export;
mod tui;
mod wastewater_amr_tab;

use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let db_path = std::env::var("EPILEDGER_DB")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("epiledger.sqlite3"));

    if args.len() > 1 {
        if args[1] == "tui" {
            tui::run_tui(&db_path)?;
            return Ok(());
        }
        if args[1] != "gui" {
            cli::run_cli(&args[1..], &db_path)?;
            return Ok(());
        }
    }

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1480.0, 920.0])
            .with_min_inner_size([1060.0, 660.0])
            .with_title("EpiLedger — Native Rust Public Health Informatics Workbench"),
        ..Default::default()
    };

    eframe::run_native(
        "EpiLedger Workbench",
        native_options,
        Box::new(move |cc| Ok(Box::new(gui::EpiLedgerGuiApp::new(cc, db_path)))),
    )
    .map_err(|e| format!("eframe launch error: {e}").into())
}
