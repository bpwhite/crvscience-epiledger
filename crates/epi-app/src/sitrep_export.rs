use epi_analytics::{
    AberrationSeries, CoriRtSeries, EnsembleBacktestResult, MultiZoneCluster, NowcastSeries,
};
use epi_core::{CellValue, SuppressionPolicy, TimeSeriesDataset};

/// Zero-JS, print-ready, self-contained HTML5 + inline vector `<svg>` Epidemiological Situation Report (`SitRep`)
/// generator embedding provenance `SHA-256` checksums, active small-cell disclosure guards (`NCHS <= 9`),
/// aberration alerts, nowcasted incidence, `K = 11` `WIS` decomposition bars, and multi-model stacked forecasts.

pub struct SitRepContext<'a> {
    pub dataset: &'a TimeSeriesDataset,
    pub suppression_policy: SuppressionPolicy,
    pub rt_series: &'a CoriRtSeries,
    pub ab_series: &'a AberrationSeries,
    pub nowcast: &'a NowcastSeries,
    pub backtest: &'a EnsembleBacktestResult,
    pub spatial_clusters: &'a [MultiZoneCluster],
}

fn html_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(ch),
        }
    }
    out
}

fn build_incidence_and_aberration_svg(ctx: &SitRepContext<'_>) -> String {
    let n_total = ctx.dataset.points.len();
    if n_total == 0 {
        return String::new();
    }
    let win = n_total.min(52);
    let start = n_total - win;

    let w_svg = 760.0_f64;
    let h_svg = 220.0_f64;
    let pad_l = 52.0_f64;
    let pad_r = 18.0_f64;
    let pad_t = 18.0_f64;
    let pad_b = 34.0_f64;
    let pw = w_svg - pad_l - pad_r;
    let ph = h_svg - pad_t - pad_b;

    let mut max_y = 5.0_f64;
    for i in start..n_total {
        max_y = max_y.max(ctx.dataset.points[i].count);
        if let Some(u) = ctx.ab_series.points[i].upper_threshold {
            max_y = max_y.max(u);
        }
    }
    max_y *= 1.12;

    let x_of = |offset: usize| -> f64 {
        if win <= 1 {
            pad_l + 0.5 * pw
        } else {
            pad_l + (offset as f64) * pw / ((win - 1) as f64)
        }
    };
    let y_of = |val: f64| -> f64 { pad_t + ph - (val.clamp(0.0, max_y) / max_y) * ph };

    let mut svg = format!(
        "<svg viewBox=\"0 0 {w_svg:.0} {h_svg:.0}\" width=\"100%\" xmlns=\"http://www.w3.org/2000/svg\">\n"
    );
    // Horizontal gridlines
    for g in 0..=4 {
        let v = max_y * (g as f64) / 4.0;
        let yg = y_of(v);
        svg.push_str(&format!(
            "  <line x1=\"{pad_l:.1}\" y1=\"{yg:.1}\" x2=\"{:.1}\" y2=\"{yg:.1}\" stroke=\"#e2e8f0\" stroke-width=\"1\"/>\n",
            w_svg - pad_r
        ));
        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"ui-monospace, monospace\" font-size=\"10\" fill=\"#475569\" text-anchor=\"end\">{v:.0}</text>\n",
            pad_l - 6.0,
            yg + 3.5
        ));
    }

    // Bars for observed weekly counts
    let bar_w = (pw / (win as f64) * 0.66).clamp(2.5, 11.0);
    for (off, i) in (start..n_total).enumerate() {
        let raw_c = ctx.dataset.points[i].count;
        let disp_c = if ctx.suppression_policy == SuppressionPolicy::NchsWonderLe9
            && raw_c > 0.0
            && raw_c <= 9.0
        {
            0.0
        } else {
            raw_c
        };
        let xc = x_of(off);
        let y_top = y_of(disp_c);
        let bh = (pad_t + ph - y_top).max(0.5);
        let fill = if ctx.ab_series.points[i].is_alert {
            "#e11d48"
        } else {
            "#3b82f6"
        };
        svg.push_str(&format!(
            "  <rect x=\"{:.1}\" y=\"{y_top:.1}\" width=\"{bar_w:.1}\" height=\"{bh:.1}\" fill=\"{fill}\" rx=\"1.5\"/>\n",
            xc - 0.5 * bar_w
        ));
    }

    // Polyline for upper aberration threshold
    let mut thr_pts = Vec::new();
    for (off, i) in (start..n_total).enumerate() {
        if let Some(u) = ctx.ab_series.points[i].upper_threshold {
            thr_pts.push(format!("{:.1},{:.1}", x_of(off), y_of(u)));
        }
    }
    if thr_pts.len() >= 2 {
        svg.push_str(&format!(
            "  <polyline fill=\"none\" stroke=\"#d97706\" stroke-width=\"1.8\" stroke-dasharray=\"4 3\" points=\"{}\"/>\n",
            thr_pts.join(" ")
        ));
    }

    // X-axis first/middle/last date labels
    let first_lbl = html_escape(&ctx.dataset.points[start].date_label);
    let mid_lbl = html_escape(&ctx.dataset.points[start + win / 2].date_label);
    let last_lbl = html_escape(&ctx.dataset.points[n_total - 1].date_label);
    let y_lbl = h_svg - 10.0;
    svg.push_str(&format!(
        "  <text x=\"{pad_l:.1}\" y=\"{y_lbl:.1}\" font-family=\"ui-monospace, monospace\" font-size=\"10\" fill=\"#475569\">{first_lbl}</text>\n"
    ));
    svg.push_str(&format!(
        "  <text x=\"{:.1}\" y=\"{y_lbl:.1}\" font-family=\"ui-monospace, monospace\" font-size=\"10\" fill=\"#475569\" text-anchor=\"middle\">{mid_lbl}</text>\n",
        pad_l + 0.5 * pw
    ));
    svg.push_str(&format!(
        "  <text x=\"{:.1}\" y=\"{y_lbl:.1}\" font-family=\"ui-monospace, monospace\" font-size=\"10\" fill=\"#475569\" text-anchor=\"end\">{last_lbl}</text>\n",
        w_svg - pad_r
    ));

    svg.push_str("</svg>\n");
    svg
}

fn build_ensemble_fan_svg(ctx: &SitRepContext<'_>) -> String {
    let n_obs = ctx.dataset.points.len();
    if n_obs == 0 || ctx.backtest.forward_steps.is_empty() {
        return String::new();
    }
    let tail = n_obs.min(20);
    let start = n_obs - tail;
    let h_len = ctx.backtest.forward_steps.len();
    let total_x = tail + h_len;

    let w_svg = 760.0_f64;
    let h_svg = 215.0_f64;
    let pad_l = 52.0_f64;
    let pad_r = 18.0_f64;
    let pad_t = 18.0_f64;
    let pad_b = 32.0_f64;
    let pw = w_svg - pad_l - pad_r;
    let ph = h_svg - pad_t - pad_b;

    let mut max_y = 5.0_f64;
    for i in start..n_obs {
        max_y = max_y.max(ctx.nowcast.points[i].nowcast_mean);
    }
    for st in &ctx.backtest.forward_steps {
        max_y = max_y.max(st.inverse_wis_ensemble.pi95_high());
    }
    max_y *= 1.12;

    let x_of = |step: usize| -> f64 {
        pad_l + (step as f64) * pw / ((total_x.saturating_sub(1).max(1)) as f64)
    };
    let y_of = |val: f64| -> f64 { pad_t + ph - (val.clamp(0.0, max_y) / max_y) * ph };

    let mut svg = format!(
        "<svg viewBox=\"0 0 {w_svg:.0} {h_svg:.0}\" width=\"100%\" xmlns=\"http://www.w3.org/2000/svg\">\n"
    );
    for g in 0..=4 {
        let v = max_y * (g as f64) / 4.0;
        let yg = y_of(v);
        svg.push_str(&format!(
            "  <line x1=\"{pad_l:.1}\" y1=\"{yg:.1}\" x2=\"{:.1}\" y2=\"{yg:.1}\" stroke=\"#e2e8f0\" stroke-width=\"1\"/>\n",
            w_svg - pad_r
        ));
        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"ui-monospace, monospace\" font-size=\"10\" fill=\"#475569\" text-anchor=\"end\">{v:.0}</text>\n",
            pad_l - 6.0,
            yg + 3.5
        ));
    }

    // Forecast origin divider
    let x_orig = x_of(tail - 1);
    svg.push_str(&format!(
        "  <line x1=\"{x_orig:.1}\" y1=\"{pad_t:.1}\" x2=\"{x_orig:.1}\" y2=\"{:.1}\" stroke=\"#94a3b8\" stroke-dasharray=\"3 3\" stroke-width=\"1.2\"/>\n",
        pad_t + ph
    ));

    // Build 95% and 50% prediction ribbon polygons anchored at last nowcast point
    let anchor_y = ctx.nowcast.points[n_obs - 1].nowcast_mean;
    let mut p95_upper = vec![(x_orig, y_of(anchor_y))];
    let mut p95_lower = vec![(x_orig, y_of(anchor_y))];
    let mut p50_upper = vec![(x_orig, y_of(anchor_y))];
    let mut p50_lower = vec![(x_orig, y_of(anchor_y))];
    let mut ens_med_pts = vec![format!("{x_orig:.1},{:.1}", y_of(anchor_y))];

    for (idx, st) in ctx.backtest.forward_steps.iter().enumerate() {
        let xs = x_of(tail + idx);
        p95_upper.push((xs, y_of(st.inverse_wis_ensemble.pi95_high())));
        p95_lower.push((xs, y_of(st.inverse_wis_ensemble.pi95_low())));
        p50_upper.push((xs, y_of(st.inverse_wis_ensemble.pi50_high())));
        p50_lower.push((xs, y_of(st.inverse_wis_ensemble.pi50_low())));
        ens_med_pts.push(format!(
            "{xs:.1},{:.1}",
            y_of(st.inverse_wis_ensemble.median())
        ));
    }

    let poly_str = |up: &[(f64, f64)], lo: &[(f64, f64)]| -> String {
        let mut pts = Vec::with_capacity(up.len() + lo.len());
        for &(x, y) in up {
            pts.push(format!("{x:.1},{y:.1}"));
        }
        for &(x, y) in lo.iter().rev() {
            pts.push(format!("{x:.1},{y:.1}"));
        }
        pts.join(" ")
    };

    svg.push_str(&format!(
        "  <polygon points=\"{}\" fill=\"#bae6fd\" fill-opacity=\"0.55\"/>\n",
        poly_str(&p95_upper, &p95_lower)
    ));
    svg.push_str(&format!(
        "  <polygon points=\"{}\" fill=\"#38bdf8\" fill-opacity=\"0.45\"/>\n",
        poly_str(&p50_upper, &p50_lower)
    ));

    // Observed + Nowcast polylines
    let obs_pts: Vec<String> = (0..tail)
        .map(|k| {
            format!(
                "{:.1},{:.1}",
                x_of(k),
                y_of(ctx.dataset.points[start + k].count)
            )
        })
        .collect();
    let nc_pts: Vec<String> = (0..tail)
        .map(|k| {
            format!(
                "{:.1},{:.1}",
                x_of(k),
                y_of(ctx.nowcast.points[start + k].nowcast_mean)
            )
        })
        .collect();

    svg.push_str(&format!(
        "  <polyline fill=\"none\" stroke=\"#334155\" stroke-width=\"2.0\" points=\"{}\"/>\n",
        obs_pts.join(" ")
    ));
    svg.push_str(&format!(
        "  <polyline fill=\"none\" stroke=\"#10b981\" stroke-width=\"1.8\" stroke-dasharray=\"4 2\" points=\"{}\"/>\n",
        nc_pts.join(" ")
    ));
    svg.push_str(&format!(
        "  <polyline fill=\"none\" stroke=\"#0284c7\" stroke-width=\"2.4\" points=\"{}\"/>\n",
        ens_med_pts.join(" ")
    ));

    svg.push_str(&format!(
        "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#0369a1\">Inverse-WIS Ensemble (50% &amp; 95% PI)</text>\n",
        x_orig + 6.0,
        pad_t + 12.0
    ));
    svg.push_str("</svg>\n");
    svg
}

fn build_wis_decomposition_svg(ctx: &SitRepContext<'_>) -> String {
    let cards = &ctx.backtest.scorecards;
    if cards.is_empty() {
        return String::new();
    }
    let w_svg = 760.0_f64;
    let row_h = 28.0_f64;
    let pad_t = 24.0_f64;
    let pad_b = 18.0_f64;
    let pad_l = 215.0_f64;
    let pad_r = 70.0_f64;
    let h_svg = pad_t + pad_b + (cards.len() as f64) * row_h;
    let bar_max_w = w_svg - pad_l - pad_r;

    let max_wis = cards
        .iter()
        .map(|c| c.mean_wis)
        .fold(1.0_f64, f64::max)
        * 1.08;

    let mut svg = format!(
        "<svg viewBox=\"0 0 {w_svg:.0} {h_svg:.0}\" width=\"100%\" xmlns=\"http://www.w3.org/2000/svg\">\n"
    );
    // Legend header
    svg.push_str("  <rect x=\"215\" y=\"4\" width=\"10\" height=\"10\" fill=\"#0284c7\"/><text x=\"229\" y=\"13\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#334155\">Sharpness (Width)</text>\n");
    svg.push_str("  <rect x=\"335\" y=\"4\" width=\"10\" height=\"10\" fill=\"#f59e0b\"/><text x=\"349\" y=\"13\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#334155\">Overprediction</text>\n");
    svg.push_str("  <rect x=\"445\" y=\"4\" width=\"10\" height=\"10\" fill=\"#e11d48\"/><text x=\"459\" y=\"13\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#334155\">Underprediction</text>\n");

    for (idx, sc) in cards.iter().enumerate() {
        let y0 = pad_t + (idx as f64) * row_h;
        let w_sh = (sc.mean_sharpness / max_wis) * bar_max_w;
        let w_ov = (sc.mean_overprediction / max_wis) * bar_max_w;
        let w_un = (sc.mean_underprediction / max_wis) * bar_max_w;

        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#0f172a\" text-anchor=\"end\">{}</text>\n",
            pad_l - 10.0,
            y0 + 14.0,
            html_escape(&sc.model_label)
        ));
        svg.push_str(&format!(
            "  <rect x=\"{pad_l:.1}\" y=\"{:.1}\" width=\"{w_sh:.1}\" height=\"16\" fill=\"#0284c7\"/>\n",
            y0 + 2.0
        ));
        svg.push_str(&format!(
            "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{w_ov:.1}\" height=\"16\" fill=\"#f59e0b\"/>\n",
            pad_l + w_sh,
            y0 + 2.0
        ));
        svg.push_str(&format!(
            "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{w_un:.1}\" height=\"16\" fill=\"#e11d48\"/>\n",
            pad_l + w_sh + w_ov,
            y0 + 2.0
        ));
        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"ui-monospace, monospace\" font-size=\"10.5\" fill=\"#334155\">{:.2} ({:.2}x)</text>\n",
            pad_l + w_sh + w_ov + w_un + 6.0,
            y0 + 14.0,
            sc.mean_wis,
            sc.relative_wis_vs_naive
        ));
    }

    svg.push_str("</svg>\n");
    svg
}

/// Render a standalone, zero-JS HTML5 + inline vector `<svg>` Epidemiological Situation Report (`SitRep`).
pub fn render_standalone_sitrep_html(ctx: &SitRepContext<'_>) -> String {
    let latest_pt = ctx.dataset.points.last();
    let latest_disp = latest_pt
        .map(|p| {
            CellValue::from_count(Some(p.count), p.suppressed, ctx.suppression_policy)
                .display_string(0)
        })
        .unwrap_or_else(|| "--".to_string());
    let latest_period = latest_pt
        .map(|p| html_escape(&p.date_label))
        .unwrap_or_else(|| "--".to_string());

    let latest_rt = ctx
        .rt_series
        .estimates
        .iter()
        .rev()
        .flatten()
        .next()
        .map(|r| format!("{:.2} [{:.2}, {:.2}]", r.mean, r.cri_low_95, r.cri_high_95))
        .unwrap_or_else(|| "--".to_string());

    let inv_card = ctx
        .backtest
        .scorecards
        .iter()
        .find(|s| s.model_id == "inv_wis_ens");
    let rel_skill = inv_card
        .map(|s| format!("{:.2}x vs Naive", s.relative_wis_vs_naive))
        .unwrap_or_else(|| "--".to_string());

    let next_wk_fc = ctx
        .backtest
        .forward_steps
        .first()
        .map(|st| {
            format!(
                "{:.1} [{:.0}, {:.0}]",
                st.inverse_wis_ensemble.median(),
                st.inverse_wis_ensemble.pi95_low(),
                st.inverse_wis_ensemble.pi95_high()
            )
        })
        .unwrap_or_else(|| "--".to_string());

    let fig1 = build_incidence_and_aberration_svg(ctx);
    let fig2 = build_ensemble_fan_svg(ctx);
    let fig3 = build_wis_decomposition_svg(ctx);

    let mut score_rows = String::new();
    for sc in &ctx.backtest.scorecards {
        let w_txt = if sc.final_ensemble_weight > 0.0 {
            format!("{:.1}%", sc.final_ensemble_weight * 100.0)
        } else {
            "--".to_string()
        };
        score_rows.push_str(&format!(
            "<tr><td>{}</td><td class=\"num\">{:.2}</td><td class=\"num\">{:.2}</td><td class=\"num\">{:.2}</td><td class=\"num\">{:.2}</td><td class=\"num\">{:.3}</td><td class=\"num\">{:.2}x</td><td class=\"num\">{:.1}%</td><td class=\"num\">{}</td></tr>\n",
            html_escape(&sc.model_label),
            sc.mean_wis,
            sc.mean_sharpness,
            sc.mean_overprediction,
            sc.mean_underprediction,
            sc.mean_log_wis,
            sc.relative_wis_vs_naive,
            sc.coverage_95_pct,
            w_txt
        ));
    }

    let mut fc_rows = String::new();
    for st in &ctx.backtest.forward_steps {
        let ens = &st.inverse_wis_ensemble;
        fc_rows.push_str(&format!(
            "<tr><td>Week +{}</td><td class=\"num\">{:.1}</td><td class=\"num\">[{:.1}, {:.1}]</td><td class=\"num\">[{:.1}, {:.1}]</td><td class=\"num\">{:.1}</td><td class=\"num\">{:.1}</td><td class=\"num\">{:.1}</td></tr>\n",
            st.horizon_week,
            ens.median(),
            ens.pi50_low(),
            ens.pi50_high(),
            ens.pi95_low(),
            ens.pi95_high(),
            st.seir_median,
            st.renewal_median,
            st.trend_median
        ));
    }

    let mut cluster_rows = String::new();
    for (idx, cl) in ctx.spatial_clusters.iter().take(5).enumerate() {
        cluster_rows.push_str(&format!(
            "<tr><td>#{}</td><td>{}</td><td class=\"num\">{:.0} km</td><td class=\"num\">{:.0}</td><td class=\"num\">{:.1}</td><td class=\"num\">{:.2}x</td><td class=\"num\">{:.1}</td><td class=\"num\">{:.3}</td></tr>\n",
            idx + 1,
            html_escape(&cl.center_name),
            cl.radius_km,
            cl.observed_cases,
            cl.expected_cases,
            cl.relative_risk,
            cl.kulldorff_llr,
            cl.monte_carlo_p_value
        ));
    }

    let mask_badge = match ctx.suppression_policy {
        SuppressionPolicy::Unmasked => "Internal Unmasked View",
        SuppressionPolicy::NchsWonderLe9 => "NCHS CDC WONDER Small-Cell Guard (Counts <= 9 Suppressed)",
    };

    format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>EpiLedger Epidemiological Situation Report — {title}</title>
  <style>
    @page {{ size: letter; margin: 0.5in; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; color: #0f172a; background: #f8fafc; margin: 0; padding: 24px; line-height: 1.45; }}
    .page {{ max-width: 920px; margin: 0 auto; background: #ffffff; border: 1px solid #cbd5e1; border-radius: 10px; padding: 28px 32px; }}
    h1 {{ font-size: 20px; margin: 0 0 4px 0; color: #0f172a; }}
    h2 {{ font-size: 14px; text-transform: uppercase; letter-spacing: 0.04em; color: #334155; border-bottom: 1px solid #e2e8f0; padding-bottom: 5px; margin: 22px 0 10px 0; }}
    .meta {{ font-family: ui-monospace, monospace; font-size: 11.5px; color: #475569; margin-bottom: 16px; }}
    .kpi-grid {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px; margin-bottom: 16px; }}
    .kpi {{ border: 1px solid #e2e8f0; border-radius: 8px; padding: 10px 12px; background: #f8fafc; }}
    .kpi-lbl {{ font-size: 10.5px; text-transform: uppercase; color: #64748b; font-weight: 600; }}
    .kpi-val {{ font-family: ui-monospace, monospace; font-size: 17px; font-weight: 700; color: #0f172a; margin-top: 3px; }}
    .kpi-sub {{ font-family: ui-monospace, monospace; font-size: 10.5px; color: #475569; }}
    .chart-box {{ border: 1px solid #e2e8f0; border-radius: 8px; padding: 10px 12px; margin-bottom: 14px; background: #ffffff; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 12px; margin-bottom: 12px; }}
    th {{ background: #f1f5f9; color: #334155; text-align: left; padding: 6px 8px; border-bottom: 2px solid #cbd5e1; }}
    td {{ padding: 5px 8px; border-bottom: 1px solid #e2e8f0; }}
    td.num, th.num {{ text-align: right; font-family: ui-monospace, monospace; }}
    .refs {{ font-size: 11px; color: #475569; margin-top: 16px; }}
    .refs a {{ color: #0284c7; text-decoration: none; }}
    @media print {{ body {{ background: #ffffff; padding: 0; }} .page {{ border: none; padding: 0; }} }}
  </style>
</head>
<body>
<div class="page">
  <h1>EpiLedger Epidemiological Situation Report (SitRep)</h1>
  <div class="meta">
    Active Vintage: <strong>{vid}</strong> | Title: {title} | SHA-256: <strong>{hash}</strong> | Disclosure Policy: {mask_badge}
  </div>

  <div class="kpi-grid">
    <div class="kpi">
      <div class="kpi-lbl">Latest Observed</div>
      <div class="kpi-val">{latest_disp}</div>
      <div class="kpi-sub">Period {latest_period} (Nowcast +{nc_imp:.1})</div>
    </div>
    <div class="kpi">
      <div class="kpi-lbl">Instantaneous Cori Rt</div>
      <div class="kpi-val">{latest_rt}</div>
      <div class="kpi-sub">Aberration alerts: {total_alerts}</div>
    </div>
    <div class="kpi">
      <div class="kpi-lbl">Inverse-WIS Skill</div>
      <div class="kpi-val">{rel_skill}</div>
      <div class="kpi-sub">{origins} rolling origins ({horizon}w horizon)</div>
    </div>
    <div class="kpi">
      <div class="kpi-lbl">Wk +1 Stacked Forecast</div>
      <div class="kpi-val">{next_wk_fc}</div>
      <div class="kpi-sub">Median [95% Prediction Interval]</div>
    </div>
  </div>

  <h2>1. Surveillance Incidence &amp; Aberration Thresholds (Trailing 52 Steps)</h2>
  <div class="chart-box">{fig1}</div>

  <h2>2. Nowcast-Adjusted Trajectory &amp; Trained Inverse-WIS Probabilistic Ensemble</h2>
  <div class="chart-box">{fig2}</div>

  <h2>3. Out-of-Sample Weighted Interval Score (K = 11 WIS) Decomposition</h2>
  <div class="chart-box">{fig3}</div>

  <table>
    <thead>
      <tr>
        <th>Forecasting Model</th>
        <th class="num">Mean WIS</th>
        <th class="num">Sharpness</th>
        <th class="num">Overpred</th>
        <th class="num">Underpred</th>
        <th class="num">Log-WIS</th>
        <th class="num">Rel. WIS</th>
        <th class="num">95% Cov</th>
        <th class="num">Stack Wt</th>
      </tr>
    </thead>
    <tbody>
      {score_rows}
    </tbody>
  </table>

  <h2>4. Multi-Model Forward Forecast Horizon (Trained Inverse-WIS Linear Pool)</h2>
  <table>
    <thead>
      <tr>
        <th>Horizon</th>
        <th class="num">Ensemble Median</th>
        <th class="num">50% PI</th>
        <th class="num">95% PI</th>
        <th class="num">SEIR Med</th>
        <th class="num">Renewal Med</th>
        <th class="num">Trend Med</th>
      </tr>
    </thead>
    <tbody>
      {fc_rows}
    </tbody>
  </table>

  <h2>5. Top Multi-Zone Circular Kulldorff Spatial Clusters</h2>
  <table>
    <thead>
      <tr>
        <th>Rank</th>
        <th>Cluster Center</th>
        <th class="num">Radius</th>
        <th class="num">Observed</th>
        <th class="num">Expected</th>
        <th class="num">Rel. Risk</th>
        <th class="num">LLR</th>
        <th class="num">MC p-val</th>
      </tr>
    </thead>
    <tbody>
      {cluster_rows}
    </tbody>
  </table>

  <div class="refs">
    <strong>Methodological References:</strong>
    <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC3816335/">Cori Rt Renewal (PMC3816335)</a> |
    <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC3692796/">Farrington Outbreak GLM (PMC3692796)</a> |
    <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7162546/">NobBS Delay Nowcasting (PMC7162546)</a> |
    <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7880475/">Bracher K=11 WIS (PMC7880475)</a> |
    <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10495027/">Transformed Log-WIS (PMC10495027)</a> |
    <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8963571/">Taylor Inverse-Score Ensembles (PMC8963571)</a> |
    <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5816564/">Kulldorff Spatial Scan (PMC5816564)</a>
  </div>
</div>
</body>
</html>
"#,
        title = html_escape(&ctx.dataset.title),
        vid = html_escape(&ctx.dataset.vintage_id),
        hash = html_escape(&ctx.dataset.checksum_hex),
        mask_badge = mask_badge,
        latest_disp = html_escape(&latest_disp),
        latest_period = latest_period,
        nc_imp = ctx.nowcast.total_imputed_unreported,
        latest_rt = html_escape(&latest_rt),
        total_alerts = ctx.ab_series.total_alerts,
        rel_skill = html_escape(&rel_skill),
        origins = ctx.backtest.origins_evaluated,
        horizon = ctx.backtest.horizon_weeks,
        next_wk_fc = html_escape(&next_wk_fc),
        fig1 = fig1,
        fig2 = fig2,
        fig3 = fig3,
        score_rows = score_rows,
        fc_rows = fc_rows,
        cluster_rows = cluster_rows,
    )
}
