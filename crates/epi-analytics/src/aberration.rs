use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AberrationMethod {
    FarringtonQuasiPoisson,
    EarsC1,
    EarsC2,
    EarsC3Cusum,
}

impl AberrationMethod {
    pub fn label(&self) -> &'static str {
        match self {
            Self::FarringtonQuasiPoisson => "Farrington Quasi-Poisson",
            Self::EarsC1 => "CDC EARS C1 (0-Wk Guard)",
            Self::EarsC2 => "CDC EARS C2 (2-Wk Guard)",
            Self::EarsC3Cusum => "CDC EARS C3 (3-Wk CUSUM)",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AberrationConfig {
    pub method: AberrationMethod,
    pub z_alpha: f64,
    pub history_weeks: usize,
    pub guard_weeks: usize,
}

impl Default for AberrationConfig {
    fn default() -> Self {
        Self {
            method: AberrationMethod::FarringtonQuasiPoisson,
            z_alpha: 2.58,
            history_weeks: 52,
            guard_weeks: 3,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AberrationPoint {
    pub expected_baseline: Option<f64>,
    pub upper_threshold: Option<f64>,
    pub score_or_ratio: Option<f64>,
    pub is_alert: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AberrationSeries {
    pub points: Vec<AberrationPoint>,
    pub total_alerts: usize,
}

/// Evaluate Farrington Flexible-inspired Quasi-Poisson GLM (PMC3692796) or CDC EARS C1/C2/C3 (PMC10576884).
pub fn compute_aberrations(counts: &[f64], cfg: &AberrationConfig) -> AberrationSeries {
    let n = counts.len();
    let mut points = vec![
        AberrationPoint {
            expected_baseline: None,
            upper_threshold: None,
            score_or_ratio: None,
            is_alert: false,
        };
        n
    ];

    match cfg.method {
        AberrationMethod::FarringtonQuasiPoisson => {
            let guard = cfg.guard_weeks.max(1);
            let hist = cfg.history_weeks.max(8);

            for t in (guard + 5)..n {
                let end_idx = t - guard;
                let start_idx = end_idx.saturating_sub(hist - 1);
                let indices: Vec<usize> = (start_idx..=end_idx).collect();
                if indices.len() < 5 {
                    continue;
                }

                // Boost seasonal phase matches (+/- 4 weeks mod 52)
                let mut weights: Vec<f64> = indices
                    .iter()
                    .map(|&i| {
                        let lag = (t - i) % 52;
                        let dist = lag.min(52 - lag);
                        if dist <= 4 { 2.2 } else { 1.0 }
                    })
                    .collect();

                let mut beta0 = 0.0;
                let mut beta1 = 0.0;

                // Two-step weighted log-linear fit with past-aberration residual downweighting
                for iter in 0..2 {
                    let mut sw = 0.0;
                    let mut sx = 0.0;
                    let mut sy = 0.0;
                    let mut sxx = 0.0;
                    let mut sxy = 0.0;

                    for (&i, &w) in indices.iter().zip(weights.iter()) {
                        let x = (i - start_idx) as f64;
                        let y_log = (counts[i].max(0.0) + 0.5).ln();
                        sw += w;
                        sx += w * x;
                        sy += w * y_log;
                        sxx += w * x * x;
                        sxy += w * x * y_log;
                    }

                    let denom = sw * sxx - sx * sx;
                    if denom.abs() > 1e-8 {
                        beta1 = ((sw * sxy - sx * sy) / denom).clamp(-0.06, 0.06);
                        beta0 = (sy - beta1 * sx) / sw;
                    } else {
                        beta1 = 0.0;
                        beta0 = sy / sw.max(1e-6);
                    }

                    if iter == 0 {
                        for (&i, w) in indices.iter().zip(weights.iter_mut()) {
                            let mu_i = (beta0 + beta1 * ((i - start_idx) as f64)).exp().max(0.25);
                            let pearson_r = (counts[i].max(0.0) - mu_i) / mu_i.sqrt();
                            if pearson_r > 2.58 {
                                *w *= (2.58 / pearson_r).powi(2);
                            }
                        }
                    }
                }

                let mu_hat = (beta0 + beta1 * ((t - start_idx) as f64)).exp().max(0.30);
                let mut disp_sum = 0.0;
                let mut w_sum = 0.0;
                for (&i, &w) in indices.iter().zip(weights.iter()) {
                    let mu_i = (beta0 + beta1 * ((i - start_idx) as f64)).exp().max(0.30);
                    disp_sum += w * (counts[i].max(0.0) - mu_i).powi(2) / mu_i;
                    w_sum += w;
                }
                let phi = (disp_sum / (w_sum - 2.0).max(1.0)).max(1.0);

                // Farrington 2/3-power skewness transformation upper prediction bound
                let inner = (1.0 + (2.0 / 3.0) * cfg.z_alpha * (phi / mu_hat).sqrt()).max(0.0);
                let u_bound = (mu_hat * inner.powf(1.5)).max(1.0);
                let obs = counts[t].max(0.0);
                let is_alert = obs > u_bound && obs >= 2.0;

                points[t] = AberrationPoint {
                    expected_baseline: Some(mu_hat),
                    upper_threshold: Some(u_bound),
                    score_or_ratio: Some(obs / u_bound),
                    is_alert,
                };
            }
        }
        AberrationMethod::EarsC1 | AberrationMethod::EarsC2 | AberrationMethod::EarsC3Cusum => {
            let guard = if cfg.method == AberrationMethod::EarsC1 {
                0
            } else {
                2
            };
            let win = 7usize;
            let mut z_scores = vec![0.0; n];

            for t in (win + guard)..n {
                let slice = &counts[t - guard - win..t - guard];
                let mean = slice.iter().copied().sum::<f64>() / (win as f64);
                let var = slice.iter().map(|&v| (v - mean).powi(2)).sum::<f64>()
                    / ((win - 1).max(1) as f64);
                let sd = var.sqrt().max(0.5);
                let obs = counts[t].max(0.0);
                let z = (obs - mean) / sd;
                z_scores[t] = z;

                if cfg.method == AberrationMethod::EarsC3Cusum {
                    let start_c = t.saturating_sub(2).max(win + guard);
                    let cusum: f64 = z_scores[start_c..=t]
                        .iter()
                        .map(|&zv| (zv - 1.0).max(0.0))
                        .sum();
                    let prior_cusum: f64 = if t > start_c {
                        z_scores[start_c..t]
                            .iter()
                            .map(|&zv| (zv - 1.0).max(0.0))
                            .sum()
                    } else {
                        0.0
                    };
                    let needed_z = (cfg.z_alpha - prior_cusum).max(0.0) + 1.0;
                    let u_bound = mean + needed_z * sd;
                    points[t] = AberrationPoint {
                        expected_baseline: Some(mean),
                        upper_threshold: Some(u_bound),
                        score_or_ratio: Some(cusum),
                        is_alert: cusum > cfg.z_alpha && obs >= 2.0,
                    };
                } else {
                    let u_bound = mean + cfg.z_alpha * sd;
                    points[t] = AberrationPoint {
                        expected_baseline: Some(mean),
                        upper_threshold: Some(u_bound),
                        score_or_ratio: Some(z),
                        is_alert: z > cfg.z_alpha && obs >= 2.0,
                    };
                }
            }
        }
    }

    let total_alerts = points.iter().filter(|p| p.is_alert).count();
    AberrationSeries {
        points,
        total_alerts,
    }
}
