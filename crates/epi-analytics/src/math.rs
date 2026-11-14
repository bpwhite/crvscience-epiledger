/// Pure-Rust special functions and non-parametric rank correlation helpers.

pub fn ln_gamma(z: f64) -> f64 {
    const COEFFS: [f64; 6] = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.001208650973866179,
        -0.000005395239384953,
    ];
    let x = z;
    let mut y = z;
    let mut tmp = x + 5.5;
    tmp -= (x + 0.5) * tmp.ln();
    let mut ser = 1.000000000190015;
    for &c in &COEFFS {
        y += 1.0;
        ser += c / y;
    }
    -tmp + (2.5066282746310005 * ser / x).ln()
}

pub fn gamma_pdf(x: f64, shape: f64, rate: f64) -> f64 {
    if x <= 0.0 || shape <= 0.0 || rate <= 0.0 {
        return 0.0;
    }
    let ln_p = shape * rate.ln() + (shape - 1.0) * x.ln() - rate * x - ln_gamma(shape);
    ln_p.exp()
}

/// Lower regularized incomplete gamma function P(a, x) via series (x < a + 1) or Lentz continued fraction.
pub fn regularized_gamma_p(a: f64, x: f64) -> f64 {
    if x <= 0.0 || a <= 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        let mut ap = a;
        let mut sum = 1.0 / a;
        let mut del = sum;
        for _ in 0..200 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-12 {
                break;
            }
        }
        (sum * (-x + a * x.ln() - ln_gamma(a)).exp()).clamp(0.0, 1.0)
    } else {
        (1.0 - regularized_gamma_q_cf(a, x)).clamp(0.0, 1.0)
    }
}

fn regularized_gamma_q_cf(a: f64, x: f64) -> f64 {
    let fp_min = 1e-30;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / fp_min;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..=200 {
        let i_f = f64::from(i);
        let an = -i_f * (i_f - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < fp_min {
            d = fp_min;
        }
        c = b + an / c;
        if c.abs() < fp_min {
            c = fp_min;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-12 {
            break;
        }
    }
    ((-x + a * x.ln() - ln_gamma(a)).exp() * h).clamp(0.0, 1.0)
}

/// Exact Gamma quantile refining Wilson-Hilferty initial guess with 4 Newton-Raphson iterations.
pub fn gamma_quantile(p: f64, shape: f64, rate: f64) -> f64 {
    if p <= 0.0 || shape <= 0.0 || rate <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    let z = normal_quantile(p);
    let v = 1.0 / (9.0 * shape);
    let wh_root = (1.0 - v + z * v.sqrt()).max(0.02);
    let mut x = shape * wh_root.powi(3);

    // Refine on unit-rate Gamma(shape, 1)
    for _ in 0..4 {
        let cdf = regularized_gamma_p(shape, x);
        let pdf = gamma_pdf(x, shape, 1.0);
        if pdf > 1e-18 {
            let step = (cdf - p) / pdf;
            x = (x - step).max(1e-8);
            if step.abs() < 1e-8 * x {
                break;
            }
        }
    }
    x / rate
}

/// Acklam rational approximation for the standard normal quantile Phi^-1(p).
pub fn normal_quantile(p: f64) -> f64 {
    let p_clamped = p.clamp(1e-12, 1.0 - 1e-12);
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];

    let p_low = 0.02425;
    let p_high = 1.0 - p_low;

    if p_clamped < p_low {
        let q = (-2.0 * p_clamped.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p_clamped <= p_high {
        let q = p_clamped - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p_clamped).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

/// Spearman rank correlation rho in [-1, 1] with average-rank tie handling.
pub fn spearman_rho(xs: &[f64], ys: &[f64]) -> Option<f64> {
    let n = xs.len().min(ys.len());
    if n < 3 {
        return None;
    }
    let rx = average_ranks(&xs[..n]);
    let ry = average_ranks(&ys[..n]);
    pearson(&rx, &ry)
}

fn average_ranks(vals: &[f64]) -> Vec<f64> {
    let mut indexed: Vec<(usize, f64)> = vals.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![0.0; vals.len()];
    let mut i = 0;
    while i < indexed.len() {
        let mut j = i + 1;
        while j < indexed.len() && (indexed[j].1 - indexed[i].1).abs() < 1e-12 {
            j += 1;
        }
        let avg_rank = (i + 1 + j) as f64 / 2.0;
        for k in i..j {
            ranks[indexed[k].0] = avg_rank;
        }
        i = j;
    }
    ranks
}

fn pearson(xs: &[f64], ys: &[f64]) -> Option<f64> {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut dx2 = 0.0;
    let mut dy2 = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let dx = x - mx;
        let dy = y - my;
        num += dx * dy;
        dx2 += dx * dx;
        dy2 += dy * dy;
    }
    let denom = (dx2 * dy2).sqrt();
    if denom < 1e-12 {
        None
    } else {
        Some((num / denom).clamp(-1.0, 1.0))
    }
}
