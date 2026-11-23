# EpiLedger Scientific Methods & Privacy Guarantees (`epiledger-rs` v2.0.0)

This document summarizes the statistical, cryptographic, and compartmental models implemented across `epi-core`, `epi-store`, and `epi-analytics`.

---

## 1. CSPRNG-Keyed Differential Privacy & zCDP Accounting (`crates/epi-store/src/dp_and_audit.rs`)

### Secret Key Derivation (`build_dp_hmac_key`)
Counter-mode `HMAC-SHA256` samplers (`sample_discrete_laplace_ctr`, `sample_discrete_gaussian_ctr`) derive their key from:
- **Production / non-deterministic mode (`deterministic_seed: None`)**: 32 bytes (256 bits) drawn directly from the OS cryptographic random number generator (`getrandom::getrandom`), never persisted to SQLite or `ingested_at_utc`. Two releases executed within the same wall-clock second produce independent, cryptographically unpredictable noise vectors.
- **Reproducible test mode (`deterministic_seed: Some(seed)`)**: Explicit 64-bit seed appended for repeatable validation.

### Discrete Laplace ($\varepsilon$-DP) & Exact $\rho$-zCDP Conversion
Given L1 sensitivity $\Delta_1$ and privacy parameter $\varepsilon > 0$, each coordinate noise draw $Z \in \mathbb{Z}$ follows the two-sided geometric distribution $\mathbb{P}(Z = k) = \frac{1 - \alpha}{1 + \alpha}\alpha^{|k|}$ with $\alpha = \exp(-\varepsilon / \Delta_1)$. By Bun & Steinke (2016, Proposition 1.4), pure $\varepsilon$-DP satisfies $\rho$-zCDP with:
$$\rho_{\text{Laplace}} = \frac{\varepsilon^2}{2}.$$
Every Discrete Laplace release records `zcdp_rho_spent = 0.5 * eps * eps`, updates `cumulative_rho_after`, and preserves `zcdp_delta = 0.0` in `dp_budget_ledger`.

### Discrete Gaussian ($\rho$-zCDP) & Per-Row `delta` Preservation
`release_discrete_gaussian_zcdp_series` implements Canonne–Kamath–Steinke (2020) Algorithm 3 exact rejection sampling over a Discrete Laplace proposal with $\sigma = \Delta_2 / \sqrt{2\rho}$, recording `zcdp_rho_spent = rho` and preserving each release's target `zcdp_delta` ($\delta \in (0, 1)$) per row. When all prior zCDP rows share the same $\delta$, cumulative $(\varepsilon, \delta)$-DP uses Bun & Steinke Proposition 1.3:
$$\varepsilon(\delta) = \rho + 2\sqrt{\rho \ln(1/\delta)},$$
and falls back to additive basic composition when distinct $\delta$ targets are mixed across releases.

---

## 2. Cori Renewal $R_t$ & Shifted Triangular Serial Interval (`crates/epi-analytics/src/renewal_rt.rs`)

`build_discrete_gamma_si` computes discrete serial interval weights $w_s$ ($s = 1, \dots, S_{\max}$, $w_0 = 0$) from a continuous $\mathrm{Gamma}(\text{shape} = \mu^2/\sigma^2,\, \text{rate} = \mu/\sigma^2)$ distribution using Cori et al. (2013, *Am J Epidemiol*, Web Appendix 1) shifted triangular-tent integration:
$$w_s = \int_{s-1}^{s} (u - (s - 1)) f_\Gamma(u)\,du + \int_{s}^{s+1} (s + 1 - u) f_\Gamma(u)\,du,$$
normalized over $s = 1 \dots S_{\max}$. Integrating over $[s - 1, s + 1]$ properly captures $[0, 0.5]$ sub-half-step generation mass in $w_1$ rather than truncating $[0, 0.5]$. When $\mu_{\text{SI}} < 7$ days on weekly data, `compute_cori_rt` performs mass-conserving daily disaggregation before evaluating the Bayesian Gamma–Poisson posterior with $N_{\min} = 12$ incident cases per window $\tau$ ($\mathrm{CV} \le 0.30$).

---

## 3. Stochastic SEIR / SEIRS & Binomial Case Observation (`crates/epi-analytics/src/seir.rs`)

Deterministic forward projections integrate the SEIR/SEIRS ODE system with 4th-order Runge–Kutta (`RK4`, sub-step $\Delta t = 0.1$ days). Parallel stochastic ensemble trajectories (`std::thread::scope`) integrate half-day ($\tau = 0.5$ d) Poisson tau-leaping compartment transitions $(S \to E \to I \to R \to S)$ across $14$ sub-steps per MMWR week to accumulate weekly true incident infections $N_{\text{inf}}$.

Weekly observed reported cases are sampled strictly via an exact Binomial observation process:
$$Y_w \sim \mathrm{Binomial}(N_{\text{inf}, w},\; p_{\text{rep}}),$$
using exact Bernoulli summation for $N_{\text{inf}, w} \le 128$, Devroye's first waiting-time geometric sampler for $N_{\text{inf}, w} p_{\text{eff}} < 25$, and clamped normal approximation for large counts. Consequently, simulated reported cases strictly satisfy $0 \le Y_w \le N_{\text{inf}, w}$ on every trajectory.

---

## 4. Disclosure Control & Spatial Cluster Detection

- **Small-Cell Confidentiality (`epi-core/src/schema.rs`)**: `SuppressionPolicy::NchsWonderLe9` suppresses all positive counts $c \in (0, 9]$ and derived crude rates (`SUPP (<=9) (SUPP/100k)`) across cartogram tiles, hover tooltips, state inspector cards, `epi_curve_plot`, and Kulldorff LLR tables while leaving true `0`-count jurisdictions unmasked.
- **Spatial Cluster Scan (`epi-analytics/src/spatial_scan_zones.rs`)**: Multi-zone Kulldorff circular Poisson likelihood-ratio scan with Monte Carlo $p$-values bounded by `max_zone_pop_frac`.
