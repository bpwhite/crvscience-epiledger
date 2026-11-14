use epi_core::JurisdictionRecord;
use serde::{Deserialize, Serialize};

/// Multi-jurisdiction Kulldorff Circular Spatial Scan Statistic with Multi-Core `std::thread::scope`
/// Parallel Multinomial Monte Carlo Hypothesis Testing (`p = rank / (M + 1)`, up to `M = 9,999`,
/// deterministic per-replication SplitMix64 seeding so outputs are invariant to thread count).
/// Reference: Kulldorff (1997) & PMC5816564 / PMC2691403.

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatialScanInputZone {
    pub id: String,
    pub name: String,
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub observed_cases: f64,
    pub population: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiZoneCluster {
    pub center_id: String,
    pub center_name: String,
    pub member_ids: Vec<String>,
    pub member_names: Vec<String>,
    pub radius_km: f64,
    pub observed_cases: f64,
    pub expected_cases: f64,
    pub population: f64,
    pub relative_risk: f64,
    pub kulldorff_llr: f64,
    pub monte_carlo_p_value: f64,
}

/// Great-circle Haversine distance in kilometers on WGS-84 mean spherical radius (`R = 6371.0 km`).
pub fn haversine_km(lat1_deg: f64, lon1_deg: f64, lat2_deg: f64, lon2_deg: f64) -> f64 {
    let r_earth = 6371.0_f64;
    let dlat = (lat2_deg - lat1_deg).to_radians();
    let dlon = (lon2_deg - lon1_deg).to_radians();
    let a = (dlat * 0.5).sin().powi(2)
        + lat1_deg.to_radians().cos()
            * lat2_deg.to_radians().cos()
            * (dlon * 0.5).sin().powi(2);
    2.0 * r_earth * a.clamp(0.0, 1.0).sqrt().asin()
}

/// Build `SpatialScanInputZone` slice for all 50 US states + DC using verified geographic centroids.
pub fn build_us_state_scan_zones(
    records: &[JurisdictionRecord],
    use_2025: bool,
) -> Vec<SpatialScanInputZone> {
    records
        .iter()
        .filter_map(|r| {
            let (lat, lon) = us_state_centroid_lat_lon(&r.fips)?;
            let obs = if let Some(custom) = r.active_custom_value {
                custom.max(0.0)
            } else if use_2025 {
                f64::from(r.cases_2025)
            } else {
                f64::from(r.cases_2024)
            };
            Some(SpatialScanInputZone {
                id: r.fips.clone(),
                name: r.name.clone(),
                lat_deg: lat,
                lon_deg: lon,
                observed_cases: obs,
                population: r.census_pop_2025 as f64,
            })
        })
        .collect()
}

fn splitmix64_seed(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    let out = z ^ (z >> 31);
    if out == 0 {
        0xD1B5_4A32_D192_ED03
    } else {
        out
    }
}

/// Evaluate all concentric circular candidate windows up to `max_pop_fraction` (default `0.50`)
/// and assess statistical significance via multi-threaded `std::thread::scope` Monte Carlo simulations
/// (`monte_carlo_reps`, up to `9,999`), returning non-overlapping top clusters ordered by `kulldorff_llr`.
pub fn run_kulldorff_circular_scan(
    zones: &[SpatialScanInputZone],
    max_pop_fraction: f64,
    monte_carlo_reps: usize,
    max_reported_clusters: usize,
) -> Vec<MultiZoneCluster> {
    let n = zones.len();
    if n < 2 {
        return Vec::new();
    }

    let total_pop: f64 = zones.iter().map(|z| z.population.max(0.0)).sum();
    let total_cases: f64 = zones.iter().map(|z| z.observed_cases.max(0.0)).sum();
    if total_pop <= 0.0 || total_cases < 1.0 {
        return Vec::new();
    }

    let pop_cap = total_pop * max_pop_fraction.clamp(0.05, 0.50);

    let mut candidate_windows: Vec<(usize, f64, Vec<usize>, f64)> = Vec::new();
    for i in 0..n {
        let mut neighbors: Vec<(usize, f64)> = (0..n)
            .map(|j| {
                let d = if i == j {
                    0.0
                } else {
                    haversine_km(
                        zones[i].lat_deg,
                        zones[i].lon_deg,
                        zones[j].lat_deg,
                        zones[j].lon_deg,
                    )
                };
                (j, d)
            })
            .collect();
        neighbors.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let mut members = Vec::new();
        let mut cum_pop = 0.0;
        for (j, dist_km) in neighbors {
            if cum_pop + zones[j].population > pop_cap && !members.is_empty() {
                break;
            }
            cum_pop += zones[j].population;
            members.push(j);
            if cum_pop <= pop_cap {
                candidate_windows.push((i, dist_km, members.clone(), cum_pop));
            }
        }
    }

    if candidate_windows.is_empty() {
        return Vec::new();
    }

    let mut scored_windows: Vec<(usize, f64, f64, f64)> = candidate_windows
        .iter()
        .enumerate()
        .map(|(w_idx, (_, _, mem, win_pop))| {
            let obs_z: f64 = mem.iter().map(|&idx| zones[idx].observed_cases.max(0.0)).sum();
            let exp_z = total_cases * (win_pop / total_pop);
            let llr = poisson_llr(obs_z, exp_z, total_cases);
            (w_idx, llr, obs_z, exp_z)
        })
        .collect();

    scored_windows.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let reps = monte_carlo_reps.clamp(99, 9_999);
    let total_cases_int = total_cases.round() as usize;
    let mut cdf_probs = Vec::with_capacity(n);
    let mut acc_p = 0.0;
    for z in zones {
        acc_p += z.population.max(0.0) / total_pop;
        cdf_probs.push(acc_p.min(1.0));
    }
    if let Some(last) = cdf_probs.last_mut() {
        *last = 1.0;
    }

    let base_seed = 0xD1B5_4A32_D192_ED03_u64 ^ (total_cases_int as u64);
    let num_workers = std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(4)
        .clamp(1, 16)
        .min(reps);
    let chunk_size = reps.div_ceil(num_workers);

    // Zero-lock scoped multi-threading across OS cores; deterministic per-rep seed ensures thread-count invariance
    let null_max_llrs: Vec<f64> = std::thread::scope(|s| {
        let mut handles = Vec::with_capacity(num_workers);
        let cdf_ref = &cdf_probs;
        let win_ref = &candidate_windows;

        for worker_id in 0..num_workers {
            let r_start = worker_id * chunk_size;
            let r_end = ((worker_id + 1) * chunk_size).min(reps);
            if r_start >= r_end {
                continue;
            }
            handles.push(s.spawn(move || {
                let mut local_maxes = Vec::with_capacity(r_end - r_start);
                let mut sim_counts = vec![0.0_f64; n];

                for rep_idx in r_start..r_end {
                    let mut rng_state = splitmix64_seed(base_seed ^ (rep_idx as u64));
                    sim_counts.fill(0.0);

                    if total_cases_int <= 5_000 {
                        for _ in 0..total_cases_int {
                            let u = next_u01(&mut rng_state);
                            let cell = cdf_ref.partition_point(|&p| p < u).min(n - 1);
                            sim_counts[cell] += 1.0;
                        }
                    } else {
                        let mut sim_total = 0.0;
                        for (idx, z) in zones.iter().enumerate() {
                            let lam = total_cases * (z.population.max(0.0) / total_pop);
                            let c_draw = sample_poisson(&mut rng_state, lam);
                            sim_counts[idx] = c_draw;
                            sim_total += c_draw;
                        }
                        if sim_total > 0.0 {
                            let scale = total_cases / sim_total;
                            for c in &mut sim_counts {
                                *c *= scale;
                            }
                        }
                    }

                    let mut max_sim_llr = 0.0_f64;
                    for (_, _, mem, win_pop) in win_ref {
                        let obs_z: f64 = mem.iter().map(|&idx| sim_counts[idx]).sum();
                        let exp_z = total_cases * (win_pop / total_pop);
                        let llr = poisson_llr(obs_z, exp_z, total_cases);
                        if llr > max_sim_llr {
                            max_sim_llr = llr;
                        }
                    }
                    local_maxes.push(max_sim_llr);
                }
                local_maxes
            }));
        }

        let mut merged = Vec::with_capacity(reps);
        for h in handles {
            if let Ok(part) = h.join() {
                merged.extend(part);
            }
        }
        merged
    });

    let mut used_zones = vec![false; n];
    let mut results = Vec::new();

    for (w_idx, llr, obs_z, exp_z) in scored_windows {
        if llr <= 0.0 || results.len() >= max_reported_clusters {
            break;
        }
        let (seed_i, radius_km, ref mem, win_pop) = candidate_windows[w_idx];
        if mem.iter().any(|&idx| used_zones[idx]) {
            continue;
        }
        for &idx in mem {
            used_zones[idx] = true;
        }

        let exceed_count = null_max_llrs.iter().filter(|&&v| v >= llr).count();
        let p_val = ((exceed_count + 1) as f64) / ((reps + 1) as f64);

        let out_obs = (total_cases - obs_z).max(1e-9);
        let out_pop = (total_pop - win_pop).max(1.0);
        let rr = ((obs_z / win_pop.max(1.0)) / (out_obs / out_pop)).max(0.0);

        results.push(MultiZoneCluster {
            center_id: zones[seed_i].id.clone(),
            center_name: zones[seed_i].name.clone(),
            member_ids: mem.iter().map(|&idx| zones[idx].id.clone()).collect(),
            member_names: mem.iter().map(|&idx| zones[idx].name.clone()).collect(),
            radius_km,
            observed_cases: obs_z,
            expected_cases: exp_z,
            population: win_pop,
            relative_risk: rr,
            kulldorff_llr: llr,
            monte_carlo_p_value: p_val,
        });
    }

    results
}

fn poisson_llr(c: f64, e: f64, total_c: f64) -> f64 {
    if c <= e || e <= 0.0 || (total_c - e) <= 0.0 {
        return 0.0;
    }
    let out_c = (total_c - c).max(0.0);
    let out_e = total_c - e;
    let term1 = c * (c / e).ln();
    let term2 = if out_c > 0.0 {
        out_c * (out_c / out_e).ln()
    } else {
        0.0
    };
    (term1 + term2).max(0.0)
}

fn next_u01(state: &mut u64) -> f64 {
    let mut x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    let r = x.wrapping_mul(0x2545_F491_4F6C_DD1D);
    ((r >> 11) as f64 / ((1u64 << 53) as f64)).clamp(1e-12, 1.0 - 1e-12)
}

fn sample_poisson(state: &mut u64, lambda: f64) -> f64 {
    if lambda <= 0.0 {
        return 0.0;
    }
    if lambda < 25.0 {
        let l = (-lambda).exp();
        let mut k = 0u32;
        let mut p = 1.0;
        loop {
            k += 1;
            p *= next_u01(state);
            if p <= l || k > 400 {
                break;
            }
        }
        f64::from(k - 1)
    } else {
        let u1 = next_u01(state);
        let u2 = next_u01(state);
        let z = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
        (lambda + lambda.sqrt() * z).round().max(0.0)
    }
}

fn us_state_centroid_lat_lon(fips: &str) -> Option<(f64, f64)> {
    match fips {
        "01" => Some((32.8067, -86.7911)),
        "02" => Some((61.3707, -152.4044)),
        "04" => Some((33.7298, -111.4312)),
        "05" => Some((34.9697, -92.3731)),
        "06" => Some((36.1162, -119.6816)),
        "08" => Some((39.0598, -105.3111)),
        "09" => Some((41.5978, -72.7554)),
        "10" => Some((39.3185, -75.5071)),
        "11" => Some((38.8974, -77.0268)),
        "12" => Some((27.7663, -81.6868)),
        "13" => Some((33.0406, -83.6431)),
        "15" => Some((21.0943, -157.4983)),
        "16" => Some((44.2405, -114.4788)),
        "17" => Some((40.3495, -88.9861)),
        "18" => Some((39.8494, -86.2583)),
        "19" => Some((42.0115, -93.2105)),
        "20" => Some((38.5266, -96.7265)),
        "21" => Some((37.6681, -84.6701)),
        "22" => Some((31.1695, -91.8678)),
        "23" => Some((44.6939, -69.3819)),
        "24" => Some((39.0639, -76.8021)),
        "25" => Some((42.2302, -71.5301)),
        "26" => Some((43.3266, -84.5361)),
        "27" => Some((45.6945, -93.9002)),
        "28" => Some((32.7416, -89.6787)),
        "29" => Some((38.4561, -92.2884)),
        "30" => Some((46.9219, -110.4544)),
        "31" => Some((41.1254, -98.2681)),
        "32" => Some((38.3135, -117.0554)),
        "33" => Some((43.4525, -71.5639)),
        "34" => Some((40.2989, -74.5210)),
        "35" => Some((34.8405, -106.2485)),
        "36" => Some((42.1657, -74.9481)),
        "37" => Some((35.6301, -79.8064)),
        "38" => Some((47.5289, -99.7840)),
        "39" => Some((40.3888, -82.7649)),
        "40" => Some((35.5653, -96.9289)),
        "41" => Some((44.5720, -122.0709)),
        "42" => Some((40.5908, -77.2098)),
        "44" => Some((41.6809, -71.5118)),
        "45" => Some((33.8569, -80.9450)),
        "46" => Some((44.2998, -99.4388)),
        "47" => Some((35.7478, -86.6923)),
        "48" => Some((31.0545, -97.5635)),
        "49" => Some((40.1500, -111.8624)),
        "50" => Some((44.0459, -72.7107)),
        "51" => Some((37.7693, -78.1700)),
        "53" => Some((47.4009, -121.4905)),
        "54" => Some((38.4912, -80.9545)),
        "55" => Some((44.2685, -89.6165)),
        "56" => Some((42.7560, -107.3025)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use epi_core::builtin_pertussis_spatial_records;

    #[test]
    fn test_kulldorff_multizone_parallel_reproducibility() {
        let recs = builtin_pertussis_spatial_records();
        let zones = build_us_state_scan_zones(&recs, true);
        let c1 = run_kulldorff_circular_scan(&zones, 0.25, 299, 5);
        let c2 = run_kulldorff_circular_scan(&zones, 0.25, 299, 5);
        assert!(!c1.is_empty());
        assert_eq!(c1[0].center_id, c2[0].center_id);
        assert!((c1[0].monte_carlo_p_value - c2[0].monte_carlo_p_value).abs() < 1e-12);
        assert!(c1[0].kulldorff_llr > 500.0);
    }
}
