use epi_core::EpiError;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::Path;

/// Pure-Rust GeoJSON (`FeatureCollection` of `Polygon` / `MultiPolygon` / `Point`) spatial parser
/// extracting county, census-tract, or state centroids `(lat_deg, lon_deg)` for Kulldorff spatial scans.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeoFeatureCentroid {
    pub geoid: String,
    pub name: String,
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub cases: Option<f64>,
    pub population: Option<f64>,
}

pub fn parse_geojson_centroids_str(geojson_str: &str) -> Result<Vec<GeoFeatureCentroid>, EpiError> {
    let root: Value = serde_json::from_str(geojson_str)?;
    let features = root
        .get("features")
        .and_then(Value::as_array)
        .ok_or_else(|| {
            EpiError::Schema("GeoJSON missing top-level `features` array".to_string())
        })?;

    let mut out = Vec::with_capacity(features.len());
    for (idx, feat) in features.iter().enumerate() {
        let props = feat.get("properties");
        let geoid = extract_prop_str(props, &["GEOID", "FIPS", "fips", "STATEFP", "id"])
            .or_else(|| feat.get("id").and_then(Value::as_str).map(str::to_string))
            .unwrap_or_else(|| format!("zone_{idx:04}"));
        let name = extract_prop_str(props, &["NAME", "name", "NAMELSAD", "county"])
            .unwrap_or_else(|| geoid.clone());
        let cases = extract_prop_f64(props, &["cases", "count", "observed"]);
        let population = extract_prop_f64(props, &["population", "pop", "POP2025", "ACS_POP"]);

        if let Some(geom) = feat.get("geometry") {
            if let Some((lat_deg, lon_deg)) = compute_geometry_centroid(geom) {
                out.push(GeoFeatureCentroid {
                    geoid,
                    name,
                    lat_deg,
                    lon_deg,
                    cases,
                    population,
                });
            }
        }
    }

    if out.is_empty() {
        return Err(EpiError::Schema(
            "No valid Polygon, MultiPolygon, or Point coordinates found in GeoJSON".to_string(),
        ));
    }
    Ok(out)
}

pub fn parse_geojson_centroids_file<P: AsRef<Path>>(
    path: P,
) -> Result<Vec<GeoFeatureCentroid>, EpiError> {
    let txt = std::fs::read_to_string(path)?;
    parse_geojson_centroids_str(&txt)
}

fn compute_geometry_centroid(geom: &Value) -> Option<(f64, f64)> {
    let gtype = geom.get("type").and_then(Value::as_str)?;
    let coords = geom.get("coordinates")?;

    let mut sum_lon = 0.0;
    let mut sum_lat = 0.0;
    let mut count = 0usize;
    accumulate_lon_lat(gtype, coords, &mut sum_lon, &mut sum_lat, &mut count);

    if count == 0 {
        None
    } else {
        Some((sum_lat / (count as f64), sum_lon / (count as f64)))
    }
}

fn accumulate_lon_lat(
    gtype: &str,
    coords: &Value,
    sum_lon: &mut f64,
    sum_lat: &mut f64,
    count: &mut usize,
) {
    match gtype {
        "Point" => {
            if let Some(arr) = coords.as_array() {
                if arr.len() >= 2 {
                    if let (Some(lon), Some(lat)) = (arr[0].as_f64(), arr[1].as_f64()) {
                        *sum_lon += lon;
                        *sum_lat += lat;
                        *count += 1;
                    }
                }
            }
        }
        "Polygon" => {
            // Outer ring is coordinates[0]
            if let Some(ring) = coords.as_array().and_then(| rings| rings.first()).and_then(Value::as_array) {
                for pt in ring {
                    if let Some(pair) = pt.as_array() {
                        if pair.len() >= 2 {
                            if let (Some(lon), Some(lat)) = (pair[0].as_f64(), pair[1].as_f64()) {
                                *sum_lon += lon;
                                *sum_lat += lat;
                                *count += 1;
                            }
                        }
                    }
                }
            }
        }
        "MultiPolygon" => {
            if let Some(polys) = coords.as_array() {
                for poly in polys {
                    accumulate_lon_lat("Polygon", poly, sum_lon, sum_lat, count);
                }
            }
        }
        _ => {}
    }
}

fn extract_prop_str(props: Option<&Value>, keys: &[&str]) -> Option<String> {
    let obj = props?.as_object()?;
    for &k in keys {
        if let Some(v) = obj.get(k) {
            if let Some(s) = v.as_str() {
                if !s.trim().is_empty() {
                    return Some(s.trim().to_string());
                }
            } else if let Some(n) = v.as_i64() {
                return Some(n.to_string());
            }
        }
    }
    None
}

fn extract_prop_f64(props: Option<&Value>, keys: &[&str]) -> Option<f64> {
    let obj = props?.as_object()?;
    for &k in keys {
        if let Some(v) = obj.get(k) {
            if let Some(n) = v.as_f64() {
                return Some(n);
            }
            if let Some(s) = v.as_str() {
                if let Ok(p) = s.trim().parse::<f64>() {
                    return Some(p);
                }
            }
        }
    }
    None
}
