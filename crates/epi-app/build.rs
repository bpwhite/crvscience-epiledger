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

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn blend_pixel(dst: &mut [u8], r: f32, g: f32, b: f32, alpha: f32) {
    let a = alpha.clamp(0.0, 1.0);
    if a <= 0.0 {
        return;
    }
    let da = dst[3] as f32 / 255.0;
    let out_a = a + da * (1.0 - a);
    if out_a <= 1e-6 {
        return;
    }
    let dr = dst[0] as f32;
    let dg = dst[1] as f32;
    let db = dst[2] as f32;
    dst[0] = ((r * a + dr * da * (1.0 - a)) / out_a)
        .round()
        .clamp(0.0, 255.0) as u8;
    dst[1] = ((g * a + dg * da * (1.0 - a)) / out_a)
        .round()
        .clamp(0.0, 255.0) as u8;
    dst[2] = ((b * a + db * da * (1.0 - a)) / out_a)
        .round()
        .clamp(0.0, 255.0) as u8;
    dst[3] = (out_a * 255.0).round().clamp(0.0, 255.0) as u8;
}

fn sd_rounded_box(px: f32, py: f32, hx: f32, hy: f32, r: f32) -> f32 {
    let qx = px.abs() - hx + r;
    let qy = py.abs() - hy + r;
    let out_d = (qx.max(0.0).powi(2) + qy.max(0.0).powi(2)).sqrt();
    let in_d = qx.max(qy).min(0.0);
    out_d + in_d - r
}

fn dist_segment(px: f32, py: f32, ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
    let vx = bx - ax;
    let vy = by - ay;
    let wx = px - ax;
    let wy = py - ay;
    let len2 = vx * vx + vy * vy;
    let t = if len2 > 1e-8 {
        ((wx * vx + wy * vy) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let dx = px - (ax + t * vx);
    let dy = py - (ay + t * vy);
    (dx * dx + dy * dy).sqrt()
}

/// Deterministically synthesize the EpiLedger icon (deep clinical navy badge, cyan border,
/// ledger bars, emerald Rt=1.0 equilibrium line, and bright cyan/white epidemic pulse curve).
fn render_epiledger_icon_rgba(size: u32) -> Vec<u8> {
    let s = size as f32;
    let mut rgba = vec![0u8; (size * size * 4) as usize];
    let aa = 1.15 / s;

    // Epidemic curve polyline points in normalized [0, 1]^2 coordinates
    let curve_pts: [(f32, f32); 9] = [
        (0.14, 0.72),
        (0.24, 0.70),
        (0.34, 0.62),
        (0.43, 0.44),
        (0.52, 0.24),
        (0.60, 0.34),
        (0.70, 0.56),
        (0.79, 0.64),
        (0.87, 0.66),
    ];

    // Ledger bars below the epidemic curve (normalized x_center, top_y)
    let bars: [(f32, f32, bool); 6] = [
        (0.21, 0.74, false),
        (0.32, 0.66, false),
        (0.43, 0.49, true),
        (0.54, 0.31, true),
        (0.65, 0.52, true),
        (0.76, 0.66, false),
    ];

    for y in 0..size {
        for x in 0..size {
            let u = (x as f32 + 0.5) / s;
            let v = (y as f32 + 0.5) / s;
            let idx = ((y * size + x) * 4) as usize;
            let px = &mut rgba[idx..idx + 4];

            // Outer rounded badge
            let d_badge = sd_rounded_box(u - 0.5, v - 0.5, 0.44, 0.44, 0.13);
            let badge_cov = (0.5 - d_badge / aa).clamp(0.0, 1.0);
            if badge_cov <= 0.0 {
                continue;
            }

            // Vertical gradient from deep slate (#1E293B) to clinical navy (#0B1324)
            let bg_r = 28.0 * (1.0 - v) + 11.0 * v;
            let bg_g = 41.0 * (1.0 - v) + 19.0 * v;
            let bg_b = 62.0 * (1.0 - v) + 36.0 * v;
            blend_pixel(px, bg_r, bg_g, bg_b, badge_cov);

            // Subtle ledger horizontal lines inside badge
            let inner_cov = (0.5 - (d_badge + 0.04) / aa).clamp(0.0, 1.0);
            if inner_cov > 0.0 {
                for gy in [0.28f32, 0.44, 0.60, 0.76] {
                    let gd = (v - gy).abs();
                    if gd < 0.014 && (0.14..=0.86).contains(&u) {
                        let line_a = (1.0 - gd / 0.014) * 0.25 * inner_cov;
                        blend_pixel(px, 56.0, 189.0, 248.0, line_a);
                    }
                }

                // Ledger incidence bars
                for &(bx, top_y, is_alert) in &bars {
                    if (u - bx).abs() <= 0.036 && v >= top_y && v <= 0.82 {
                        let edge_x = ((0.036 - (u - bx).abs()) / aa).clamp(0.0, 1.0);
                        let edge_y = ((v - top_y) / aa).min((0.82 - v) / aa).clamp(0.0, 1.0);
                        let bar_a = edge_x * edge_y * 0.68 * inner_cov;
                        if is_alert {
                            blend_pixel(px, 244.0, 63.0, 94.0, bar_a);
                        } else {
                            blend_pixel(px, 14.0, 165.0, 233.0, bar_a);
                        }
                    }
                }

                // Emerald Rt = 1.0 equilibrium baseline
                let eq_d = (v - 0.56).abs();
                if (0.14..=0.86).contains(&u) && eq_d < 0.022 {
                    let eq_a = ((0.016 - eq_d) / aa + 0.5).clamp(0.0, 1.0) * 0.85 * inner_cov;
                    blend_pixel(px, 16.0, 185.0, 129.0, eq_a);
                }

                // Epidemic curve polyline
                let mut min_d = f32::MAX;
                for w in curve_pts.windows(2) {
                    let d = dist_segment(u, v, w[0].0, w[0].1, w[1].0, w[1].1);
                    if d < min_d {
                        min_d = d;
                    }
                }
                let stroke_r = 0.034;
                let glow_r = 0.068;
                if min_d < glow_r {
                    let glow_a = (1.0 - min_d / glow_r) * 0.35 * inner_cov;
                    blend_pixel(px, 56.0, 189.0, 248.0, glow_a);
                }
                let line_cov = ((stroke_r - min_d) / aa + 0.5).clamp(0.0, 1.0) * inner_cov;
                if line_cov > 0.0 {
                    blend_pixel(px, 224.0, 247.0, 255.0, line_cov);
                }

                // Peak outbreak highlight dot at curve_pts[4]
                let peak_d = ((u - curve_pts[4].0).powi(2) + (v - curve_pts[4].1).powi(2)).sqrt();
                let dot_cov = ((0.062 - peak_d) / aa + 0.5).clamp(0.0, 1.0) * inner_cov;
                if dot_cov > 0.0 {
                    blend_pixel(px, 52.0, 211.0, 153.0, dot_cov);
                }

                // Crisp small medical/shield cross accent in top-left corner
                let cx = 0.25;
                let cy = 0.26;
                let h_bar = (u - cx).abs() <= 0.065 && (v - cy).abs() <= 0.021;
                let v_bar = (u - cx).abs() <= 0.021 && (v - cy).abs() <= 0.065;
                if h_bar || v_bar {
                    blend_pixel(px, 56.0, 189.0, 248.0, 0.92 * inner_cov);
                }
            }

            // Glowing cyan border ring
            let ring_d = (d_badge + 0.018).abs();
            let ring_cov = ((0.020 - ring_d) / aa + 0.5).clamp(0.0, 1.0);
            if ring_cov > 0.0 {
                blend_pixel(px, 56.0, 189.0, 248.0, ring_cov * 0.92);
            }
        }
    }

    rgba
}

fn write_ico_file(path: &Path, sizes: &[u32]) -> std::io::Result<()> {
    let mut images: Vec<(u32, Vec<u8>)> = Vec::with_capacity(sizes.len());
    for &sz in sizes {
        images.push((sz, render_epiledger_icon_rgba(sz)));
    }

    let count = images.len() as u16;
    let mut out = Vec::new();
    // ICONDIR
    out.extend_from_slice(&0u16.to_le_bytes()); // Reserved
    out.extend_from_slice(&1u16.to_le_bytes()); // Type = 1 (ICON)
    out.extend_from_slice(&count.to_le_bytes()); // Count

    let mut payload_offset = 6u32 + 16u32 * (images.len() as u32);
    let mut blobs: Vec<Vec<u8>> = Vec::with_capacity(images.len());

    for (sz, rgba) in &images {
        let w = *sz;
        let h = *sz;
        let and_stride = ((w + 31) / 32) * 4;
        let xor_bytes = w * h * 4;
        let and_bytes = and_stride * h;
        let dib_size = 40 + xor_bytes + and_bytes;

        let mut dib = Vec::with_capacity(dib_size as usize);
        // BITMAPINFOHEADER
        dib.extend_from_slice(&40u32.to_le_bytes()); // biSize
        dib.extend_from_slice(&(w as i32).to_le_bytes()); // biWidth
        dib.extend_from_slice(&((h * 2) as i32).to_le_bytes()); // biHeight (XOR + AND)
        dib.extend_from_slice(&1u16.to_le_bytes()); // biPlanes
        dib.extend_from_slice(&32u16.to_le_bytes()); // biBitCount
        dib.extend_from_slice(&0u32.to_le_bytes()); // biCompression = BI_RGB
        dib.extend_from_slice(&(xor_bytes + and_bytes).to_le_bytes()); // biSizeImage
        dib.extend_from_slice(&0i32.to_le_bytes()); // biXPelsPerMeter
        dib.extend_from_slice(&0i32.to_le_bytes()); // biYPelsPerMeter
        dib.extend_from_slice(&0u32.to_le_bytes()); // biClrUsed
        dib.extend_from_slice(&0u32.to_le_bytes()); // biClrImportant

        // Bottom-up BGRA XOR mask
        for row in (0..h).rev() {
            for col in 0..w {
                let idx = ((row * w + col) * 4) as usize;
                dib.push(rgba[idx + 2]); // B
                dib.push(rgba[idx + 1]); // G
                dib.push(rgba[idx]);     // R
                dib.push(rgba[idx + 3]); // A
            }
        }
        // 1-bpp AND mask (0 = use alpha channel, 1 = transparent when alpha == 0)
        for row in (0..h).rev() {
            let mut row_buf = vec![0u8; and_stride as usize];
            for col in 0..w {
                let idx = ((row * w + col) * 4) as usize;
                if rgba[idx + 3] == 0 {
                    row_buf[(col / 8) as usize] |= 1 << (7 - (col % 8));
                }
            }
            dib.extend_from_slice(&row_buf);
        }

        // ICONDIRENTRY
        out.push(if w >= 256 { 0 } else { w as u8 });
        out.push(if h >= 256 { 0 } else { h as u8 });
        out.push(0); // bColorCount
        out.push(0); // bReserved
        out.extend_from_slice(&1u16.to_le_bytes()); // wPlanes
        out.extend_from_slice(&32u16.to_le_bytes()); // wBitCount
        out.extend_from_slice(&(dib.len() as u32).to_le_bytes());
        out.extend_from_slice(&payload_offset.to_le_bytes());

        payload_offset += dib.len() as u32;
        blobs.push(dib);
    }

    for blob in blobs {
        out.extend_from_slice(&blob);
    }
    fs::write(path, out)
}

fn find_rc_exe() -> Option<PathBuf> {
    if let Ok(path_var) = env::var("PATH") {
        for dir in env::split_paths(&path_var) {
            let cand = dir.join("rc.exe");
            if cand.is_file() {
                return Some(cand);
            }
        }
    }
    let kits_root = Path::new(r"C:\Program Files (x86)\Windows Kits\10\bin");
    if kits_root.is_dir() {
        if let Ok(entries) = fs::read_dir(kits_root) {
            let mut versions: Vec<PathBuf> = entries.flatten().map(|e| e.path()).collect();
            versions.sort();
            versions.reverse();
            for ver_dir in versions {
                let cand = ver_dir.join("x64").join("rc.exe");
                if cand.is_file() {
                    return Some(cand);
                }
            }
        }
    }
    None
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=CARGO_PKG_VERSION");

    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR missing"));

    // Always emit 64x64 RGBA icon for ViewportBuilder::with_icon on all platforms
    let icon_rgba = render_epiledger_icon_rgba(64);
    fs::write(out_dir.join("icon_64x64.rgba"), icon_rgba)
        .expect("failed to write icon_64x64.rgba");

    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    if target_os != "windows" {
        return;
    }

    let ico_path = out_dir.join("epiledger.ico");
    write_ico_file(&ico_path, &[16, 32, 48, 64]).expect("failed to write epiledger.ico");

    let maj: u16 = env::var("CARGO_PKG_VERSION_MAJOR")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2);
    let min: u16 = env::var("CARGO_PKG_VERSION_MINOR")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let pat: u16 = env::var("CARGO_PKG_VERSION_PATCH")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let ver_str = format!("{maj}.{min}.{pat}.0");

    let rc_content = format!(
        r#"1 ICON "epiledger.ico"

1 VERSIONINFO
FILEVERSION {maj},{min},{pat},0
PRODUCTVERSION {maj},{min},{pat},0
FILEFLAGSMASK 0x3fL
FILEFLAGS 0x0L
FILEOS 0x40004L
FILETYPE 0x1L
FILESUBTYPE 0x0L
BEGIN
  BLOCK "StringFileInfo"
  BEGIN
    BLOCK "040904B0"
    BEGIN
      VALUE "CompanyName", "CRV Science Network\0"
      VALUE "FileDescription", "EpiLedger — Native Rust Public Health Informatics Workbench\0"
      VALUE "FileVersion", "{ver_str}\0"
      VALUE "InternalName", "epiledger\0"
      VALUE "LegalCopyright", "Copyright CRV Science Network (GPL-3.0-or-later)\0"
      VALUE "OriginalFilename", "epiledger.exe\0"
      VALUE "ProductName", "EpiLedger\0"
      VALUE "ProductVersion", "{ver_str}\0"
    END
  END
  BLOCK "VarFileInfo"
  BEGIN
    VALUE "Translation", 0x0409, 1200
  END
END
"#
    );

    let rc_path = out_dir.join("epiledger.rc");
    fs::write(&rc_path, rc_content).expect("failed to write epiledger.rc");

    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
    if target_env == "gnu" {
        let obj_path = out_dir.join("app_res.o");
        let compilers = ["windres", "x86_64-w64-mingw32-windres", "i686-w64-mingw32-windres"];
        let mut compiled = false;
        for tool in compilers {
            if let Ok(status) = Command::new(tool)
                .arg("-I")
                .arg(&out_dir)
                .arg("-i")
                .arg(&rc_path)
                .arg("-O")
                .arg("coff")
                .arg("-o")
                .arg(&obj_path)
                .status()
            {
                if status.success() {
                    println!("cargo:rustc-link-arg-bin=epiledger={}", obj_path.display());
                    compiled = true;
                    break;
                }
            }
        }
        if !compiled {
            println!("cargo:warning=windres not found; Windows PE icon/VERSIONINFO resource skipped");
        }
    } else if target_env == "msvc" {
        let res_path = out_dir.join("app_res.res");
        if let Some(rc_exe) = find_rc_exe() {
            if let Ok(status) = Command::new(rc_exe)
                .arg("/nologo")
                .arg("/fo")
                .arg(&res_path)
                .arg(&rc_path)
                .status()
            {
                if status.success() {
                    println!("cargo:rustc-link-arg-bin=epiledger={}", res_path.display());
                    return;
                }
            }
        }
        println!("cargo:warning=rc.exe not found; Windows PE icon/VERSIONINFO resource skipped");
    }
}
