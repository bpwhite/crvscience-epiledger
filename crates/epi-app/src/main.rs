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

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
#![deny(unsafe_code)]

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

use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;

const APP_ICON_RGBA: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/icon_64x64.rgba"));

#[cfg(windows)]
#[allow(unsafe_code)]
fn attach_parent_console_if_cli() {
    if std::env::args_os().len() > 1 {
        const ATTACH_PARENT_PROCESS: u32 = u32::MAX;
        #[link(name = "kernel32")]
        extern "system" {
            fn AttachConsole(dw_process_id: u32) -> i32;
        }
        unsafe {
            let _ = AttachConsole(ATTACH_PARENT_PROCESS);
        }
    }
}

#[cfg(not(windows))]
fn attach_parent_console_if_cli() {}

fn is_help_flag(arg: &str) -> bool {
    matches!(arg, "--help" | "-h" | "help")
}

fn is_version_flag(arg: &str) -> bool {
    matches!(arg, "--version" | "-V")
}

fn parse_screenshot_arg(
    args: &[String],
) -> Result<Option<(PathBuf, u8)>, Box<dyn std::error::Error>> {
    let mut out_path: Option<PathBuf> = None;
    let mut view: u8 = 1;
    let mut i = 1;
    while i < args.len() {
        let arg = &args[i];
        if arg == "--screenshot" {
            let next = args
                .get(i + 1)
                .ok_or("--screenshot requires an output PNG path argument")?;
            out_path = Some(PathBuf::from(next));
            i += 2;
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--screenshot=") {
            if rest.is_empty() {
                return Err("--screenshot= requires a non-empty output PNG path".into());
            }
            out_path = Some(PathBuf::from(rest));
            i += 1;
            continue;
        }
        if arg == "--screenshot-view" {
            let next = args
                .get(i + 1)
                .ok_or("--screenshot-view requires a view number (1, 2, or 3)")?;
            view = next
                .parse::<u8>()
                .map_err(|e| format!("invalid --screenshot-view '{next}': {e}"))?;
            i += 2;
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--screenshot-view=") {
            view = rest
                .parse::<u8>()
                .map_err(|e| format!("invalid --screenshot-view '{rest}': {e}"))?;
            i += 1;
            continue;
        }
        i += 1;
    }
    Ok(out_path.map(|p| (p, view)))
}

fn block_on<F: std::future::Future>(future: F) -> F::Output {
    use std::sync::Arc;
    use std::task::{Context, Poll, Wake, Waker};

    struct ThreadWaker(std::thread::Thread);
    impl Wake for ThreadWaker {
        fn wake(self: Arc<Self>) {
            self.0.unpark();
        }
        fn wake_by_ref(self: &Arc<Self>) {
            self.0.unpark();
        }
    }

    let waker = Waker::from(Arc::new(ThreadWaker(std::thread::current())));
    let mut cx = Context::from_waker(&waker);
    let mut future = std::pin::pin!(future);
    loop {
        match future.as_mut().poll(&mut cx) {
            Poll::Ready(out) => return out,
            Poll::Pending => std::thread::park(),
        }
    }
}

fn write_rgba_png(
    path: &std::path::Path,
    width: u32,
    height: u32,
    rgba: &[u8],
) -> std::io::Result<()> {
    fn crc32(data: &[u8]) -> u32 {
        let mut crc = 0xFFFF_FFFFu32;
        for &b in data {
            crc ^= b as u32;
            for _ in 0..8 {
                crc = if (crc & 1) != 0 { (crc >> 1) ^ 0xEDB8_8320 } else { crc >> 1 };
            }
        }
        !crc
    }
    fn push_chunk(out: &mut Vec<u8>, tag: &[u8; 4], payload: &[u8]) {
        out.extend_from_slice(&(payload.len() as u32).to_be_bytes());
        let start = out.len();
        out.extend_from_slice(tag);
        out.extend_from_slice(payload);
        out.extend_from_slice(&crc32(&out[start..]).to_be_bytes());
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let row_bytes = width as usize * 4;
    let mut raw = Vec::with_capacity(height as usize * (1 + row_bytes));
    for row in rgba.chunks_exact(row_bytes) {
        // PNG Sub filter (type 1): subtract left RGBA pixel (bpp = 4) for high DEFLATE compression
        raw.push(1);
        for (x, &b) in row.iter().enumerate() {
            let left = if x >= 4 { row[x - 4] } else { 0 };
            raw.push(b.wrapping_sub(left));
        }
    }
    let mut encoder = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(&raw)?;
    let zlib = encoder.finish()?;
    let mut ihdr = Vec::with_capacity(13);
    ihdr.extend_from_slice(&width.to_be_bytes());
    ihdr.extend_from_slice(&height.to_be_bytes());
    ihdr.extend_from_slice(&[8, 6, 0, 0, 0]);
    let mut png = Vec::with_capacity(zlib.len() + 64);
    png.extend_from_slice(b"\x89PNG\r\n\x1a\n");
    push_chunk(&mut png, b"IHDR", &ihdr);
    push_chunk(&mut png, b"IDAT", &zlib);
    push_chunk(&mut png, b"IEND", &[]);
    std::fs::write(path, png)
}

fn capture_headless_screenshot(
    out_png: &std::path::Path,
    view: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    use eframe::{egui_wgpu, wgpu};

    const W: u32 = 1600;
    const H: u32 = 1000;

    let temp_db = std::env::temp_dir().join(format!(
        "epiledger_screenshot_{}.sqlite3",
        std::process::id()
    ));

    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends: wgpu::Backends::all(),
        ..Default::default()
    });
    let adapter = block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        compatible_surface: None,
        force_fallback_adapter: false,
    }))
    .ok_or("No suitable wgpu adapter found for offscreen capture")?;

    let (device, queue) = block_on(adapter.request_device(
        &wgpu::DeviceDescriptor {
            label: Some("epiledger_headless_device"),
            required_features: wgpu::Features::empty(),
            required_limits: wgpu::Limits::default(),
        },
        None,
    ))?;

    let format = wgpu::TextureFormat::Rgba8Unorm;
    let mut renderer = egui_wgpu::Renderer::new(&device, format, None, 1);

    let tex_extent = wgpu::Extent3d {
        width: W,
        height: H,
        depth_or_array_layers: 1,
    };
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("epiledger_offscreen_target"),
        size: tex_extent,
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let tex_view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    let bytes_per_row = W * 4;
    let readback_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("epiledger_offscreen_readback"),
        size: (bytes_per_row * H) as u64,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let ctx = egui::Context::default();
    ctx.set_pixels_per_point(1.0);
    let mut app = gui::EpiLedgerGuiApp::new_headless(&ctx, temp_db.clone());
    app.configure_screenshot_view(view);

    for frame_idx in 0..4u32 {
        let mut raw_input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(W as f32, H as f32),
            )),
            time: Some(frame_idx as f64 * 0.016),
            predicted_dt: 1.0 / 60.0,
            ..Default::default()
        };
        raw_input
            .viewports
            .entry(egui::ViewportId::ROOT)
            .or_default()
            .native_pixels_per_point = Some(1.0);

        let full_output = ctx.run(raw_input, |ctx| app.update_ui(ctx));
        let clipped_primitives =
            ctx.tessellate(full_output.shapes, full_output.pixels_per_point);
        let screen_descriptor = egui_wgpu::ScreenDescriptor {
            size_in_pixels: [W, H],
            pixels_per_point: 1.0,
        };

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("epiledger_offscreen_encoder"),
        });
        for (id, image_delta) in &full_output.textures_delta.set {
            renderer.update_texture(&device, &queue, *id, image_delta);
        }
        let user_cmds = renderer.update_buffers(
            &device,
            &queue,
            &mut encoder,
            &clipped_primitives,
            &screen_descriptor,
        );

        let bg = ctx.style().visuals.panel_fill.to_normalized_gamma_f32();
        {
            let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("epiledger_offscreen_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &tex_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: bg[0] as f64,
                            g: bg[1] as f64,
                            b: bg[2] as f64,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                occlusion_query_set: None,
                timestamp_writes: None,
            });
            renderer.render(&mut rpass, &clipped_primitives, &screen_descriptor);
        }

        if frame_idx == 3 {
            encoder.copy_texture_to_buffer(
                texture.as_image_copy(),
                wgpu::ImageCopyBuffer {
                    buffer: &readback_buf,
                    layout: wgpu::ImageDataLayout {
                        offset: 0,
                        bytes_per_row: Some(bytes_per_row),
                        rows_per_image: Some(H),
                    },
                },
                tex_extent,
            );
        }

        queue.submit(user_cmds.into_iter().chain(std::iter::once(encoder.finish())));
        for id in &full_output.textures_delta.free {
            renderer.free_texture(id);
        }
    }

    let slice = readback_buf.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |res| {
        let _ = tx.send(res);
    });
    let _ = device.poll(wgpu::Maintain::Wait);
    rx.recv()??;
    let rgba = slice.get_mapped_range().to_vec();
    readback_buf.unmap();

    drop(app);
    for suffix in ["", "-wal", "-shm"] {
        let p = PathBuf::from(format!("{}{}", temp_db.display(), suffix));
        let _ = std::fs::remove_file(p);
    }

    write_rgba_png(out_png, W, H, &rgba)?;
    println!(
        "Saved offscreen GPU screenshot ({}x{} RGBA) to {}",
        W,
        H,
        out_png.display()
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    attach_parent_console_if_cli();
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).is_some_and(|a| is_help_flag(a)) {
        cli::print_cli_help();
        return Ok(());
    }
    if args.get(1).is_some_and(|a| is_version_flag(a)) {
        println!("epiledger {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    if let Some((out_png, view)) = parse_screenshot_arg(&args)? {
        return capture_headless_screenshot(&out_png, view);
    }
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
            .with_title("EpiLedger — Native Rust Public Health Informatics Workbench")
            .with_icon(Arc::new(egui::IconData {
                rgba: APP_ICON_RGBA.to_vec(),
                width: 64,
                height: 64,
            })),
        ..Default::default()
    };

    eframe::run_native(
        "EpiLedger Workbench",
        native_options,
        Box::new(move |cc| Ok(Box::new(gui::EpiLedgerGuiApp::new(cc, db_path)))),
    )
    .map_err(|e| format!("eframe launch error: {e}").into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_rgba_png_deflate_under_1mb_and_valid_ihdr() {
        let width = 1600u32;
        let height = 1000u32;
        let mut rgba = vec![0u8; (width * height * 4) as usize];
        for y in 0..height {
            for x in 0..width {
                let idx = ((y * width + x) * 4) as usize;
                // UI-like dark background with periodic chart bands
                if (200..600).contains(&y) && (x % 12 < 8) {
                    rgba[idx] = 34;
                    rgba[idx + 1] = 211;
                    rgba[idx + 2] = 238;
                } else {
                    rgba[idx] = 18;
                    rgba[idx + 1] = 22;
                    rgba[idx + 2] = 30;
                }
                rgba[idx + 3] = 255;
            }
        }
        let tmp = std::env::temp_dir().join(format!(
            "epiledger_png_test_{}.png",
            std::process::id()
        ));
        write_rgba_png(&tmp, width, height, &rgba).expect("write PNG");
        let bytes = std::fs::read(&tmp).expect("read PNG");
        let _ = std::fs::remove_file(&tmp);

        assert!(bytes.len() < 1_000_000, "PNG size {} >= 1 MB", bytes.len());
        assert_eq!(&bytes[0..8], b"\x89PNG\r\n\x1a\n");
        assert_eq!(&bytes[12..16], b"IHDR");
        let w = u32::from_be_bytes(bytes[16..20].try_into().unwrap());
        let h = u32::from_be_bytes(bytes[20..24].try_into().unwrap());
        assert_eq!((w, h), (width, height));
        assert_eq!(bytes[24], 8); // 8-bit depth
        assert_eq!(bytes[25], 6); // RGBA color type
    }

    #[test]
    fn test_cli_help_and_screenshot_arg_parsing() {
        assert!(is_help_flag("--help"));
        assert!(is_help_flag("-h"));
        assert!(is_help_flag("help"));
        assert!(!is_help_flag("--screenshot"));
        assert!(is_version_flag("--version"));
        assert!(is_version_flag("-V"));

        let parsed = parse_screenshot_arg(&[
            "epiledger".to_string(),
            "--screenshot=out.png".to_string(),
            "--screenshot-view".to_string(),
            "3".to_string(),
        ])
        .unwrap();
        assert_eq!(parsed, Some((PathBuf::from("out.png"), 3)));

        let tmp_db = std::env::temp_dir().join(format!(
            "epiledger_no_db_on_help_{}.sqlite3",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&tmp_db);
        cli::run_cli(&["--help".to_string()], &tmp_db).expect("--help exits Ok");
        assert!(
            !tmp_db.exists(),
            "--help must not create a SQLite database file"
        );
    }
}
