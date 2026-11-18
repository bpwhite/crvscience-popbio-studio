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

fn blend_rgba(dst: &mut [u8], idx: usize, sr: f32, sg: f32, sb: f32, sa: f32) {
    let sa = sa.clamp(0.0, 1.0);
    if sa <= 0.0 {
        return;
    }
    let dr = dst[idx] as f32 / 255.0;
    let dg = dst[idx + 1] as f32 / 255.0;
    let db = dst[idx + 2] as f32 / 255.0;
    let da = dst[idx + 3] as f32 / 255.0;

    let out_a = sa + da * (1.0 - sa);
    if out_a <= 1e-6 {
        dst[idx] = 0;
        dst[idx + 1] = 0;
        dst[idx + 2] = 0;
        dst[idx + 3] = 0;
        return;
    }
    let out_r = (sr * sa + dr * da * (1.0 - sa)) / out_a;
    let out_g = (sg * sa + dg * da * (1.0 - sa)) / out_a;
    let out_b = (sb * sa + db * da * (1.0 - sa)) / out_a;

    dst[idx] = (out_r.clamp(0.0, 1.0) * 255.0).round() as u8;
    dst[idx + 1] = (out_g.clamp(0.0, 1.0) * 255.0).round() as u8;
    dst[idx + 2] = (out_b.clamp(0.0, 1.0) * 255.0).round() as u8;
    dst[idx + 3] = (out_a.clamp(0.0, 1.0) * 255.0).round() as u8;
}

fn dist_to_segment(px: f32, py: f32, ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
    let vx = bx - ax;
    let vy = by - ay;
    let wx = px - ax;
    let wy = py - ay;
    let len2 = vx * vx + vy * vy;
    if len2 <= 1e-8 {
        return ((px - ax).powi(2) + (py - ay).powi(2)).sqrt();
    }
    let t = (wx * vx + wy * vy) / len2;
    let t = t.clamp(0.0, 1.0);
    let cx = ax + t * vx;
    let cy = ay + t * vy;
    ((px - cx).powi(2) + (py - cy).powi(2)).sqrt()
}

fn smoothstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Render a distinct PopBio-Studio RGBA icon (Sewall Wright dual-peak fitness landscape +
/// connected haplotype network triad on a deep navy rounded tile).
fn render_popbio_icon(size: u32) -> Vec<u8> {
    let mut rgba = vec![0u8; (size * size * 4) as usize];
    let s = size as f32;
    let aa = 1.2 / s;

    for y in 0..size {
        for x in 0..size {
            let u = (x as f32 + 0.5) / s;
            let v = (y as f32 + 0.5) / s;
            let idx = ((y * size + x) * 4) as usize;

            // Rounded rectangle SDF
            let half = 0.44;
            let rad = 0.14;
            let qx = (u - 0.5).abs() - (half - rad);
            let qy = (v - 0.5).abs() - (half - rad);
            let out_d = qx.max(0.0).hypot(qy.max(0.0)) + qx.max(qy).min(0.0) - rad;
            let tile_alpha = 1.0 - smoothstep(-aa, aa, out_d);
            if tile_alpha <= 0.0 {
                continue;
            }

            // Deep navy gradient tile background
            let bg_r = 0.05 + 0.04 * (1.0 - v);
            let bg_g = 0.08 + 0.06 * (1.0 - v);
            let bg_b = 0.13 + 0.08 * (1.0 - u * 0.4);
            blend_rgba(&mut rgba, idx, bg_r, bg_g, bg_b, tile_alpha);

            // Subtle cyan/slate border ring
            let border_dist = (out_d + 0.025).abs();
            let border_a = (1.0 - smoothstep(0.0, 0.022 + aa, border_dist)) * 0.75 * tile_alpha;
            blend_rgba(&mut rgba, idx, 0.22, 0.74, 0.97, border_a);

            // Dual-peak Sewall Wright adaptive landscape silhouette in lower half
            let peak1 = 0.20 * (-((u - 0.31) / 0.11).powi(2)).exp();
            let peak2 = 0.34 * (-((u - 0.70) / 0.13).powi(2)).exp();
            let ridge_y = 0.78 - peak1 - peak2;
            if v >= ridge_y - aa && out_d < -0.02 {
                let depth_t = ((v - ridge_y) / 0.36).clamp(0.0, 1.0);
                let ridge_alpha = (1.0 - smoothstep(-aa, aa, ridge_y - v)) * tile_alpha;
                let lr = 0.08 + 0.15 * (1.0 - depth_t);
                let lg = 0.52 + 0.32 * (1.0 - depth_t);
                let lb = 0.65 - 0.18 * depth_t;
                blend_rgba(&mut rgba, idx, lr, lg, lb, ridge_alpha * 0.88);
            }
            let crest_dist = (v - ridge_y).abs();
            if out_d < -0.025 {
                let crest_a = (1.0 - smoothstep(0.012, 0.025 + aa, crest_dist)) * tile_alpha;
                blend_rgba(&mut rgba, idx, 0.20, 0.92, 0.65, crest_a);
            }

            // Haplotype network triad in upper region
            let n0 = (0.28_f32, 0.31_f32);
            let n1 = (0.52_f32, 0.22_f32);
            let n2 = (0.72_f32, 0.36_f32);
            let edge_d = dist_to_segment(u, v, n0.0, n0.1, n1.0, n1.1)
                .min(dist_to_segment(u, v, n1.0, n1.1, n2.0, n2.1))
                .min(dist_to_segment(u, v, n0.0, n0.1, n2.0, n2.1));
            let edge_a = (1.0 - smoothstep(0.014, 0.024 + aa, edge_d)) * 0.95 * tile_alpha;
            blend_rgba(&mut rgba, idx, 0.89, 0.91, 0.94, edge_a);

            let nodes = [
                (n0.0, n0.1, 0.068_f32, (0.22_f32, 0.74_f32, 0.97_f32)),
                (n1.0, n1.1, 0.078_f32, (0.98_f32, 0.80_f32, 0.08_f32)),
                (n2.0, n2.1, 0.064_f32, (0.98_f32, 0.57_f32, 0.24_f32)),
            ];
            for &(nx, ny, nr, (cr, cg, cb)) in &nodes {
                let d = ((u - nx).powi(2) + (v - ny).powi(2)).sqrt();
                let ring_a = (1.0 - smoothstep(nr, nr + 0.014 + aa, d)) * tile_alpha;
                blend_rgba(&mut rgba, idx, 0.97, 0.98, 0.99, ring_a);
                let fill_a = (1.0 - smoothstep(nr - 0.015, nr - 0.015 + aa, d)) * tile_alpha;
                blend_rgba(&mut rgba, idx, cr, cg, cb, fill_a);
            }
        }
    }
    rgba
}

/// Encode multiple RGBA images as a Windows multi-size `.ico` (32-bit BGRA BMP DIBs + 1bpp AND-masks).
fn encode_multi_ico(sizes: &[u32]) -> Vec<u8> {
    let mut images: Vec<(u32, Vec<u8>)> = Vec::with_capacity(sizes.len());
    for &sz in sizes {
        let rgba = render_popbio_icon(sz);
        let mut dib = Vec::new();
        // BITMAPINFOHEADER (40 bytes)
        dib.extend_from_slice(&40u32.to_le_bytes()); // biSize
        dib.extend_from_slice(&(sz as i32).to_le_bytes()); // biWidth
        dib.extend_from_slice(&((sz * 2) as i32).to_le_bytes()); // biHeight (XOR + AND)
        dib.extend_from_slice(&1u16.to_le_bytes()); // biPlanes
        dib.extend_from_slice(&32u16.to_le_bytes()); // biBitCount
        dib.extend_from_slice(&0u32.to_le_bytes()); // biCompression = BI_RGB
        dib.extend_from_slice(&(sz * sz * 4).to_le_bytes()); // biSizeImage
        dib.extend_from_slice(&0i32.to_le_bytes()); // biXPelsPerMeter
        dib.extend_from_slice(&0i32.to_le_bytes()); // biYPelsPerMeter
        dib.extend_from_slice(&0u32.to_le_bytes()); // biClrUsed
        dib.extend_from_slice(&0u32.to_le_bytes()); // biClrImportant

        // Bottom-up BGRA XOR mask
        for y in (0..sz).rev() {
            for x in 0..sz {
                let idx = ((y * sz + x) * 4) as usize;
                dib.push(rgba[idx + 2]); // B
                dib.push(rgba[idx + 1]); // G
                dib.push(rgba[idx]);     // R
                dib.push(rgba[idx + 3]); // A
            }
        }

        // Bottom-up 1bpp AND mask, rows padded to 4-byte boundary
        let and_row_bytes = (((sz + 31) / 32) * 4) as usize;
        for y in (0..sz).rev() {
            let mut row = vec![0u8; and_row_bytes];
            for x in 0..sz {
                let idx = ((y * sz + x) * 4) as usize;
                if rgba[idx + 3] < 128 {
                    row[(x / 8) as usize] |= 1 << (7 - (x % 8));
                }
            }
            dib.extend_from_slice(&row);
        }
        images.push((sz, dib));
    }

    let count = images.len() as u16;
    let mut ico = Vec::new();
    // ICONDIR (6 bytes)
    ico.extend_from_slice(&0u16.to_le_bytes()); // idReserved
    ico.extend_from_slice(&1u16.to_le_bytes()); // idType = 1 (ICO)
    ico.extend_from_slice(&count.to_le_bytes()); // idCount

    let mut offset = 6u32 + (count as u32) * 16u32;
    for (sz, dib) in &images {
        let dim_byte = if *sz >= 256 { 0u8 } else { *sz as u8 };
        ico.push(dim_byte); // bWidth
        ico.push(dim_byte); // bHeight
        ico.push(0u8); // bColorCount
        ico.push(0u8); // bReserved
        ico.extend_from_slice(&1u16.to_le_bytes()); // wPlanes
        ico.extend_from_slice(&32u16.to_le_bytes()); // wBitCount
        ico.extend_from_slice(&(dib.len() as u32).to_le_bytes()); // dwBytesInRes
        ico.extend_from_slice(&offset.to_le_bytes()); // dwImageOffset
        offset += dib.len() as u32;
    }

    for (_, dib) in images {
        ico.extend_from_slice(&dib);
    }
    ico
}

fn find_rc_exe() -> Option<PathBuf> {
    if let Ok(status) = Command::new("rc.exe").arg("/?").output() {
        let _ = status;
        return Some(PathBuf::from("rc.exe"));
    }
    let kits_roots = [
        env::var_os("ProgramFiles(x86)").map(PathBuf::from),
        env::var_os("ProgramFiles").map(PathBuf::from),
    ];
    for root_opt in kits_roots.into_iter().flatten() {
        let bin_dir = root_opt.join("Windows Kits").join("10").join("bin");
        if let Ok(entries) = fs::read_dir(&bin_dir) {
            let mut versions: Vec<PathBuf> = entries.flatten().map(|e| e.path()).collect();
            versions.sort();
            versions.reverse();
            for ver_dir in versions {
                let candidate = ver_dir.join("x64").join("rc.exe");
                if candidate.is_file() {
                    return Some(candidate);
                }
            }
        }
    }
    None
}

fn find_windres() -> Option<&'static str> {
    for tool in ["x86_64-w64-mingw32-windres", "windres"] {
        if Command::new(tool).arg("--version").output().is_ok() {
            return Some(tool);
        }
    }
    None
}

fn compile_windows_resources(out_dir: &Path) {
    let rc_path = out_dir.join("app_res.rc");
    let rc_content = r#"1 ICON "popbio_studio.ico"
1 VERSIONINFO
FILEVERSION 1,2,0,0
PRODUCTVERSION 1,2,0,0
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
      VALUE "CompanyName", "CRV Science Network"
      VALUE "FileDescription", "PopBio-Studio: Population Biology, Adaptive Landscapes, Matrix Models, Bayesian Skyline & Haplotype Network Simulator"
      VALUE "FileVersion", "1.2.0"
      VALUE "InternalName", "popbio-studio"
      VALUE "LegalCopyright", "Copyright CRV Science Network (GPL-3.0-or-later)"
      VALUE "OriginalFilename", "popbio-studio.exe"
      VALUE "ProductName", "PopBio Studio"
      VALUE "ProductVersion", "1.2.0"
    END
  END
  BLOCK "VarFileInfo"
  BEGIN
    VALUE "Translation", 0x409, 1200
  END
END
"#;
    fs::write(&rc_path, rc_content).expect("write app_res.rc");

    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
    if target_env == "gnu" {
        if let Some(windres) = find_windres() {
            let obj_path = out_dir.join("app_res.o");
            let status = Command::new(windres)
                .current_dir(out_dir)
                .args(["-i", "app_res.rc", "-o", "app_res.o", "-O", "coff"])
                .status();
            match status {
                Ok(s) if s.success() => {
                    println!("cargo:rustc-link-arg-bin=popbio-studio={}", obj_path.display());
                }
                Ok(s) => {
                    println!("cargo:warning={windres} exited with status {s}; skipping PE resource embed");
                }
                Err(e) => {
                    println!("cargo:warning=Failed to run {windres}: {e}; skipping PE resource embed");
                }
            }
        } else {
            println!("cargo:warning=windres not found on PATH; skipping PE icon/VERSIONINFO embed");
        }
    } else if target_env == "msvc" {
        if let Some(rc_exe) = find_rc_exe() {
            let res_path = out_dir.join("app_res.res");
            let status = Command::new(&rc_exe)
                .current_dir(out_dir)
                .args(["/nologo", "/fo", "app_res.res", "app_res.rc"])
                .status();
            match status {
                Ok(s) if s.success() => {
                    println!("cargo:rustc-link-arg-bin=popbio-studio={}", res_path.display());
                }
                Ok(s) => {
                    println!("cargo:warning=rc.exe exited with status {s}; skipping PE resource embed");
                }
                Err(e) => {
                    println!("cargo:warning=Failed to run rc.exe: {e}; skipping PE resource embed");
                }
            }
        } else {
            println!("cargo:warning=rc.exe not found; skipping PE icon/VERSIONINFO embed");
        }
    }
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=Cargo.toml");

    let out_dir = PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR"));
    let rgba_64 = render_popbio_icon(64);
    fs::write(out_dir.join("app_icon_64.rgba"), &rgba_64).expect("write app_icon_64.rgba");

    let ico_bytes = encode_multi_ico(&[16, 32, 48, 64]);
    fs::write(out_dir.join("popbio_studio.ico"), &ico_bytes).expect("write popbio_studio.ico");

    if env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("windows") {
        compile_windows_resources(&out_dir);
    }
}
