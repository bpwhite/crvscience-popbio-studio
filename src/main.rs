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

use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use clap::{Parser, Subcommand};

use popbio_studio::gui::{render_headless_screenshot, PopBioApp};
use popbio_studio::haplotype_network::{build_haplotype_network, NetworkAlgorithm};
use popbio_studio::model::ModuleId;
use popbio_studio::popgen_alignment::{
    analyze_alignment, parse_sequence_alignment, GapHandling,
};
use popbio_studio::popgen_skyline::{run_skyline_and_expansion, run_sliding_window_scan};
use popbio_studio::workset::WorksetState;

#[derive(Parser, Debug)]
#[command(
    name = "popbio-studio",
    version = "1.2.0",
    about = "Native Rust GUI & CLI educational population biology, theoretical genetics, Bayesian Skyline, and haplotype network simulator"
)]
struct Cli {
    /// Jump directly into a specific module slug (e.g. `v1-hare-lynx`, `m1-ecoli`, `e3-bayesian-skyline`, `e4-sliding-window`)
    #[arg(long)]
    module: Option<String>,

    /// Restore a saved `.popworkset.json` Workset file on startup
    #[arg(long)]
    load_workset: Option<PathBuf>,

    /// Render a 1600x1000 offscreen GPU screenshot to PNG and exit without opening a window
    #[arg(long, value_name = "OUT_PATH")]
    headless_screenshot: Option<PathBuf>,

    /// Select screenshot view (1 = Sewall Wright 3D Landscape, 2 = Bayesian Skyline, 3 = Haplotype Network)
    #[arg(long, default_value_t = 1)]
    screenshot_view: u8,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Run a headless dynamical simulation and write its time-series + growth phases to CSV
    Simulate {
        #[arg(long, default_value = "v1-hare-lynx")]
        module: String,

        #[arg(long, default_value_t = 400)]
        steps: usize,

        #[arg(long, default_value = "simulation_output.csv")]
        out: PathBuf,
    },

    /// Analyze an empirical FASTA / PHYLIP / NEXUS DNA alignment for PopGen neutrality, Bayesian Skyline demographic expansion, sliding-window scans, and Median-Joining Haplotype Networks
    AnalyzeFasta {
        /// Input FASTA, PHYLIP, or NEXUS alignment path
        input: PathBuf,

        /// Bandelt Median-Joining homoplasy tolerance epsilon (`0..=3`)
        #[arg(long, default_value_t = 0)]
        epsilon: usize,

        /// Number of piecewise-constant Bayesian Skyline groups `m`
        #[arg(long, default_value_t = 4)]
        bsp_groups: usize,

        /// Number of MCMC iterations for Bayesian Skyline posterior sampling
        #[arg(long, default_value_t = 2500)]
        bsp_mcmc_steps: usize,

        /// Molecular clock rate (substitutions / site / generation)
        #[arg(long, default_value_t = 1.0e-7)]
        clock_rate: f64,

        /// Sliding-window width in base pairs
        #[arg(long, default_value_t = 24)]
        window_bp: usize,

        /// Sliding-window step size in base pairs
        #[arg(long, default_value_t = 6)]
        step_bp: usize,

        /// Optional JSON output path for PopGen summary statistics
        #[arg(long)]
        out_stats: Option<PathBuf>,

        /// Optional CSV output path for Bayesian Skyline & Generalized Skyline trajectory
        #[arg(long)]
        out_skyline_csv: Option<PathBuf>,

        /// Optional CSV output path for Sliding-Window Neutrality & F_ST scan
        #[arg(long)]
        out_windows_csv: Option<PathBuf>,

        /// Optional SVG output path for rendered Haplotype Network
        #[arg(long)]
        out_svg: Option<PathBuf>,

        /// Optional Gephi GEXF 1.3 output path
        #[arg(long)]
        out_gexf: Option<PathBuf>,
    },
}

#[cfg(windows)]
#[allow(unsafe_code)]
fn attach_parent_console_if_cli() {
    if std::env::args_os().len() > 1 {
        const ATTACH_PARENT_PROCESS: u32 = u32::MAX;
        #[link(name = "kernel32")]
        extern "system" {
            fn AttachConsole(dw_process_id: u32) -> i32;
        }
        // SAFETY: Attaching to the parent console when CLI arguments are passed is optional;
        // Win32 ignores or returns 0 if no parent console exists.
        unsafe {
            let _ = AttachConsole(ATTACH_PARENT_PROCESS);
        }
    }
}

#[cfg(not(windows))]
fn attach_parent_console_if_cli() {}

fn load_window_icon() -> Arc<eframe::egui::IconData> {
    Arc::new(eframe::egui::IconData {
        rgba: include_bytes!(concat!(env!("OUT_DIR"), "/app_icon_64.rgba")).to_vec(),
        width: 64,
        height: 64,
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    attach_parent_console_if_cli();
    let cli = Cli::parse();

    if let Some(out_path) = &cli.headless_screenshot {
        render_headless_screenshot(out_path, cli.screenshot_view)?;
        println!(
            "Saved {} (1600x1000 RGBA offscreen screenshot, view {})",
            out_path.display(),
            cli.screenshot_view
        );
        return Ok(());
    }

    if let Some(cmd) = cli.command {
        match cmd {
            Commands::Simulate { module, steps, out } => {
                let mod_id = ModuleId::from_slug(&module).ok_or_else(|| {
                    format!("Unknown module slug '{module}'. Try e.g. v1-hare-lynx or m1-ecoli.")
                })?;
                let mut ws = WorksetState::default();
                for _ in 0..steps {
                    ws.step_module(mod_id);
                }
                let csv = ws.export_module_csv(mod_id);
                fs::write(&out, csv)?;
                println!(
                    "Simulated {} for {} steps -> {} (final phase: {})",
                    mod_id.slug(),
                    steps,
                    out.display(),
                    ws.current_phase(mod_id).label()
                );
                return Ok(());
            }

            Commands::AnalyzeFasta {
                input,
                epsilon,
                bsp_groups,
                bsp_mcmc_steps,
                clock_rate,
                window_bp,
                step_bp,
                out_stats,
                out_skyline_csv,
                out_windows_csv,
                out_svg,
                out_gexf,
            } => {
                let raw = fs::read_to_string(&input)?;
                let seqs = parse_sequence_alignment(&raw)?;
                let name = input
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("Empirical Alignment");
                let summary = analyze_alignment(name, &seqs, GapHandling::PairwiseDeletion)?;
                let net = build_haplotype_network(
                    &seqs,
                    NetworkAlgorithm::BandeltMedianJoining,
                    epsilon,
                )?;
                let sky = run_skyline_and_expansion(
                    &seqs,
                    &summary.mismatch_histogram,
                    bsp_groups,
                    bsp_mcmc_steps,
                    clock_rate,
                    1.0,
                )?;
                let win_scan = run_sliding_window_scan(
                    &seqs,
                    window_bp,
                    step_bp,
                    GapHandling::PairwiseDeletion,
                )?;

                println!(
                    "Alignment: {} (n = {}, L = {} bp, S = {}, K = {} haplotypes)",
                    summary.dataset_name,
                    summary.num_sequences,
                    summary.alignment_length,
                    summary.segregating_sites_s,
                    summary.num_haplotypes_k
                );
                println!(
                    "Hd = {:.5} | pi = {:.6} | theta_W/site = {:.6} | Tajima's D = {:?}",
                    summary.haplotype_diversity_hd,
                    summary.nucleotide_diversity_pi,
                    summary.watterson_theta_per_site,
                    summary.tajimas_d
                );
                println!(
                    "Skyline: {} | Rogers-Harpending tau = {:.4} (SSD = {:.5}, p_boot = {:.3})",
                    sky.demographic_trend,
                    sky.expansion.tau,
                    sky.expansion.ssd_expansion,
                    sky.expansion.bootstrap_p_ssd
                );

                if let Some(p) = out_stats {
                    fs::write(&p, serde_json::to_string_pretty(&summary)?)?;
                    println!("Wrote PopGen JSON summary -> {}", p.display());
                }
                if let Some(p) = out_skyline_csv {
                    fs::write(&p, sky.to_csv())?;
                    println!("Wrote Bayesian Skyline CSV -> {}", p.display());
                }
                if let Some(p) = out_windows_csv {
                    fs::write(&p, win_scan.to_csv())?;
                    println!("Wrote Sliding-Window Neutrality CSV -> {}", p.display());
                }
                if let Some(p) = out_svg {
                    fs::write(&p, net.to_svg(1024, 720))?;
                    println!("Wrote Median-Joining SVG network -> {}", p.display());
                }
                if let Some(p) = out_gexf {
                    fs::write(&p, net.to_gexf())?;
                    println!("Wrote GEXF 1.3 network -> {}", p.display());
                }
                return Ok(());
            }
        }
    }

    let mut initial_ws = if let Some(path) = &cli.load_workset {
        WorksetState::load_from_file(path)?
    } else {
        WorksetState::default()
    };

    let mut start_in_hub = cli.load_workset.is_none();
    if let Some(slug) = &cli.module {
        let m = ModuleId::from_slug(slug).ok_or_else(|| {
            format!("Unknown module slug '{slug}'. Try e.g. v1-hare-lynx, m1-ecoli, or e2-haplotype-network.")
        })?;
        initial_ws.active_module = m;
        start_in_hub = false;
    }

    let native_options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 880.0])
            .with_min_inner_size([960.0, 640.0])
            .with_icon(load_window_icon())
            .with_title("PopBio Studio v1.2.0 — Population Biology, Demographic Skyline & Haplotype Workstation"),
        ..Default::default()
    };

    eframe::run_native(
        "PopBio Studio v1.2.0",
        native_options,
        Box::new(move |cc| Ok(Box::new(PopBioApp::new(cc, initial_ws, start_in_hub)))),
    )
    .map_err(|e| format!("GUI runtime error: {e}").into())
}

#[cfg(test)]
mod tests {
    use super::load_window_icon;

    #[test]
    fn embedded_window_icon_is_valid_64x64_rgba() {
        let icon = load_window_icon();
        assert_eq!(icon.width, 64);
        assert_eq!(icon.height, 64);
        assert_eq!(icon.rgba.len(), 64 * 64 * 4);
        assert!(icon.rgba.iter().any(|&b| b > 0));
    }
}
