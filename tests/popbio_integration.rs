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

#![forbid(unsafe_code)]

use popbio_studio::gui::{AppScreen, PopBioApp, WorkspaceTab};
use popbio_studio::haplotype_network::{build_haplotype_network, NetworkAlgorithm};
use popbio_studio::model::{GrowthPhase, ModuleId};
use popbio_studio::popgen_alignment::{
    analyze_alignment, builtin_alignment_presets, parse_sequence_alignment, GapHandling,
};
use popbio_studio::popgen_skyline::{run_skyline_and_expansion, run_sliding_window_scan};
use popbio_studio::sim_landscapes_matrices::{LeslieLefkovitchSim, MatrixPreset, WrightLandscapeSim};
use popbio_studio::sim_vertebrates::HareLynxSim;
use popbio_studio::workset::{WorksetDb, WorksetState};

#[test]
fn all_dynamical_modules_step_and_export_csv() {
    let mut ws = WorksetState::default();
    for mod_id in ModuleId::ALL {
        for _ in 0..25 {
            ws.step_module(mod_id);
        }
        if !WorksetState::is_empirical_module(mod_id) {
            let csv = ws.export_module_csv(mod_id);
            assert!(
                csv.lines().count() >= 25,
                "Expected >= 25 CSV lines for {}",
                mod_id.slug()
            );
        }
    }
}

#[test]
fn ecoli_batch_transitions_across_growth_phases() {
    let mut ws = WorksetState::default();
    assert_eq!(
        ws.current_phase(ModuleId::M1EcoliBatchChemostat),
        GrowthPhase::LagPhase
    );
    let steps1 = ws.advance_to_next_growth_phase(ModuleId::M1EcoliBatchChemostat, 300);
    assert!(steps1 > 0);
    assert_eq!(
        ws.current_phase(ModuleId::M1EcoliBatchChemostat),
        GrowthPhase::ExponentialLog
    );
}

#[test]
fn empirical_alignment_popgen_and_median_joining_network() {
    for (title, raw_fasta) in builtin_alignment_presets() {
        let seqs = parse_sequence_alignment(raw_fasta).expect("parse builtin FASTA");
        let sum = analyze_alignment(title, &seqs, GapHandling::CompleteDeletion)
            .expect("compute PopGen summary");
        assert!(sum.segregating_sites_s >= 4);
        assert!(sum.haplotype_diversity_hd > 0.3 && sum.haplotype_diversity_hd <= 1.0);
        assert!(sum.nucleotide_diversity_pi > 0.0);
        assert!(sum.tajimas_d.is_some());

        let net = build_haplotype_network(&seqs, NetworkAlgorithm::BandeltMedianJoining, 0)
            .expect("build Median-Joining network");
        assert!(net.nodes.len() >= sum.num_haplotypes_k);
        assert!(!net.edges.is_empty());
        let svg = net.to_svg(800, 600);
        assert!(svg.contains("<svg") && svg.contains("</svg>"));
    }

    let star_fasta = "\
        >T1|PopA\nCAAAAA\n\
        >T2|PopB\nACAAAA\n\
        >T3|PopC\nAACAAA\n";
    let star_seqs = parse_sequence_alignment(star_fasta).expect("star parse");
    let mj_net = build_haplotype_network(&star_seqs, NetworkAlgorithm::BandeltMedianJoining, 0)
        .expect("star MJ");
    assert!(
        mj_net.nodes.iter().any(|n| n.is_median_vector),
        "Expected inferred ancestral median vector AAAAAA connecting T1, T2, T3"
    );
}

#[test]
fn bayesian_skyline_and_sliding_window_detect_expansion_and_selection() {
    let presets = builtin_alignment_presets();

    // Preset 2: E. coli star-phylogeny expansion -> check Bayesian Skyline & Rogers-Harpending fit
    let (_, star_fasta) = presets[2];
    let star_seqs = parse_sequence_alignment(star_fasta).expect("parse star alignment");
    let star_sum =
        analyze_alignment("star", &star_seqs, GapHandling::CompleteDeletion).expect("star sum");
    let sky = run_skyline_and_expansion(
        &star_seqs,
        &star_sum.mismatch_histogram,
        3,
        1500,
        1.0e-7,
        1.0,
    )
    .expect("run Bayesian Skyline");
    assert_eq!(sky.grid.len(), 56);
    for pt in &sky.grid {
        assert!(pt.ne_lower_95 <= pt.ne_median + 1e-9);
        assert!(pt.ne_median <= pt.ne_upper_95 + 1e-9);
    }
    assert!(sky.expansion.tau > 0.0);
    assert!(sky.to_csv().lines().count() >= 56);

    // Preset 3: 180-bp Chromosomal Selective Sweep (bp 1-60) & Balancing F_ST Peak (bp 121-180)
    let (_, chrom_fasta) = presets[3];
    let chrom_seqs = parse_sequence_alignment(chrom_fasta).expect("parse 180bp locus");
    let scan = run_sliding_window_scan(&chrom_seqs, 50, 20, GapHandling::CompleteDeletion)
        .expect("sliding window scan");
    assert!(scan.windows.len() >= 6);
    let first_w = &scan.windows[0];
    let last_w = scan.windows.last().expect("last window");
    assert!(
        first_w.tajimas_d < -1.2,
        "Expected negative Tajima's D in sweep window, got {}",
        first_w.tajimas_d
    );
    assert!(
        last_w.tajimas_d > 1.2 && last_w.mean_hudson_fst > 0.65,
        "Expected positive Tajima's D and elevated F_ST in balancing window, got D={} F_ST={}",
        last_w.tajimas_d,
        last_w.mean_hudson_fst
    );
}

#[test]
fn workset_envelope_and_sqlite_roundtrip() {
    let mut ws = WorksetState::default();
    ws.workset_title = "v1.1 Skyline Checkpoint".to_string();
    ws.active_module = ModuleId::E3BayesianSkyline;
    ws.empirical.load_preset(3);

    let json = ws.to_envelope_json().expect("serialize workset");
    let restored = WorksetState::from_envelope_json(&json).expect("deserialize workset");
    assert_eq!(restored.active_module, ModuleId::E3BayesianSkyline);
    assert!(restored.empirical.skyline_report.is_some());
    assert!(restored.empirical.sliding_report.is_some());

    let tmp_db = std::env::temp_dir().join("popbio_v1_1_slots.sqlite3");
    let _ = std::fs::remove_file(&tmp_db);
    let db = WorksetDb::open_path(&tmp_db).expect("open sqlite");
    db.save_slot("Slot_Skyline", &mut ws).expect("save slot");
    let loaded = db.load_slot("Slot_Skyline").expect("load slot");
    assert_eq!(loaded.workset_title, "v1.1 Skyline Checkpoint");
    let _ = std::fs::remove_file(&tmp_db);
}

#[test]
fn wright_3d_landscape_and_leslie_lefkovitch_eigen_elasticity() {
    let mut wl = WrightLandscapeSim::default();
    let w_local = wl.fitness_at(0.24, 0.26);
    let w_global = wl.fitness_at(0.76, 0.74);
    assert!(w_global > w_local, "Global Peak 2 must exceed Local Peak 1");
    for _ in 0..40 {
        wl.step();
    }
    assert_eq!(wl.demes.len(), wl.num_demes);

    let mut turtle = LeslieLefkovitchSim::default();
    // Crouse et al. (1987) Caretta caretta baseline lambda_1 ~ 0.945 (< 1.0 declining without TEDs)
    assert!(
        (turtle.eigen.dominant_lambda - 0.945).abs() < 0.015,
        "Expected loggerhead lambda_1 ~ 0.945, got {}",
        turtle.eigen.dominant_lambda
    );
    let e_sum: f64 = turtle.eigen.stage_elasticity_totals.iter().sum();
    assert!((e_sum - 1.0).abs() < 1e-6);
    // Large juvenile / subadult elasticity exceeds egg/hatchling elasticity
    assert!(turtle.eigen.stage_elasticity_totals[1] + turtle.eigen.stage_elasticity_totals[2] > turtle.eigen.stage_elasticity_totals[0]);

    // Switching to strict 4-age Songbird Leslie matrix converges realized_lambda -> dominant_lambda
    turtle.load_preset(MatrixPreset::SongbirdLeslie4);
    for _ in 0..35 {
        turtle.step();
    }
    assert!((turtle.realized_lambda - turtle.eigen.dominant_lambda).abs() < 0.01);
}

#[test]
fn hare_lynx_default_exhibits_coupled_limit_cycles_with_lagged_peaks() {
    let mut sim = HareLynxSim::default();
    // Pre-warmed 500 steps (25.0 years) immediately on default()
    assert!(sim.history.len() >= 500);
    assert!((sim.t - 25.0).abs() < 1e-6);

    // Advance additional 100 steps (5.0 years) for a 30-year total horizon
    for _ in 0..100 {
        sim.step();
    }
    assert!((sim.t - 30.0).abs() < 1e-6);

    let mut hare_peaks: Vec<(f64, f64)> = Vec::new();
    let mut lynx_peaks: Vec<(f64, f64)> = Vec::new();
    for w in sim.history.windows(3) {
        let (h0, h1, h2) = (w[0].values[0], w[1].values[0], w[2].values[0]);
        if h1 > h0 && h1 >= h2 && h1 > 70.0 {
            hare_peaks.push((w[1].t, h1));
        }
        let (l0, l1, l2) = (w[0].values[1], w[1].values[1], w[2].values[1]);
        if l1 > l0 && l1 >= l2 && l1 > 12.0 {
            lynx_peaks.push((w[1].t, l1));
        }
    }

    assert!(
        hare_peaks.len() >= 2,
        "Expected >= 2 distinct Snowshoe Hare peaks over 30 years, got {:?}",
        hare_peaks
    );
    assert!(
        lynx_peaks.len() >= 2,
        "Expected >= 2 distinct Canada Lynx peaks over 30 years, got {:?}",
        lynx_peaks
    );

    // Every Hare peak (prior to the end boundary) must be followed by a quarter-phase-lagged Lynx peak
    let mut paired_lags = 0usize;
    for &(ht, _) in &hare_peaks {
        if let Some(&(lt, _)) = lynx_peaks.iter().find(|&&(lt, _)| lt > ht) {
            let lag = lt - ht;
            assert!(
                lag > 0.8 && lag < 5.0,
                "Expected quarter-phase Lynx lag after Hare peak at t={ht:.2}, got lag={lag:.2} (Lynx peak t={lt:.2})"
            );
            paired_lags += 1;
        }
    }
    assert!(
        paired_lags >= 2,
        "Expected >= 2 Hare->Lynx phase-lagged cycle pairs, got {paired_lags}"
    );

    // Verify scaled plot helper returns Hares, 4x-scaled Lynx, and kill rate
    let plot_series = sim.time_series_for_plot();
    assert_eq!(plot_series[0].1.len(), sim.history.len());
    assert!(plot_series[1].0.contains("×4"));
    let lynx_plot_max = plot_series[1]
        .1
        .iter()
        .map(|p| p[1])
        .fold(0.0_f64, f64::max);
    assert!(
        lynx_plot_max > 60.0,
        "Expected 4x-scaled Lynx peaks to span alongside Hares (>60), got {lynx_plot_max:.2}"
    );
}

#[test]
fn custom_matrix_and_skyline_log10_roundtrip_in_workset_envelope() {
    let mut ws = WorksetState::default();
    ws.t6_matrix
        .load_from_csv("StageA,StageB,StageC\n0.0,1.2,3.5\n0.4,0.1,0.0\n0.0,0.6,0.75\n")
        .expect("load custom CSV");
    ws.empirical.skyline_log10_y = true;

    let json = ws.to_envelope_json().expect("serialize custom workset");
    let restored = WorksetState::from_envelope_json(&json).expect("deserialize custom workset");
    assert_eq!(restored.t6_matrix.preset, MatrixPreset::Custom);
    assert_eq!(restored.t6_matrix.stage_names, vec!["StageA", "StageB", "StageC"]);
    assert!(restored.empirical.skyline_log10_y);

    // Ensure legacy payloads lacking `skyline_log10_y` still deserialize cleanly via #[serde(default)]
    let raw_payload = serde_json::to_string(&ws).expect("serialize payload");
    let stripped_payload = raw_payload.replace(",\"skyline_log10_y\":true", "");
    assert!(!stripped_payload.contains("skyline_log10_y"));
    let legacy_restored: WorksetState =
        serde_json::from_str(&stripped_payload).expect("deserialize legacy payload");
    assert!(!legacy_restored.empirical.skyline_log10_y);
}

#[test]
fn headless_f1_escape_priority_and_space_e2_haplotype_playback() {
    let ctx = egui::Context::default();
    let mut ws = WorksetState::default();
    ws.active_module = ModuleId::E2HaplotypeNetwork;

    let mut app = PopBioApp {
        screen: AppScreen::ModuleWorkspace,
        active_tab: WorkspaceTab::VisualizationAndDistribution,
        workset: ws,
        playing: false,
        steps_per_frame: 1,
        filter_query: String::new(),
        new_slot_name: "Lab_Checkpoint_1".to_string(),
        saved_slots: Vec::new(),
        db: None,
        db_error: None,
        dragged_hap_node: None,
        show_help_modal: false,
        show_about_modal: false,
        status_message: "Ready".to_string(),
    };

    let run_frame = |ctx: &egui::Context, app: &mut PopBioApp, events: Vec<egui::Event>| {
        let raw = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1400.0, 880.0),
            )),
            events,
            ..Default::default()
        };
        let _ = ctx.run(raw, |c| app.ui_step(c));
    };

    let key = |k: egui::Key| egui::Event::Key {
        key: k,
        physical_key: None,
        pressed: true,
        repeat: false,
        modifiers: egui::Modifiers::NONE,
    };

    // Warm up frame
    run_frame(&ctx, &mut app, vec![]);
    assert!(!app.show_help_modal);
    assert!(!app.workset.empirical.hap_sim_playing);

    // F1 toggles Help modal open and closed
    run_frame(&ctx, &mut app, vec![key(egui::Key::F1)]);
    assert!(app.show_help_modal);
    run_frame(&ctx, &mut app, vec![key(egui::Key::F1)]);
    assert!(!app.show_help_modal);

    // Space toggles E2 Haplotype Network playback without flipping dynamical playing flag
    run_frame(&ctx, &mut app, vec![key(egui::Key::Space)]);
    assert!(app.workset.empirical.hap_sim_playing);
    assert!(!app.playing);
    run_frame(&ctx, &mut app, vec![key(egui::Key::Space)]);
    assert!(!app.workset.empirical.hap_sim_playing);

    // Re-enable playback and open both modals to verify layered Escape priority
    app.workset.empirical.hap_sim_playing = true;
    app.show_about_modal = true;
    app.show_help_modal = true;

    run_frame(&ctx, &mut app, vec![key(egui::Key::Escape)]);
    assert!(!app.show_about_modal);
    assert!(app.show_help_modal);
    assert_eq!(app.screen, AppScreen::ModuleWorkspace);
    assert!(app.workset.empirical.hap_sim_playing);

    run_frame(&ctx, &mut app, vec![key(egui::Key::Escape)]);
    assert!(!app.show_help_modal);
    assert_eq!(app.screen, AppScreen::ModuleWorkspace);
    assert!(app.workset.empirical.hap_sim_playing);

    // Final Escape returns to FrontMenu and stops playback
    run_frame(&ctx, &mut app, vec![key(egui::Key::Escape)]);
    assert_eq!(app.screen, AppScreen::FrontMenu);
    assert!(!app.playing);
    assert!(!app.workset.empirical.hap_sim_playing);
}
