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

use popbio_studio::haplotype_network::{build_haplotype_network, NetworkAlgorithm};
use popbio_studio::model::{GrowthPhase, ModuleId};
use popbio_studio::popgen_alignment::{
    analyze_alignment, builtin_alignment_presets, parse_sequence_alignment, GapHandling,
};
use popbio_studio::popgen_skyline::{run_skyline_and_expansion, run_sliding_window_scan};
use popbio_studio::sim_landscapes_matrices::{LeslieLefkovitchSim, MatrixPreset, WrightLandscapeSim};
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
