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

use std::fs;
use std::path::Path;

use chrono::Utc;
use rusqlite::{params, Connection};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::haplotype_network::{build_haplotype_network, HaplotypeNetwork, NetworkAlgorithm};
use crate::model::{GrowthPhase, ModuleId, PhaseTransitionEvent, TimeRecord};
use crate::popgen_alignment::{
    analyze_alignment, builtin_alignment_presets, parse_sequence_alignment, AlignedSequence,
    GapHandling, PopGenSummary,
};
use crate::popgen_skyline::{
    run_skyline_and_expansion, run_sliding_window_scan, SkylineAndExpansionReport,
    SlidingWindowScanReport,
};
use crate::sim_invertebrates::{AphidParasitoidSim, HuffakerMiteSim, TriboliumSim};
use crate::sim_microbiology::{EcoliGrowthSim, GauseProtistSim, PhagePlaqueSim};
use crate::sim_landscapes_matrices::{LeslieLefkovitchSim, WrightLandscapeSim};
use crate::sim_theoretical::{
    KingmanCoalescentSim, LevinsMetapopSim, MayBifurcationSim, WrightFisherSim,
};
use crate::sim_vertebrates::{HareLynxSim, IsleRoyaleSim, SoaySheepSim};

fn default_hap_mu() -> f64 {
    0.004
}

fn default_hap_seed() -> u64 {
    0x5EED_CAFE_2026
}

/// Shared Empirical FASTA Workspace powering Modules E1 (Global PopGen), E2 (Haplotype Network),
/// E3 (Bayesian Skyline & Rogers-Harpending Expansion), and E4 (Sliding-Window Neutrality Scan).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmpiricalWorkspaceState {
    pub dataset_name: String,
    pub raw_alignment_text: String,
    #[serde(default)]
    pub baseline_alignment_text: String,
    #[serde(default)]
    pub hap_sim_playing: bool,
    #[serde(default)]
    pub hap_sim_generation: usize,
    #[serde(default = "default_hap_mu")]
    pub hap_mutation_rate_per_bp: f64,
    #[serde(default = "default_hap_seed")]
    pub hap_sim_seed: u64,
    pub gap_handling: GapHandling,
    pub network_algorithm: NetworkAlgorithm,
    pub mj_epsilon: usize,
    pub show_mutation_labels: bool,
    pub bsp_groups: usize,
    pub bsp_mcmc_steps: usize,
    pub molecular_clock_rate: f64,
    pub generation_time_years: f64,
    pub sliding_window_bp: usize,
    pub sliding_step_bp: usize,
    pub sequences: Vec<AlignedSequence>,
    pub summary: Option<PopGenSummary>,
    pub network: Option<HaplotypeNetwork>,
    pub skyline_report: Option<SkylineAndExpansionReport>,
    pub sliding_report: Option<SlidingWindowScanReport>,
    pub last_error: Option<String>,
    #[serde(default)]
    pub skyline_log10_y: bool,
    #[serde(skip)]
    pub pending_recompute_at: Option<std::time::Instant>,
}

impl Default for EmpiricalWorkspaceState {
    fn default() -> Self {
        let presets = builtin_alignment_presets();
        let (title, fasta) = presets[0];
        let mut s = Self {
            dataset_name: title.to_string(),
            raw_alignment_text: fasta.to_string(),
            baseline_alignment_text: fasta.to_string(),
            hap_sim_playing: false,
            hap_sim_generation: 0,
            hap_mutation_rate_per_bp: default_hap_mu(),
            hap_sim_seed: default_hap_seed(),
            gap_handling: GapHandling::PairwiseDeletion,
            network_algorithm: NetworkAlgorithm::BandeltMedianJoining,
            mj_epsilon: 0,
            show_mutation_labels: true,
            bsp_groups: 4,
            bsp_mcmc_steps: 2400,
            molecular_clock_rate: 1.0e-7,
            generation_time_years: 1.0,
            sliding_window_bp: 20,
            sliding_step_bp: 5,
            sequences: Vec::new(),
            summary: None,
            network: None,
            skyline_report: None,
            sliding_report: None,
            last_error: None,
            skyline_log10_y: false,
            pending_recompute_at: None,
        };
        s.recompute();
        s
    }
}

impl EmpiricalWorkspaceState {
    pub fn load_preset(&mut self, index: usize) {
        let presets = builtin_alignment_presets();
        let (title, fasta) = presets[index % presets.len()];
        self.dataset_name = title.to_string();
        self.raw_alignment_text = fasta.to_string();
        self.baseline_alignment_text = fasta.to_string();
        self.hap_sim_playing = false;
        self.hap_sim_generation = 0;
        // Auto-scale sliding window nicely for longer 240 bp sweep/balancing locus vs 60 bp presets
        if index == 3 {
            self.sliding_window_bp = 48;
            self.sliding_step_bp = 12;
        } else if self.sliding_window_bp > 30 {
            self.sliding_window_bp = 20;
            self.sliding_step_bp = 5;
        }
        self.pending_recompute_at = None;
        self.recompute();
    }

    pub fn reset_haplotype_simulation(&mut self) {
        if !self.baseline_alignment_text.trim().is_empty() {
            self.raw_alignment_text = self.baseline_alignment_text.clone();
        }
        self.hap_sim_playing = false;
        self.hap_sim_generation = 0;
        self.pending_recompute_at = None;
        self.recompute();
    }

    pub fn schedule_debounced_recompute(&mut self, now: std::time::Instant, delay_ms: u64) {
        self.pending_recompute_at = Some(now + std::time::Duration::from_millis(delay_ms));
    }

    pub fn flush_debounced_recompute(&mut self, now: std::time::Instant) -> bool {
        if let Some(deadline) = self.pending_recompute_at {
            if now >= deadline {
                self.pending_recompute_at = None;
                self.recompute();
                return true;
            }
        }
        false
    }

    pub fn recompute_immediately(&mut self) {
        self.pending_recompute_at = None;
        self.recompute();
    }

    pub fn step_haplotype_evolution(&mut self) {
        if self.sequences.len() < 2 {
            return;
        }
        if self.baseline_alignment_text.trim().is_empty() {
            self.baseline_alignment_text = self.raw_alignment_text.clone();
        }

        let mut prev_pos: std::collections::BTreeMap<Vec<u8>, (f64, f64)> =
            std::collections::BTreeMap::new();
        if let Some(net) = &self.network {
            for node in &net.nodes {
                if !node.is_median_vector {
                    prev_pos.insert(node.sequence.clone(), (node.x, node.y));
                }
            }
        }

        let mut by_pop: std::collections::BTreeMap<String, Vec<&AlignedSequence>> =
            std::collections::BTreeMap::new();
        for s in &self.sequences {
            by_pop.entry(s.population.clone()).or_default().push(s);
        }

        let mut rng = crate::rng::SimRng::new(
            self.hap_sim_seed ^ ((self.hap_sim_generation as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15)),
        );
        let mu = self.hap_mutation_rate_per_bp.clamp(0.0, 0.05);
        let mut next_fasta = String::new();
        let mut seq_idx = 0usize;

        for (pop, pool) in by_pop {
            let k = pool.len();
            for _ in 0..k {
                seq_idx += 1;
                let parent = pool[rng.gen_index(k)];
                let mut mutated = Vec::with_capacity(parent.sequence.len());
                for ch in parent.sequence.chars() {
                    let up = ch.to_ascii_uppercase();
                    if matches!(up, 'A' | 'C' | 'G' | 'T') && rng.next_f64() < mu {
                        let alts = match up {
                            'A' => ['C', 'G', 'T'],
                            'C' => ['A', 'G', 'T'],
                            'G' => ['A', 'C', 'T'],
                            _ => ['A', 'C', 'G'],
                        };
                        mutated.push(alts[rng.gen_index(3)]);
                    } else {
                        mutated.push(up);
                    }
                }
                let seq_str: String = mutated.into_iter().collect();
                next_fasta.push_str(&format!(">Sim{seq_idx}|{pop}\n{seq_str}\n"));
            }
        }

        self.hap_sim_generation += 1;
        self.raw_alignment_text = next_fasta;
        self.recompute();

        if let Some(net) = &mut self.network {
            if !prev_pos.is_empty() {
                for node in &mut net.nodes {
                    if let Some(&(px, py)) = prev_pos.get(&node.sequence) {
                        node.x = px;
                        node.y = py;
                    }
                }
                net.relax_layout_from_current(32);
            }
        }
    }

    pub fn recompute(&mut self) {
        match parse_sequence_alignment(&self.raw_alignment_text) {
            Ok(seqs) => {
                let sum_res = analyze_alignment(&self.dataset_name, &seqs, self.gap_handling);
                let net_res =
                    build_haplotype_network(&seqs, self.network_algorithm, self.mj_epsilon);
                match (sum_res, net_res) {
                    (Ok(summary), Ok(net)) => {
                        let sky_res = run_skyline_and_expansion(
                            &seqs,
                            &summary.mismatch_histogram,
                            self.bsp_groups,
                            self.bsp_mcmc_steps,
                            self.molecular_clock_rate,
                            self.generation_time_years,
                        )
                        .ok();
                        let win_res = run_sliding_window_scan(
                            &seqs,
                            self.sliding_window_bp,
                            self.sliding_step_bp,
                            self.gap_handling,
                        )
                        .ok();
                        self.sequences = seqs;
                        self.summary = Some(summary);
                        self.network = Some(net);
                        self.skyline_report = sky_res;
                        self.sliding_report = win_res;
                        self.last_error = None;
                    }
                    (Err(e), _) | (_, Err(e)) => {
                        self.last_error = Some(e);
                    }
                }
            }
            Err(e) => {
                self.last_error = Some(e);
            }
        }
    }
}

/// Complete multi-module "Save Game" / Workset snapshot preserving all 19 modules,
/// active simulation histories, RNG seeds, custom parameters, study notes, and empirical analyses.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorksetState {
    pub workset_title: String,
    pub study_notes: String,
    pub saved_at_utc: String,
    pub active_module: ModuleId,
    pub v1_hare_lynx: HareLynxSim,
    pub v2_isle_royale: IsleRoyaleSim,
    pub v3_soay_sheep: SoaySheepSim,
    pub i1_tribolium: TriboliumSim,
    pub i2_huffaker: HuffakerMiteSim,
    pub i3_aphid: AphidParasitoidSim,
    pub m1_ecoli: EcoliGrowthSim,
    pub m2_gause: GauseProtistSim,
    pub m3_phage: PhagePlaqueSim,
    pub t1_levins: LevinsMetapopSim,
    pub t2_may: MayBifurcationSim,
    pub t3_wright_fisher: WrightFisherSim,
    pub t4_coalescent: KingmanCoalescentSim,
    pub t5_landscape: WrightLandscapeSim,
    pub t6_matrix: LeslieLefkovitchSim,
    pub empirical: EmpiricalWorkspaceState,
}

impl Default for WorksetState {
    fn default() -> Self {
        Self {
            workset_title: "Default Population Biology & Demographic Skyline Workset".to_string(),
            study_notes: "Compare boreal vertebrate cycles, stage-structured Tribolium cannibalism, E. coli growth phases, Bayesian Skyline expansion, and sliding-window neutrality scans.".to_string(),
            saved_at_utc: Utc::now().to_rfc3339(),
            active_module: ModuleId::V1HareLynx,
            v1_hare_lynx: HareLynxSim::default(),
            v2_isle_royale: IsleRoyaleSim::default(),
            v3_soay_sheep: SoaySheepSim::default(),
            i1_tribolium: TriboliumSim::default(),
            i2_huffaker: HuffakerMiteSim::default(),
            i3_aphid: AphidParasitoidSim::default(),
            m1_ecoli: EcoliGrowthSim::default(),
            m2_gause: GauseProtistSim::default(),
            m3_phage: PhagePlaqueSim::default(),
            t1_levins: LevinsMetapopSim::default(),
            t2_may: MayBifurcationSim::default(),
            t3_wright_fisher: WrightFisherSim::default(),
            t4_coalescent: KingmanCoalescentSim::default(),
            t5_landscape: WrightLandscapeSim::default(),
            t6_matrix: LeslieLefkovitchSim::default(),
            empirical: EmpiricalWorkspaceState::default(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorksetEnvelope {
    pub format: String,
    pub schema_version: u32,
    pub sha256_hex: String,
    pub payload_json: String,
}

impl WorksetState {
    pub fn is_empirical_module(id: ModuleId) -> bool {
        matches!(
            id,
            ModuleId::E1SequencePopGen
                | ModuleId::E2HaplotypeNetwork
                | ModuleId::E3BayesianSkyline
                | ModuleId::E4SlidingWindowScan
        )
    }

    pub fn step_module(&mut self, id: ModuleId) {
        match id {
            ModuleId::V1HareLynx => self.v1_hare_lynx.step(),
            ModuleId::V2IsleRoyale => self.v2_isle_royale.step(),
            ModuleId::V3SoaySheep => self.v3_soay_sheep.step(),
            ModuleId::I1TriboliumLpa => self.i1_tribolium.step(),
            ModuleId::I2HuffakerMites => self.i2_huffaker.step(),
            ModuleId::I3AphidParasitoid => self.i3_aphid.step(),
            ModuleId::M1EcoliBatchChemostat => self.m1_ecoli.step(),
            ModuleId::M2GauseProtists => self.m2_gause.step(),
            ModuleId::M3PhagePlaque => self.m3_phage.step(),
            ModuleId::T1LevinsMetapop => self.t1_levins.step(),
            ModuleId::T2MayBifurcation => self.t2_may.step(),
            ModuleId::T3WrightFisher => self.t3_wright_fisher.step(),
            ModuleId::T4KingmanCoalescent => self.t4_coalescent.step(),
            ModuleId::T5WrightLandscape => self.t5_landscape.step(),
            ModuleId::T6LeslieLefkovitch => self.t6_matrix.step(),
            ModuleId::E1SequencePopGen
            | ModuleId::E2HaplotypeNetwork
            | ModuleId::E3BayesianSkyline
            | ModuleId::E4SlidingWindowScan => {
                self.empirical.recompute();
            }
        }
    }

    pub fn current_phase(&self, id: ModuleId) -> GrowthPhase {
        match id {
            ModuleId::V1HareLynx => self.v1_hare_lynx.current_phase,
            ModuleId::V2IsleRoyale => self.v2_isle_royale.current_phase,
            ModuleId::V3SoaySheep => self.v3_soay_sheep.current_phase,
            ModuleId::I1TriboliumLpa => self.i1_tribolium.current_phase,
            ModuleId::I2HuffakerMites => self.i2_huffaker.current_phase,
            ModuleId::I3AphidParasitoid => self.i3_aphid.current_phase,
            ModuleId::M1EcoliBatchChemostat => self.m1_ecoli.current_phase,
            ModuleId::M2GauseProtists => self.m2_gause.current_phase,
            ModuleId::M3PhagePlaque => self.m3_phage.current_phase,
            ModuleId::T1LevinsMetapop => self.t1_levins.current_phase,
            ModuleId::T2MayBifurcation => self.t2_may.current_phase,
            ModuleId::T3WrightFisher => self.t3_wright_fisher.current_phase,
            ModuleId::T4KingmanCoalescent => self.t4_coalescent.current_phase,
            ModuleId::T5WrightLandscape => self.t5_landscape.current_phase,
            ModuleId::T6LeslieLefkovitch => self.t6_matrix.current_phase,
            ModuleId::E1SequencePopGen
            | ModuleId::E2HaplotypeNetwork
            | ModuleId::E3BayesianSkyline
            | ModuleId::E4SlidingWindowScan => GrowthPhase::EmpiricalAnalysis,
        }
    }

    pub fn advance_to_next_growth_phase(&mut self, id: ModuleId, max_steps: usize) -> usize {
        let initial = self.current_phase(id);
        if initial == GrowthPhase::EmpiricalAnalysis {
            self.empirical.recompute();
            return 0;
        }
        let mut taken = 0usize;
        for _ in 0..max_steps.max(1) {
            self.step_module(id);
            taken += 1;
            if self.current_phase(id) != initial {
                break;
            }
        }
        taken
    }

    pub fn reset_module(&mut self, id: ModuleId) {
        match id {
            ModuleId::V1HareLynx => self.v1_hare_lynx.reset(),
            ModuleId::V2IsleRoyale => self.v2_isle_royale.reset(),
            ModuleId::V3SoaySheep => self.v3_soay_sheep.reset(),
            ModuleId::I1TriboliumLpa => self.i1_tribolium.reset(),
            ModuleId::I2HuffakerMites => self.i2_huffaker.reset(),
            ModuleId::I3AphidParasitoid => self.i3_aphid.reset(),
            ModuleId::M1EcoliBatchChemostat => self.m1_ecoli.reset(),
            ModuleId::M2GauseProtists => self.m2_gause.reset(),
            ModuleId::M3PhagePlaque => self.m3_phage.reset(),
            ModuleId::T1LevinsMetapop => self.t1_levins.reset(),
            ModuleId::T2MayBifurcation => self.t2_may.reset(),
            ModuleId::T3WrightFisher => self.t3_wright_fisher.reset(),
            ModuleId::T4KingmanCoalescent => self.t4_coalescent.reset(),
            ModuleId::T5WrightLandscape => self.t5_landscape.reset(),
            ModuleId::T6LeslieLefkovitch => self.t6_matrix.reset(),
            ModuleId::E1SequencePopGen
            | ModuleId::E2HaplotypeNetwork
            | ModuleId::E3BayesianSkyline
            | ModuleId::E4SlidingWindowScan => {
                self.empirical = EmpiricalWorkspaceState::default();
            }
        }
    }

    pub fn inject_perturbation(&mut self, id: ModuleId) {
        match id {
            ModuleId::V1HareLynx => self.v1_hare_lynx.inject_perturbation(),
            ModuleId::V2IsleRoyale => self.v2_isle_royale.inject_perturbation(),
            ModuleId::V3SoaySheep => self.v3_soay_sheep.inject_perturbation(),
            ModuleId::I1TriboliumLpa => self.i1_tribolium.inject_perturbation(),
            ModuleId::I2HuffakerMites => self.i2_huffaker.inject_perturbation(),
            ModuleId::I3AphidParasitoid => self.i3_aphid.inject_perturbation(),
            ModuleId::M1EcoliBatchChemostat => self.m1_ecoli.inject_perturbation(),
            ModuleId::M2GauseProtists => self.m2_gause.inject_perturbation(),
            ModuleId::M3PhagePlaque => self.m3_phage.inject_perturbation(),
            ModuleId::T1LevinsMetapop => self.t1_levins.inject_perturbation(),
            ModuleId::T2MayBifurcation => self.t2_may.inject_perturbation(),
            ModuleId::T3WrightFisher => self.t3_wright_fisher.inject_perturbation(),
            ModuleId::T4KingmanCoalescent => self.t4_coalescent.inject_perturbation(),
            ModuleId::T5WrightLandscape => self.t5_landscape.inject_perturbation(),
            ModuleId::T6LeslieLefkovitch => self.t6_matrix.inject_perturbation(),
            ModuleId::E1SequencePopGen
            | ModuleId::E2HaplotypeNetwork
            | ModuleId::E3BayesianSkyline
            | ModuleId::E4SlidingWindowScan => {
                self.empirical.recompute();
            }
        }
    }

    pub fn history_and_names(
        &self,
        id: ModuleId,
    ) -> (&[TimeRecord], &[PhaseTransitionEvent], [&'static str; 5]) {
        match id {
            ModuleId::V1HareLynx => (
                &self.v1_hare_lynx.history,
                &self.v1_hare_lynx.transitions,
                self.v1_hare_lynx.series_names(),
            ),
            ModuleId::V2IsleRoyale => (
                &self.v2_isle_royale.history,
                &self.v2_isle_royale.transitions,
                self.v2_isle_royale.series_names(),
            ),
            ModuleId::V3SoaySheep => (
                &self.v3_soay_sheep.history,
                &self.v3_soay_sheep.transitions,
                self.v3_soay_sheep.series_names(),
            ),
            ModuleId::I1TriboliumLpa => (
                &self.i1_tribolium.history,
                &self.i1_tribolium.transitions,
                self.i1_tribolium.series_names(),
            ),
            ModuleId::I2HuffakerMites => (
                &self.i2_huffaker.history,
                &self.i2_huffaker.transitions,
                self.i2_huffaker.series_names(),
            ),
            ModuleId::I3AphidParasitoid => (
                &self.i3_aphid.history,
                &self.i3_aphid.transitions,
                self.i3_aphid.series_names(),
            ),
            ModuleId::M1EcoliBatchChemostat => (
                &self.m1_ecoli.history,
                &self.m1_ecoli.transitions,
                self.m1_ecoli.series_names(),
            ),
            ModuleId::M2GauseProtists => (
                &self.m2_gause.history,
                &self.m2_gause.transitions,
                self.m2_gause.series_names(),
            ),
            ModuleId::M3PhagePlaque => (
                &self.m3_phage.history,
                &self.m3_phage.transitions,
                self.m3_phage.series_names(),
            ),
            ModuleId::T1LevinsMetapop => (
                &self.t1_levins.history,
                &self.t1_levins.transitions,
                self.t1_levins.series_names(),
            ),
            ModuleId::T2MayBifurcation => (
                &self.t2_may.history,
                &self.t2_may.transitions,
                self.t2_may.series_names(),
            ),
            ModuleId::T3WrightFisher => (
                &self.t3_wright_fisher.history,
                &self.t3_wright_fisher.transitions,
                self.t3_wright_fisher.series_names(),
            ),
            ModuleId::T4KingmanCoalescent => (
                &self.t4_coalescent.history,
                &self.t4_coalescent.transitions,
                self.t4_coalescent.series_names(),
            ),
            ModuleId::T5WrightLandscape => (
                &self.t5_landscape.history,
                &self.t5_landscape.transitions,
                self.t5_landscape.series_names(),
            ),
            ModuleId::T6LeslieLefkovitch => (
                &self.t6_matrix.history,
                &self.t6_matrix.transitions,
                self.t6_matrix.series_names(),
            ),
            ModuleId::E1SequencePopGen
            | ModuleId::E2HaplotypeNetwork
            | ModuleId::E3BayesianSkyline
            | ModuleId::E4SlidingWindowScan => (&[], &[], ["", "", "", "", ""]),
        }
    }

    pub fn export_module_csv(&self, id: ModuleId) -> String {
        if id == ModuleId::E3BayesianSkyline {
            return self
                .empirical
                .skyline_report
                .as_ref()
                .map(|r| r.to_csv())
                .unwrap_or_default();
        }
        if id == ModuleId::E4SlidingWindowScan {
            return self
                .empirical
                .sliding_report
                .as_ref()
                .map(|r| r.to_csv())
                .unwrap_or_default();
        }
        let (records, _, names) = self.history_and_names(id);
        let mut csv = format!(
            "time,{},{},{},{},{},per_capita_growth,growth_phase\n",
            names[0], names[1], names[2], names[3], names[4]
        );
        for rec in records {
            csv.push_str(&format!(
                "{:.4},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},\"{}\"\n",
                rec.t,
                rec.values[0],
                rec.values[1],
                rec.values[2],
                rec.values[3],
                rec.values[4],
                rec.per_capita_growth,
                rec.phase.label()
            ));
        }
        csv
    }

    pub fn to_envelope_json(&mut self) -> Result<String, String> {
        self.saved_at_utc = Utc::now().to_rfc3339();
        let payload_json =
            serde_json::to_string_pretty(self).map_err(|e| format!("Serialize error: {e}"))?;
        let mut hasher = Sha256::new();
        hasher.update(payload_json.as_bytes());
        let sha256_hex = format!("{:x}", hasher.finalize());
        let env = WorksetEnvelope {
            format: "popbio-studio-workset".to_string(),
            schema_version: 2,
            sha256_hex,
            payload_json,
        };
        serde_json::to_string_pretty(&env).map_err(|e| format!("Envelope JSON error: {e}"))
    }

    pub fn from_envelope_json(raw: &str) -> Result<Self, String> {
        let env: WorksetEnvelope =
            serde_json::from_str(raw).map_err(|e| format!("Invalid workset JSON: {e}"))?;
        let mut hasher = Sha256::new();
        hasher.update(env.payload_json.as_bytes());
        let computed = format!("{:x}", hasher.finalize());
        if computed != env.sha256_hex {
            return Err(format!(
                "Workset SHA-256 checksum mismatch (expected {}, got {}).",
                env.sha256_hex, computed
            ));
        }
        serde_json::from_str(&env.payload_json).map_err(|e| format!("Payload parse error: {e}"))
    }

    pub fn save_to_file(&mut self, path: &Path) -> Result<(), String> {
        let txt = self.to_envelope_json()?;
        fs::write(path, txt).map_err(|e| format!("Failed writing {}: {e}", path.display()))
    }

    pub fn load_from_file(path: &Path) -> Result<Self, String> {
        let txt =
            fs::read_to_string(path).map_err(|e| format!("Failed reading {}: {e}", path.display()))?;
        Self::from_envelope_json(&txt)
    }
}

#[derive(Debug, Clone)]
pub struct SaveSlotInfo {
    pub slot_name: String,
    pub workset_title: String,
    pub active_module: String,
    pub saved_at_utc: String,
}

pub struct WorksetDb {
    conn: Connection,
}

impl WorksetDb {
    pub fn open_default() -> Result<Self, String> {
        Self::open_path(Path::new("popbio_worksets.sqlite3"))
    }

    pub fn open_path(path: &Path) -> Result<Self, String> {
        let conn =
            Connection::open(path).map_err(|e| format!("SQLite open error ({}): {e}", path.display()))?;
        conn.execute_batch(
            "PRAGMA journal_mode = WAL;\n\
             CREATE TABLE IF NOT EXISTS workset_slots (\n\
                 slot_name TEXT PRIMARY KEY,\n\
                 workset_title TEXT NOT NULL,\n\
                 active_module TEXT NOT NULL,\n\
                 saved_at_utc TEXT NOT NULL,\n\
                 envelope_json TEXT NOT NULL\n\
             );",
        )
        .map_err(|e| format!("SQLite schema init error: {e}"))?;
        Ok(Self { conn })
    }

    pub fn save_slot(&self, slot_name: &str, state: &mut WorksetState) -> Result<(), String> {
        let env_json = state.to_envelope_json()?;
        self.conn
            .execute(
                "INSERT INTO workset_slots (slot_name, workset_title, active_module, saved_at_utc, envelope_json)\n\
                 VALUES (?1, ?2, ?3, ?4, ?5)\n\
                 ON CONFLICT(slot_name) DO UPDATE SET\n\
                     workset_title = excluded.workset_title,\n\
                     active_module = excluded.active_module,\n\
                     saved_at_utc = excluded.saved_at_utc,\n\
                     envelope_json = excluded.envelope_json;",
                params![
                    slot_name,
                    state.workset_title,
                    state.active_module.slug(),
                    state.saved_at_utc,
                    env_json,
                ],
            )
            .map_err(|e| format!("SQLite save_slot error: {e}"))?;
        Ok(())
    }

    pub fn load_slot(&self, slot_name: &str) -> Result<WorksetState, String> {
        let env_json: String = self
            .conn
            .query_row(
                "SELECT envelope_json FROM workset_slots WHERE slot_name = ?1",
                params![slot_name],
                |row| row.get(0),
            )
            .map_err(|e| format!("Slot '{slot_name}' not found: {e}"))?;
        WorksetState::from_envelope_json(&env_json)
    }

    pub fn delete_slot(&self, slot_name: &str) -> Result<(), String> {
        self.conn
            .execute(
                "DELETE FROM workset_slots WHERE slot_name = ?1",
                params![slot_name],
            )
            .map_err(|e| format!("Delete slot error: {e}"))?;
        Ok(())
    }

    pub fn list_slots(&self) -> Vec<SaveSlotInfo> {
        let mut out = Vec::new();
        if let Ok(mut stmt) = self.conn.prepare(
            "SELECT slot_name, workset_title, active_module, saved_at_utc \
             FROM workset_slots ORDER BY saved_at_utc DESC",
        ) {
            if let Ok(rows) = stmt.query_map([], |row| {
                Ok(SaveSlotInfo {
                    slot_name: row.get(0)?,
                    workset_title: row.get(1)?,
                    active_module: row.get(2)?,
                    saved_at_utc: row.get(3)?,
                })
            }) {
                for item in rows.flatten() {
                    out.push(item);
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{Duration, Instant};

    #[test]
    fn empirical_mcmc_debounce_coalesces_until_deadline_or_immediate_flush() {
        let mut emp = EmpiricalWorkspaceState::default();
        let t0 = Instant::now();

        emp.bsp_mcmc_steps = 4000;
        emp.schedule_debounced_recompute(t0, 220);
        assert!(emp.pending_recompute_at.is_some());
        assert!(!emp.flush_debounced_recompute(t0 + Duration::from_millis(100)));
        assert!(emp.pending_recompute_at.is_some());

        // Subsequent drag pushes deadline forward
        let t1 = t0 + Duration::from_millis(120);
        emp.bsp_mcmc_steps = 5200;
        emp.schedule_debounced_recompute(t1, 220);
        assert!(!emp.flush_debounced_recompute(t0 + Duration::from_millis(250)));

        // After full delay elapses, flush executes recompute once
        assert!(emp.flush_debounced_recompute(t1 + Duration::from_millis(230)));
        assert!(emp.pending_recompute_at.is_none());
        assert!(emp.skyline_report.is_some());

        // Immediate flush clears pending deadline
        emp.schedule_debounced_recompute(t1, 500);
        emp.recompute_immediately();
        assert!(emp.pending_recompute_at.is_none());
    }
}
