use serde::{Deserialize, Serialize};

use crate::model::{GrowthPhase, PhaseTransitionEvent, TimeRecord};
use crate::rng::SimRng;

const MAX_HISTORY: usize = 2400;

fn push_record(
    history: &mut Vec<TimeRecord>,
    transitions: &mut Vec<PhaseTransitionEvent>,
    rec: TimeRecord,
) {
    let changed = history
        .last()
        .map(|prev| prev.phase != rec.phase)
        .unwrap_or(true);
    if changed {
        transitions.push(PhaseTransitionEvent {
            t: rec.t,
            phase: rec.phase,
        });
    }
    if history.len() >= MAX_HISTORY {
        history.drain(0..200);
    }
    history.push(rec);
}

// ============================================================================
// T1: Levins (1969) Metapopulation & Habitat Destruction Network (64 Patches)
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PatchNode {
    pub x: f64,
    pub y: f64,
    pub destroyed: bool,
    pub occupied: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LevinsMetapopSim {
    pub t: f64,
    pub dt: f64,
    pub det_occupancy_p: f64,
    pub colonization_c: f64,
    pub extinction_e: f64,
    pub destruction_h: f64,
    pub patches: Vec<PatchNode>,
    pub current_phase: GrowthPhase,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for LevinsMetapopSim {
    fn default() -> Self {
        let mut rng = SimRng::new(1969);
        let mut patches = Vec::with_capacity(64);
        for i in 0..64 {
            let r = i / 8;
            let c = i % 8;
            patches.push(PatchNode {
                x: ((c as f64 + 0.5) / 8.0 + rng.range_f64(-0.03, 0.03)).clamp(0.05, 0.95),
                y: ((r as f64 + 0.5) / 8.0 + rng.range_f64(-0.03, 0.03)).clamp(0.05, 0.95),
                destroyed: i < 8, // initial h = 0.125
                occupied: i >= 8 && i % 3 != 0,
            });
        }
        let mut s = Self {
            t: 0.0,
            dt: 0.1,
            det_occupancy_p: 0.55,
            colonization_c: 0.62,
            extinction_e: 0.20,
            destruction_h: 0.125,
            patches,
            current_phase: GrowthPhase::DeceleratingTransition,
            rng,
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.apply_destruction_fraction();
        s.record_state(0.0);
        s
    }
}

impl LevinsMetapopSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Stochastic Patch Occupancy (%)",
            "Deterministic Levins p(t) (%)",
            "Analytical Equilibrium p* (%)",
            "Habitat Destroyed h (%)",
            "",
        ]
    }

    pub fn critical_destruction_threshold(&self) -> f64 {
        (1.0 - self.extinction_e / self.colonization_c.max(0.01)).clamp(0.0, 1.0)
    }

    pub fn analytical_equilibrium(&self) -> f64 {
        ((1.0 - self.destruction_h) - (self.extinction_e / self.colonization_c.max(0.01))).max(0.0)
    }

    pub fn apply_destruction_fraction(&mut self) {
        let n_destroy = ((self.destruction_h * self.patches.len() as f64).round() as usize)
            .min(self.patches.len());
        for (i, patch) in self.patches.iter_mut().enumerate() {
            patch.destroyed = i < n_destroy;
            if patch.destroyed {
                patch.occupied = false;
            }
        }
    }

    fn stochastic_occupancy(&self) -> f64 {
        let occ = self
            .patches
            .iter()
            .filter(|p| !p.destroyed && p.occupied)
            .count();
        occ as f64 / (self.patches.len().max(1) as f64)
    }

    fn record_state(&mut self, dp: f64) {
        let stoch_p = self.stochastic_occupancy();
        let p_star = self.analytical_equilibrium();
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.t,
                values: [
                    stoch_p * 100.0,
                    self.det_occupancy_p * 100.0,
                    p_star * 100.0,
                    self.destruction_h * 100.0,
                    0.0,
                ],
                per_capita_growth: dp,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        self.apply_destruction_fraction();
        let dt = self.dt;
        let dp = self.colonization_c
            * self.det_occupancy_p
            * (1.0 - self.destruction_h - self.det_occupancy_p)
            - self.extinction_e * self.det_occupancy_p;
        self.det_occupancy_p =
            (self.det_occupancy_p + dt * dp).clamp(0.0, (1.0 - self.destruction_h).max(0.0));

        let stoch_p = self.stochastic_occupancy();
        let p_col = (1.0 - (-self.colonization_c * stoch_p * dt).exp()).clamp(0.0, 1.0);
        let p_ext = (1.0 - (-self.extinction_e * dt).exp()).clamp(0.0, 1.0);

        for patch in &mut self.patches {
            if patch.destroyed {
                patch.occupied = false;
            } else if patch.occupied {
                if self.rng.next_f64() < p_ext {
                    patch.occupied = false;
                }
            } else if self.rng.next_f64() < p_col {
                patch.occupied = true;
            }
        }
        self.t += dt;

        self.current_phase = if self.destruction_h >= self.critical_destruction_threshold() {
            GrowthPhase::DeathSenescence
        } else if self.stochastic_occupancy() < 0.05 {
            GrowthPhase::PopulationCrashBottleneck
        } else if (self.det_occupancy_p - self.analytical_equilibrium()).abs() < 0.02 {
            GrowthPhase::StationaryCarryingCapacity
        } else if dp > 0.02 {
            GrowthPhase::ExponentialLog
        } else {
            GrowthPhase::DeceleratingTransition
        };

        self.record_state(dp);
    }

    pub fn inject_perturbation(&mut self) {
        // Restore 10% destroyed habitat and recolonize 8 suitable patches
        self.destruction_h = (self.destruction_h - 0.10).max(0.0);
        self.apply_destruction_fraction();
        let mut seeded = 0usize;
        for p in &mut self.patches {
            if !p.destroyed && !p.occupied && seeded < 8 {
                p.occupied = true;
                seeded += 1;
            }
        }
        self.det_occupancy_p = self.stochastic_occupancy().max(0.08);
    }
}

// ============================================================================
// T2: May (1976) Discrete Logistic & Ricker Bifurcation / Chaos Lab
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MayBifurcationSim {
    pub generation: usize,
    pub ricker_map: bool,
    pub r_param: f64,
    pub carrying_k: f64,
    pub n_current: f64,
    pub lyapunov_sum: f64,
    pub lyapunov_count: usize,
    pub cobweb_points: Vec<(f64, f64)>,
    pub current_phase: GrowthPhase,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for MayBifurcationSim {
    fn default() -> Self {
        let mut s = Self {
            generation: 0,
            ricker_map: false,
            r_param: 2.62,
            carrying_k: 100.0,
            n_current: 18.0,
            lyapunov_sum: 0.0,
            lyapunov_count: 0,
            cobweb_points: vec![(18.0, 0.0)],
            current_phase: GrowthPhase::LimitCycleOscillation,
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.record_state(0.0);
        s
    }
}

impl MayBifurcationSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Population N_t",
            "Carrying Capacity K",
            "Lyapunov Exponent lambda (x 50)",
            "",
            "",
        ]
    }

    pub fn map_eval(&self, n: f64, r: f64) -> f64 {
        let k = self.carrying_k.max(1.0);
        if self.ricker_map {
            n * (r * (1.0 - n / k)).exp()
        } else {
            (n * (1.0 + r * (1.0 - n / k))).max(0.0)
        }
    }

    pub fn map_derivative_abs(&self, n: f64, r: f64) -> f64 {
        let k = self.carrying_k.max(1.0);
        if self.ricker_map {
            ((1.0 - r * n / k) * (r * (1.0 - n / k)).exp()).abs()
        } else {
            (1.0 + r - (2.0 * r * n) / k).abs()
        }
    }

    pub fn lyapunov_exponent(&self) -> f64 {
        if self.lyapunov_count == 0 {
            0.0
        } else {
            self.lyapunov_sum / (self.lyapunov_count as f64)
        }
    }

    fn record_state(&mut self, pcg: f64) {
        let lyap = self.lyapunov_exponent();
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.generation as f64,
                values: [self.n_current, self.carrying_k, lyap * 50.0, 0.0, 0.0],
                per_capita_growth: pcg,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        let prev = self.n_current;
        let next = self.map_eval(prev, self.r_param);
        let deriv = self.map_derivative_abs(prev, self.r_param).max(1e-9);
        self.lyapunov_sum += deriv.ln();
        self.lyapunov_count += 1;

        self.cobweb_points.push((prev, next));
        self.cobweb_points.push((next, next));
        if self.cobweb_points.len() > 160 {
            self.cobweb_points.drain(0..20);
        }

        self.n_current = next.clamp(0.0, self.carrying_k * 4.0);
        self.generation += 1;

        let lyap = self.lyapunov_exponent();
        let pcg = if prev > 1e-6 && self.n_current > 1e-6 {
            (self.n_current / prev).ln()
        } else {
            0.0
        };

        self.current_phase = if lyap > 0.01 || self.r_param > 2.57 {
            GrowthPhase::ChaoticRegime
        } else if self.r_param > 2.0 {
            GrowthPhase::LimitCycleOscillation
        } else if (self.n_current - self.carrying_k).abs() < 0.8 {
            GrowthPhase::StationaryCarryingCapacity
        } else if self.n_current < 0.5 * self.carrying_k {
            GrowthPhase::ExponentialLog
        } else {
            GrowthPhase::DeceleratingTransition
        };

        self.record_state(pcg);
    }

    pub fn inject_perturbation(&mut self) {
        self.n_current = 0.12 * self.carrying_k;
        self.lyapunov_sum = 0.0;
        self.lyapunov_count = 0;
        self.cobweb_points.clear();
        self.cobweb_points.push((self.n_current, 0.0));
    }
}

// ============================================================================
// T3: Wright-Fisher Genetic Drift, Selection, Mutation & Island Gene Flow
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WrightFisherSim {
    pub generation: usize,
    pub n_diploid: usize,
    pub num_demes: usize,
    pub initial_p0: f64,
    pub w_aa_dom: f64,
    pub w_aa_het: f64,
    pub w_aa_rec: f64,
    pub mu_forward: f64,
    pub nu_backward: f64,
    pub migration_m: f64,
    pub deme_freqs: Vec<f64>,
    pub deme_histories: Vec<Vec<f64>>,
    pub current_phase: GrowthPhase,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for WrightFisherSim {
    fn default() -> Self {
        let num_demes = 6;
        let p0 = 0.50;
        let mut s = Self {
            generation: 0,
            n_diploid: 60,
            num_demes,
            initial_p0: p0,
            w_aa_dom: 1.00,
            w_aa_het: 1.00,
            w_aa_rec: 0.94,
            mu_forward: 0.0005,
            nu_backward: 0.0005,
            migration_m: 0.01,
            deme_freqs: vec![p0; num_demes],
            deme_histories: vec![vec![p0]; num_demes],
            current_phase: GrowthPhase::DriftSelectionTrajectory,
            rng: SimRng::new(1931),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.record_state(0.0);
        s
    }
}

impl WrightFisherSim {
    pub fn reset(&mut self) {
        let p0 = self.initial_p0.clamp(0.01, 0.99);
        self.generation = 0;
        self.deme_freqs = vec![p0; self.num_demes];
        self.deme_histories = vec![vec![p0]; self.num_demes];
        self.history.clear();
        self.transitions.clear();
        self.current_phase = GrowthPhase::DriftSelectionTrajectory;
        self.record_state(0.0);
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Mean Allele Frequency p_bar",
            "Within-Deme Heterozygosity H_S",
            "Total Heterozygosity H_T",
            "Wright's Fixation Index F_ST",
            "Deme 1 Allele Frequency p_1",
        ]
    }

    pub fn summary_popgen(&self) -> (f64, f64, f64, f64) {
        let k = self.deme_freqs.len().max(1) as f64;
        let p_bar = self.deme_freqs.iter().sum::<f64>() / k;
        let h_s = self
            .deme_freqs
            .iter()
            .map(|&p| 2.0 * p * (1.0 - p))
            .sum::<f64>()
            / k;
        let h_t = 2.0 * p_bar * (1.0 - p_bar);
        let fst = if h_t > 1e-6 {
            ((h_t - h_s) / h_t).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (p_bar, h_s, h_t, fst)
    }

    fn record_state(&mut self, dp: f64) {
        let (p_bar, h_s, h_t, fst) = self.summary_popgen();
        let p1 = self.deme_freqs.first().copied().unwrap_or(p_bar);
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.generation as f64,
                values: [p_bar, h_s, h_t, fst, p1],
                per_capita_growth: dp,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        if self.deme_freqs.len() != self.num_demes {
            self.deme_freqs.resize(self.num_demes, self.initial_p0);
            self.deme_histories
                .resize(self.num_demes, vec![self.initial_p0]);
        }
        let (prev_p_bar, _, _, _) = self.summary_popgen();
        let migrant_pool = prev_p_bar;
        let two_n = (2 * self.n_diploid).max(4);

        for i in 0..self.num_demes {
            let p = self.deme_freqs[i];
            let q = 1.0 - p;
            // Viability selection
            let w_bar = (p * p * self.w_aa_dom
                + 2.0 * p * q * self.w_aa_het
                + q * q * self.w_aa_rec)
                .max(1e-9);
            let p_sel = (p * p * self.w_aa_dom + p * q * self.w_aa_het) / w_bar;
            // Bidirectional mutation
            let p_mut = p_sel * (1.0 - self.mu_forward) + (1.0 - p_sel) * self.nu_backward;
            // Island migration
            let p_mig = (1.0 - self.migration_m) * p_mut + self.migration_m * migrant_pool;
            // Wright-Fisher binomial gamete sampling
            let copies = self.rng.binomial(two_n, p_mig);
            let p_next = (copies as f64) / (two_n as f64);
            self.deme_freqs[i] = p_next;
            if self.deme_histories[i].len() >= 1200 {
                self.deme_histories[i].drain(0..100);
            }
            self.deme_histories[i].push(p_next);
        }
        self.generation += 1;

        let (p_bar, h_s, _, _) = self.summary_popgen();
        let all_fixed_or_lost = self
            .deme_freqs
            .iter()
            .all(|&p| p <= 0.001 || p >= 0.999);
        self.current_phase = if all_fixed_or_lost {
            GrowthPhase::FixationOrLoss
        } else if h_s < 0.05 {
            GrowthPhase::DeceleratingTransition
        } else {
            GrowthPhase::DriftSelectionTrajectory
        };

        self.record_state(p_bar - prev_p_bar);
    }

    pub fn inject_perturbation(&mut self) {
        // Severe founder bottleneck: sample 8 alleles (2N = 8) once across all demes
        for i in 0..self.deme_freqs.len() {
            let k = self.rng.binomial(8, self.deme_freqs[i]);
            self.deme_freqs[i] = k as f64 / 8.0;
        }
    }
}

// ============================================================================
// T4: Kingman (1982) Coalescent Genealogy & Neutral Site Frequency Spectrum
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoalescentNode {
    pub id: usize,
    pub time_gen: f64,
    pub x_pos: f64,
    pub left_child: Option<usize>,
    pub right_child: Option<usize>,
    pub descendant_tips: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoalescentBranch {
    pub child_id: usize,
    pub parent_id: usize,
    pub length_gen: f64,
    pub mutations: usize,
    pub descendant_tips: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KingmanCoalescentSim {
    pub replicate_index: usize,
    pub sample_size_n: usize,
    pub effective_ne: f64,
    pub theta_locus: f64,
    pub nodes: Vec<CoalescentNode>,
    pub branches: Vec<CoalescentBranch>,
    pub t_mrca: f64,
    pub total_tree_length: f64,
    pub observed_sfs: Vec<usize>,
    pub cumulative_sfs: Vec<f64>,
    pub current_phase: GrowthPhase,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for KingmanCoalescentSim {
    fn default() -> Self {
        let mut s = Self {
            replicate_index: 0,
            sample_size_n: 14,
            effective_ne: 1000.0,
            theta_locus: 12.0,
            nodes: Vec::new(),
            branches: Vec::new(),
            t_mrca: 0.0,
            total_tree_length: 0.0,
            observed_sfs: Vec::new(),
            cumulative_sfs: Vec::new(),
            current_phase: GrowthPhase::DriftSelectionTrajectory,
            rng: SimRng::new(1982),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.step();
        s
    }
}

impl KingmanCoalescentSim {
    pub fn reset(&mut self) {
        self.replicate_index = 0;
        self.cumulative_sfs.clear();
        self.history.clear();
        self.transitions.clear();
        self.step();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Simulated T_MRCA (generations)",
            "Expected E[T_MRCA] = 4*Ne*(1 - 1/n)",
            "Segregating Mutations S on Tree",
            "",
            "",
        ]
    }

    pub fn expected_tmrca(&self) -> f64 {
        let n = self.sample_size_n.max(2) as f64;
        4.0 * self.effective_ne * (1.0 - 1.0 / n)
    }

    /// Generates one complete Kingman coalescent genealogy with Poisson infinite-sites mutations.
    pub fn step(&mut self) {
        let n = self.sample_size_n.clamp(3, 40);
        self.sample_size_n = n;
        if self.cumulative_sfs.len() != n - 1 {
            self.cumulative_sfs = vec![0.0; n - 1];
        }

        let mut nodes = Vec::with_capacity(2 * n - 1);
        let mut active: Vec<usize> = Vec::with_capacity(n);
        for i in 0..n {
            nodes.push(CoalescentNode {
                id: i,
                time_gen: 0.0,
                x_pos: (i as f64 + 0.5) / (n as f64),
                left_child: None,
                right_child: None,
                descendant_tips: 1,
            });
            active.push(i);
        }

        let mut branches = Vec::with_capacity(2 * n - 2);
        let mut t_curr = 0.0;

        while active.len() > 1 {
            let k = active.len() as f64;
            // Coalescence rate lambda_k = k*(k - 1) / (4 * Ne)
            let rate = (k * (k - 1.0)) / (4.0 * self.effective_ne.max(10.0));
            let wait = self.rng.exponential(rate);
            t_curr += wait;

            let idx_a = self.rng.gen_index(active.len());
            let child_a = active.remove(idx_a);
            let idx_b = self.rng.gen_index(active.len());
            let child_b = active.remove(idx_b);

            let parent_id = nodes.len();
            let x_mid = 0.5 * (nodes[child_a].x_pos + nodes[child_b].x_pos);
            let tips_sum = nodes[child_a].descendant_tips + nodes[child_b].descendant_tips;

            for &cid in &[child_a, child_b] {
                let len = (t_curr - nodes[cid].time_gen).max(0.0);
                // Expected mutations along branch of length `len` is `(theta / 2) * (len / (2 * Ne))` = `theta * len / (4 * Ne)`
                let mean_mut = self.theta_locus * len / (4.0 * self.effective_ne.max(10.0));
                // Poisson sample via Knuth product for modest `mean_mut`
                let mut mut_count = 0usize;
                let l_lim = (-mean_mut.min(30.0)).exp();
                let mut prod = self.rng.next_f64();
                while prod > l_lim && mut_count < 60 {
                    mut_count += 1;
                    prod *= self.rng.next_f64();
                }

                branches.push(CoalescentBranch {
                    child_id: cid,
                    parent_id,
                    length_gen: len,
                    mutations: mut_count,
                    descendant_tips: nodes[cid].descendant_tips,
                });
            }

            nodes.push(CoalescentNode {
                id: parent_id,
                time_gen: t_curr,
                x_pos: x_mid,
                left_child: Some(child_a),
                right_child: Some(child_b),
                descendant_tips: tips_sum,
            });
            active.push(parent_id);
        }

        let mut sfs = vec![0usize; n - 1];
        let mut total_muts = 0usize;
        let mut total_len = 0.0;
        for br in &branches {
            total_len += br.length_gen;
            total_muts += br.mutations;
            if br.descendant_tips >= 1 && br.descendant_tips < n {
                sfs[br.descendant_tips - 1] += br.mutations;
            }
        }
        for i in 0..(n - 1) {
            self.cumulative_sfs[i] += sfs[i] as f64;
        }

        self.nodes = nodes;
        self.branches = branches;
        self.t_mrca = t_curr;
        self.total_tree_length = total_len;
        self.observed_sfs = sfs;
        self.replicate_index += 1;

        let exp_tmrca = self.expected_tmrca();
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.replicate_index as f64,
                values: [
                    self.t_mrca,
                    exp_tmrca,
                    total_muts as f64,
                    0.0,
                    0.0,
                ],
                per_capita_growth: 0.0,
                phase: self.current_phase,
            },
        );
    }

    pub fn inject_perturbation(&mut self) {
        // Simulate 25 rapid coalescent replicates to converge empirical SFS toward Watterson's theta / i
        for _ in 0..25 {
            self.step();
        }
    }
}
