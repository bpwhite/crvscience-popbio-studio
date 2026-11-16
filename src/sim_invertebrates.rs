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

use serde::{Deserialize, Serialize};

use crate::model::{Agent2D, GrowthPhase, PhaseTransitionEvent, TimeRecord};
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
// I1: Flour Beetles (Tribolium castaneum & Tribolium confusum) LPA + Competition
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriboliumSim {
    pub step_index: usize,
    pub larvae: f64,
    pub pupae: f64,
    pub adults: f64,
    pub confusum_adults: f64,
    pub park_competition_mode: bool,
    pub b_recruitment: f64,
    pub c_el: f64,
    pub c_ea: f64,
    pub c_pa: f64,
    pub mu_l: f64,
    pub mu_a: f64,
    pub interspecific_alpha: f64,
    pub current_phase: GrowthPhase,
    pub agents: Vec<Agent2D>,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for TriboliumSim {
    fn default() -> Self {
        // Classic Costantino et al. (1997) Science parameters capable of equilibrium, 2-cycle, or chaos as mu_a / c_pa vary
        let mut s = Self {
            step_index: 0,
            larvae: 120.0,
            pupae: 35.0,
            adults: 65.0,
            confusum_adults: 60.0,
            park_competition_mode: false,
            b_recruitment: 10.45,
            c_el: 0.0173,
            c_ea: 0.0131,
            c_pa: 0.35,
            mu_l: 0.20,
            mu_a: 0.96,
            interspecific_alpha: 1.12,
            current_phase: GrowthPhase::ExponentialLog,
            agents: Vec::new(),
            rng: SimRng::new(2001),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.sync_agents();
        s.record_state(0.0);
        s
    }
}

impl TriboliumSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Larvae L_t (T. castaneum)",
            "Pupae P_t (T. castaneum)",
            "Adults A_t (T. castaneum)",
            "Total Biomass Count",
            "Competitor T. confusum Adults",
        ]
    }

    fn record_state(&mut self, pcg: f64) {
        let total = self.larvae + self.pupae + self.adults;
        let comp = if self.park_competition_mode {
            self.confusum_adults
        } else {
            0.0
        };
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.step_index as f64,
                values: [self.larvae, self.pupae, self.adults, total, comp],
                per_capita_growth: pcg,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        let prev_total = (self.larvae + self.pupae + self.adults).max(1.0);
        let comp_effect = if self.park_competition_mode {
            self.interspecific_alpha * self.confusum_adults
        } else {
            0.0
        };

        let next_l = self.b_recruitment
            * self.adults
            * (-self.c_el * self.larvae - self.c_ea * (self.adults + comp_effect)).exp();
        let next_p = (1.0 - self.mu_l) * self.larvae;
        let next_a = self.pupae * (-self.c_pa * (self.adults + 0.5 * comp_effect)).exp()
            + (1.0 - self.mu_a) * self.adults;

        if self.park_competition_mode {
            let r_conf = 1.48;
            let k_conf: f64 = 190.0;
            let next_conf = self.confusum_adults
                * (r_conf
                    * (1.0
                        - (self.confusum_adults + 0.94 * self.adults) / k_conf.max(10.0))
                    .exp()
                    .clamp(0.15, 2.8));
            self.confusum_adults = next_conf.clamp(0.0, 600.0);
        }

        self.larvae = next_l.clamp(0.1, 1500.0);
        self.pupae = next_p.clamp(0.1, 1200.0);
        self.adults = next_a.clamp(0.5, 1000.0);
        self.step_index += 1;

        let total = self.larvae + self.pupae + self.adults;
        let pcg = (total / prev_total).ln();

        self.current_phase = if self.c_pa >= 0.25 && self.mu_a >= 0.90 {
            GrowthPhase::ChaoticRegime
        } else if pcg.abs() < 0.015 {
            GrowthPhase::StationaryCarryingCapacity
        } else if pcg > 0.15 {
            GrowthPhase::ExponentialLog
        } else if pcg < -0.20 {
            GrowthPhase::PopulationCrashBottleneck
        } else {
            GrowthPhase::LimitCycleOscillation
        };

        self.record_state(pcg);
        self.sync_agents();
    }

    pub fn inject_perturbation(&mut self) {
        // Sieving flour medium: harvest 60% of adults, triggering larval cannibalism release!
        self.adults = (self.adults * 0.40).max(2.0);
        self.sync_agents();
    }

    fn sync_agents(&mut self) {
        // 3-depth flour vial column: Surface/Upper = Adults (kind 2), Mid = Larvae (kind 0), Lower = Pupae (kind 1)
        let n_l = (self.larvae * 0.35).round().clamp(1.0, 120.0) as usize;
        let n_p = (self.pupae * 0.35).round().clamp(1.0, 90.0) as usize;
        let n_a = (self.adults * 0.50).round().clamp(1.0, 100.0) as usize;

        self.agents.clear();
        for _ in 0..n_l {
            self.agents.push(Agent2D {
                x: self.rng.range_f64(0.06, 0.94),
                y: self.rng.range_f64(0.32, 0.72),
                vx: self.rng.range_f64(-0.01, 0.01),
                vy: self.rng.range_f64(-0.006, 0.006),
                kind: 0,
                aux: 0.5,
            });
        }
        for _ in 0..n_p {
            self.agents.push(Agent2D {
                x: self.rng.range_f64(0.06, 0.94),
                y: self.rng.range_f64(0.68, 0.95),
                vx: 0.0,
                vy: 0.0,
                kind: 1,
                aux: 0.8,
            });
        }
        for _ in 0..n_a {
            self.agents.push(Agent2D {
                x: self.rng.range_f64(0.06, 0.94),
                y: self.rng.range_f64(0.06, 0.42),
                vx: self.rng.range_f64(-0.015, 0.015),
                vy: self.rng.range_f64(-0.008, 0.008),
                kind: 2,
                aux: 1.0,
            });
        }
    }
}

// ============================================================================
// I2: Huffaker (1958) Spatial Citrus-Tray Mites (6x4 = 24 Orange Patches)
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HuffakerMiteSim {
    pub t: f64,
    pub dt: f64,
    pub cols: usize,
    pub rows: usize,
    pub prey_patches: Vec<f64>,
    pub pred_patches: Vec<f64>,
    pub r_prey: f64,
    pub k_patch: f64,
    pub attack_rate: f64,
    pub handling_time: f64,
    pub pred_conversion: f64,
    pub pred_mortality: f64,
    pub prey_dispersal: f64,
    pub pred_dispersal: f64,
    /// `0.0` = open tray (rapid global crash), `0.78` = Vaseline maze + silk launch posts (persistent waves!)
    pub vaseline_barrier: f64,
    pub current_phase: GrowthPhase,
    pub agents: Vec<Agent2D>,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for HuffakerMiteSim {
    fn default() -> Self {
        let cols = 6;
        let rows = 4;
        let n = cols * rows;
        let mut rng = SimRng::new(1958);
        let mut prey_patches = vec![0.0; n];
        let mut pred_patches = vec![0.0; n];
        for i in 0..n {
            if i % 3 == 0 {
                prey_patches[i] = rng.range_f64(25.0, 90.0);
            } else if i % 5 == 0 {
                prey_patches[i] = 15.0;
                pred_patches[i] = 6.0;
            }
        }
        pred_patches[0] = 8.0;
        let mut s = Self {
            t: 0.0,
            dt: 0.1,
            cols,
            rows,
            prey_patches,
            pred_patches,
            r_prey: 0.85,
            k_patch: 220.0,
            attack_rate: 0.045,
            handling_time: 0.12,
            pred_conversion: 0.28,
            pred_mortality: 0.36,
            prey_dispersal: 0.22,
            pred_dispersal: 0.09,
            vaseline_barrier: 0.72,
            current_phase: GrowthPhase::ExponentialLog,
            agents: Vec::new(),
            rng,
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.sync_agents();
        s.record_state(0.0);
        s
    }
}

impl HuffakerMiteSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Total Six-Spotted Mites (Prey)",
            "Total Predatory Mites (x10)",
            "Occupied Prey Oranges (of 24)",
            "Occupied Predator Oranges",
            "",
        ]
    }

    fn totals(&self) -> (f64, f64, f64, f64) {
        let tot_prey: f64 = self.prey_patches.iter().sum();
        let tot_pred: f64 = self.pred_patches.iter().sum();
        let occ_prey = self.prey_patches.iter().filter(|&&x| x >= 1.0).count() as f64;
        let occ_pred = self.pred_patches.iter().filter(|&&x| x >= 0.5).count() as f64;
        (tot_prey, tot_pred, occ_prey, occ_pred)
    }

    fn record_state(&mut self, pcg: f64) {
        let (tot_prey, tot_pred, occ_prey, occ_pred) = self.totals();
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.t,
                values: [tot_prey, tot_pred * 10.0, occ_prey, occ_pred, 0.0],
                per_capita_growth: pcg,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        let (prev_prey, _, _, _) = self.totals();
        let dt = self.dt;
        let n_cells = self.cols * self.rows;
        let mut next_prey = self.prey_patches.clone();
        let mut next_pred = self.pred_patches.clone();

        // Prey mites balloon across Vaseline barriers via silk threads better than walking predatory mites!
        let eff_dn = self.prey_dispersal * (1.0 - 0.55 * self.vaseline_barrier);
        let eff_dp = self.pred_dispersal * (1.0 - 0.92 * self.vaseline_barrier);

        for r in 0..self.rows {
            for c in 0..self.cols {
                let idx = r * self.cols + c;
                let n = self.prey_patches[idx];
                let p = self.pred_patches[idx];

                let mut lap_n = 0.0;
                let mut lap_p = 0.0;
                let neighbors = [
                    (r.wrapping_sub(1), c, r > 0),
                    (r + 1, c, r + 1 < self.rows),
                    (r, c.wrapping_sub(1), c > 0),
                    (r, c + 1, c + 1 < self.cols),
                ];
                for &(nr, nc, valid) in &neighbors {
                    if valid {
                        let nidx = nr * self.cols + nc;
                        lap_n += self.prey_patches[nidx] - n;
                        lap_p += self.pred_patches[nidx] - p;
                    }
                }

                let fn_rate = (self.attack_rate * n) / (1.0 + self.attack_rate * self.handling_time * n);
                let dn = self.r_prey * n * (1.0 - n / self.k_patch.max(10.0)) - fn_rate * p + eff_dn * lap_n;
                let dp = self.pred_conversion * fn_rate * p - self.pred_mortality * p + eff_dp * lap_p;

                let mut nn = (n + dt * dn).max(0.0);
                let mut np = (p + dt * dp).max(0.0);
                // Local extinction threshold on depleted single orange
                if nn < 0.08 && self.rng.next_f64() < 0.15 {
                    nn = 0.0;
                }
                if np < 0.05 && self.rng.next_f64() < 0.18 {
                    np = 0.0;
                }
                // Occasional aerial silk colonization if prey present elsewhere
                if nn == 0.0 && prev_prey > 30.0 && self.rng.next_f64() < 0.015 * (1.0 - 0.3 * self.vaseline_barrier) {
                    nn = 2.5;
                }
                next_prey[idx] = nn.min(self.k_patch * 1.5);
                next_pred[idx] = np.min(self.k_patch);
            }
        }
        let _ = n_cells;
        self.prey_patches = next_prey;
        self.pred_patches = next_pred;
        self.t += dt;

        let (tot_prey, tot_pred, occ_prey, _) = self.totals();
        let pcg = if prev_prey > 1e-3 && tot_prey > 1e-3 {
            (tot_prey / prev_prey).ln() / dt
        } else {
            0.0
        };

        self.current_phase = if tot_prey < 5.0 {
            GrowthPhase::PopulationCrashBottleneck
        } else if tot_pred > tot_prey * 0.25 && pcg < 0.0 {
            GrowthPhase::PredatorOvershoot
        } else if occ_prey >= 6.0 && tot_pred > 2.0 {
            GrowthPhase::LimitCycleOscillation
        } else {
            GrowthPhase::ExponentialLog
        };

        self.record_state(pcg);
        self.sync_agents();
    }

    pub fn inject_perturbation(&mut self) {
        // Reseed 3 corner oranges with fresh six-spotted mite colonies
        let last = self.cols * self.rows - 1;
        self.prey_patches[0] += 45.0;
        self.prey_patches[self.cols - 1] += 45.0;
        self.prey_patches[last] += 45.0;
        self.sync_agents();
    }

    fn sync_agents(&mut self) {
        self.agents.clear();
        for r in 0..self.rows {
            for c in 0..self.cols {
                let idx = r * self.cols + c;
                let cx = (c as f64 + 0.5) / (self.cols as f64);
                let cy = (r as f64 + 0.5) / (self.rows as f64);
                let n_dots = (self.prey_patches[idx] / 22.0).round().clamp(0.0, 9.0) as usize;
                let p_dots = (self.pred_patches[idx] / 4.0).round().clamp(0.0, 6.0) as usize;
                for _ in 0..n_dots {
                    self.agents.push(Agent2D {
                        x: cx + self.rng.range_f64(-0.055, 0.055),
                        y: cy + self.rng.range_f64(-0.075, 0.075),
                        vx: 0.0,
                        vy: 0.0,
                        kind: 0,
                        aux: idx as f64,
                    });
                }
                for _ in 0..p_dots {
                    self.agents.push(Agent2D {
                        x: cx + self.rng.range_f64(-0.05, 0.05),
                        y: cy + self.rng.range_f64(-0.07, 0.07),
                        vx: 0.0,
                        vy: 0.0,
                        kind: 1,
                        aux: idx as f64,
                    });
                }
            }
        }
    }
}

// ============================================================================
// I3: Pea Aphid (Acyrthosiphon pisum) & Parasitoid Wasp (Aphidius ervi)
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AphidParasitoidSim {
    pub generation: usize,
    pub aphids: f64,
    pub wasps: f64,
    pub mummies: f64,
    pub lambda_host: f64,
    pub k_host: f64,
    pub search_area_a: f64,
    /// Negative-binomial aggregation `k`: `k < 1` stabilizes Nicholson-Bailey via spatial clumping!
    pub clumping_k: f64,
    pub stem_refuge_fraction: f64,
    pub wasp_emergence_c: f64,
    pub current_phase: GrowthPhase,
    pub agents: Vec<Agent2D>,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for AphidParasitoidSim {
    fn default() -> Self {
        let mut s = Self {
            generation: 0,
            aphids: 140.0,
            wasps: 28.0,
            mummies: 15.0,
            lambda_host: 2.2,
            k_host: 1200.0,
            search_area_a: 0.032,
            clumping_k: 0.78,
            stem_refuge_fraction: 0.12,
            wasp_emergence_c: 0.92,
            current_phase: GrowthPhase::ExponentialLog,
            agents: Vec::new(),
            rng: SimRng::new(3003),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.sync_agents();
        push_record(
            &mut s.history,
            &mut s.transitions,
            TimeRecord {
                t: 0.0,
                values: [s.aphids, s.wasps, s.mummies, 0.0, 0.0],
                per_capita_growth: 0.0,
                phase: s.current_phase,
            },
        );
        s
    }
}

impl AphidParasitoidSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Pea Aphids H_t",
            "Adult Aphidius Wasps W_t",
            "Parasitized Mummies",
            "",
            "",
        ]
    }

    pub fn step(&mut self) {
        let prev_h = self.aphids.max(1e-3);
        let k = self.clumping_k.max(0.08);
        let nb_escape = (1.0 + (self.search_area_a * self.wasps) / k).powf(-k);
        let escape_frac =
            (self.stem_refuge_fraction + (1.0 - self.stem_refuge_fraction) * nb_escape).clamp(0.0, 1.0);
        let parasitized = self.aphids * (1.0 - escape_frac);

        let density_factor = (-(self.aphids / self.k_host.max(50.0)) * (self.lambda_host.ln().max(0.1))).exp();
        let next_h = self.lambda_host * self.aphids * density_factor * escape_frac;
        let next_w = self.wasp_emergence_c * parasitized;

        self.mummies = parasitized;
        self.aphids = next_h.clamp(0.5, 4000.0);
        self.wasps = next_w.clamp(0.1, 3000.0);
        self.generation += 1;

        let pcg = (self.aphids / prev_h).ln();
        self.current_phase = if self.aphids < 15.0 {
            GrowthPhase::PopulationCrashBottleneck
        } else if self.wasps > self.aphids * 0.65 {
            GrowthPhase::PredatorOvershoot
        } else if pcg.abs() < 0.03 {
            GrowthPhase::StationaryCarryingCapacity
        } else if pcg > 0.18 {
            GrowthPhase::ExponentialLog
        } else {
            GrowthPhase::LimitCycleOscillation
        };

        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.generation as f64,
                values: [self.aphids, self.wasps, self.mummies, 0.0, 0.0],
                per_capita_growth: pcg,
                phase: self.current_phase,
            },
        );
        self.sync_agents();
    }

    pub fn inject_perturbation(&mut self) {
        // Biological control release of 60 Aphidius ervi parasitoid wasps
        self.wasps = (self.wasps + 60.0).min(2000.0);
        self.sync_agents();
    }

    fn sync_agents(&mut self) {
        // 5 vertical legume stems at x = 0.16, 0.33, 0.50, 0.67, 0.84
        let n_h = (self.aphids / 8.0).round().clamp(1.0, 150.0) as usize;
        let n_w = (self.wasps / 6.0).round().clamp(1.0, 65.0) as usize;
        let n_m = (self.mummies / 10.0).round().clamp(0.0, 50.0) as usize;
        self.agents.clear();
        for i in 0..n_h {
            let stem_x = 0.16 + 0.17 * ((i % 5) as f64);
            self.agents.push(Agent2D {
                x: stem_x + self.rng.range_f64(-0.03, 0.03),
                y: self.rng.range_f64(0.12, 0.88),
                vx: 0.0,
                vy: 0.0,
                kind: 0,
                aux: if (i as f64 / n_h as f64) < self.stem_refuge_fraction {
                    1.0
                } else {
                    0.0
                },
            });
        }
        for _ in 0..n_w {
            self.agents.push(Agent2D {
                x: self.rng.range_f64(0.06, 0.94),
                y: self.rng.range_f64(0.08, 0.90),
                vx: self.rng.range_f64(-0.02, 0.02),
                vy: self.rng.range_f64(-0.02, 0.02),
                kind: 1,
                aux: 1.0,
            });
        }
        for i in 0..n_m {
            let stem_x = 0.16 + 0.17 * ((i % 5) as f64);
            self.agents.push(Agent2D {
                x: stem_x + self.rng.range_f64(-0.025, 0.025),
                y: self.rng.range_f64(0.18, 0.86),
                vx: 0.0,
                vy: 0.0,
                kind: 2,
                aux: 0.5,
            });
        }
    }
}
