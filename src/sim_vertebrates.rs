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

use crate::model::{Agent2D, FunctionalResponse, GrowthPhase, PhaseTransitionEvent, TimeRecord};
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
// V1: Snowshoe Hare (Lepus americanus) & Canada Lynx (Lynx canadensis)
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HareLynxSim {
    pub t: f64,
    pub dt: f64,
    pub hares: f64,
    pub lynx: f64,
    pub r_hare: f64,
    pub k_hare: f64,
    pub attack_rate: f64,
    pub handling_time: f64,
    pub conversion_eff: f64,
    pub lynx_mortality: f64,
    pub functional_response: FunctionalResponse,
    pub demographic_noise: bool,
    pub noise_sigma: f64,
    pub spatial_dispersal: f64,
    pub current_phase: GrowthPhase,
    pub agents: Vec<Agent2D>,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for HareLynxSim {
    fn default() -> Self {
        let mut s = Self {
            t: 0.0,
            dt: 0.05,
            hares: 45.0,
            lynx: 9.0,
            r_hare: 0.92,
            k_hare: 180.0,
            attack_rate: 0.028,
            handling_time: 0.32,
            conversion_eff: 0.24,
            lynx_mortality: 0.42,
            functional_response: FunctionalResponse::TypeIIHyperbolic,
            demographic_noise: false,
            noise_sigma: 0.04,
            spatial_dispersal: 0.35,
            current_phase: GrowthPhase::ExponentialLog,
            agents: Vec::new(),
            rng: SimRng::new(1001),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.sync_agents();
        s.record_initial();
        s
    }
}

impl HareLynxSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Snowshoe Hares (density/km^2)",
            "Canada Lynx (density/100 km^2)",
            "Predation Rate f(H)*L",
            "",
            "",
        ]
    }

    fn derivatives(&self, h: f64, l: f64) -> (f64, f64) {
        let fh = self
            .functional_response
            .evaluate(h, self.attack_rate, self.handling_time);
        let dh = self.r_hare * h * (1.0 - h / self.k_hare.max(1.0)) - fh * l;
        let dl = self.conversion_eff * fh * l - self.lynx_mortality * l;
        (dh, dl)
    }

    fn classify_phase(&self, dh: f64, dl: f64) -> GrowthPhase {
        if self.hares < 0.08 * self.k_hare {
            GrowthPhase::PopulationCrashBottleneck
        } else if dl > 0.02 && dh < 0.0 {
            GrowthPhase::PredatorOvershoot
        } else if dh > 0.0 && self.hares < 0.45 * self.k_hare {
            GrowthPhase::ExponentialLog
        } else if dh > 0.0 && self.hares >= 0.45 * self.k_hare {
            GrowthPhase::DeceleratingTransition
        } else if dh.abs() < 0.05 && dl.abs() < 0.02 {
            GrowthPhase::StationaryCarryingCapacity
        } else {
            GrowthPhase::LimitCycleOscillation
        }
    }

    fn record_initial(&mut self) {
        let (dh, dl) = self.derivatives(self.hares, self.lynx);
        self.current_phase = self.classify_phase(dh, dl);
        let pred_flux = self
            .functional_response
            .evaluate(self.hares, self.attack_rate, self.handling_time)
            * self.lynx;
        let pcg = if self.hares > 1e-6 {
            dh / self.hares
        } else {
            0.0
        };
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.t,
                values: [self.hares, self.lynx, pred_flux, 0.0, 0.0],
                per_capita_growth: pcg,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        let dt = self.dt;
        let (h0, l0) = (self.hares, self.lynx);
        let (k1h, k1l) = self.derivatives(h0, l0);
        let (k2h, k2l) = self.derivatives(h0 + 0.5 * dt * k1h, l0 + 0.5 * dt * k1l);
        let (k3h, k3l) = self.derivatives(h0 + 0.5 * dt * k2h, l0 + 0.5 * dt * k2l);
        let (k4h, k4l) = self.derivatives(h0 + dt * k3h, l0 + dt * k3l);

        let mut h_next = h0 + (dt / 6.0) * (k1h + 2.0 * k2h + 2.0 * k3h + k4h);
        let mut l_next = l0 + (dt / 6.0) * (k1l + 2.0 * k2l + 2.0 * k3l + k4l);

        if self.demographic_noise {
            let sq_dt = dt.sqrt();
            h_next += self.noise_sigma * h_next.max(0.0).sqrt() * sq_dt * self.rng.standard_normal();
            l_next += self.noise_sigma * l_next.max(0.0).sqrt() * sq_dt * self.rng.standard_normal();
        }

        self.hares = h_next.clamp(0.05, self.k_hare * 3.0);
        self.lynx = l_next.clamp(0.02, self.k_hare * 1.5);
        self.t += dt;

        let (dh, dl) = self.derivatives(self.hares, self.lynx);
        self.current_phase = self.classify_phase(dh, dl);
        let pred_flux = self
            .functional_response
            .evaluate(self.hares, self.attack_rate, self.handling_time)
            * self.lynx;
        let pcg = dh / self.hares.max(1e-6);

        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.t,
                values: [self.hares, self.lynx, pred_flux, 0.0, 0.0],
                per_capita_growth: pcg,
                phase: self.current_phase,
            },
        );
        self.update_spatial_agents();
    }

    pub fn inject_perturbation(&mut self) {
        // Boreal willow browse dieback / severe icing winter reduces hare abundance by 45%
        self.hares = (self.hares * 0.55).max(1.0);
        self.sync_agents();
    }

    fn sync_agents(&mut self) {
        let target_hares = (self.hares * 1.1).round().clamp(1.0, 220.0) as usize;
        let target_lynx = (self.lynx * 1.8).round().clamp(1.0, 75.0) as usize;

        let mut h_agents: Vec<Agent2D> = self
            .agents
            .iter()
            .filter(|a| a.kind == 0)
            .cloned()
            .collect();
        let mut l_agents: Vec<Agent2D> = self
            .agents
            .iter()
            .filter(|a| a.kind == 1)
            .cloned()
            .collect();

        while h_agents.len() < target_hares {
            h_agents.push(Agent2D {
                x: self.rng.next_f64(),
                y: self.rng.next_f64(),
                vx: self.rng.range_f64(-0.015, 0.015),
                vy: self.rng.range_f64(-0.015, 0.015),
                kind: 0,
                aux: 1.0,
            });
        }
        h_agents.truncate(target_hares);

        while l_agents.len() < target_lynx {
            l_agents.push(Agent2D {
                x: self.rng.next_f64(),
                y: self.rng.next_f64(),
                vx: self.rng.range_f64(-0.022, 0.022),
                vy: self.rng.range_f64(-0.022, 0.022),
                kind: 1,
                aux: 1.0,
            });
        }
        l_agents.truncate(target_lynx);

        self.agents.clear();
        self.agents.extend(h_agents);
        self.agents.extend(l_agents);
    }

    fn update_spatial_agents(&mut self) {
        self.sync_agents();
        let disp = self.spatial_dispersal;
        let snapshot = self.agents.clone();
        for agent in &mut self.agents {
            let mut steer_x = 0.0;
            let mut steer_y = 0.0;
            if agent.kind == 0 {
                // Hares evade nearby lynx
                for other in &snapshot {
                    if other.kind == 1 {
                        let dx = agent.x - other.x;
                        let dy = agent.y - other.y;
                        let d2 = dx * dx + dy * dy;
                        if d2 > 1e-6 && d2 < 0.025 {
                            let inv = 0.0025 / d2.sqrt();
                            steer_x += dx * inv;
                            steer_y += dy * inv;
                        }
                    }
                }
            } else {
                // Lynx pursue nearest hare
                let mut best_d2 = 0.06;
                for other in &snapshot {
                    if other.kind == 0 {
                        let dx = other.x - agent.x;
                        let dy = other.y - agent.y;
                        let d2 = dx * dx + dy * dy;
                        if d2 < best_d2 && d2 > 1e-6 {
                            best_d2 = d2;
                            let inv = 0.004 / d2.sqrt();
                            steer_x = dx * inv;
                            steer_y = dy * inv;
                        }
                    }
                }
            }
            agent.vx = 0.86 * agent.vx + steer_x + self.rng.range_f64(-0.006, 0.006) * disp;
            agent.vy = 0.86 * agent.vy + steer_y + self.rng.range_f64(-0.006, 0.006) * disp;
            let max_v = if agent.kind == 1 { 0.022 } else { 0.016 };
            agent.vx = agent.vx.clamp(-max_v, max_v);
            agent.vy = agent.vy.clamp(-max_v, max_v);
            agent.x = (agent.x + agent.vx).rem_euclid(1.0);
            agent.y = (agent.y + agent.vy).rem_euclid(1.0);
        }
    }
}

// ============================================================================
// V2: Isle Royale Moose (Alces alces) & Gray Wolf (Canis lupus)
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsleRoyaleSim {
    pub t: f64,
    pub dt: f64,
    pub moose: f64,
    pub wolves: f64,
    pub fir_browse: f64,
    pub r_moose: f64,
    pub k_moose_base: f64,
    pub kill_rate: f64,
    pub handling_time: f64,
    pub wolf_conversion: f64,
    pub wolf_mortality: f64,
    pub allee_threshold: f64,
    pub winter_severity: f64,
    pub ice_bridge_rescue: bool,
    pub current_phase: GrowthPhase,
    pub agents: Vec<Agent2D>,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for IsleRoyaleSim {
    fn default() -> Self {
        let mut s = Self {
            t: 1959.0,
            dt: 0.25,
            moose: 650.0,
            wolves: 22.0,
            fir_browse: 85.0,
            r_moose: 0.24,
            k_moose_base: 1850.0,
            kill_rate: 14.5,
            handling_time: 0.042,
            wolf_conversion: 0.026,
            wolf_mortality: 0.21,
            allee_threshold: 6.0,
            winter_severity: 0.35,
            ice_bridge_rescue: false,
            current_phase: GrowthPhase::ExponentialLog,
            agents: Vec::new(),
            rng: SimRng::new(1959),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.sync_agents();
        s.record_initial();
        s
    }
}

impl IsleRoyaleSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Moose Abundance (individuals)",
            "Gray Wolves (individuals x 20 scaled)",
            "Actual Wolf Count",
            "Balsam Fir Browse Index (%)",
            "",
        ]
    }

    fn derivatives(&self, m: f64, w: f64, fir: f64) -> (f64, f64, f64) {
        let k_eff = (self.k_moose_base * (fir / 85.0).clamp(0.35, 1.35)).max(150.0);
        let r_eff = self.r_moose * (1.0 - 0.38 * self.winter_severity);
        // Ratio-dependent predation g(M/W)
        let denom = (w + self.kill_rate * self.handling_time * m).max(1e-6);
        let total_kills = (self.kill_rate * m * w) / denom;
        let dm = r_eff * m * (1.0 - m / k_eff) - total_kills;

        // Pack Allee / inbreeding factor when W drops near allee_threshold
        let allee = if w <= 0.1 {
            0.0
        } else {
            (w / (w + self.allee_threshold)).clamp(0.0, 1.0)
        };
        let rescue_inflow = if self.ice_bridge_rescue && w < 8.0 {
            0.45
        } else {
            0.0
        };
        let dw = self.wolf_conversion * total_kills * allee - self.wolf_mortality * w + rescue_inflow;
        let dfir = 0.14 * fir * (1.0 - fir / 100.0) - 0.0085 * m * (fir / 100.0);
        (dm, dw, dfir)
    }

    fn classify_phase(&self, dm: f64, dw: f64) -> GrowthPhase {
        if self.wolves < self.allee_threshold {
            GrowthPhase::PopulationCrashBottleneck
        } else if dw > 0.3 && dm < 0.0 {
            GrowthPhase::PredatorOvershoot
        } else if dm > 15.0 && self.moose < 0.55 * self.k_moose_base {
            GrowthPhase::ExponentialLog
        } else if dm > 0.0 {
            GrowthPhase::DeceleratingTransition
        } else {
            GrowthPhase::LimitCycleOscillation
        }
    }

    fn record_initial(&mut self) {
        let (dm, dw, _) = self.derivatives(self.moose, self.wolves, self.fir_browse);
        self.current_phase = self.classify_phase(dm, dw);
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.t,
                values: [
                    self.moose,
                    self.wolves * 20.0,
                    self.wolves,
                    self.fir_browse,
                    0.0,
                ],
                per_capita_growth: dm / self.moose.max(1.0),
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        let dt = self.dt;
        let (m0, w0, f0) = (self.moose, self.wolves, self.fir_browse);
        let (k1m, k1w, k1f) = self.derivatives(m0, w0, f0);
        let (k2m, k2w, k2f) = self.derivatives(
            m0 + 0.5 * dt * k1m,
            w0 + 0.5 * dt * k1w,
            f0 + 0.5 * dt * k1f,
        );
        let (k3m, k3w, k3f) = self.derivatives(
            m0 + 0.5 * dt * k2m,
            w0 + 0.5 * dt * k2w,
            f0 + 0.5 * dt * k2f,
        );
        let (k4m, k4w, k4f) = self.derivatives(m0 + dt * k3m, w0 + dt * k3w, f0 + dt * k3f);

        self.moose = (m0 + (dt / 6.0) * (k1m + 2.0 * k2m + 2.0 * k3m + k4m)).clamp(10.0, 3500.0);
        self.wolves = (w0 + (dt / 6.0) * (k1w + 2.0 * k2w + 2.0 * k3w + k4w)).clamp(0.0, 90.0);
        self.fir_browse =
            (f0 + (dt / 6.0) * (k1f + 2.0 * k2f + 2.0 * k3f + k4f)).clamp(10.0, 100.0);
        self.t += dt;

        let (dm, dw, _) = self.derivatives(self.moose, self.wolves, self.fir_browse);
        self.current_phase = self.classify_phase(dm, dw);
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.t,
                values: [
                    self.moose,
                    self.wolves * 20.0,
                    self.wolves,
                    self.fir_browse,
                    0.0,
                ],
                per_capita_growth: dm / self.moose.max(1.0),
                phase: self.current_phase,
            },
        );
        self.sync_agents();
    }

    pub fn inject_perturbation(&mut self) {
        // Genetic rescue translocation (+8 wolves introduced across packs)
        self.wolves = (self.wolves + 8.0).min(80.0);
        self.sync_agents();
    }

    fn sync_agents(&mut self) {
        // Elongated Isle Royale island silhouette inside [0.1..0.9] x [0.28..0.72]
        let target_moose = (self.moose / 12.0).round().clamp(2.0, 200.0) as usize;
        let target_wolves = self.wolves.round().clamp(0.0, 65.0) as usize;

        self.agents.retain(|a| a.kind <= 1);
        let mut m_list: Vec<Agent2D> = self
            .agents
            .iter()
            .filter(|a| a.kind == 0)
            .cloned()
            .collect();
        let mut w_list: Vec<Agent2D> = self
            .agents
            .iter()
            .filter(|a| a.kind == 1)
            .cloned()
            .collect();

        while m_list.len() < target_moose {
            let u = self.rng.range_f64(0.10, 0.90);
            let v = 0.50 + 0.16 * (u - 0.5) + self.rng.range_f64(-0.14, 0.14);
            m_list.push(Agent2D {
                x: u,
                y: v.clamp(0.22, 0.78),
                vx: self.rng.range_f64(-0.008, 0.008),
                vy: self.rng.range_f64(-0.005, 0.005),
                kind: 0,
                aux: self.fir_browse / 100.0,
            });
        }
        m_list.truncate(target_moose);

        while w_list.len() < target_wolves {
            let pack = (w_list.len() % 3) as f64;
            let cx = 0.25 + 0.25 * pack;
            let cy = 0.42 + 0.08 * pack;
            w_list.push(Agent2D {
                x: (cx + self.rng.range_f64(-0.06, 0.06)).clamp(0.12, 0.88),
                y: (cy + self.rng.range_f64(-0.05, 0.05)).clamp(0.25, 0.75),
                vx: self.rng.range_f64(-0.012, 0.012),
                vy: self.rng.range_f64(-0.008, 0.008),
                kind: 1,
                aux: pack,
            });
        }
        w_list.truncate(target_wolves);

        for m in &mut m_list {
            m.vx = (0.85 * m.vx + self.rng.range_f64(-0.005, 0.005)).clamp(-0.012, 0.012);
            m.vy = (0.85 * m.vy + self.rng.range_f64(-0.004, 0.004)).clamp(-0.008, 0.008);
            m.x = (m.x + m.vx).clamp(0.10, 0.90);
            let center_y = 0.50 + 0.16 * (m.x - 0.5);
            m.y = (m.y + m.vy).clamp(center_y - 0.15, center_y + 0.15);
        }
        for w in &mut w_list {
            let pack_cx = 0.25 + 0.25 * w.aux;
            let pack_cy = 0.42 + 0.08 * w.aux;
            w.vx = (0.82 * w.vx + 0.04 * (pack_cx - w.x) + self.rng.range_f64(-0.008, 0.008))
                .clamp(-0.018, 0.018);
            w.vy = (0.82 * w.vy + 0.04 * (pack_cy - w.y) + self.rng.range_f64(-0.006, 0.006))
                .clamp(-0.014, 0.014);
            w.x = (w.x + w.vx).clamp(0.11, 0.89);
            w.y = (w.y + w.vy).clamp(0.24, 0.76);
        }

        self.agents.clear();
        self.agents.extend(m_list);
        self.agents.extend(w_list);
    }
}

// ============================================================================
// V3: St. Kilda Soay Sheep (Ovis aries) Overcompensatory Dynamics
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SoaySheepSim {
    pub year: usize,
    pub population: f64,
    pub prev_population: f64,
    pub lambda_fecundity: f64,
    pub carrying_scale_k: f64,
    pub overcompensation_beta: f64,
    pub gale_probability: f64,
    pub gale_severity: f64,
    pub last_storm_mortality: f64,
    pub current_phase: GrowthPhase,
    pub agents: Vec<Agent2D>,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for SoaySheepSim {
    fn default() -> Self {
        let mut s = Self {
            year: 1985,
            population: 620.0,
            prev_population: 620.0,
            lambda_fecundity: 2.35,
            carrying_scale_k: 1150.0,
            overcompensation_beta: 3.4,
            gale_probability: 0.28,
            gale_severity: 0.42,
            last_storm_mortality: 0.0,
            current_phase: GrowthPhase::ExponentialLog,
            agents: Vec::new(),
            rng: SimRng::new(1985),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.sync_agents();
        push_record(
            &mut s.history,
            &mut s.transitions,
            TimeRecord {
                t: s.year as f64,
                values: [s.population, s.carrying_scale_k, 0.0, 0.0, 0.0],
                per_capita_growth: 0.0,
                phase: s.current_phase,
            },
        );
        s
    }
}

impl SoaySheepSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Hirta Soay Sheep Census",
            "Forage Carrying Scale K",
            "Winter Gale Severity (%)",
            "",
            "",
        ]
    }

    pub fn step(&mut self) {
        self.prev_population = self.population;
        let storm_hit = self.rng.next_f64() < self.gale_probability;
        let w_storm = if storm_hit {
            self.gale_severity * self.rng.range_f64(0.75, 1.25)
        } else {
            0.03 * self.rng.next_f64()
        };
        self.last_storm_mortality = (1.0 - (-w_storm).exp()) * 100.0;

        let density_denom =
            1.0 + (self.population / self.carrying_scale_k.max(50.0)).powf(self.overcompensation_beta);
        let next_n = (self.population * self.lambda_fecundity * (-w_storm).exp()) / density_denom;
        self.population = next_n.clamp(25.0, 3500.0);
        self.year += 1;

        let pcg = (self.population / self.prev_population.max(1.0)).ln();
        self.current_phase = if pcg < -0.28 {
            GrowthPhase::PopulationCrashBottleneck
        } else if pcg > 0.18 && self.population < 0.85 * self.carrying_scale_k {
            GrowthPhase::ExponentialLog
        } else if self.population > 1.15 * self.carrying_scale_k {
            GrowthPhase::PredatorOvershoot
        } else if pcg.abs() < 0.05 {
            GrowthPhase::StationaryCarryingCapacity
        } else if self.overcompensation_beta > 3.8 {
            GrowthPhase::ChaoticRegime
        } else {
            GrowthPhase::LimitCycleOscillation
        };

        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.year as f64,
                values: [
                    self.population,
                    self.carrying_scale_k,
                    self.last_storm_mortality,
                    0.0,
                    0.0,
                ],
                per_capita_growth: pcg,
                phase: self.current_phase,
            },
        );
        self.sync_agents();
    }

    pub fn inject_perturbation(&mut self) {
        // Trigger extreme North Atlantic winter gale immediately
        self.prev_population = self.population;
        self.population = (self.population * 0.48).max(40.0);
        self.last_storm_mortality = 52.0;
        self.current_phase = GrowthPhase::PopulationCrashBottleneck;
        self.sync_agents();
    }

    fn sync_agents(&mut self) {
        let count = (self.population / 9.0).round().clamp(4.0, 240.0) as usize;
        while self.agents.len() < count {
            self.agents.push(Agent2D {
                x: self.rng.range_f64(0.08, 0.92),
                y: self.rng.range_f64(0.08, 0.92),
                vx: self.rng.range_f64(-0.012, 0.012),
                vy: self.rng.range_f64(-0.012, 0.012),
                kind: 0,
                aux: 1.0,
            });
        }
        self.agents.truncate(count);
        for a in &mut self.agents {
            a.vx = (0.75 * a.vx + self.rng.range_f64(-0.015, 0.015)).clamp(-0.025, 0.025);
            a.vy = (0.75 * a.vy + self.rng.range_f64(-0.015, 0.015)).clamp(-0.025, 0.025);
            a.x = (a.x + a.vx).clamp(0.06, 0.94);
            a.y = (a.y + a.vy).clamp(0.06, 0.94);
        }
    }
}
