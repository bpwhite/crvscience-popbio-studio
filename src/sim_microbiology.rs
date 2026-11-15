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
// M1: Escherichia coli K-12 Batch Growth Phases & Monod Chemostat
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EcoliGrowthSim {
    pub t_hours: f64,
    pub dt: f64,
    /// Viable cell density (`10^7 CFU/mL` units)
    pub viable_cells: f64,
    /// Dead/lysed cell density (`10^7 cells/mL` units, contributes to OD600 turbidity)
    pub dead_cells: f64,
    /// Primary carbon source (Glucose, `g/L`)
    pub glucose_g_l: f64,
    /// Secondary carbon source (Lactose, `g/L`, active when `diauxic_mode = true`)
    pub lactose_g_l: f64,
    /// Intracellular ribosomal/enzymatic acclimation state `0..1` governing the Lag Phase
    pub lag_enzyme_state: f64,
    /// Lac-operon induction state `0..1` governing the Diauxic Lag between Glucose and Lactose
    pub lac_operon_state: f64,
    /// Accumulated toxic organic acid / waste inhibitor
    pub waste_inhibitor: f64,
    pub mu_max_glucose: f64,
    pub ks_glucose: f64,
    pub mu_max_lactose: f64,
    pub ks_lactose: f64,
    pub yield_coeff: f64,
    pub lag_rate_h: f64,
    pub death_rate_h: f64,
    pub chemostat_mode: bool,
    pub dilution_rate_d: f64,
    pub feed_glucose_s0: f64,
    pub diauxic_mode: bool,
    pub current_phase: GrowthPhase,
    pub agents: Vec<Agent2D>,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for EcoliGrowthSim {
    fn default() -> Self {
        let mut s = Self {
            t_hours: 0.0,
            dt: 0.08,
            viable_cells: 1.5,
            dead_cells: 0.05,
            glucose_g_l: 4.0,
            lactose_g_l: 3.0,
            lag_enzyme_state: 0.06,
            lac_operon_state: 0.02,
            waste_inhibitor: 0.0,
            mu_max_glucose: 0.95,
            ks_glucose: 0.12,
            mu_max_lactose: 0.58,
            ks_lactose: 0.25,
            yield_coeff: 28.0,
            lag_rate_h: 0.85,
            death_rate_h: 0.09,
            chemostat_mode: false,
            dilution_rate_d: 0.25,
            feed_glucose_s0: 4.0,
            diauxic_mode: false,
            current_phase: GrowthPhase::LagPhase,
            agents: Vec::new(),
            rng: SimRng::new(4001),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.sync_agents();
        s.record_state(0.0);
        s
    }
}

impl EcoliGrowthSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Viable Cells (10^7 CFU/mL)",
            "log10(Viable CFU/mL) x 20",
            "Total Turbidity (Viable + Dead)",
            "Glucose S1 (g/L x 25)",
            "Lactose S2 (g/L x 25)",
        ]
    }

    fn record_state(&mut self, net_mu: f64) {
        let log10_scaled = (self.viable_cells.max(0.01) * 1e7).log10() * 20.0;
        let total_od = self.viable_cells + self.dead_cells;
        let lac_val = if self.diauxic_mode {
            self.lactose_g_l * 25.0
        } else {
            0.0
        };
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.t_hours,
                values: [
                    self.viable_cells,
                    log10_scaled,
                    total_od,
                    self.glucose_g_l * 25.0,
                    lac_val,
                ],
                per_capita_growth: net_mu,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        let dt = self.dt;
        let dil = if self.chemostat_mode {
            self.dilution_rate_d
        } else {
            0.0
        };

        // Physiological acclimation out of Lag Phase
        self.lag_enzyme_state += dt * self.lag_rate_h * (1.0 - self.lag_enzyme_state);
        self.lag_enzyme_state = self.lag_enzyme_state.clamp(0.0, 1.0);

        // Catabolite repression: cAMP-CRP induces lac operon only when glucose drops low
        let glucose_repression = 0.08 / (0.08 + self.glucose_g_l);
        if self.diauxic_mode && self.lactose_g_l > 0.02 {
            self.lac_operon_state +=
                dt * 1.15 * glucose_repression * (1.0 - self.lac_operon_state);
            self.lac_operon_state = self.lac_operon_state.clamp(0.0, 1.0);
        }

        let tox_factor = 1.0 / (1.0 + 0.45 * self.waste_inhibitor);
        let mu_glc = self.lag_enzyme_state
            * self.mu_max_glucose
            * (self.glucose_g_l / (self.ks_glucose + self.glucose_g_l).max(1e-6))
            * tox_factor;
        let mu_lac = if self.diauxic_mode {
            self.lag_enzyme_state
                * self.lac_operon_state
                * self.mu_max_lactose
                * (self.lactose_g_l / (self.ks_lactose + self.lactose_g_l).max(1e-6))
                * glucose_repression
                * tox_factor
        } else {
            0.0
        };

        let gross_mu = mu_glc + mu_lac;
        let starvation_stress = if self.glucose_g_l + self.lactose_g_l * (self.diauxic_mode as u8 as f64) < 0.05 {
            1.0 + 1.8 * self.waste_inhibitor
        } else {
            0.15 * self.waste_inhibitor
        };
        let death_rate = if self.chemostat_mode {
            0.015
        } else {
            self.death_rate_h * starvation_stress
        };

        let net_mu = gross_mu - death_rate - dil;
        let d_viable = net_mu * self.viable_cells;
        let d_dead = death_rate * self.viable_cells - dil * self.dead_cells;

        let d_glc = -(mu_glc * self.viable_cells) / self.yield_coeff.max(1.0)
            + dil * (self.feed_glucose_s0 - self.glucose_g_l);
        let d_lac = if self.diauxic_mode {
            -(mu_lac * self.viable_cells) / self.yield_coeff.max(1.0) - dil * self.lactose_g_l
        } else {
            0.0
        };
        let d_waste = if self.chemostat_mode {
            0.002 * gross_mu * self.viable_cells - dil * self.waste_inhibitor
        } else {
            0.0035 * gross_mu * self.viable_cells + 0.02 * starvation_stress * dt
        };

        self.viable_cells = (self.viable_cells + dt * d_viable).clamp(0.01, 500.0);
        self.dead_cells = (self.dead_cells + dt * d_dead).clamp(0.0, 500.0);
        self.glucose_g_l = (self.glucose_g_l + dt * d_glc).clamp(0.0, 20.0);
        self.lactose_g_l = (self.lactose_g_l + dt * d_lac).clamp(0.0, 20.0);
        self.waste_inhibitor = (self.waste_inhibitor + dt * d_waste).clamp(0.0, 25.0);
        self.t_hours += dt;

        // Classify classic microbiology growth phase
        self.current_phase = if self.lag_enzyme_state < 0.55 {
            GrowthPhase::LagPhase
        } else if self.diauxic_mode
            && self.glucose_g_l < 0.08
            && self.lactose_g_l > 0.15
            && self.lac_operon_state < 0.65
        {
            GrowthPhase::DiauxicShift
        } else if net_mu > 0.35 * self.mu_max_glucose {
            GrowthPhase::ExponentialLog
        } else if net_mu > 0.04 {
            GrowthPhase::DeceleratingTransition
        } else if net_mu >= -0.025 {
            GrowthPhase::StationaryCarryingCapacity
        } else {
            GrowthPhase::DeathSenescence
        };

        self.record_state(net_mu);
        self.sync_agents();
    }

    pub fn inject_perturbation(&mut self) {
        // Fresh glucose nutrient bolus (+2.5 g/L) & partial medium refresh
        self.glucose_g_l = (self.glucose_g_l + 2.5).min(12.0);
        self.waste_inhibitor *= 0.4;
        self.sync_agents();
    }

    fn sync_agents(&mut self) {
        let n_live = (self.viable_cells * 1.15).round().clamp(2.0, 220.0) as usize;
        let n_dead = (self.dead_cells * 0.85).round().clamp(0.0, 90.0) as usize;

        let mut live_list: Vec<Agent2D> = self
            .agents
            .iter()
            .filter(|a| a.kind == 0)
            .cloned()
            .collect();
        let mut dead_list: Vec<Agent2D> = self
            .agents
            .iter()
            .filter(|a| a.kind == 2)
            .cloned()
            .collect();

        while live_list.len() < n_live {
            // Spawn near an existing dividing cell to form realistic microcolonies
            let (bx, by) = if !live_list.is_empty() {
                let parent = &live_list[self.rng.gen_index(live_list.len())];
                (
                    (parent.x + self.rng.range_f64(-0.035, 0.035)).clamp(0.05, 0.95),
                    (parent.y + self.rng.range_f64(-0.035, 0.035)).clamp(0.05, 0.95),
                )
            } else {
                (self.rng.range_f64(0.2, 0.8), self.rng.range_f64(0.2, 0.8))
            };
            live_list.push(Agent2D {
                x: bx,
                y: by,
                vx: self.rng.range_f64(-0.006, 0.006),
                vy: self.rng.range_f64(-0.006, 0.006),
                kind: 0,
                aux: self.rng.next_f64(), // rod orientation angle [0..1)
            });
        }
        live_list.truncate(n_live);

        while dead_list.len() < n_dead {
            dead_list.push(Agent2D {
                x: self.rng.range_f64(0.05, 0.95),
                y: self.rng.range_f64(0.05, 0.95),
                vx: 0.0,
                vy: 0.0,
                kind: 2,
                aux: self.rng.next_f64(),
            });
        }
        dead_list.truncate(n_dead);

        for c in &mut live_list {
            c.vx = (0.75 * c.vx + self.rng.range_f64(-0.005, 0.005)).clamp(-0.01, 0.01);
            c.vy = (0.75 * c.vy + self.rng.range_f64(-0.005, 0.005)).clamp(-0.01, 0.01);
            c.x = (c.x + c.vx).clamp(0.04, 0.96);
            c.y = (c.y + c.vy).clamp(0.04, 0.96);
        }

        self.agents.clear();
        self.agents.extend(live_list);
        self.agents.extend(dead_list);
    }
}

// ============================================================================
// M2: Gause (1934) Paramecium Competition & Didinium Predation
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GauseProtistSim {
    pub t_days: f64,
    pub dt: f64,
    /// `true` = P. aurelia vs P. caudatum competition; `false` = P. caudatum + Didinium nasutum predation
    pub competition_mode: bool,
    pub p_caudatum: f64,
    pub p_aurelia: f64,
    pub didinium: f64,
    pub r_caudatum: f64,
    pub k_caudatum: f64,
    pub r_aurelia: f64,
    pub k_aurelia: f64,
    pub alpha_caud_on_aur: f64,
    pub alpha_aur_on_caud: f64,
    pub didinium_attack: f64,
    pub didinium_handling: f64,
    pub sediment_refuge: f64,
    pub periodic_immigration: bool,
    pub current_phase: GrowthPhase,
    pub agents: Vec<Agent2D>,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for GauseProtistSim {
    fn default() -> Self {
        let mut s = Self {
            t_days: 0.0,
            dt: 0.1,
            competition_mode: false,
            p_caudatum: 32.0,
            p_aurelia: 28.0,
            didinium: 8.0,
            r_caudatum: 0.78,
            k_caudatum: 195.0,
            r_aurelia: 1.05,
            k_aurelia: 240.0,
            alpha_caud_on_aur: 0.62,
            alpha_aur_on_caud: 1.38,
            didinium_attack: 0.045,
            didinium_handling: 0.18,
            sediment_refuge: 0.22,
            periodic_immigration: true,
            current_phase: GrowthPhase::ExponentialLog,
            agents: Vec::new(),
            rng: SimRng::new(1934),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.sync_agents();
        s.record_state(0.0);
        s
    }
}

impl GauseProtistSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Paramecium caudatum (cells/0.5mL)",
            "P. aurelia (comp) OR Didinium (pred)",
            "Sediment Refuge Fraction (%)",
            "",
            "",
        ]
    }

    fn record_state(&mut self, pcg: f64) {
        let second = if self.competition_mode {
            self.p_aurelia
        } else {
            self.didinium
        };
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.t_days,
                values: [
                    self.p_caudatum,
                    second,
                    self.sediment_refuge * 100.0,
                    0.0,
                    0.0,
                ],
                per_capita_growth: pcg,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        let dt = self.dt;
        let prev_c = self.p_caudatum.max(0.01);

        if self.competition_mode {
            let dc = self.r_caudatum
                * self.p_caudatum
                * (1.0
                    - (self.p_caudatum + self.alpha_aur_on_caud * self.p_aurelia)
                        / self.k_caudatum.max(10.0));
            let da = self.r_aurelia
                * self.p_aurelia
                * (1.0
                    - (self.p_aurelia + self.alpha_caud_on_aur * self.p_caudatum)
                        / self.k_aurelia.max(10.0));
            self.p_caudatum = (self.p_caudatum + dt * dc).clamp(0.0, 400.0);
            self.p_aurelia = (self.p_aurelia + dt * da).clamp(0.0, 400.0);
            self.t_days += dt;

            let pcg = dc / prev_c;
            self.current_phase = if self.p_caudatum < 2.0 && self.p_aurelia > 150.0 {
                GrowthPhase::StationaryCarryingCapacity
            } else if dc > 0.0 && da > 0.0 {
                GrowthPhase::ExponentialLog
            } else {
                GrowthPhase::DeceleratingTransition
            };
            self.record_state(pcg);
        } else {
            let exposed = (1.0 - self.sediment_refuge) * self.p_caudatum;
            let pred_rate = (self.didinium_attack * exposed)
                / (1.0 + self.didinium_attack * self.didinium_handling * exposed);
            let dc = self.r_caudatum * self.p_caudatum * (1.0 - self.p_caudatum / self.k_caudatum)
                - pred_rate * self.didinium;
            let dd = 0.36 * pred_rate * self.didinium - 0.42 * self.didinium;

            self.p_caudatum = (self.p_caudatum + dt * dc).clamp(0.0, 400.0);
            self.didinium = (self.didinium + dt * dd).clamp(0.0, 250.0);

            // Gause's periodic 3-day 1-cell immigration experiment
            let prev_mod = (self.t_days % 3.0).abs();
            self.t_days += dt;
            if self.periodic_immigration && (self.t_days % 3.0) < prev_mod {
                self.p_caudatum += 1.5;
                self.didinium += 1.0;
            }

            let pcg = dc / prev_c;
            self.current_phase = if self.p_caudatum < 3.0 {
                GrowthPhase::PopulationCrashBottleneck
            } else if dd > 0.2 && dc < 0.0 {
                GrowthPhase::PredatorOvershoot
            } else if dc > 1.5 {
                GrowthPhase::ExponentialLog
            } else {
                GrowthPhase::LimitCycleOscillation
            };
            self.record_state(pcg);
        }
        self.sync_agents();
    }

    pub fn inject_perturbation(&mut self) {
        // Add oat-sediment refuge + pipette inoculation of 15 Paramecium caudatum
        self.p_caudatum = (self.p_caudatum + 15.0).min(300.0);
        if !self.competition_mode {
            self.didinium = (self.didinium + 4.0).min(150.0);
        }
        self.sync_agents();
    }

    fn sync_agents(&mut self) {
        let n0 = (self.p_caudatum * 0.85).round().clamp(0.0, 140.0) as usize;
        let n1 = if self.competition_mode {
            (self.p_aurelia * 0.85).round().clamp(0.0, 140.0) as usize
        } else {
            (self.didinium * 1.4).round().clamp(0.0, 90.0) as usize
        };
        self.agents.clear();
        for i in 0..n0 {
            let in_refuge = !self.competition_mode && (i as f64 / n0.max(1) as f64) < self.sediment_refuge;
            let y_range = if in_refuge { (0.78, 0.96) } else { (0.06, 0.77) };
            self.agents.push(Agent2D {
                x: self.rng.range_f64(0.06, 0.94),
                y: self.rng.range_f64(y_range.0, y_range.1),
                vx: self.rng.range_f64(-0.015, 0.015),
                vy: self.rng.range_f64(-0.012, 0.012),
                kind: 0,
                aux: if in_refuge { 1.0 } else { 0.0 },
            });
        }
        for _ in 0..n1 {
            self.agents.push(Agent2D {
                x: self.rng.range_f64(0.06, 0.94),
                y: self.rng.range_f64(0.06, 0.78),
                vx: self.rng.range_f64(-0.018, 0.018),
                vy: self.rng.range_f64(-0.015, 0.015),
                kind: 1,
                aux: 1.0,
            });
        }
    }
}

// ============================================================================
// M3: T4 Bacteriophage & Escherichia coli Lytic Plaque Wave
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhagePlaqueSim {
    pub t_hours: f64,
    pub dt: f64,
    pub susceptible_bact: f64,
    pub infected_bact: f64,
    pub resistant_bact: f64,
    pub free_phage: f64,
    pub nutrient_r: f64,
    pub mu_bact: f64,
    pub adsorption_phi: f64,
    pub lysis_rate: f64,
    pub burst_size_beta: f64,
    pub phage_decay: f64,
    pub resistance_mut_rate: f64,
    pub grid_side: usize,
    /// `0.0..1.0` bacterial lawn density on a `grid_side x grid_side` agar plate (plaques appear as low-density rings)
    pub agar_lawn: Vec<f64>,
    pub agar_phage: Vec<f64>,
    pub current_phase: GrowthPhase,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for PhagePlaqueSim {
    fn default() -> Self {
        let side = 24;
        let mut lawn = vec![0.65; side * side];
        let mut phage_grid = vec![0.0; side * side];
        // Seed 3 initial plaque-forming units (PFUs) on the agar lawn
        for &(r, c) in &[(7usize, 8usize), (15, 16), (16, 7)] {
            let idx = r * side + c;
            phage_grid[idx] = 12.0;
            lawn[idx] = 0.25;
        }

        let mut s = Self {
            t_hours: 0.0,
            dt: 0.05,
            susceptible_bact: 85.0,
            infected_bact: 1.2,
            resistant_bact: 0.05,
            free_phage: 14.0,
            nutrient_r: 100.0,
            mu_bact: 0.82,
            adsorption_phi: 0.0085,
            lysis_rate: 1.45,
            burst_size_beta: 45.0,
            phage_decay: 0.12,
            resistance_mut_rate: 0.0008,
            grid_side: side,
            agar_lawn: lawn,
            agar_phage: phage_grid,
            current_phase: GrowthPhase::ExponentialLog,
            rng: SimRng::new(5005),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.record_state(0.0);
        s
    }
}

impl PhagePlaqueSim {
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Susceptible E. coli S",
            "Infected Eclipse Cells I",
            "Resistant Mutants R",
            "Free T4 Virions (x 0.1)",
            "",
        ]
    }

    fn record_state(&mut self, pcg: f64) {
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.t_hours,
                values: [
                    self.susceptible_bact,
                    self.infected_bact,
                    self.resistant_bact,
                    self.free_phage * 0.1,
                    0.0,
                ],
                per_capita_growth: pcg,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        let dt = self.dt;
        let prev_s = self.susceptible_bact.max(0.01);
        let monod = self.nutrient_r / (18.0 + self.nutrient_r);
        let mu_s = self.mu_bact * monod;
        let mu_r = 0.84 * mu_s; // Fitness cost of phage-receptor mutation

        let adsorp = self.adsorption_phi * self.susceptible_bact * self.free_phage;
        let lysis = self.lysis_rate * self.infected_bact;
        let mut_flux = self.resistance_mut_rate * self.susceptible_bact;

        let ds = mu_s * self.susceptible_bact - adsorp - mut_flux;
        let di = adsorp - lysis;
        let dr = mu_r * self.resistant_bact + mut_flux;
        let dv = self.burst_size_beta * lysis - adsorp - self.phage_decay * self.free_phage;
        let d_nut = -0.18 * (mu_s * self.susceptible_bact + mu_r * self.resistant_bact)
            + 0.08 * (100.0 - self.nutrient_r);

        self.susceptible_bact = (self.susceptible_bact + dt * ds).clamp(0.01, 500.0);
        self.infected_bact = (self.infected_bact + dt * di).clamp(0.0, 500.0);
        self.resistant_bact = (self.resistant_bact + dt * dr).clamp(0.0, 500.0);
        self.free_phage = (self.free_phage + dt * dv).clamp(0.0, 8000.0);
        self.nutrient_r = (self.nutrient_r + dt * d_nut).clamp(0.5, 150.0);
        self.t_hours += dt;

        // Update 2D agar-lawn radial plaque diffusion
        let side = self.grid_side;
        let mut next_lawn = self.agar_lawn.clone();
        let mut next_phage = self.agar_phage.clone();
        for r in 0..side {
            for c in 0..side {
                let idx = r * side + c;
                let b = self.agar_lawn[idx];
                let v = self.agar_phage[idx];
                let mut lap_v = 0.0;
                if r > 0 {
                    lap_v += self.agar_phage[(r - 1) * side + c] - v;
                }
                if r + 1 < side {
                    lap_v += self.agar_phage[(r + 1) * side + c] - v;
                }
                if c > 0 {
                    lap_v += self.agar_phage[r * side + (c - 1)] - v;
                }
                if c + 1 < side {
                    lap_v += self.agar_phage[r * side + (c + 1)] - v;
                }
                let inf = 0.055 * b * v;
                next_lawn[idx] = (b + dt * (0.32 * b * (1.0 - b) - inf)).clamp(0.02, 1.0);
                next_phage[idx] =
                    (v + dt * (14.0 * inf - 0.18 * v + 0.55 * lap_v)).clamp(0.0, 60.0);
            }
        }
        self.agar_lawn = next_lawn;
        self.agar_phage = next_phage;

        let pcg = ds / prev_s;
        self.current_phase = if lysis > 4.0 && ds < 0.0 {
            GrowthPhase::PredatorOvershoot
        } else if self.susceptible_bact < 4.0 && self.resistant_bact > 20.0 {
            GrowthPhase::DiauxicShift
        } else if self.susceptible_bact < 4.0 {
            GrowthPhase::PopulationCrashBottleneck
        } else {
            GrowthPhase::ExponentialLog
        };
        self.record_state(pcg);
    }

    pub fn inject_perturbation(&mut self) {
        // Inoculate a new PFU droplet onto a random spot on the agar lawn
        let idx = self.rng.gen_index(self.grid_side * self.grid_side);
        self.agar_phage[idx] += 25.0;
        self.free_phage = (self.free_phage + 80.0).min(8000.0);
    }
}
