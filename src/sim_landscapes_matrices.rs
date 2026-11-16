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
// T5: Sewall Wright (1932) 3D Adaptive Fitness Landscape & Shifting Balance
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LandscapeDeme {
    pub x: f64,
    pub y: f64,
    pub fitness: f64,
    pub trail: Vec<[f64; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WrightLandscapeSim {
    pub generation: usize,
    pub num_demes: usize,
    pub effective_ne: f64,
    pub selection_strength: f64,
    pub interdeme_migration_m: f64,
    /// Height bonus of Peak 1 (local suboptimal peak at `(0.24, 0.26)`)
    pub peak1_height: f64,
    /// Height bonus of Peak 2 (global adaptive peak at `(0.76, 0.74)`)
    pub peak2_height: f64,
    /// Epistatic ruggedness / maladaptive valley depth separating the adaptive peaks
    pub valley_ruggedness: f64,
    /// Interactive 3D camera angles (degrees) for rotating the 3D surface in the GUI
    pub camera_azimuth_deg: f64,
    pub camera_elevation_deg: f64,
    pub demes: Vec<LandscapeDeme>,
    pub current_phase: GrowthPhase,
    pub rng: SimRng,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for WrightLandscapeSim {
    fn default() -> Self {
        let mut s = Self {
            generation: 0,
            num_demes: 8,
            effective_ne: 38.0,
            selection_strength: 0.16,
            interdeme_migration_m: 0.018,
            peak1_height: 0.55,
            peak2_height: 1.05,
            valley_ruggedness: 0.32,
            camera_azimuth_deg: -38.0,
            camera_elevation_deg: 30.0,
            demes: Vec::new(),
            current_phase: GrowthPhase::DriftSelectionTrajectory,
            rng: SimRng::new(1932),
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.seed_demes_near_suboptimal_peak();
        s.record_state(0.0);
        s
    }
}

impl WrightLandscapeSim {
    pub fn reset(&mut self) {
        let az = self.camera_azimuth_deg;
        let el = self.camera_elevation_deg;
        *self = Self::default();
        self.camera_azimuth_deg = az;
        self.camera_elevation_deg = el;
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Mean Metapopulation Fitness W_bar",
            "Highest Deme Fitness max(W)",
            "Fraction of Demes on Global Peak (%)",
            "Mean Coordinate X (Trait/Allele 1)",
            "Mean Coordinate Y (Trait/Allele 2)",
        ]
    }

    /// Evaluate mean population fitness `W_bar(x, y)` on `[0, 1]^2`.
    pub fn fitness_at(&self, x: f64, y: f64) -> f64 {
        let xc = x.clamp(0.0, 1.0);
        let yc = y.clamp(0.0, 1.0);
        let d1_sq = (xc - 0.24).powi(2) + (yc - 0.26).powi(2);
        let d2_sq = (xc - 0.76).powi(2) + (yc - 0.74).powi(2);
        let d3_sq = (xc - 0.25).powi(2) + (yc - 0.78).powi(2);

        let p1 = self.peak1_height * (-d1_sq / 0.036).exp();
        let p2 = self.peak2_height * (-d2_sq / 0.042).exp();
        let p3 = 0.38 * self.valley_ruggedness * (-d3_sq / 0.028).exp();
        // Epistatic saddle trench between (0.24,0.26) and (0.76,0.74) centered near (0.48, 0.48)
        let d_val = (xc - 0.48).powi(2) + (yc - 0.48).powi(2);
        let valley_dip = 0.35 * self.valley_ruggedness * (-d_val / 0.032).exp();
        let ripple = 0.06
            * self.valley_ruggedness
            * (std::f64::consts::TAU * 2.0 * xc).sin()
            * (std::f64::consts::TAU * 2.0 * yc).cos();

        (0.50 + p1 + p2 + p3 - valley_dip + ripple).max(0.08)
    }

    /// Numerical fitness gradient `(dW/dx, dW/dy)` via central finite differences.
    pub fn fitness_gradient(&self, x: f64, y: f64) -> (f64, f64) {
        let eps = 0.004;
        let dw_dx = (self.fitness_at((x + eps).min(1.0), y)
            - self.fitness_at((x - eps).max(0.0), y))
            / (2.0 * eps);
        let dw_dy = (self.fitness_at(x, (y + eps).min(1.0))
            - self.fitness_at(x, (y - eps).max(0.0)))
            / (2.0 * eps);
        (dw_dx, dw_dy)
    }

    fn seed_demes_near_suboptimal_peak(&mut self) {
        self.demes.clear();
        let k = self.num_demes.clamp(1, 16);
        for _ in 0..k {
            let x = self.rng.range_f64(0.14, 0.34);
            let y = self.rng.range_f64(0.15, 0.36);
            let w = self.fitness_at(x, y);
            self.demes.push(LandscapeDeme {
                x,
                y,
                fitness: w,
                trail: vec![[x, y, w]],
            });
        }
    }

    fn summary_metrics(&self) -> (f64, f64, f64, f64, f64) {
        let k = self.demes.len().max(1) as f64;
        let mean_w = self.demes.iter().map(|d| d.fitness).sum::<f64>() / k;
        let max_w = self
            .demes
            .iter()
            .map(|d| d.fitness)
            .fold(0.0f64, f64::max);
        let global_peak_frac = 100.0
            * (self
                .demes
                .iter()
                .filter(|d| (d.x - 0.76).hypot(d.y - 0.74) < 0.22)
                .count() as f64)
            / k;
        let mean_x = self.demes.iter().map(|d| d.x).sum::<f64>() / k;
        let mean_y = self.demes.iter().map(|d| d.y).sum::<f64>() / k;
        (mean_w, max_w, global_peak_frac, mean_x, mean_y)
    }

    fn record_state(&mut self, dw_dt: f64) {
        let (mean_w, max_w, frac_global, mx, my) = self.summary_metrics();
        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.generation as f64,
                values: [mean_w, max_w, frac_global, mx, my],
                per_capita_growth: dw_dt,
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        let target_k = self.num_demes.clamp(1, 16);
        while self.demes.len() < target_k {
            let x = self.rng.range_f64(0.16, 0.34);
            let y = self.rng.range_f64(0.16, 0.34);
            let w = self.fitness_at(x, y);
            self.demes.push(LandscapeDeme {
                x,
                y,
                fitness: w,
                trail: vec![[x, y, w]],
            });
        }
        self.demes.truncate(target_k);

        let (prev_mean_w, _, _, _, _) = self.summary_metrics();

        // Wright Shifting Balance Phase III: fitness-weighted migrant pool (higher-fitness peaks export more migrants!)
        let mut weight_sum = 0.0;
        let mut mig_x = 0.0;
        let mut mig_y = 0.0;
        for d in &self.demes {
            let w_sq = d.fitness.powi(3);
            weight_sum += w_sq;
            mig_x += w_sq * d.x;
            mig_y += w_sq * d.y;
        }
        if weight_sum > 1e-9 {
            mig_x /= weight_sum;
            mig_y /= weight_sum;
        }

        // Genetic drift standard deviation scales as 1 / sqrt(2 * N_e)
        let drift_sigma = 0.22 / (2.0 * self.effective_ne.max(4.0)).sqrt();

        // Collect updated coordinates first so `&self` gradient calls don't conflict with `&mut self.demes`
        let mut updates = Vec::with_capacity(self.demes.len());
        for d in &self.demes {
            let (gx, gy) = self.fitness_gradient(d.x, d.y);
            let w_cur = d.fitness.max(0.1);
            // Lande / Wright steepest-ascent step `(sel / W_bar) * grad(W)`
            let sel_dx = self.selection_strength * 0.035 * (gx / w_cur);
            let sel_dy = self.selection_strength * 0.035 * (gy / w_cur);
            let mig_dx = self.interdeme_migration_m * (mig_x - d.x);
            let mig_dy = self.interdeme_migration_m * (mig_y - d.y);
            updates.push((d.x + sel_dx + mig_dx, d.y + sel_dy + mig_dy));
        }

        for (i, (bx, by)) in updates.into_iter().enumerate() {
            let nx = (bx + drift_sigma * self.rng.standard_normal()).clamp(0.03, 0.97);
            let ny = (by + drift_sigma * self.rng.standard_normal()).clamp(0.03, 0.97);
            let nw = self.fitness_at(nx, ny);
            let d = &mut self.demes[i];
            d.x = nx;
            d.y = ny;
            d.fitness = nw;
            if d.trail.len() >= 90 {
                d.trail.drain(0..15);
            }
            d.trail.push([nx, ny, nw]);
        }

        self.generation += 1;
        let (mean_w, _, frac_global, _, _) = self.summary_metrics();

        // Classify Wright's 3-Phase Shifting Balance stage!
        self.current_phase = if frac_global >= 85.0 {
            GrowthPhase::StationaryCarryingCapacity
        } else if frac_global >= 12.0 {
            // Phase III: At least one deme crossed the valley and is pulling others via interdeme dispersion!
            GrowthPhase::ExponentialLog
        } else if self.effective_ne <= 45.0 {
            // Phase I: Exploratory genetic drift across the epistatic saddle
            GrowthPhase::DriftSelectionTrajectory
        } else {
            // Trapped on local suboptimal Peak 1 under strong selection / large N_e
            GrowthPhase::DeceleratingTransition
        };

        self.record_state(mean_w - prev_mean_w);
    }

    pub fn inject_perturbation(&mut self) {
        // Environmental shift / founder drift pulse: scatter demes across the epistatic saddle
        for i in 0..self.demes.len() {
            let nx = (self.demes[i].x + self.rng.range_f64(-0.22, 0.28)).clamp(0.05, 0.95);
            let ny = (self.demes[i].y + self.rng.range_f64(-0.22, 0.28)).clamp(0.05, 0.95);
            let nw = self.fitness_at(nx, ny);
            self.demes[i].x = nx;
            self.demes[i].y = ny;
            self.demes[i].fitness = nw;
            self.demes[i].trail.push([nx, ny, nw]);
        }
    }
}

// ============================================================================
// T6: Age- & Stage-Structured Leslie (1945) & Lefkovitch (1965) Matrix Model
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MatrixPreset {
    LoggerheadTurtleLefkovitch7,
    SongbirdLeslie4,
    WoodlandPerennialLefkovitch5,
}

impl MatrixPreset {
    pub fn label(&self) -> &'static str {
        match self {
            Self::LoggerheadTurtleLefkovitch7 => {
                "Loggerhead Sea Turtle (Caretta caretta, 7-Stage Lefkovitch — Crouse et al. 1987)"
            }
            Self::SongbirdLeslie4 => {
                "Passerine Songbird (4-Age Strict Leslie Matrix — Leslie 1945)"
            }
            Self::WoodlandPerennialLefkovitch5 => {
                "Perennial Woodland Herb (5-Stage Seed/Vegetative/Flowering Lefkovitch)"
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatrixEigenSummary {
    pub dominant_lambda: f64,
    pub intrinsic_r: f64,
    pub net_reproductive_r0: f64,
    pub damping_ratio: f64,
    pub stable_distribution_w: Vec<f64>,
    pub reproductive_value_v: Vec<f64>,
    pub sensitivity: Vec<Vec<f64>>,
    pub elasticity: Vec<Vec<f64>>,
    /// Per-stage summed elasticity `(fecundity_sum, growth_sum, stasis_sum)`
    pub stage_elasticity_totals: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LeslieLefkovitchSim {
    pub step_year: usize,
    pub preset: MatrixPreset,
    pub stage_names: Vec<String>,
    pub base_matrix: Vec<Vec<f64>>,
    pub fecundity_mult: f64,
    pub juvenile_survival_mult: f64,
    pub adult_survival_mult: f64,
    /// Optional carrying-capacity Beverton-Holt density regulation on recruitment (`false` = pure linear projection)
    pub density_regulated: bool,
    pub recruitment_k: f64,
    pub abundances: Vec<f64>,
    pub realized_lambda: f64,
    pub eigen: MatrixEigenSummary,
    pub current_phase: GrowthPhase,
    pub history: Vec<TimeRecord>,
    pub transitions: Vec<PhaseTransitionEvent>,
}

impl Default for LeslieLefkovitchSim {
    fn default() -> Self {
        let mut s = Self {
            step_year: 0,
            preset: MatrixPreset::LoggerheadTurtleLefkovitch7,
            stage_names: Vec::new(),
            base_matrix: Vec::new(),
            fecundity_mult: 1.0,
            juvenile_survival_mult: 1.0,
            adult_survival_mult: 1.0,
            density_regulated: false,
            recruitment_k: 5000.0,
            abundances: Vec::new(),
            realized_lambda: 0.945,
            eigen: MatrixEigenSummary {
                dominant_lambda: 1.0,
                intrinsic_r: 0.0,
                net_reproductive_r0: 1.0,
                damping_ratio: 1.1,
                stable_distribution_w: Vec::new(),
                reproductive_value_v: Vec::new(),
                sensitivity: Vec::new(),
                elasticity: Vec::new(),
                stage_elasticity_totals: Vec::new(),
            },
            current_phase: GrowthPhase::DeathSenescence,
            history: Vec::new(),
            transitions: Vec::new(),
        };
        s.load_preset(MatrixPreset::LoggerheadTurtleLefkovitch7);
        s
    }
}

impl LeslieLefkovitchSim {
    pub fn reset(&mut self) {
        let p = self.preset;
        self.fecundity_mult = 1.0;
        self.juvenile_survival_mult = 1.0;
        self.adult_survival_mult = 1.0;
        self.density_regulated = false;
        self.load_preset(p);
    }

    pub fn series_names(&self) -> [&'static str; 5] {
        [
            "Total Population N(t)",
            "Stage 1 (Eggs / Juveniles / Seeds)",
            "Mid-Stages (Juveniles / Vegetative)",
            "Reproductive Adults",
            "Asymptotic Eigenvalue lambda_1 (x 100)",
        ]
    }

    pub fn load_preset(&mut self, preset: MatrixPreset) {
        self.preset = preset;
        self.step_year = 0;
        self.history.clear();
        self.transitions.clear();

        match preset {
            MatrixPreset::LoggerheadTurtleLefkovitch7 => {
                // Classic Crouse, Crowder & Caswell (1987) Ecology 68:1412-1423 Table 4 Caretta caretta matrix
                self.stage_names = vec![
                    "1. Eggs/Hatchlings".to_string(),
                    "2. Small Juveniles".to_string(),
                    "3. Large Juveniles".to_string(),
                    "4. Subadults".to_string(),
                    "5. Novice Breeders".to_string(),
                    "6. 1st-Yr Remigrants".to_string(),
                    "7. Mature Breeders".to_string(),
                ];
                self.base_matrix = vec![
                    vec![0.0, 0.0, 0.0, 0.0, 127.0, 4.0, 80.0],
                    vec![0.6747, 0.7370, 0.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0486, 0.6610, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0147, 0.6907, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0518, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.8091, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.8091, 0.8089],
                ];
                self.abundances = vec![900.0, 420.0, 140.0, 35.0, 6.0, 5.0, 18.0];
            }
            MatrixPreset::SongbirdLeslie4 => {
                // Strict age-classified Leslie (1945) matrix (zeros on diagonal, survivals on subdiagonal)
                self.stage_names = vec![
                    "Age 0-1 (Fledglings)".to_string(),
                    "Age 1-2 (Yearlings)".to_string(),
                    "Age 2-3 (Prime Adults)".to_string(),
                    "Age 3-4 (Senescent)".to_string(),
                ];
                self.base_matrix = vec![
                    vec![0.45, 1.65, 2.10, 0.90],
                    vec![0.38, 0.0, 0.0, 0.0],
                    vec![0.0, 0.58, 0.0, 0.0],
                    vec![0.0, 0.0, 0.44, 0.0],
                ];
                self.abundances = vec![500.0, 60.0, 20.0, 5.0];
            }
            MatrixPreset::WoodlandPerennialLefkovitch5 => {
                self.stage_names = vec![
                    "1. Dormant Seed Bank".to_string(),
                    "2. Seedlings".to_string(),
                    "3. Small Vegetative".to_string(),
                    "4. Single-Flower Adult".to_string(),
                    "5. Multi-Flower Large".to_string(),
                ];
                self.base_matrix = vec![
                    vec![0.35, 0.0, 0.0, 4.2, 11.5],
                    vec![0.18, 0.0, 0.0, 0.35, 0.85],
                    vec![0.0, 0.42, 0.64, 0.08, 0.0],
                    vec![0.0, 0.0, 0.22, 0.72, 0.12],
                    vec![0.0, 0.0, 0.0, 0.14, 0.84],
                ];
                self.abundances = vec![400.0, 120.0, 95.0, 45.0, 22.0];
            }
        }

        self.recompute_eigen();
        self.realized_lambda = self.eigen.dominant_lambda;
        self.update_phase_and_record();
    }

    /// Construct the active `k x k` projection matrix `A` incorporating user conservation/vital-rate multipliers.
    pub fn active_matrix(&self) -> Vec<Vec<f64>> {
        let k = self.base_matrix.len();
        let mut mat = self.base_matrix.clone();
        let adult_start = (k * 3) / 5;

        for r in 0..k {
            for c in 0..k {
                if r == 0 && (c > 0 || self.preset == MatrixPreset::SongbirdLeslie4) && mat[r][c] > 0.0 {
                    // Top-row fecundity (keep seed-bank stasis [0][0] separate on Lefkovitch plant model)
                    if !(self.preset == MatrixPreset::WoodlandPerennialLefkovitch5 && c == 0) {
                        mat[r][c] = (mat[r][c] * self.fecundity_mult).max(0.0);
                    }
                } else if r > 0 && mat[r][c] > 0.0 {
                    let mult = if c < adult_start {
                        self.juvenile_survival_mult
                    } else {
                        self.adult_survival_mult
                    };
                    mat[r][c] = (mat[r][c] * mult).clamp(0.0, 0.995);
                }
            }
        }
        // Ensure column survival sums below row 0 never exceed 0.998
        for c in 0..k {
            let surv_sum: f64 = (1..k).map(|r| mat[r][c]).sum();
            let diag_zero = if r_Zero_is_stasis(self.preset) && c == 0 {
                mat[0][0]
            } else {
                0.0
            };
            let tot = surv_sum + diag_zero;
            if tot > 0.995 {
                let scale = 0.995 / tot;
                for r in 1..k {
                    mat[r][c] *= scale;
                }
                if diag_zero > 0.0 {
                    mat[0][0] *= scale;
                }
            }
        }
        mat
    }

    /// Recompute exact dominant eigenvalue `lambda_1`, stable stage distribution `w`,
    /// Fisher reproductive value `v`, sensitivity `S_ij`, and elasticity `E_ij` via power iteration.
    pub fn recompute_eigen(&mut self) {
        let mat = self.active_matrix();
        let k = mat.len();
        if k == 0 {
            return;
        }

        // Right eigenvector w (A * w = lambda * w)
        let mut w = vec![1.0 / (k as f64); k];
        let mut lambda = 1.0;
        for _ in 0..160 {
            let mut next_w = vec![0.0; k];
            for r in 0..k {
                for c in 0..k {
                    next_w[r] += mat[r][c] * w[c];
                }
            }
            let sum: f64 = next_w.iter().sum();
            if sum <= 1e-15 {
                lambda = 0.0;
                break;
            }
            lambda = sum;
            for x in &mut next_w {
                *x /= sum;
            }
            w = next_w;
        }

        // Left eigenvector v (A^T * v = lambda * v -> Fisher's Reproductive Value)
        let mut v = vec![1.0; k];
        for _ in 0..160 {
            let mut next_v = vec![0.0; k];
            for c in 0..k {
                for r in 0..k {
                    next_v[c] += mat[r][c] * v[r];
                }
            }
            let norm = next_v[0].max(next_v.iter().cloned().fold(0.0f64, f64::max)).max(1e-12);
            for x in &mut next_v {
                *x /= norm;
            }
            v = next_v;
        }
        if v[0] > 1e-9 {
            let v0 = v[0];
            for x in &mut v {
                *x /= v0;
            }
        }

        // Estimate subdominant magnitude |lambda_2| via deflated orthogonal iteration for damping ratio
        let vw_dot: f64 = v.iter().zip(w.iter()).map(|(vi, wi)| vi * wi).sum::<f64>().max(1e-12);
        let mut u = (0..k)
            .map(|i| if i % 2 == 0 { 1.0 } else { -0.7 })
            .collect::<Vec<_>>();
        let mut lam2_mag = 0.5 * lambda;
        for _ in 0..60 {
            // Deflate dominant component: u <- u - (<v, u> / <v, w>) * w
            let vu: f64 = v.iter().zip(u.iter()).map(|(vi, ui)| vi * ui).sum();
            for i in 0..k {
                u[i] -= (vu / vw_dot) * w[i];
            }
            let mut nu = vec![0.0; k];
            for r in 0..k {
                for c in 0..k {
                    nu[r] += mat[r][c] * u[c];
                }
            }
            let norm_u = u.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-12);
            let norm_nu = nu.iter().map(|x| x * x).sum::<f64>().sqrt();
            lam2_mag = norm_nu / norm_u;
            if norm_nu > 1e-12 {
                for i in 0..k {
                    u[i] = nu[i] / norm_nu;
                }
            }
        }
        let damping_ratio = if lam2_mag > 1e-6 {
            (lambda / lam2_mag).max(1.0)
        } else {
            9.99
        };

        // Sensitivity S_ij = (v_i * w_j) / <v, w> and Elasticity E_ij = (a_ij / lambda) * S_ij
        let mut sensitivity = vec![vec![0.0; k]; k];
        let mut elasticity = vec![vec![0.0; k]; k];
        let mut stage_elasticity_totals = vec![0.0; k];

        for r in 0..k {
            for c in 0..k {
                let s_ij = (v[r] * w[c]) / vw_dot;
                sensitivity[r][c] = s_ij;
                if mat[r][c] > 0.0 && lambda > 1e-12 {
                    let e_ij = (mat[r][c] / lambda) * s_ij;
                    elasticity[r][c] = e_ij;
                    stage_elasticity_totals[c] += e_ij;
                }
            }
        }
        // Normalize elasticities strictly to sum = 1.0
        let e_sum: f64 = elasticity.iter().flatten().sum();
        if e_sum > 1e-9 {
            for r in 0..k {
                for c in 0..k {
                    elasticity[r][c] /= e_sum;
                }
            }
            for c in 0..k {
                stage_elasticity_totals[c] /= e_sum;
            }
        }

        // Net reproductive rate R_0 approximation via cumulative cohort survivorship * fecundity
        let mut r0 = 0.0;
        let mut cohort = vec![0.0; k];
        cohort[0] = 1.0;
        for _ in 0..80 {
            for c in 1..k {
                r0 += mat[0][c] * cohort[c];
            }
            if self.preset == MatrixPreset::SongbirdLeslie4 {
                r0 += mat[0][0] * cohort[0];
            }
            let mut next_c = vec![0.0; k];
            for r in 1..k {
                for c in 0..k {
                    next_c[r] += mat[r][c] * cohort[c];
                }
            }
            cohort = next_c;
            if cohort.iter().sum::<f64>() < 1e-6 {
                break;
            }
        }

        self.eigen = MatrixEigenSummary {
            dominant_lambda: lambda,
            intrinsic_r: lambda.max(1e-9).ln(),
            net_reproductive_r0: r0,
            damping_ratio,
            stable_distribution_w: w,
            reproductive_value_v: v,
            sensitivity,
            elasticity,
            stage_elasticity_totals,
        };
    }

    fn update_phase_and_record(&mut self) {
        let tot: f64 = self.abundances.iter().sum::<f64>().max(1e-6);
        // L1 distance from stable stage distribution w
        let l1_dev: f64 = self
            .abundances
            .iter()
            .zip(self.eigen.stable_distribution_w.iter())
            .map(|(&n_i, &w_i)| ((n_i / tot) - w_i).abs())
            .sum();

        self.current_phase = if l1_dev > 0.14 {
            GrowthPhase::DeceleratingTransition // Transient stage-wave convergence
        } else if self.density_regulated && (self.realized_lambda - 1.0).abs() < 0.015 {
            GrowthPhase::StationaryCarryingCapacity
        } else if self.eigen.dominant_lambda > 1.005 {
            GrowthPhase::ExponentialLog
        } else if self.eigen.dominant_lambda < 0.995 {
            GrowthPhase::DeathSenescence
        } else {
            GrowthPhase::StationaryCarryingCapacity
        };

        let k = self.abundances.len();
        let s1 = self.abundances.first().copied().unwrap_or(0.0);
        let adult_start = (k * 3) / 5;
        let mid: f64 = self.abundances[1..adult_start.max(1)].iter().sum();
        let adults: f64 = self.abundances[adult_start..].iter().sum();

        push_record(
            &mut self.history,
            &mut self.transitions,
            TimeRecord {
                t: self.step_year as f64,
                values: [
                    tot,
                    s1,
                    mid,
                    adults,
                    self.eigen.dominant_lambda * 100.0,
                ],
                per_capita_growth: self.realized_lambda.max(1e-9).ln(),
                phase: self.current_phase,
            },
        );
    }

    pub fn step(&mut self) {
        self.recompute_eigen();
        let mat = self.active_matrix();
        let k = mat.len();
        let prev_tot: f64 = self.abundances.iter().sum::<f64>().max(1e-6);

        let mut next_n = vec![0.0; k];
        let bh_factor = if self.density_regulated {
            1.0 / (1.0 + prev_tot / self.recruitment_k.max(50.0))
        } else {
            1.0
        };

        for r in 0..k {
            for c in 0..k {
                let rate = if r == 0 && !(r_Zero_is_stasis(self.preset) && c == 0) {
                    mat[r][c] * bh_factor
                } else {
                    mat[r][c]
                };
                next_n[r] += rate * self.abundances[c];
            }
        }

        // Prevent numerical overflow or underflow over centuries of exploration
        let new_tot: f64 = next_n.iter().sum();
        self.realized_lambda = if prev_tot > 1e-9 {
            new_tot / prev_tot
        } else {
            self.eigen.dominant_lambda
        };
        if new_tot > 1.0e7 {
            let scale = 1.0e7 / new_tot;
            for x in &mut next_n {
                *x *= scale;
            }
        }
        self.abundances = next_n;
        self.step_year += 1;
        self.update_phase_and_record();
    }

    pub fn inject_perturbation(&mut self) {
        // Skew cohort 100% into Stage 1 (founder hatchling/seed pulse) to visualize transient oscillations converging to w!
        let tot: f64 = self.abundances.iter().sum::<f64>().clamp(200.0, 5000.0);
        for x in &mut self.abundances {
            *x = 0.0;
        }
        if let Some(first) = self.abundances.first_mut() {
            *first = tot;
        }
        self.update_phase_and_record();
    }
}

#[allow(non_snake_case)]
fn r_Zero_is_stasis(preset: MatrixPreset) -> bool {
    matches!(preset, MatrixPreset::WoodlandPerennialLefkovitch5)
}
