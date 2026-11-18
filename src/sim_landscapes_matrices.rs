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
    Custom,
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
            Self::Custom => "Custom / CSV Imported Matrix",
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
    #[serde(skip)]
    pub last_matrix_error: Option<String>,
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
            last_matrix_error: None,
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
        self.last_matrix_error = None;

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
            MatrixPreset::Custom => {
                if self.base_matrix.is_empty() || self.stage_names.len() != self.base_matrix.len() {
                    self.stage_names = (1..=4).map(|i| format!("Stage {i}")).collect();
                    self.base_matrix = vec![
                        vec![0.0, 1.2, 1.8, 0.8],
                        vec![0.45, 0.0, 0.0, 0.0],
                        vec![0.0, 0.60, 0.0, 0.0],
                        vec![0.0, 0.0, 0.50, 0.0],
                    ];
                    self.abundances = vec![100.0; 4];
                } else {
                    let k = self.base_matrix.len();
                    if self.abundances.len() != k {
                        self.abundances = vec![100.0; k];
                    }
                }
            }
        }

        self.recompute_eigen();
        self.realized_lambda = self.eigen.dominant_lambda;
        self.update_phase_and_record();
    }

    /// Resize the projection matrix to `k x k` (`2..=12`), preserving overlapping entries and
    /// keeping `stage_names`, `base_matrix`, `abundances`, and `eigen` dimensions synchronized.
    pub fn resize_stages(&mut self, new_k: usize) {
        let new_k = new_k.clamp(2, 12);
        let old_k = self.base_matrix.len();

        let mut next_names = Vec::with_capacity(new_k);
        let mut next_mat = vec![vec![0.0; new_k]; new_k];
        let mut next_abund = Vec::with_capacity(new_k);

        for i in 0..new_k {
            if let Some(existing) = self.stage_names.get(i) {
                next_names.push(existing.clone());
            } else {
                next_names.push(format!("Stage {}", i + 1));
            }
            if let Some(&n_i) = self.abundances.get(i) {
                next_abund.push(n_i.max(0.0));
            } else {
                next_abund.push(50.0);
            }
        }

        for r in 0..old_k.min(new_k) {
            for c in 0..old_k.min(new_k) {
                next_mat[r][c] = self.base_matrix[r].get(c).copied().unwrap_or(0.0).max(0.0);
            }
        }
        if new_k > old_k {
            for i in old_k.max(1)..new_k {
                next_mat[i][i - 1] = 0.5;
            }
        }

        self.preset = MatrixPreset::Custom;
        self.stage_names = next_names;
        self.base_matrix = next_mat;
        self.abundances = next_abund;
        self.step_year = 0;
        self.history.clear();
        self.transitions.clear();
        self.last_matrix_error = None;
        self.recompute_eigen();
        self.realized_lambda = self.eigen.dominant_lambda;
        self.update_phase_and_record();
    }

    /// Parse a single CSV line with RFC-4180 quote handling (`"a, b"`, escaped `""`).
    fn parse_csv_row(line: &str, row_num_1based: usize) -> Result<Vec<String>, String> {
        let mut cells = Vec::new();
        let mut cur = String::new();
        let mut in_quotes = false;
        let mut just_closed_quote = false;
        let mut chars = line.chars().peekable();

        while let Some(c) = chars.next() {
            if in_quotes {
                if c == '"' {
                    if chars.peek() == Some(&'"') {
                        chars.next();
                        cur.push('"');
                    } else {
                        in_quotes = false;
                        just_closed_quote = true;
                    }
                } else {
                    cur.push(c);
                }
            } else if c == ',' {
                cells.push(cur.trim().to_string());
                cur.clear();
                just_closed_quote = false;
            } else if c == '"' {
                if cur.trim().is_empty() && !just_closed_quote {
                    cur.clear();
                    in_quotes = true;
                } else {
                    return Err(format!(
                        "CSV row {row_num_1based} contains an unexpected quote character."
                    ));
                }
            } else if just_closed_quote {
                if !c.is_whitespace() {
                    return Err(format!(
                        "CSV row {row_num_1based} contains trailing characters after a closing quote."
                    ));
                }
            } else {
                cur.push(c);
            }
        }

        if in_quotes {
            return Err(format!(
                "CSV row {row_num_1based} contains an unclosed quoted field."
            ));
        }
        cells.push(cur.trim().to_string());
        Ok(cells)
    }

    /// Parse a square Leslie/Lefkovitch projection matrix CSV supporting:
    /// 1) Plain `k x k` numeric CSV
    /// 2) Header row (`k` stage labels) + `k` numeric rows
    /// 3) Header row (`k+1` cells) + `k` rows with leading row label (`(k+1) x (k+1)`)
    pub fn parse_matrix_csv(csv_text: &str) -> Result<(Vec<String>, Vec<Vec<f64>>), String> {
        let stripped = csv_text.trim_start_matches('\u{feff}');
        let mut rows: Vec<Vec<String>> = Vec::new();
        for (line_idx, line) in stripped.lines().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            rows.push(Self::parse_csv_row(trimmed, line_idx + 1)?);
        }

        if rows.len() < 2 {
            return Err("CSV must contain at least a 2x2 matrix.".to_string());
        }
        let cols = rows[0].len();
        if cols < 2 {
            return Err("CSV must have at least 2 columns.".to_string());
        }
        for (idx, row) in rows.iter().enumerate() {
            if row.len() != cols {
                return Err(format!(
                    "CSV row {} has {} columns, expected {} columns.",
                    idx + 1,
                    row.len(),
                    cols
                ));
            }
        }

        let parse_val = |r_idx: usize, c_idx: usize, s: &str| -> Result<f64, String> {
            let v = s.parse::<f64>().map_err(|_| {
                format!(
                    "Invalid numeric entry at row {}, col {}: '{}'",
                    r_idx + 1,
                    c_idx + 1,
                    s
                )
            })?;
            if !v.is_finite() || v < 0.0 {
                return Err(format!(
                    "Matrix entry at row {}, col {} must be finite and >= 0 (got {}).",
                    r_idx + 1,
                    c_idx + 1,
                    s
                ));
            }
            Ok(v)
        };

        let (names, matrix) = if rows.len() == cols {
            // Classify plain k x k numeric matrices by float-token syntax (including negatives,
            // NaN, and inf) so malformed numeric rows never fall through into the labeled branch.
            let all_f64_tokens = rows
                .iter()
                .all(|row| row.iter().all(|s| s.parse::<f64>().is_ok()));

            if all_f64_tokens {
                let k = cols;
                if !(2..=16).contains(&k) {
                    return Err(format!("Matrix stage count k={k} must be in 2..=16."));
                }
                let mut mat = vec![vec![0.0; k]; k];
                for r in 0..k {
                    for c in 0..k {
                        mat[r][c] = parse_val(r, c, &rows[r][c])?;
                    }
                }
                let names = (1..=k).map(|i| format!("Stage {i}")).collect();
                (names, mat)
            } else {
                let k = cols - 1;
                if !(2..=16).contains(&k) {
                    return Err(format!("Labeled matrix stage count k={k} must be in 2..=16."));
                }
                for (c, cell) in rows[0].iter().enumerate() {
                    if cell.trim().parse::<f64>().is_ok() {
                        return Err(format!(
                            "Invalid numeric stage header '{}' at row 1, col {}.",
                            cell,
                            c + 1
                        ));
                    }
                }
                for (r, row) in rows.iter().enumerate().skip(1) {
                    if row[0].trim().parse::<f64>().is_ok() {
                        return Err(format!(
                            "Invalid numeric row label '{}' at row {}, col 1.",
                            row[0],
                            r + 1
                        ));
                    }
                }
                let mut mat = vec![vec![0.0; k]; k];
                let mut names = Vec::with_capacity(k);
                for c in 0..k {
                    let hdr = rows[0][c + 1].trim();
                    let row_lbl = rows[c + 1][0].trim();
                    let label = if !hdr.is_empty() {
                        hdr.to_string()
                    } else if !row_lbl.is_empty() {
                        row_lbl.to_string()
                    } else {
                        format!("Stage {}", c + 1)
                    };
                    names.push(label);
                }
                for r in 0..k {
                    for c in 0..k {
                        mat[r][c] = parse_val(r + 1, c + 1, &rows[r + 1][c + 1])?;
                    }
                }
                (names, mat)
            }
        } else if rows.len() == cols + 1 {
            let k = cols;
            if !(2..=16).contains(&k) {
                return Err(format!("Matrix stage count k={k} must be in 2..=16."));
            }
            for (c, cell) in rows[0].iter().enumerate() {
                if cell.trim().parse::<f64>().is_ok() {
                    return Err(format!(
                        "Invalid numeric stage header '{}' at row 1, col {}.",
                        cell,
                        c + 1
                    ));
                }
            }
            let mut mat = vec![vec![0.0; k]; k];
            for r in 0..k {
                for c in 0..k {
                    mat[r][c] = parse_val(r + 1, c, &rows[r + 1][c])?;
                }
            }
            let names = rows[0]
                .iter()
                .enumerate()
                .map(|(i, s)| {
                    if s.trim().is_empty() {
                        format!("Stage {}", i + 1)
                    } else {
                        s.trim().to_string()
                    }
                })
                .collect();
            (names, mat)
        } else {
            return Err(format!(
                "CSV shape {}x{} is not a square k x k, (k+1) x k, or (k+1) x (k+1) matrix.",
                rows.len(),
                cols
            ));
        };

        let mat_sum: f64 = matrix.iter().flatten().sum();
        if mat_sum <= 0.0 {
            return Err("Projection matrix entries cannot all be zero.".to_string());
        }
        Ok((names, matrix))
    }

    pub fn load_from_csv(&mut self, csv_text: &str) -> Result<(), String> {
        let (names, mat) = match Self::parse_matrix_csv(csv_text) {
            Ok(parsed) => parsed,
            Err(err) => {
                self.last_matrix_error = Some(err.clone());
                return Err(err);
            }
        };
        let k = mat.len();
        self.preset = MatrixPreset::Custom;
        self.stage_names = names;
        self.base_matrix = mat;
        self.abundances = vec![100.0; k];
        self.fecundity_mult = 1.0;
        self.juvenile_survival_mult = 1.0;
        self.adult_survival_mult = 1.0;
        self.step_year = 0;
        self.history.clear();
        self.transitions.clear();
        self.last_matrix_error = None;
        self.recompute_eigen();
        self.realized_lambda = self.eigen.dominant_lambda;
        self.update_phase_and_record();
        Ok(())
    }

    pub fn to_matrix_csv(&self) -> String {
        let k = self.base_matrix.len();
        let mut out = String::from("Stage");
        for name in &self.stage_names {
            out.push(',');
            out.push_str(&name.replace(',', ";"));
        }
        out.push('\n');
        for r in 0..k {
            let row_lbl = self
                .stage_names
                .get(r)
                .map(|s| s.replace(',', ";"))
                .unwrap_or_else(|| format!("Stage {}", r + 1));
            out.push_str(&row_lbl);
            for c in 0..k {
                let val = self.base_matrix[r].get(c).copied().unwrap_or(0.0);
                out.push_str(&format!(",{val:.6}"));
            }
            out.push('\n');
        }
        out
    }

    /// Construct the active `k x k` projection matrix `A` incorporating user conservation/vital-rate multipliers.
    pub fn active_matrix(&self) -> Vec<Vec<f64>> {
        let k = self.base_matrix.len();
        let mut mat = self.base_matrix.clone();
        let adult_start = (k * 3) / 5;
        let max_surv = if self.preset == MatrixPreset::Custom
            && (self.juvenile_survival_mult - 1.0).abs() < 1e-12
            && (self.adult_survival_mult - 1.0).abs() < 1e-12
        {
            1.0
        } else {
            0.995
        };

        for r in 0..k {
            for c in 0..k {
                if r == 0
                    && (c > 0
                        || matches!(
                            self.preset,
                            MatrixPreset::SongbirdLeslie4 | MatrixPreset::Custom
                        ))
                    && mat[r][c] > 0.0
                {
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
                    mat[r][c] = (mat[r][c] * mult).clamp(0.0, max_surv);
                }
            }
        }
        // Ensure column survival sums below row 0 never exceed max_surv
        for c in 0..k {
            let surv_sum: f64 = (1..k).map(|r| mat[r][c]).sum();
            let diag_zero = if r_Zero_is_stasis(self.preset) && c == 0 {
                mat[0][0]
            } else {
                0.0
            };
            let tot = surv_sum + diag_zero;
            if tot > max_surv {
                let scale = max_surv / tot;
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
    /// Fisher reproductive value `v`, sensitivity `S_ij`, and elasticity `E_ij` via
    /// Tarjan SCC decomposition and Osborne-balanced Perron-shifted squaring + M-matrix
    /// polishing so periodic/imprimitive Leslie matrices (`k = 2..=16`) and reducible /
    /// defective Jordan-block Lefkovitch matrices converge to machine precision.
    pub fn recompute_eigen(&mut self) {
        let mat = self.active_matrix();
        let k = mat.len();
        if k == 0 {
            return;
        }

        let max_entry = mat.iter().flatten().copied().fold(0.0_f64, f64::max);
        if max_entry <= 0.0 {
            self.eigen = MatrixEigenSummary {
                dominant_lambda: 0.0,
                intrinsic_r: 1e-9_f64.ln(),
                net_reproductive_r0: 0.0,
                damping_ratio: 1.0,
                stable_distribution_w: vec![1.0 / (k as f64); k],
                reproductive_value_v: vec![1.0; k],
                sensitivity: vec![vec![1.0 / (k as f64); k]; k],
                elasticity: vec![vec![0.0; k]; k],
                stage_elasticity_totals: vec![0.0; k],
            };
            return;
        }

        // Decompose into Strongly Connected Components (edges r -> c where mat[r][c] > 0)
        // in reverse topological order (if C_a -> C_b then b < a).
        let sccs = tarjan_scc(&mat);
        let p = sccs.len();
        let mut scc_of = vec![0usize; k];
        for (s_idx, comp) in sccs.iter().enumerate() {
            for &u in comp {
                scc_of[u] = s_idx;
            }
        }

        let mut rho_scc = vec![0.0_f64; p];
        let mut w_scc: Vec<Vec<f64>> = Vec::with_capacity(p);
        let mut v_scc: Vec<Vec<f64>> = Vec::with_capacity(p);

        for (s_idx, comp) in sccs.iter().enumerate() {
            if comp.len() == 1 {
                let u = comp[0];
                rho_scc[s_idx] = mat[u][u];
                w_scc.push(vec![1.0]);
                v_scc.push(vec![1.0]);
            } else {
                let m = comp.len();
                let mut sub = vec![vec![0.0; m]; m];
                for (r_loc, &r_glob) in comp.iter().enumerate() {
                    for (c_loc, &c_glob) in comp.iter().enumerate() {
                        sub[r_loc][c_loc] = mat[r_glob][c_glob];
                    }
                }
                let (rho_c, w_c, v_c) = solve_irreducible_perron(&sub);
                rho_scc[s_idx] = rho_c;
                w_scc.push(w_c);
                v_scc.push(v_c);
            }
        }

        let lambda = rho_scc.iter().copied().fold(0.0_f64, f64::max);
        let dom_tol = 1e-12 * lambda.max(1e-300);
        let dom_count = rho_scc
            .iter()
            .filter(|&&r| r >= lambda - 1e-9 * lambda.max(1e-12))
            .count();

        let (lambda, w, v) = if lambda <= 1e-300 {
            // Nilpotent nonnegative matrix (lambda_1 = 0): any nonnegative combination of
            // zero columns satisfies A*w = 0, and zero rows satisfies A^T*v = 0.
            let mut zero_cols: Vec<usize> = (0..k)
                .filter(|&c| (0..k).all(|r| mat[r][c] <= 0.0))
                .collect();
            if zero_cols.is_empty() {
                zero_cols = (0..k)
                    .filter(|&c| (0..k).all(|r| mat[r][c] <= 1e-15))
                    .collect();
            }
            let mut w = vec![0.0; k];
            if zero_cols.is_empty() {
                w.fill(1.0 / (k as f64));
            } else {
                let val = 1.0 / (zero_cols.len() as f64);
                for c in zero_cols {
                    w[c] = val;
                }
            }

            let mut zero_rows: Vec<usize> = (0..k)
                .filter(|&r| (0..k).all(|c| mat[r][c] <= 0.0))
                .collect();
            if zero_rows.is_empty() {
                zero_rows = (0..k)
                    .filter(|&r| (0..k).all(|c| mat[r][c] <= 1e-15))
                    .collect();
            }
            let mut v = vec![0.0; k];
            if zero_rows.is_empty() {
                v.fill(1.0);
            } else {
                let val = 1.0 / (zero_rows.len() as f64);
                for r in zero_rows {
                    v[r] = val;
                }
                if v[0] > 1e-9 {
                    let v0 = v[0];
                    for x in &mut v {
                        *x /= v0;
                    }
                }
            }
            (0.0, w, v)
        } else if p == 1 {
            let w = w_scc.into_iter().next().unwrap_or_else(|| vec![1.0 / (k as f64); k]);
            let mut v = v_scc.into_iter().next().unwrap_or_else(|| vec![1.0 / (k as f64); k]);
            if v[0] > 1e-9 {
                let v0 = v[0];
                for x in &mut v {
                    *x /= v0;
                }
            }
            (lambda, w, v)
        } else {
            // Reducible nonnegative matrix: Perron-Frobenius decomposition on SCC condensation DAG
            let is_dom: Vec<bool> = rho_scc.iter().map(|&r| r >= lambda - dom_tol).collect();
            let mut reaches = vec![vec![false; p]; p];
            for r in 0..k {
                for c in 0..k {
                    if mat[r][c] > 0.0 && scc_of[r] != scc_of[c] {
                        reaches[scc_of[r]][scc_of[c]] = true;
                    }
                }
            }
            for mid in 0..p {
                for i in 0..p {
                    if reaches[i][mid] {
                        for j in 0..p {
                            if reaches[mid][j] {
                                reaches[i][j] = true;
                            }
                        }
                    }
                }
            }

            let valid_w: Vec<bool> = (0..p)
                .map(|s| is_dom[s] && !(0..p).any(|d| d != s && is_dom[d] && reaches[d][s]))
                .collect();
            let valid_v: Vec<bool> = (0..p)
                .map(|s| is_dom[s] && !(0..p).any(|d| d != s && is_dom[d] && reaches[s][d]))
                .collect();

            // Prefer an isolated dominant SCC valid for both w and v so <v, w> > 0
            let best_scc = |mask: &[bool]| -> usize {
                (0..p)
                    .filter(|&s| mask[s])
                    .max_by(|&a, &b| {
                        rho_scc[a]
                            .partial_cmp(&rho_scc[b])
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .unwrap_or(0)
            };
            let both_valid: Vec<bool> = (0..p).map(|s| valid_w[s] && valid_v[s]).collect();
            let (s_w, s_v) = if both_valid.iter().any(|&x| x) {
                let s_iso = best_scc(&both_valid);
                (s_iso, s_iso)
            } else {
                (best_scc(&valid_w), best_scc(&valid_v))
            };

            // Solve (lambda*I - A) w = 0 in Tarjan order (s = 0..p; C_s -> C_t implies t < s)
            let mut w = vec![0.0_f64; k];
            for s in 0..p {
                if s == s_w {
                    for (idx, &u) in sccs[s].iter().enumerate() {
                        w[u] = w_scc[s][idx];
                    }
                } else if is_dom[s] {
                    // Other dominant SCCs remain 0
                } else {
                    let comp = &sccs[s];
                    let m = comp.len();
                    let mut rhs = vec![0.0_f64; m];
                    for (r_loc, &u) in comp.iter().enumerate() {
                        let mut acc = 0.0_f64;
                        for j in 0..k {
                            if scc_of[j] != s {
                                acc += mat[u][j] * w[j];
                            }
                        }
                        rhs[r_loc] = acc;
                    }
                    if rhs.iter().any(|&x| x > 0.0) {
                        if m == 1 {
                            let u = comp[0];
                            let denom = lambda - mat[u][u];
                            if denom > 0.0 {
                                w[u] = (rhs[0] / denom).max(0.0);
                            }
                        } else {
                            let mut a_loc = vec![vec![0.0_f64; m]; m];
                            for (r_loc, &r_glob) in comp.iter().enumerate() {
                                for (c_loc, &c_glob) in comp.iter().enumerate() {
                                    a_loc[r_loc][c_loc] = if r_loc == c_loc {
                                        lambda - mat[r_glob][c_glob]
                                    } else {
                                        -mat[r_glob][c_glob]
                                    };
                                }
                            }
                            if let Some(sol) = solve_linear_system(a_loc, rhs) {
                                for (idx, &u) in comp.iter().enumerate() {
                                    w[u] = sol[idx].max(0.0);
                                }
                            }
                        }
                    }
                }
            }
            let w_sum: f64 = w.iter().sum();
            if w_sum > 1e-300 {
                for x in &mut w {
                    *x /= w_sum;
                }
            } else {
                w.fill(1.0 / (k as f64));
            }

            // Solve (lambda*I - A^T) v = 0 in reverse Tarjan order (s = p-1..=0)
            let mut v = vec![0.0_f64; k];
            for s in (0..p).rev() {
                if s == s_v {
                    for (idx, &u) in sccs[s].iter().enumerate() {
                        v[u] = v_scc[s][idx];
                    }
                } else if is_dom[s] {
                    // Other dominant SCCs remain 0
                } else {
                    let comp = &sccs[s];
                    let m = comp.len();
                    let mut rhs = vec![0.0_f64; m];
                    for (r_loc, &u) in comp.iter().enumerate() {
                        let mut acc = 0.0_f64;
                        for i in 0..k {
                            if scc_of[i] != s {
                                acc += mat[i][u] * v[i];
                            }
                        }
                        rhs[r_loc] = acc;
                    }
                    if rhs.iter().any(|&x| x > 0.0) {
                        if m == 1 {
                            let u = comp[0];
                            let denom = lambda - mat[u][u];
                            if denom > 0.0 {
                                v[u] = (rhs[0] / denom).max(0.0);
                            }
                        } else {
                            let mut a_loc = vec![vec![0.0_f64; m]; m];
                            for (r_loc, &r_glob) in comp.iter().enumerate() {
                                for (c_loc, &c_glob) in comp.iter().enumerate() {
                                    a_loc[r_loc][c_loc] = if r_loc == c_loc {
                                        lambda - mat[c_glob][r_glob]
                                    } else {
                                        -mat[c_glob][r_glob]
                                    };
                                }
                            }
                            if let Some(sol) = solve_linear_system(a_loc, rhs) {
                                for (idx, &u) in comp.iter().enumerate() {
                                    v[u] = sol[idx].max(0.0);
                                }
                            }
                        }
                    }
                }
            }
            let v_sum: f64 = v.iter().sum();
            if v_sum > 1e-300 {
                for x in &mut v {
                    *x /= v_sum;
                }
            } else {
                v.fill(1.0 / (k as f64));
            }
            if v[0] > 1e-9 {
                let v0 = v[0];
                for x in &mut v {
                    *x /= v0;
                }
            }
            (lambda, w, v)
        };

        // Estimate subdominant magnitude |lambda_2| via deflated orthogonal iteration for damping ratio
        let vw_dot_raw: f64 = v.iter().zip(w.iter()).map(|(vi, wi)| vi * wi).sum();
        let vw_dot: f64 = if vw_dot_raw > 1e-300 {
            vw_dot_raw
        } else {
            1e-12
        };
        let mut lam2_mag = 0.0_f64;
        if lambda > 1e-12 && (dom_count > 1 || vw_dot_raw <= 1e-200) {
            lam2_mag = lambda;
        } else if lambda > 1e-12 {
            let mut u = (0..k)
                .map(|i| if i % 2 == 0 { 1.0 } else { -0.7 })
                .collect::<Vec<_>>();
            let mut log_ratio_sum = 0.0_f64;
            let mut log_ratio_cnt = 0usize;
            for iter in 0..100 {
                // Deflate dominant component: u <- u - (<v, u> / <v, w>) * w
                let vu: f64 = v.iter().zip(u.iter()).map(|(vi, ui)| vi * ui).sum();
                for i in 0..k {
                    u[i] -= (vu / vw_dot_raw) * w[i];
                }
                let norm_u = u.iter().map(|x| x * x).sum::<f64>().sqrt();
                if norm_u <= 1e-12 {
                    lam2_mag = 0.0;
                    log_ratio_cnt = 0;
                    break;
                }
                let mut nu = vec![0.0; k];
                for r in 0..k {
                    for c in 0..k {
                        nu[r] += mat[r][c] * u[c];
                    }
                }
                let norm_nu = nu.iter().map(|x| x * x).sum::<f64>().sqrt();
                if norm_nu <= 1e-12 {
                    lam2_mag = 0.0;
                    log_ratio_cnt = 0;
                    break;
                }
                let ratio = norm_nu / norm_u;
                lam2_mag = ratio;
                if iter >= 40 {
                    log_ratio_sum += ratio.ln();
                    log_ratio_cnt += 1;
                }
                for i in 0..k {
                    u[i] = nu[i] / norm_nu;
                }
            }
            if log_ratio_cnt > 0 {
                lam2_mag = (log_ratio_sum / (log_ratio_cnt as f64)).exp();
            }
        }
        let damping_ratio = if lambda <= 1e-12 {
            1.0
        } else if lam2_mag > 1e-6 {
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
                if mat[r][c] > 0.0 && lambda > 1e-15 {
                    let e_ij = (mat[r][c] / lambda) * s_ij;
                    elasticity[r][c] = e_ij;
                    stage_elasticity_totals[c] += e_ij;
                }
            }
        }
        // Normalize elasticities strictly to sum = 1.0
        let e_sum: f64 = elasticity.iter().flatten().sum();
        if e_sum > 1e-12 {
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
            if matches!(
                self.preset,
                MatrixPreset::SongbirdLeslie4 | MatrixPreset::Custom
            ) {
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

/// Decompose directed graph `r -> c` (`mat[r][c] > 0.0`) into Strongly Connected Components
/// in reverse topological order (if `C_a -> C_b` then `b < a`).
fn tarjan_scc(mat: &[Vec<f64>]) -> Vec<Vec<usize>> {
    struct State<'a> {
        mat: &'a [Vec<f64>],
        index: usize,
        indices: Vec<usize>,
        lowlink: Vec<usize>,
        on_stack: Vec<bool>,
        stack: Vec<usize>,
        sccs: Vec<Vec<usize>>,
    }

    fn strongconnect(v: usize, st: &mut State<'_>) {
        let idx = st.index;
        st.indices[v] = idx;
        st.lowlink[v] = idx;
        st.index += 1;
        st.stack.push(v);
        st.on_stack[v] = true;

        let k = st.mat.len();
        for w in 0..k {
            if st.mat[v][w] > 0.0 {
                if st.indices[w] == usize::MAX {
                    strongconnect(w, st);
                    if st.lowlink[w] < st.lowlink[v] {
                        st.lowlink[v] = st.lowlink[w];
                    }
                } else if st.on_stack[w] && st.indices[w] < st.lowlink[v] {
                    st.lowlink[v] = st.indices[w];
                }
            }
        }

        if st.lowlink[v] == st.indices[v] {
            let mut comp = Vec::new();
            while let Some(w) = st.stack.pop() {
                st.on_stack[w] = false;
                comp.push(w);
                if w == v {
                    break;
                }
            }
            comp.sort_unstable();
            st.sccs.push(comp);
        }
    }

    let k = mat.len();
    let mut st = State {
        mat,
        index: 0,
        indices: vec![usize::MAX; k],
        lowlink: vec![0; k],
        on_stack: vec![false; k],
        stack: Vec::with_capacity(k),
        sccs: Vec::new(),
    };
    for v in 0..k {
        if st.indices[v] == usize::MAX {
            strongconnect(v, &mut st);
        }
    }
    st.sccs
}

/// Dense Gaussian elimination with scaled partial pivoting (`n <= 16`).
fn solve_linear_system(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let mut row_scale = vec![0.0_f64; n];
    for r in 0..n {
        let s = a[r].iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        if !s.is_finite() || s <= 0.0 {
            return None;
        }
        row_scale[r] = s;
    }

    for col in 0..n {
        let mut piv = col;
        let mut best_ratio = -1.0_f64;
        for r in col..n {
            let ratio = a[r][col].abs() / row_scale[r];
            if ratio.is_finite() && ratio > best_ratio {
                best_ratio = ratio;
                piv = r;
            }
        }
        if !a[piv][col].is_finite() || a[piv][col].abs() <= 1e-300 * row_scale[piv] {
            return None;
        }
        if piv != col {
            a.swap(col, piv);
            b.swap(col, piv);
            row_scale.swap(col, piv);
        }
        let diag = a[col][col];
        for r in (col + 1)..n {
            let factor = a[r][col] / diag;
            a[r][col] = 0.0;
            if factor != 0.0 {
                for c in (col + 1)..n {
                    a[r][c] -= factor * a[col][c];
                }
                b[r] -= factor * b[col];
            }
        }
    }

    let mut x = vec![0.0_f64; n];
    for r in (0..n).rev() {
        let mut sum = 0.0_f64;
        for c in (r + 1)..n {
            sum += a[r][c] * x[c];
        }
        let val = (b[r] - sum) / a[r][r];
        if !val.is_finite() {
            return None;
        }
        x[r] = val;
    }
    Some(x)
}

fn rayleigh_quotient(mat: &[Vec<f64>], w: &[f64], v: &[f64]) -> Option<f64> {
    let m = mat.len();
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for r in 0..m {
        let mut aw_r = 0.0_f64;
        for c in 0..m {
            aw_r += mat[r][c] * w[c];
        }
        num += v[r] * aw_r;
        den += v[r] * w[r];
    }
    if den > 1e-300 && num.is_finite() && den.is_finite() {
        Some(num / den)
    } else {
        None
    }
}

fn refine_eigenvector_m_matrix(
    mat: &[Vec<f64>],
    rho: f64,
    approx: &[f64],
    transpose: bool,
) -> Option<Vec<f64>> {
    let m = mat.len();
    if m < 2 || rho <= 0.0 || !rho.is_finite() {
        return None;
    }
    let mut i_max = 0usize;
    let mut best_val = -1.0_f64;
    for (i, &val) in approx.iter().enumerate() {
        if val > best_val {
            best_val = val;
            i_max = i;
        }
    }

    let n = m - 1;
    let mut sub = vec![vec![0.0_f64; n]; n];
    let mut rhs = vec![0.0_f64; n];
    let mut r_sub = 0usize;
    for r in 0..m {
        if r == i_max {
            continue;
        }
        rhs[r_sub] = if transpose {
            mat[i_max][r]
        } else {
            mat[r][i_max]
        };
        let mut c_sub = 0usize;
        for c in 0..m {
            if c == i_max {
                continue;
            }
            let entry = if transpose { mat[c][r] } else { mat[r][c] };
            sub[r_sub][c_sub] = if r == c { rho - entry } else { -entry };
            c_sub += 1;
        }
        r_sub += 1;
    }

    let sol = solve_linear_system(sub, rhs)?;
    let max_sol = sol.iter().copied().fold(1.0_f64, f64::max);
    if sol.iter().any(|&x| !x.is_finite() || x < -1e-12 * max_sol) {
        return None;
    }

    let mut out = vec![0.0_f64; m];
    let mut idx = 0usize;
    for i in 0..m {
        if i == i_max {
            out[i] = 1.0;
        } else {
            out[i] = sol[idx].max(0.0);
            idx += 1;
        }
    }
    let sum: f64 = out.iter().sum();
    if sum <= 1e-300 || !sum.is_finite() {
        return None;
    }
    for x in &mut out {
        *x /= sum;
    }
    Some(out)
}

fn eigen_residual_inf(mat: &[Vec<f64>], lam: f64, w: &[f64], v_l1: &[f64]) -> f64 {
    let m = mat.len();
    let v_scale = if v_l1[0] > 1e-9 { 1.0 / v_l1[0] } else { 1.0 };
    let mut err = 0.0_f64;
    for r in 0..m {
        let mut aw_r = 0.0_f64;
        let mut atv_r = 0.0_f64;
        for c in 0..m {
            aw_r += mat[r][c] * w[c];
            atv_r += mat[c][r] * (v_l1[c] * v_scale);
        }
        err = err.max((aw_r - lam * w[r]).abs());
        err = err.max((atv_r - lam * (v_l1[r] * v_scale)).abs());
    }
    err
}

/// Compute exact Perron eigenvalue `rho > 0` and L1-normalized right (`w`) and left (`v`)
/// eigenvectors for an irreducible nonnegative `m x m` matrix (`m >= 2`).
fn solve_irreducible_perron(b: &[Vec<f64>]) -> (f64, Vec<f64>, Vec<f64>) {
    let m = b.len();
    let scale = b.iter().flatten().copied().fold(0.0_f64, f64::max);
    if scale <= 0.0 || !scale.is_finite() {
        return (0.0, vec![1.0 / (m as f64); m], vec![1.0 / (m as f64); m]);
    }

    let mut b_scaled = vec![vec![0.0_f64; m]; m];
    for r in 0..m {
        for c in 0..m {
            b_scaled[r][c] = b[r][c] / scale;
        }
    }

    // Osborne diagonal similarity balancing: b_bal = D^-1 * b_scaled * D
    let mut d = vec![1.0_f64; m];
    let mut b_bal = b_scaled.clone();
    for _sweep in 0..32 {
        let mut max_step = 1.0_f64;
        for i in 0..m {
            let mut row_sum = 0.0_f64;
            let mut col_sum = 0.0_f64;
            for j in 0..m {
                if j != i {
                    row_sum += b_bal[i][j];
                    col_sum += b_bal[j][i];
                }
            }
            if row_sum > 1e-300 && col_sum > 1e-300 {
                let target_d = (d[i] * (row_sum / col_sum).sqrt().clamp(1e-8, 1e8))
                    .clamp(1e-150, 1e150);
                let f = target_d / d[i];
                if (f - 1.0).abs() > 1e-13 {
                    d[i] = target_d;
                    for j in 0..m {
                        if j != i {
                            b_bal[i][j] /= f;
                            b_bal[j][i] *= f;
                        }
                    }
                    max_step = max_step.max(f.max(1.0 / f));
                }
            }
        }
        if (max_step - 1.0) < 1e-12 {
            break;
        }
    }
    for r in 0..m {
        for c in 0..m {
            b_bal[r][c] = b_scaled[r][c] * (d[c] / d[r]);
        }
    }

    // Estimate Perron shift sigma on balanced matrix
    let mut sigma_est_mat = b_bal.clone();
    let mut log2_scale = 0.0_f64;
    let mut est_ok = true;
    for _ in 0..4 {
        let mut sq = vec![vec![0.0_f64; m]; m];
        let mut sq_max = 0.0_f64;
        for r in 0..m {
            for c in 0..m {
                let mut acc = 0.0_f64;
                for j in 0..m {
                    acc += sigma_est_mat[r][j] * sigma_est_mat[j][c];
                }
                sq[r][c] = acc;
                if acc > sq_max {
                    sq_max = acc;
                }
            }
        }
        if sq_max <= 1e-300 || !sq_max.is_finite() {
            est_ok = false;
            break;
        }
        for row in &mut sq {
            for x in row {
                *x /= sq_max;
            }
        }
        log2_scale = 2.0 * log2_scale + sq_max.log2();
        sigma_est_mat = sq;
    }
    let row_sum_max = b_bal
        .iter()
        .map(|row| row.iter().sum::<f64>())
        .fold(0.0_f64, f64::max);
    let sigma = if est_ok {
        (log2_scale / 16.0).exp2().max(1e-12)
    } else {
        (0.5 * row_sum_max).max(1e-12)
    };

    let build_shifted_norm = |shift: f64| -> Vec<Vec<f64>> {
        let mut mat = b_bal.clone();
        let mut max_e = 0.0_f64;
        for r in 0..m {
            mat[r][r] += shift;
            for c in 0..m {
                if mat[r][c] > max_e {
                    max_e = mat[r][c];
                }
            }
        }
        if max_e > 0.0 {
            for row in &mut mat {
                for x in row {
                    *x /= max_e;
                }
            }
        }
        mat
    };

    // Normalized nonnegative matrix squaring of (b_bal + sigma*I) with relative convergence
    let mut m_pow = build_shifted_norm(sigma);
    let mut w_bal = vec![1.0 / (m as f64); m];
    let mut v_bal = vec![1.0 / (m as f64); m];
    let mut steps_since_reset = 0usize;

    for step in 0..56 {
        let mut sq = vec![vec![0.0_f64; m]; m];
        let mut sq_max = 0.0_f64;
        for r in 0..m {
            for c in 0..m {
                let mut acc = 0.0_f64;
                for j in 0..m {
                    acc += m_pow[r][j] * m_pow[j][c];
                }
                sq[r][c] = acc;
                if acc > sq_max {
                    sq_max = acc;
                }
            }
        }
        if sq_max <= 1e-300 || !sq_max.is_finite() {
            break;
        }
        for row in &mut sq {
            for x in row {
                *x /= sq_max;
            }
        }
        m_pow = sq;
        steps_since_reset += 1;

        // Extract w_curr from largest-sum column of m_pow
        let mut best_c = 0usize;
        let mut best_csum = -1.0_f64;
        for c in 0..m {
            let csum: f64 = (0..m).map(|r| m_pow[r][c]).sum();
            if csum > best_csum {
                best_csum = csum;
                best_c = c;
            }
        }
        let mut w_curr = vec![0.0_f64; m];
        if best_csum > 1e-300 {
            for r in 0..m {
                w_curr[r] = m_pow[r][best_c] / best_csum;
            }
        } else {
            w_curr.fill(1.0 / (m as f64));
        }

        // Extract v_curr from largest-sum row of m_pow
        let mut best_r = 0usize;
        let mut best_rsum = -1.0_f64;
        for r in 0..m {
            let rsum: f64 = m_pow[r].iter().sum();
            if rsum > best_rsum {
                best_rsum = rsum;
                best_r = r;
            }
        }
        let mut v_curr = vec![0.0_f64; m];
        if best_rsum > 1e-300 {
            for c in 0..m {
                v_curr[c] = m_pow[best_r][c] / best_rsum;
            }
        } else {
            v_curr.fill(1.0 / (m as f64));
        }

        let mut rel_w = 0.0_f64;
        let mut rel_v = 0.0_f64;
        for i in 0..m {
            let dw = (w_curr[i] - w_bal[i]).abs() / w_curr[i].abs().max(1e-300);
            let dv = (v_curr[i] - v_bal[i]).abs() / v_curr[i].abs().max(1e-300);
            if dw > rel_w {
                rel_w = dw;
            }
            if dv > rel_v {
                rel_v = dv;
            }
        }
        w_bal = w_curr;
        v_bal = v_curr;

        if steps_since_reset >= 2 && rel_w <= 1e-14 && rel_v <= 1e-14 {
            break;
        }
        if step == 6 && (rel_w > 1e-6 || rel_v > 1e-6) {
            if let Some(rho_mid) = rayleigh_quotient(&b_bal, &w_bal, &v_bal) {
                if rho_mid > 1e-12 {
                    m_pow = build_shifted_norm(rho_mid);
                    steps_since_reset = 0;
                }
            }
        }
    }

    // Unbalance w and v back to b_scaled coordinates
    let mut w0: Vec<f64> = (0..m).map(|i| w_bal[i] * d[i]).collect();
    let w0_sum: f64 = w0.iter().sum();
    if w0_sum > 1e-300 {
        for x in &mut w0 {
            *x /= w0_sum;
        }
    } else {
        w0.fill(1.0 / (m as f64));
    }

    let mut v0: Vec<f64> = (0..m).map(|i| v_bal[i] / d[i]).collect();
    let v0_sum: f64 = v0.iter().sum();
    if v0_sum > 1e-300 {
        for x in &mut v0 {
            *x /= v0_sum;
        }
    } else {
        v0.fill(1.0 / (m as f64));
    }

    let mut best_rho = rayleigh_quotient(&b_bal, &w_bal, &v_bal)
        .or_else(|| rayleigh_quotient(&b_scaled, &w0, &v0))
        .unwrap_or_else(|| {
            (0..m)
                .map(|r| (0..m).map(|c| b_scaled[r][c] * w0[c]).sum::<f64>())
                .sum()
        });
    let mut best_w = w0;
    let mut best_v = v0;
    let mut best_err = eigen_residual_inf(&b_scaled, best_rho, &best_w, &best_v);

    // M-matrix direct refinement (up to 2 passes) to reach 0-2 ULP residual on extreme scales
    for _ in 0..2 {
        let w_cand = refine_eigenvector_m_matrix(&b_scaled, best_rho, &best_w, false);
        let v_cand = refine_eigenvector_m_matrix(&b_scaled, best_rho, &best_v, true);
        if let (Some(wc), Some(vc)) = (w_cand, v_cand) {
            let rho_c = rayleigh_quotient(&b_scaled, &wc, &vc).unwrap_or(best_rho);
            let err_c = eigen_residual_inf(&b_scaled, rho_c, &wc, &vc);
            if err_c < best_err {
                best_rho = rho_c;
                best_w = wc;
                best_v = vc;
                best_err = err_c;
            } else {
                break;
            }
        } else {
            break;
        }
    }

    (best_rho * scale, best_w, best_v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn leslie_cell_edit_and_stage_resize_preserve_invariants() {
        let mut sim = LeslieLefkovitchSim::default();
        let lam_before = sim.eigen.dominant_lambda;
        sim.preset = MatrixPreset::Custom;
        sim.base_matrix[0][6] *= 2.0;
        sim.recompute_eigen();
        assert!(sim.eigen.dominant_lambda > lam_before);
        let e_sum: f64 = sim.eigen.elasticity.iter().flatten().sum();
        assert!((e_sum - 1.0).abs() < 1e-6);

        for new_k in [4usize, 9usize, 2usize, 6usize] {
            sim.resize_stages(new_k);
            assert_eq!(sim.stage_names.len(), new_k);
            assert_eq!(sim.base_matrix.len(), new_k);
            assert!(sim.base_matrix.iter().all(|row| row.len() == new_k));
            assert_eq!(sim.abundances.len(), new_k);
            assert_eq!(sim.eigen.stable_distribution_w.len(), new_k);
            assert_eq!(sim.eigen.elasticity.len(), new_k);
            sim.step();
            assert!(sim.eigen.dominant_lambda.is_finite());
        }
    }

    #[test]
    fn leslie_csv_import_handles_valid_formats_and_rejects_invalid() {
        let mut sim = LeslieLefkovitchSim::default();

        // 1. Plain k x k numeric
        let plain = "0.0, 2.5, 4.0\n0.5, 0.0, 0.0\n0.0, 0.6, 0.0\n";
        sim.load_from_csv(plain).expect("plain 3x3 CSV");
        assert_eq!(sim.preset, MatrixPreset::Custom);
        assert_eq!(sim.stage_names.len(), 3);
        assert!(sim.eigen.dominant_lambda > 1.0);

        // 2. Header row + k numeric rows
        let hdr_csv = "Seed,Juvenile,Adult\n0.0,0.0,8.0\n0.3,0.5,0.0\n0.0,0.4,0.8\n";
        sim.load_from_csv(hdr_csv).expect("header + 3x3 CSV");
        assert_eq!(sim.stage_names, vec!["Seed", "Juvenile", "Adult"]);

        // 3. BOM + comments + (k+1) x (k+1) labeled roundtrip
        let exported = format!("\u{feff}# Test matrix\n{}", sim.to_matrix_csv());
        sim.load_from_csv(&exported).expect("labeled roundtrip CSV");
        assert_eq!(sim.stage_names, vec!["Seed", "Juvenile", "Adult"]);
        assert!(sim.last_matrix_error.is_none());

        // Reject invalid inputs while preserving existing valid state
        let saved_k = sim.stage_names.len();
        let saved_lambda = sim.eigen.dominant_lambda;
        for bad in [
            "",
            "1.0,2.0\n3.0\n",
            "0.0,-0.5\n0.4,0.0\n",
            "0.0,abc\n0.4,0.0\n",
            "0.0,0.0\n0.0,0.0\n",
        ] {
            assert!(sim.load_from_csv(bad).is_err(), "expected error for: {bad}");
            assert!(sim.last_matrix_error.is_some());
            assert_eq!(sim.stage_names.len(), saved_k);
            assert!((sim.eigen.dominant_lambda - saved_lambda).abs() < 1e-12);
        }

        // Regression: 3x3 and 4x4 negative, NaN, inf in row 0 / col 0 must return Err and preserve state
        let saved_names = sim.stage_names.clone();
        let saved_matrix = sim.base_matrix.clone();
        for bad_k34 in [
            "0.0,-1.0,2.0\n0.4,0.0,0.0\n0.0,0.5,0.0\n",
            "0.0,NaN,1.0\n0.4,0.0,0.0\n0.0,0.5,0.0\n",
            "0.0,inf,1.0\n0.4,0.0,0.0\n0.0,0.5,0.0\n",
            "-0.2,1.0,2.0\n0.4,0.0,0.0\n0.0,0.5,0.0\n",
            "0.0,1.0,2.0\n-0.4,0.0,0.0\n0.0,0.5,0.0\n",
            "0.0,1.0,2.0,3.0\nNaN,0.0,0.0,0.0\n0.0,0.5,0.0,0.0\n0.0,0.0,0.6,0.0\n",
            "0.0,1.0,-inf,3.0\n0.4,0.0,0.0,0.0\n0.0,0.5,0.0,0.0\n0.0,0.0,0.6,0.0\n",
            "0.0,abc,1.0\n0.4,0.0,0.0\n0.0,0.5,0.0\n",
            "0.0,1.0,2.0\n0.4,0.0,0.0\n0.0,0.5,0.0\n0.0,0.0,0.2\n",
            "\"Unclosed, quote\",1.0\n0.4,0.0\n",
        ] {
            assert!(
                sim.load_from_csv(bad_k34).is_err(),
                "expected error for malformed CSV:\n{bad_k34}"
            );
            assert!(sim.last_matrix_error.is_some());
            assert_eq!(sim.stage_names, saved_names);
            assert_eq!(sim.base_matrix, saved_matrix);
            assert!((sim.eigen.dominant_lambda - saved_lambda).abs() < 1e-12);
        }

        // Regression: quote-aware CSV splitting preserves embedded commas and escaped "" quotes
        let quoted_csv = "\"Stage 1, Eggs\",\"Stage \"\"2\"\", Adults\"\r\n0.0,2.5\r\n0.4,0.7\r\n";
        sim.load_from_csv(quoted_csv)
            .expect("quoted CSV with commas and escaped quotes");
        assert_eq!(
            sim.stage_names,
            vec!["Stage 1, Eggs", "Stage \"2\", Adults"]
        );
        assert_eq!(sim.base_matrix, vec![vec![0.0, 2.5], vec![0.4, 0.7]]);
        assert!(sim.last_matrix_error.is_none());
    }

    fn assert_eigen_residuals(sim: &LeslieLefkovitchSim, label: &str) {
        let mat = sim.active_matrix();
        let k = mat.len();
        let lam = sim.eigen.dominant_lambda;
        let w = &sim.eigen.stable_distribution_w;
        let v = &sim.eigen.reproductive_value_v;
        assert_eq!(w.len(), k, "{label}: w length");
        assert_eq!(v.len(), k, "{label}: v length");
        assert!(lam.is_finite(), "{label}: dominant_lambda finite");
        assert!(sim.eigen.intrinsic_r.is_finite(), "{label}: intrinsic_r finite");
        assert!(
            sim.eigen.net_reproductive_r0.is_finite(),
            "{label}: net_reproductive_r0 finite"
        );
        assert!(
            sim.eigen.damping_ratio.is_finite(),
            "{label}: damping_ratio finite"
        );

        let mut rw_err = 0.0_f64;
        let mut lv_err = 0.0_f64;
        for r in 0..k {
            let aw_r: f64 = (0..k).map(|c| mat[r][c] * w[c]).sum();
            rw_err = rw_err.max((aw_r - lam * w[r]).abs());
            let atv_r: f64 = (0..k).map(|c| mat[c][r] * v[c]).sum();
            lv_err = lv_err.max((atv_r - lam * v[r]).abs());
        }
        assert!(
            rw_err < 1e-6,
            "{label}: ||A*w - lambda*w||_inf = {rw_err} >= 1e-6 (lam={lam})"
        );
        assert!(
            lv_err < 1e-6,
            "{label}: ||A^T*v - lambda*v||_inf = {lv_err} >= 1e-6 (lam={lam})"
        );
    }

    #[test]
    fn periodic_2x2_and_3x3_leslie_matrices_converge_to_exact_eigen_and_elasticity() {
        let mut sim = LeslieLefkovitchSim::default();

        // Canonical 2x2 semelparous Leslie matrix: lambda_1 = 1.0, E_01 = E_10 = 0.5
        sim.load_from_csv("0.0,2.5\n0.4,0.0\n")
            .expect("2x2 periodic Leslie CSV");
        assert_eigen_residuals(&sim, "periodic_2x2");
        assert!(
            (sim.eigen.dominant_lambda - 1.0).abs() < 1e-6,
            "expected lambda=1.0, got {}",
            sim.eigen.dominant_lambda
        );
        assert!(
            (sim.eigen.elasticity[0][1] - 0.5).abs() < 1e-6,
            "expected E[0][1]=0.5, got {}",
            sim.eigen.elasticity[0][1]
        );
        assert!(
            (sim.eigen.elasticity[1][0] - 0.5).abs() < 1e-6,
            "expected E[1][0]=0.5, got {}",
            sim.eigen.elasticity[1][0]
        );
        let e_sum_2: f64 = sim.eigen.elasticity.iter().flatten().sum();
        assert!((e_sum_2 - 1.0).abs() < 1e-6);
        assert!((sim.eigen.net_reproductive_r0 - 1.0).abs() < 1e-6);
        assert!((sim.eigen.damping_ratio - 1.0).abs() < 1e-6);

        // Canonical 3x3 semelparous Leslie matrix: lambda_1 = 1.0, E_02 = E_10 = E_21 = 1/3
        sim.load_from_csv("0.0,0.0,4.0\n0.5,0.0,0.0\n0.0,0.5,0.0\n")
            .expect("3x3 periodic Leslie CSV");
        assert_eigen_residuals(&sim, "periodic_3x3");
        assert!(
            (sim.eigen.dominant_lambda - 1.0).abs() < 1e-6,
            "expected lambda=1.0, got {}",
            sim.eigen.dominant_lambda
        );
        let third = 1.0 / 3.0;
        assert!((sim.eigen.elasticity[0][2] - third).abs() < 1e-6);
        assert!((sim.eigen.elasticity[1][0] - third).abs() < 1e-6);
        assert!((sim.eigen.elasticity[2][1] - third).abs() < 1e-6);
        let e_sum_3: f64 = sim.eigen.elasticity.iter().flatten().sum();
        assert!((e_sum_3 - 1.0).abs() < 1e-6);
        assert!((sim.eigen.net_reproductive_r0 - 1.0).abs() < 1e-6);
        assert!((sim.eigen.damping_ratio - 1.0).abs() < 1e-6);
    }

    #[test]
    fn scaled_periodic_and_resized_reducible_matrices_preserve_eigen_invariants() {
        let mut sim = LeslieLefkovitchSim::default();

        // Scaled 2x2 periodic matrix with biological survival <= 1.0 (true lambda = 10.0)
        sim.load_from_csv("0.0,200.0\n0.5,0.0\n")
            .expect("scaled 2x2 periodic CSV");
        assert_eigen_residuals(&sim, "scaled_2x2_lambda10");
        assert!((sim.eigen.dominant_lambda - 10.0).abs() < 1e-6);
        assert!((sim.eigen.elasticity[0][1] - 0.5).abs() < 1e-6);
        assert!((sim.eigen.elasticity[1][0] - 0.5).abs() < 1e-6);

        // Scaled 2x2 where active_matrix() clamps r>0 survival 5.0 -> 1.0 (active lambda = sqrt(20))
        sim.load_from_csv("0.0,20.0\n5.0,0.0\n")
            .expect("scaled 2x2 survival-clamped CSV");
        assert_eigen_residuals(&sim, "scaled_2x2_clamped");
        assert!((sim.eigen.dominant_lambda - 20.0_f64.sqrt()).abs() < 1e-6);
        let e_sum_c: f64 = sim.eigen.elasticity.iter().flatten().sum();
        assert!((e_sum_c - 1.0).abs() < 1e-6);

        // Scaled 3x3 periodic matrix (true lambda = 10.0)
        sim.load_from_csv("0.0,0.0,8000.0\n0.5,0.0,0.0\n0.0,0.25,0.0\n")
            .expect("scaled 3x3 periodic CSV");
        assert_eigen_residuals(&sim, "scaled_3x3_lambda10");
        assert!((sim.eigen.dominant_lambda - 10.0).abs() < 1e-6);
        let third = 1.0 / 3.0;
        assert!((sim.eigen.elasticity[0][2] - third).abs() < 1e-6);
        assert!((sim.eigen.elasticity[1][0] - third).abs() < 1e-6);
        assert!((sim.eigen.elasticity[2][1] - third).abs() < 1e-6);

        // Skewed unbalanced 2x2 matrix (max_row_sum = 1000 while true lambda = 1.0)
        sim.load_from_csv("0.0,1000.0\n0.001,0.0\n")
            .expect("skewed 2x2 CSV");
        assert_eigen_residuals(&sim, "skewed_2x2");
        assert!((sim.eigen.dominant_lambda - 1.0).abs() < 1e-6);

        // resize_stages(2..=8) after CSV import (k=2 is periodic; k=4..=8 have post-reproductive sinks)
        sim.load_from_csv("Stage1,Stage2,Stage3\n0.0,1.5,3.0\n0.5,0.0,0.0\n0.0,0.6,0.2\n")
            .expect("3x3 resize base CSV");
        for k in 2..=8 {
            sim.resize_stages(k);
            assert_eigen_residuals(&sim, &format!("resize_k_{k}"));
            let e_sum: f64 = sim.eigen.elasticity.iter().flatten().sum();
            assert!(
                (e_sum - 1.0).abs() < 1e-6,
                "resize_k_{k}: sum(elasticity) = {e_sum}"
            );
            if k == 2 {
                assert!(
                    (sim.eigen.dominant_lambda - 0.75_f64.sqrt()).abs() < 1e-6,
                    "resize_k_2: expected sqrt(0.75), got {}",
                    sim.eigen.dominant_lambda
                );
            }
        }
    }

    #[test]
    fn all_zero_and_nilpotent_matrices_remain_finite_and_exact() {
        let mut sim = LeslieLefkovitchSim::default();
        sim.preset = MatrixPreset::Custom;
        sim.stage_names = vec!["S1".into(), "S2".into(), "S3".into()];
        sim.base_matrix = vec![vec![0.0; 3]; 3];
        sim.abundances = vec![100.0; 3];
        sim.recompute_eigen();
        assert_eigen_residuals(&sim, "all_zero_3x3");
        assert_eq!(sim.eigen.dominant_lambda, 0.0);
        assert_eq!(sim.eigen.net_reproductive_r0, 0.0);

        sim.load_from_csv("0.0,2.0\n0.0,0.0\n")
            .expect("nilpotent 2x2 CSV");
        assert_eigen_residuals(&sim, "nilpotent_2x2");
        assert_eq!(sim.eigen.dominant_lambda, 0.0);
        assert_eq!(sim.eigen.net_reproductive_r0, 0.0);
    }

    fn build_periodic_leslie_csv(k: usize, s: f64, target_lam: f64) -> String {
        let f = target_lam.powi(k as i32) / s.powi((k - 1) as i32);
        let mut lines = Vec::with_capacity(k);
        for r in 0..k {
            let mut row = Vec::with_capacity(k);
            for c in 0..k {
                let val = if r == 0 && c == k - 1 {
                    f
                } else if r > 0 && c == r - 1 {
                    s
                } else {
                    0.0
                };
                row.push(format!("{val:.17e}"));
            }
            lines.push(row.join(","));
        }
        lines.join("\n")
    }

    #[test]
    fn high_order_periodic_leslie_matrices_k10_to_k16_converge_to_exact_left_right_eigenvectors() {
        let mut sim = LeslieLefkovitchSim::default();

        // Exact cases flagged in Challenger 1 Iter 4 plus full k = 2..=16 sweeps
        let cases: &[(usize, f64, f64, &str)] = &[
            (10, 0.25, 1.0, "periodic_s0.25_lam1.0_k10"),
            (12, 0.25, 1.0, "periodic_s0.25_lam1.0_k12"),
            (10, 0.12, 1.0, "periodic_s0.12_lam1.0_k10"),
            (12, 0.20, 1.0, "periodic_s0.2_lam1.0_k12"),
            (12, 0.50, 1.5, "periodic_s0.5_lam1.5_k12"),
            (14, 0.50, 1.5, "periodic_s0.5_lam1.5_k14"),
            (16, 0.50, 1.5, "periodic_s0.5_lam1.5_k16"),
            (16, 0.25, 1.0, "periodic_s0.25_lam1.0_k16"),
            (16, 0.20, 1.0, "periodic_s0.2_lam1.0_k16"),
        ];

        for &(k, s, target_lam, label) in cases {
            let csv = build_periodic_leslie_csv(k, s, target_lam);
            sim.load_from_csv(&csv).expect(label);
            assert_eigen_residuals(&sim, label);
            assert!(
                (sim.eigen.dominant_lambda - target_lam).abs() < 1e-9,
                "{label}: expected lambda={target_lam}, got {}",
                sim.eigen.dominant_lambda
            );
            let e_sum: f64 = sim.eigen.elasticity.iter().flatten().sum();
            assert!(
                (e_sum - 1.0).abs() < 1e-6,
                "{label}: sum(elasticity) = {e_sum}"
            );
            let expected_e = 1.0 / (k as f64);
            assert!(
                (sim.eigen.elasticity[0][k - 1] - expected_e).abs() < 1e-6,
                "{label}: fecundity elasticity {} != {expected_e}",
                sim.eigen.elasticity[0][k - 1]
            );
            assert!(
                (sim.eigen.damping_ratio - 1.0).abs() < 1e-4,
                "{label}: expected damping_ratio=1.0, got {}",
                sim.eigen.damping_ratio
            );
        }

        // Full k = 2..=16 sweep at s = 0.25, lam = 1.0
        for k in 2..=16 {
            let label = format!("periodic_sweep_s0.25_k{k}");
            let csv = build_periodic_leslie_csv(k, 0.25, 1.0);
            sim.load_from_csv(&csv).expect(&label);
            assert_eigen_residuals(&sim, &label);
            let e_sum: f64 = sim.eigen.elasticity.iter().flatten().sum();
            assert!((e_sum - 1.0).abs() < 1e-6, "{label}: e_sum={e_sum}");
        }
    }

    #[test]
    fn reducible_and_defective_jordan_lefkovitch_matrices_converge_genuinely() {
        let mut sim = LeslieLefkovitchSim::default();

        // 1. 2x2 lower-triangular Jordan block: lambda=1.0, w=[0,1], v=[1,0]
        sim.load_from_csv("1.0,0.0\n0.5,1.0\n")
            .expect("jordan_2x2 CSV");
        assert_eigen_residuals(&sim, "jordan_2x2");
        assert!((sim.eigen.dominant_lambda - 1.0).abs() < 1e-9);
        assert!((sim.eigen.damping_ratio - 1.0).abs() < 1e-9);

        // 2. 3x3 lower-triangular Jordan block: lambda=1.0, w=[0,0,1], v=[1,0,0]
        sim.load_from_csv("1.0,0.0,0.0\n0.5,1.0,0.0\n0.0,0.5,1.0\n")
            .expect("jordan_3x3 CSV");
        assert_eigen_residuals(&sim, "jordan_3x3");
        assert!((sim.eigen.dominant_lambda - 1.0).abs() < 1e-9);
        assert!((sim.eigen.damping_ratio - 1.0).abs() < 1e-9);

        // 3. 2x2 upper-triangular Jordan block: lambda=0.8, w=[1,0], v=[0,1]
        sim.load_from_csv("0.8,1.0\n0.0,0.8\n")
            .expect("jordan_upper_2x2 CSV");
        assert_eigen_residuals(&sim, "jordan_upper_2x2");
        assert!((sim.eigen.dominant_lambda - 0.8).abs() < 1e-9);
        assert!((sim.eigen.damping_ratio - 1.0).abs() < 1e-9);

        // 4. Bottom-dominant reducible 2x2: lambda=1.0, w=[0,1], v=[1,1.25], sum(E)=1.0, damping=2.0
        sim.load_from_csv("0.5,0.0\n0.4,1.0\n")
            .expect("reducible_bottom_dom_2x2 CSV");
        assert_eigen_residuals(&sim, "reducible_bottom_dom_2x2");
        assert!((sim.eigen.dominant_lambda - 1.0).abs() < 1e-9);
        let e_sum_b: f64 = sim.eigen.elasticity.iter().flatten().sum();
        assert!((e_sum_b - 1.0).abs() < 1e-6);
        assert!((sim.eigen.damping_ratio - 2.0).abs() < 1e-4);

        // 5. Coupled Jordan 2-cycles (4x4): two 2-cycles with rho=1.0 linked 2 -> 0
        sim.load_from_csv(
            "0.0,2.5,0.0,0.0\n0.4,0.0,0.0,0.0\n0.3,0.0,0.0,2.5\n0.0,0.0,0.4,0.0\n",
        )
        .expect("coupled_jordan_2cycle CSV");
        assert_eigen_residuals(&sim, "coupled_jordan_2cycle");
        assert!((sim.eigen.dominant_lambda - 1.0).abs() < 1e-9);

        // 6. Disconnected dominant 2-cycle + subdominant sink (3x3)
        sim.load_from_csv("0.0,2.5,0.0\n0.4,0.0,0.0\n0.0,0.3,0.5\n")
            .expect("reducible_2cycle_with_sink CSV");
        assert_eigen_residuals(&sim, "reducible_2cycle_with_sink");
        assert!((sim.eigen.dominant_lambda - 1.0).abs() < 1e-9);
        let e_sum_s: f64 = sim.eigen.elasticity.iter().flatten().sum();
        assert!((e_sum_s - 1.0).abs() < 1e-6);

        // 7. Near-reducible 2x2 with 1e-10 coupling
        sim.load_from_csv("1.0,1e-10\n1e-10,1.0\n")
            .expect("near_reducible_2x2 CSV");
        assert_eigen_residuals(&sim, "near_reducible_2x2");
        assert!((sim.eigen.dominant_lambda - (1.0 + 1e-10)).abs() < 1e-12);
        let e_sum_nr: f64 = sim.eigen.elasticity.iter().flatten().sum();
        assert!((e_sum_nr - 1.0).abs() < 1e-6);

        // 8. All presets across resize_stages(2..=16)
        for preset in [
            MatrixPreset::SongbirdLeslie4,
            MatrixPreset::LoggerheadTurtleLefkovitch7,
            MatrixPreset::WoodlandPerennialLefkovitch5,
            MatrixPreset::Custom,
        ] {
            for k in 2..=16 {
                let mut psim = LeslieLefkovitchSim::default();
                psim.load_preset(preset);
                psim.resize_stages(k);
                let lbl = format!("{preset:?}_k{k}");
                assert_eigen_residuals(&psim, &lbl);
                let e_sum: f64 = psim.eigen.elasticity.iter().flatten().sum();
                assert!((e_sum - 1.0).abs() < 1e-6, "{lbl}: e_sum={e_sum}");
            }
        }
    }
}
