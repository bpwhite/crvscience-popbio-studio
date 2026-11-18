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

pub mod front_menu;
pub mod plots_panel;
pub mod popgen_view;
pub mod spatial_canvas;

use std::fs;
use std::time::Duration;

use eframe::egui::{self, Color32, Key, Margin, RichText};

use crate::model::{FunctionalResponse, ModuleGroup, ModuleId};
use crate::sim_landscapes_matrices::MatrixPreset;
use crate::workset::{SaveSlotInfo, WorksetDb, WorksetState};

use self::front_menu::{render_front_menu, FrontMenuAction};
use self::plots_panel::{render_phase_portrait_tab, render_time_series_tab};
use self::popgen_view::{
    render_empirical_controls, render_haplotype_network_view, render_popgen_summary_view,
    render_skyline_and_expansion_view, render_sliding_window_view,
};
use self::spatial_canvas::render_spatial_canvas;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AppScreen {
    FrontMenu,
    ModuleWorkspace,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkspaceTab {
    VisualizationAndDistribution,
    TimeSeriesAndGrowthPhases,
    PhasePortraitAndDiagnostics,
    EquationsAndExport,
}

pub struct PopBioApp {
    pub screen: AppScreen,
    pub active_tab: WorkspaceTab,
    pub workset: WorksetState,
    pub playing: bool,
    pub steps_per_frame: usize,
    pub filter_query: String,
    pub new_slot_name: String,
    pub saved_slots: Vec<SaveSlotInfo>,
    pub db: Option<WorksetDb>,
    pub db_error: Option<String>,
    pub dragged_hap_node: Option<usize>,
    pub show_help_modal: bool,
    pub show_about_modal: bool,
    pub status_message: String,
}

impl PopBioApp {
    pub fn new(
        cc: &eframe::CreationContext<'_>,
        initial_workset: WorksetState,
        start_in_front_menu: bool,
    ) -> Self {
        let mut app = Self::new_headless(&cc.egui_ctx, initial_workset, start_in_front_menu);
        if !start_in_front_menu && !WorksetState::is_empirical_module(app.workset.active_module) {
            app.playing = true;
            app.steps_per_frame = 2;
        }
        app
    }

    pub fn new_headless(
        ctx: &egui::Context,
        initial_workset: WorksetState,
        start_in_front_menu: bool,
    ) -> Self {
        ctx.set_visuals(egui::Visuals::dark());
        let (db, db_error) = match WorksetDb::open_default() {
            Ok(d) => (Some(d), None),
            Err(e) => (None, Some(e)),
        };
        let saved_slots = db.as_ref().map(|d| d.list_slots()).unwrap_or_default();
        let initial_status = if let Some(err) = &db_error {
            format!("Ready (SQLite warning: {err})")
        } else {
            "Ready — Choose a module on the Front Menu or press F5 to Quicksave your Workset."
                .to_string()
        };

        Self {
            screen: if start_in_front_menu {
                AppScreen::FrontMenu
            } else {
                AppScreen::ModuleWorkspace
            },
            active_tab: WorkspaceTab::VisualizationAndDistribution,
            workset: initial_workset,
            playing: false,
            steps_per_frame: 1,
            filter_query: String::new(),
            new_slot_name: "Lab_Checkpoint_1".to_string(),
            saved_slots,
            db,
            db_error,
            dragged_hap_node: None,
            show_help_modal: false,
            show_about_modal: false,
            status_message: initial_status,
        }
    }

    fn refresh_slots(&mut self) {
        if let Some(db) = &self.db {
            self.saved_slots = db.list_slots();
        }
    }

    fn save_slot(&mut self, slot_name: &str) {
        if let Some(db) = &self.db {
            match db.save_slot(slot_name, &mut self.workset) {
                Ok(()) => {
                    self.status_message = format!("Saved Workset checkpoint to slot '{slot_name}'.");
                    self.saved_slots = self.db.as_ref().map(|d| d.list_slots()).unwrap_or_default();
                }
                Err(e) => self.status_message = e,
            }
        } else if let Some(err) = &self.db_error {
            self.status_message = format!("Cannot save slot '{slot_name}': {err}");
        }
    }

    fn load_slot(&mut self, slot_name: &str) {
        if let Some(db) = &self.db {
            match db.load_slot(slot_name) {
                Ok(ws) => {
                    self.workset = ws;
                    self.playing = false;
                    self.screen = AppScreen::ModuleWorkspace;
                    self.status_message =
                        format!("Loaded Workset checkpoint '{slot_name}' ({}).", self.workset.active_module.slug());
                }
                Err(e) => self.status_message = e,
            }
        } else if let Some(err) = &self.db_error {
            self.status_message = format!("Cannot load slot '{slot_name}': {err}");
        }
    }

    fn export_workset_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("PopBio Workset", &["json", "popworkset.json"])
            .set_file_name("session.popworkset.json")
            .save_file()
        {
            match self.workset.save_to_file(&path) {
                Ok(()) => {
                    self.status_message = format!("Saved portable Workset to {}", path.display());
                }
                Err(e) => self.status_message = e,
            }
        }
    }

    fn import_workset_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("PopBio Workset", &["json"])
            .pick_file()
        {
            match WorksetState::load_from_file(&path) {
                Ok(ws) => {
                    self.workset = ws;
                    self.playing = false;
                    self.screen = AppScreen::ModuleWorkspace;
                    self.status_message =
                        format!("Restored Workset from {} (SHA-256 verified).", path.display());
                }
                Err(e) => self.status_message = e,
            }
        }
    }

    fn render_parameter_controls(&mut self, ui: &mut egui::Ui) {
        let mod_id = self.workset.active_module;
        let is_empirical = WorksetState::is_empirical_module(mod_id);

        let phase = self.workset.current_phase(mod_id);
        let [pr, pg, pb] = phase.badge_rgb();
        egui::Frame::group(ui.style())
            .fill(Color32::from_rgb(22, 30, 42))
            .inner_margin(Margin::same(8))
            .show(ui, |ui| {
                ui.label(RichText::new("Active Growth / Regime Phase:").size(11.5).color(Color32::GRAY));
                ui.label(
                    RichText::new(format!(" {} ", phase.label()))
                        .strong()
                        .color(Color32::WHITE)
                        .background_color(Color32::from_rgb(pr, pg, pb)),
                );
            });

        ui.add_space(6.0);

        if is_empirical {
            render_empirical_controls(ui, &mut self.workset.empirical, &mut self.status_message);
            return;
        }

        // Time-Advancement Transport Controls
        ui.label(RichText::new("Time-Step & Growth Phase Transport").strong());
        ui.horizontal_wrapped(|ui| {
            let play_lbl = if self.playing {
                "|| Pause (Space)"
            } else {
                "▶ Play (Space)"
            };
            if ui.button(play_lbl).clicked() {
                self.playing = !self.playing;
            }
            if ui.button("Step +1").clicked() {
                self.playing = false;
                self.workset.step_module(mod_id);
            }
            if ui.button("Step +10").clicked() {
                self.playing = false;
                for _ in 0..10 {
                    self.workset.step_module(mod_id);
                }
            }
        });
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("⏭ Next Growth Phase")
                .on_hover_text("Advance time automatically until the population transitions to a new growth phase")
                .clicked()
            {
                self.playing = false;
                let steps = self.workset.advance_to_next_growth_phase(mod_id, 600);
                self.status_message = format!(
                    "Advanced {} steps into '{}'.",
                    steps,
                    self.workset.current_phase(mod_id).label()
                );
            }
            if ui.button("⚡ Perturb").clicked() {
                self.workset.inject_perturbation(mod_id);
                self.status_message = "Applied ecological/genetic perturbation pulse.".to_string();
            }
            if ui.button("↺ Reset").clicked() {
                self.playing = false;
                self.workset.reset_module(mod_id);
                self.status_message = "Reset module to baseline empirical parameters.".to_string();
            }
        });
        ui.add(egui::Slider::new(&mut self.steps_per_frame, 1..=20).text("Speed (steps/frame)"));

        ui.separator();
        ui.label(RichText::new("Biological Parameters & Toggles").strong());

        let mut pending_status: Option<String> = None;
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| match mod_id {
                ModuleId::V1HareLynx => {
                    let s = &mut self.workset.v1_hare_lynx;
                    ui.add(egui::Slider::new(&mut s.r_hare, 0.1..=2.2).text("Hare growth r_H"));
                    ui.add(egui::Slider::new(&mut s.k_hare, 40.0..=500.0).text("Carrying capacity K"));
                    ui.add(egui::Slider::new(&mut s.attack_rate, 0.005..=0.15).text("Lynx attack rate a"));
                    ui.add(egui::Slider::new(&mut s.handling_time, 0.02..=0.80).text("Handling time T_h"));
                    ui.add(egui::Slider::new(&mut s.conversion_eff, 0.05..=0.60).text("Conversion efficiency e"));
                    ui.add(egui::Slider::new(&mut s.lynx_mortality, 0.10..=0.95).text("Lynx mortality m_L"));
                    ui.add(egui::Slider::new(&mut s.spatial_dispersal, 0.05..=1.0).text("Foraging dispersal"));
                    ui.separator();
                    ui.label("Predator Functional Response:");
                    for fr in [
                        FunctionalResponse::TypeILinear,
                        FunctionalResponse::TypeIIHyperbolic,
                        FunctionalResponse::TypeIIISigmoidal,
                    ] {
                        ui.radio_value(&mut s.functional_response, fr, fr.label());
                    }
                    ui.checkbox(&mut s.demographic_noise, "Enable Demographic Stochasticity");
                    if s.demographic_noise {
                        ui.add(egui::Slider::new(&mut s.noise_sigma, 0.01..=0.25).text("Noise sigma"));
                    }
                }

                ModuleId::V2IsleRoyale => {
                    let s = &mut self.workset.v2_isle_royale;
                    ui.add(egui::Slider::new(&mut s.r_moose, 0.08..=0.45).text("Moose intrinsic r_M"));
                    ui.add(egui::Slider::new(&mut s.k_moose_base, 600.0..=3000.0).text("Base Fir K_moose"));
                    ui.add(egui::Slider::new(&mut s.kill_rate, 4.0..=28.0).text("Wolf kill rate"));
                    ui.add(egui::Slider::new(&mut s.wolf_mortality, 0.08..=0.48).text("Wolf mortality m_W"));
                    ui.add(egui::Slider::new(&mut s.allee_threshold, 0.0..=15.0).text("Pack Allee / inbreeding threshold"));
                    ui.add(egui::Slider::new(&mut s.winter_severity, 0.0..=1.0).text("Winter NAO severity"));
                    ui.checkbox(
                        &mut s.ice_bridge_rescue,
                        "Enable Lake Superior Ice-Bridge Wolf Gene Flow",
                    );
                }

                ModuleId::V3SoaySheep => {
                    let s = &mut self.workset.v3_soay_sheep;
                    ui.add(egui::Slider::new(&mut s.lambda_fecundity, 1.2..=3.8).text("Fecundity lambda"));
                    ui.add(egui::Slider::new(&mut s.carrying_scale_k, 400.0..=2200.0).text("Pasture scale K"));
                    ui.add(
                        egui::Slider::new(&mut s.overcompensation_beta, 0.8..=5.5)
                            .text("Overcompensation beta (>2 crashes)"),
                    );
                    ui.add(egui::Slider::new(&mut s.gale_probability, 0.0..=0.80).text("Atlantic gale probability"));
                    ui.add(egui::Slider::new(&mut s.gale_severity, 0.05..=0.85).text("Winter gale mortality index"));
                }

                ModuleId::I1TriboliumLpa => {
                    let s = &mut self.workset.i1_tribolium;
                    ui.add(egui::Slider::new(&mut s.b_recruitment, 2.0..=18.0).text("Larval recruitment b"));
                    ui.add(egui::Slider::new(&mut s.c_el, 0.001..=0.04).text("Larvae-on-egg cannibalism c_el"));
                    ui.add(egui::Slider::new(&mut s.c_ea, 0.001..=0.04).text("Adult-on-egg cannibalism c_ea"));
                    ui.add(egui::Slider::new(&mut s.c_pa, 0.01..=0.80).text("Adult-on-pupa cannibalism c_pa"));
                    ui.add(egui::Slider::new(&mut s.mu_l, 0.02..=0.60).text("Larval mortality mu_l"));
                    ui.add(
                        egui::Slider::new(&mut s.mu_a, 0.05..=0.99)
                            .text("Adult mortality mu_a (bifurcation control)"),
                    );
                    ui.checkbox(
                        &mut s.park_competition_mode,
                        "Park (1948) Two-Species Competition (+ T. confusum)",
                    );
                    if s.park_competition_mode {
                        ui.add(
                            egui::Slider::new(&mut s.interspecific_alpha, 0.2..=2.0)
                                .text("Interspecific inhibition alpha"),
                        );
                    }
                }

                ModuleId::I2HuffakerMites => {
                    let s = &mut self.workset.i2_huffaker;
                    ui.add(
                        egui::Slider::new(&mut s.vaseline_barrier, 0.0..=0.98)
                            .text("Vaseline Barrier Maze (0=open crash, 0.75=waves)"),
                    );
                    ui.add(egui::Slider::new(&mut s.r_prey, 0.2..=1.6).text("Prey mite growth r"));
                    ui.add(egui::Slider::new(&mut s.k_patch, 60.0..=400.0).text("Per-orange K"));
                    ui.add(egui::Slider::new(&mut s.prey_dispersal, 0.02..=0.60).text("Prey silk dispersal D_N"));
                    ui.add(egui::Slider::new(&mut s.pred_dispersal, 0.01..=0.40).text("Predator walking D_P"));
                }

                ModuleId::I3AphidParasitoid => {
                    let s = &mut self.workset.i3_aphid;
                    ui.add(egui::Slider::new(&mut s.lambda_host, 1.2..=4.0).text("Aphid finite rate lambda"));
                    ui.add(egui::Slider::new(&mut s.search_area_a, 0.005..=0.12).text("Wasp search area a"));
                    ui.add(
                        egui::Slider::new(&mut s.clumping_k, 0.15..=5.0)
                            .text("Neg-Binomial clumping k (<1 stabilizes)"),
                    );
                    ui.add(
                        egui::Slider::new(&mut s.stem_refuge_fraction, 0.0..=0.50)
                            .text("Physical stem refuge fraction"),
                    );
                }

                ModuleId::M1EcoliBatchChemostat => {
                    let s = &mut self.workset.m1_ecoli;
                    ui.checkbox(&mut s.chemostat_mode, "Continuous Monod Chemostat Mode");
                    ui.checkbox(&mut s.diauxic_mode, "Glucose + Lactose Diauxic Shift");
                    ui.add(egui::Slider::new(&mut s.mu_max_glucose, 0.2..=1.8).text("mu_max Glucose (h^-1)"));
                    ui.add(egui::Slider::new(&mut s.ks_glucose, 0.02..=1.2).text("Monod K_s Glucose (g/L)"));
                    ui.add(egui::Slider::new(&mut s.lag_rate_h, 0.15..=2.5).text("Lag acclimation rate (h^-1)"));
                    ui.add(egui::Slider::new(&mut s.death_rate_h, 0.01..=0.35).text("Senescence death rate k_d"));
                    if s.chemostat_mode {
                        ui.add(
                            egui::Slider::new(&mut s.dilution_rate_d, 0.02..=1.5)
                                .text("Chemostat Dilution D (washout > mu_max)"),
                        );
                        ui.add(
                            egui::Slider::new(&mut s.feed_glucose_s0, 0.5..=12.0)
                                .text("Feed Reservoir Glucose S_0 (g/L)"),
                        );
                    }
                }

                ModuleId::M2GauseProtists => {
                    let s = &mut self.workset.m2_gause;
                    ui.checkbox(
                        &mut s.competition_mode,
                        "Competitive Exclusion Mode (P. aurelia vs P. caudatum)",
                    );
                    if s.competition_mode {
                        ui.add(egui::Slider::new(&mut s.alpha_aur_on_caud, 0.2..=2.2).text("alpha (aurelia -> caudatum)"));
                        ui.add(egui::Slider::new(&mut s.alpha_caud_on_aur, 0.2..=2.2).text("alpha (caudatum -> aurelia)"));
                    } else {
                        ui.add(
                            egui::Slider::new(&mut s.sediment_refuge, 0.0..=0.65)
                                .text("Oat-Sediment Refuge Fraction"),
                        );
                        ui.checkbox(&mut s.periodic_immigration, "Gause 3-Day Periodic Immigration");
                    }
                }

                ModuleId::M3PhagePlaque => {
                    let s = &mut self.workset.m3_phage;
                    ui.add(egui::Slider::new(&mut s.adsorption_phi, 0.001..=0.03).text("Phage adsorption phi"));
                    ui.add(egui::Slider::new(&mut s.lysis_rate, 0.3..=3.0).text("Lysis rate 1/tau (h^-1)"));
                    ui.add(egui::Slider::new(&mut s.burst_size_beta, 10.0..=150.0).text("Burst size beta (PFU/cell)"));
                    ui.add(
                        egui::Slider::new(&mut s.resistance_mut_rate, 0.0..=0.01)
                            .text("Receptor-mutant emergence rate"),
                    );
                }

                ModuleId::T1LevinsMetapop => {
                    let s = &mut self.workset.t1_levins;
                    ui.add(egui::Slider::new(&mut s.colonization_c, 0.05..=1.5).text("Colonization rate c"));
                    ui.add(egui::Slider::new(&mut s.extinction_e, 0.02..=0.80).text("Local extinction rate e"));
                    if ui
                        .add(egui::Slider::new(&mut s.destruction_h, 0.0..=0.95).text("Habitat destroyed h"))
                        .changed()
                    {
                        s.apply_destruction_fraction();
                    }
                    ui.label(format!(
                        "Persistence Threshold h_crit = {:.1}% | Equilibrium p* = {:.1}%",
                        s.critical_destruction_threshold() * 100.0,
                        s.analytical_equilibrium() * 100.0
                    ));
                }

                ModuleId::T2MayBifurcation => {
                    let s = &mut self.workset.t2_may;
                    ui.checkbox(&mut s.ricker_map, "Use Ricker Exponential Map (prevent negative N)");
                    let max_r = if s.ricker_map { 4.2 } else { 3.0 };
                    ui.add(egui::Slider::new(&mut s.r_param, 0.2..=max_r).text("Intrinsic growth parameter r"));
                    ui.add(egui::Slider::new(&mut s.carrying_k, 20.0..=500.0).text("Carrying capacity K"));
                }

                ModuleId::T3WrightFisher => {
                    let s = &mut self.workset.t3_wright_fisher;
                    ui.add(egui::Slider::new(&mut s.n_diploid, 5..=1000).logarithmic(true).text("Diploid Deme Size N_e"));
                    ui.add(egui::Slider::new(&mut s.num_demes, 1..=12).text("Number of Demes"));
                    ui.add(egui::Slider::new(&mut s.initial_p0, 0.02..=0.98).text("Initial Allele Freq p_0"));
                    ui.separator();
                    ui.label("Genotype Relative Fitnesses:");
                    ui.add(egui::Slider::new(&mut s.w_aa_dom, 0.2..=1.5).text("Fitness w_AA"));
                    ui.add(egui::Slider::new(&mut s.w_aa_het, 0.2..=1.5).text("Fitness w_Aa ( >1 overdominance)"));
                    ui.add(egui::Slider::new(&mut s.w_aa_rec, 0.2..=1.5).text("Fitness w_aa"));
                    ui.separator();
                    ui.add(egui::Slider::new(&mut s.mu_forward, 0.0..=0.02).text("Forward mutation A -> a"));
                    ui.add(egui::Slider::new(&mut s.nu_backward, 0.0..=0.02).text("Back mutation a -> A"));
                    ui.add(egui::Slider::new(&mut s.migration_m, 0.0..=0.25).text("Island migration rate m"));
                }

                ModuleId::T4KingmanCoalescent => {
                    let s = &mut self.workset.t4_coalescent;
                    ui.add(egui::Slider::new(&mut s.sample_size_n, 4..=36).text("Sampled Chromosomes n"));
                    ui.add(egui::Slider::new(&mut s.effective_ne, 50.0..=10_000.0).logarithmic(true).text("Effective Size N_e"));
                    ui.add(egui::Slider::new(&mut s.theta_locus, 1.0..=50.0).text("Population Mutation 4*N_e*u"));
                }

                ModuleId::T5WrightLandscape => {
                    let s = &mut self.workset.t5_landscape;
                    ui.label("Shifting Balance & Deme Parameters:");
                    ui.add(egui::Slider::new(&mut s.num_demes, 1..=14).text("Number of Demes"));
                    ui.add(
                        egui::Slider::new(&mut s.effective_ne, 6.0..=500.0)
                            .logarithmic(true)
                            .text("Deme Effective Size N_e (low = valley drift)"),
                    );
                    ui.add(egui::Slider::new(&mut s.selection_strength, 0.02..=0.60).text("Selection gradient strength"));
                    ui.add(egui::Slider::new(&mut s.interdeme_migration_m, 0.0..=0.15).text("Interdeme dispersion m (Phase III)"));
                    ui.separator();
                    ui.label("3D Adaptive Surface Topography:");
                    ui.add(egui::Slider::new(&mut s.peak1_height, 0.10..=1.20).text("Local Peak 1 height"));
                    ui.add(egui::Slider::new(&mut s.peak2_height, 0.30..=1.60).text("Global Peak 2 height"));
                    ui.add(egui::Slider::new(&mut s.valley_ruggedness, 0.0..=0.90).text("Epistatic valley ruggedness"));
                    ui.separator();
                    ui.label("3D Orbit Camera (or drag directly on canvas):");
                    ui.add(egui::Slider::new(&mut s.camera_azimuth_deg, -180.0..=180.0).text("Azimuth (deg)"));
                    ui.add(egui::Slider::new(&mut s.camera_elevation_deg, 12.0..=72.0).text("Elevation (deg)"));
                }

                ModuleId::T6LeslieLefkovitch => {
                    let s = &mut self.workset.t6_matrix;
                    ui.label("Life-History Projection Matrix Preset:");
                    let mut chosen = s.preset;
                    for pr in [
                        MatrixPreset::LoggerheadTurtleLefkovitch7,
                        MatrixPreset::SongbirdLeslie4,
                        MatrixPreset::WoodlandPerennialLefkovitch5,
                    ] {
                        if ui.radio_value(&mut chosen, pr, pr.label()).changed() {
                            s.load_preset(chosen);
                        }
                    }
                    if s.preset == MatrixPreset::Custom {
                        ui.label(
                            RichText::new(format!(
                                "Active Mode: {} ({}x{})",
                                MatrixPreset::Custom.label(),
                                s.base_matrix.len(),
                                s.base_matrix.len()
                            ))
                            .small()
                            .color(Color32::from_rgb(56, 189, 248)),
                        );
                    }
                    ui.separator();
                    ui.horizontal(|ui| {
                        let mut k_edit = s.stage_names.len().clamp(2, 12);
                        ui.label("Stage count k:");
                        if ui
                            .add(egui::DragValue::new(&mut k_edit).range(2..=12).speed(0.1))
                            .changed()
                        {
                            s.resize_stages(k_edit);
                            pending_status = Some(format!(
                                "Resized projection matrix to {k_edit}x{k_edit} stages."
                            ));
                        }
                        if ui.button("Import CSV...").clicked() {
                            if let Some(path) = rfd::FileDialog::new()
                                .add_filter("CSV Matrix", &["csv", "txt"])
                                .pick_file()
                            {
                                match fs::read_to_string(&path) {
                                    Ok(raw) => match s.load_from_csv(&raw) {
                                        Ok(()) => {
                                            pending_status = Some(format!(
                                                "Imported {}x{} projection matrix from {}",
                                                s.base_matrix.len(),
                                                s.base_matrix.len(),
                                                path.display()
                                            ));
                                        }
                                        Err(e) => {
                                            pending_status = Some(format!("Matrix CSV error: {e}"));
                                        }
                                    },
                                    Err(e) => {
                                        let err = format!("Failed to read {}: {e}", path.display());
                                        s.last_matrix_error = Some(err.clone());
                                        pending_status = Some(err);
                                    }
                                }
                            }
                        }
                        if ui.button("Export CSV...").clicked() {
                            if let Some(path) = rfd::FileDialog::new()
                                .add_filter("CSV Matrix", &["csv"])
                                .set_file_name("leslie_lefkovitch_matrix.csv")
                                .save_file()
                            {
                                match fs::write(&path, s.to_matrix_csv()) {
                                    Ok(()) => {
                                        s.last_matrix_error = None;
                                        pending_status = Some(format!(
                                            "Exported {}x{} projection matrix to {}",
                                            s.base_matrix.len(),
                                            s.base_matrix.len(),
                                            path.display()
                                        ));
                                    }
                                    Err(e) => {
                                        let err = format!("Failed to write {}: {e}", path.display());
                                        s.last_matrix_error = Some(err.clone());
                                        pending_status = Some(err);
                                    }
                                }
                            }
                        }
                    });
                    if let Some(err) = &s.last_matrix_error {
                        ui.colored_label(Color32::from_rgb(248, 113, 113), format!("Error: {err}"));
                    }

                    let mut changed = false;
                    egui::CollapsingHeader::new("Edit k x k Projection Matrix A_ij")
                        .default_open(true)
                        .show(ui, |ui| {
                            let k = s.base_matrix.len();
                            egui::ScrollArea::horizontal()
                                .id_salt("t6_matrix_scroll")
                                .show(ui, |ui| {
                                    egui::Grid::new("t6_matrix_editor")
                                        .spacing([4.0, 3.0])
                                        .show(ui, |ui| {
                                            ui.label(RichText::new("A[r,c]").monospace().small());
                                            for c in 0..k {
                                                let lbl = format!("S{}", c + 1);
                                                let tip = s
                                                    .stage_names
                                                    .get(c)
                                                    .cloned()
                                                    .unwrap_or_else(|| lbl.clone());
                                                ui.label(RichText::new(lbl).monospace().strong().small())
                                                    .on_hover_text(tip);
                                            }
                                            ui.end_row();
                                            for r in 0..k {
                                                let r_lbl = format!("S{}", r + 1);
                                                let r_tip = s
                                                    .stage_names
                                                    .get(r)
                                                    .cloned()
                                                    .unwrap_or_else(|| r_lbl.clone());
                                                ui.label(
                                                    RichText::new(r_lbl).monospace().strong().small(),
                                                )
                                                .on_hover_text(r_tip);
                                                for c in 0..k {
                                                    if ui
                                                        .add(
                                                            egui::DragValue::new(
                                                                &mut s.base_matrix[r][c],
                                                            )
                                                            .speed(0.01)
                                                            .range(0.0..=10_000.0)
                                                            .max_decimals(4),
                                                        )
                                                        .changed()
                                                    {
                                                        s.preset = MatrixPreset::Custom;
                                                        changed = true;
                                                    }
                                                }
                                                ui.end_row();
                                            }
                                        });
                                });
                        });

                    ui.separator();
                    ui.label("Conservation & Vital-Rate Multipliers:");
                    changed |= ui
                        .add(egui::Slider::new(&mut s.fecundity_mult, 0.2..=2.5).text("Fecundity F_j multiplier"))
                        .changed();
                    changed |= ui
                        .add(
                            egui::Slider::new(&mut s.juvenile_survival_mult, 0.4..=1.35)
                                .text("Juvenile / Subadult survival (e.g. TEDs)"),
                        )
                        .changed();
                    changed |= ui
                        .add(egui::Slider::new(&mut s.adult_survival_mult, 0.4..=1.20).text("Adult breeder survival"))
                        .changed();
                    ui.checkbox(&mut s.density_regulated, "Beverton-Holt Recruitment Ceiling");
                    if s.density_regulated {
                        ui.add(egui::Slider::new(&mut s.recruitment_k, 200.0..=20_000.0).logarithmic(true).text("Recruitment scale K"));
                    }
                    if changed {
                        s.recompute_eigen();
                    }
                    ui.separator();
                    ui.label(format!(
                        "Eigen Summary: lambda_1 = {:.4} | r = {:+.4} | R_0 = {:.2} | Damping rho = {:.3}",
                        s.eigen.dominant_lambda,
                        s.eigen.intrinsic_r,
                        s.eigen.net_reproductive_r0,
                        s.eigen.damping_ratio
                    ));
                }

                ModuleId::E1SequencePopGen
                | ModuleId::E2HaplotypeNetwork
                | ModuleId::E3BayesianSkyline
                | ModuleId::E4SlidingWindowScan => {}
            });
        if let Some(msg) = pending_status {
            self.status_message = msg;
        }
    }
}

impl PopBioApp {
    pub fn ui_step(&mut self, ctx: &egui::Context) {
        // Global keyboard shortcuts
        if ctx.input(|i| i.key_pressed(Key::Space))
            && self.screen == AppScreen::ModuleWorkspace
            && !ctx.wants_keyboard_input()
        {
            self.playing = !self.playing;
        }
        if ctx.input(|i| i.key_pressed(Key::F5)) {
            self.save_slot("Quicksave");
        }
        if ctx.input(|i| i.key_pressed(Key::F9)) {
            self.load_slot("Quicksave");
        }

        if self.playing && self.screen == AppScreen::ModuleWorkspace {
            let active = self.workset.active_module;
            for _ in 0..self.steps_per_frame {
                self.workset.step_module(active);
            }
            ctx.request_repaint_after(Duration::from_millis(33));
        }

        if self.screen == AppScreen::ModuleWorkspace
            && self.workset.active_module == ModuleId::E2HaplotypeNetwork
            && self.workset.empirical.hap_sim_playing
        {
            self.workset.empirical.step_haplotype_evolution();
            ctx.request_repaint_after(Duration::from_millis(140));
        }

        // Persistent Top Navigation & Workset Bar
        egui::TopBottomPanel::top("top_nav_bar").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                if ui
                    .selectable_label(self.screen == AppScreen::FrontMenu, "Front Menu (Module Hub)")
                    .clicked()
                {
                    self.playing = false;
                    self.screen = AppScreen::FrontMenu;
                    self.saved_slots = self.db.as_ref().map(|d| d.list_slots()).unwrap_or_default();
                }

                ui.separator();

                for group in ModuleGroup::ALL {
                    ui.menu_button(group.short_badge(), |ui| {
                        for m in ModuleId::ALL.into_iter().filter(|m| m.group() == group) {
                            let meta = m.metadata();
                            if ui
                                .button(format!("[{}] {}", meta.code, meta.title))
                                .clicked()
                            {
                                self.workset.active_module = m;
                                self.screen = AppScreen::ModuleWorkspace;
                                self.active_tab = WorkspaceTab::VisualizationAndDistribution;
                                self.playing = !WorksetState::is_empirical_module(m);
                                self.steps_per_frame = 2;
                                ui.close_menu();
                            }
                        }
                    });
                }

                ui.separator();

                ui.menu_button("Workset (Save Game)", |ui| {
                    if ui.button("Quicksave Slot (F5)").clicked() {
                        self.save_slot("Quicksave");
                        ui.close_menu();
                    }
                    if ui.button("Quickload Slot (F9)").clicked() {
                        self.load_slot("Quicksave");
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Save Portable .popworkset.json...").clicked() {
                        self.export_workset_dialog();
                        ui.close_menu();
                    }
                    if ui.button("Load Portable .popworkset.json...").clicked() {
                        self.import_workset_dialog();
                        ui.close_menu();
                    }
                });

                if ui.button("Help & Hotkeys").clicked() {
                    self.show_help_modal = true;
                }
                if ui.button("About").clicked() {
                    self.show_about_modal = true;
                }
            });
        });

        // Bottom Status Bar
        egui::TopBottomPanel::bottom("bottom_status_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                let meta = self.workset.active_module.metadata();
                ui.label(
                    RichText::new(format!("[{}] {}", meta.code, meta.title))
                        .monospace()
                        .strong(),
                );
                ui.separator();
                ui.label(&self.status_message);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                        RichText::new("Copyright CRV Science Network")
                            .size(11.5)
                            .color(Color32::from_rgb(160, 174, 192)),
                    );
                });
            });
        });

        if self.show_about_modal {
            egui::Window::new("About PopBio Studio")
                .open(&mut self.show_about_modal)
                .resizable(false)
                .default_width(480.0)
                .show(ctx, |ui| {
                    ui.label(
                        RichText::new(format!("PopBio Studio v{}", env!("CARGO_PKG_VERSION")))
                            .size(16.0)
                            .strong(),
                    );
                    ui.label(
                        "Native Rust population biology, Sewall Wright 3D adaptive landscapes, \
                         Leslie/Lefkovitch stage matrices, Bayesian Skyline demographic inference, \
                         and haplotype network workstation.",
                    );
                    ui.separator();
                    ui.label(RichText::new("Copyright CRV Science Network").strong());
                    ui.label(
                        "Licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later). \
                         Provided with ABSOLUTELY NO WARRANTY.",
                    );
                    ui.separator();
                    ui.hyperlink_to(
                        "CRV Science Network (https://crvscience.com)",
                        "https://crvscience.com",
                    );
                    ui.hyperlink_to(
                        "CRV Science Open-Source Apps (https://crvscience.com/apps)",
                        "https://crvscience.com/apps",
                    );
                });
        }

        if self.show_help_modal {
            egui::Window::new("PopBio Studio — Controls, Worksets & Citations")
                .open(&mut self.show_help_modal)
                .resizable(true)
                .default_width(520.0)
                .show(ctx, |ui| {
                    ui.label(RichText::new("Keyboard Shortcuts & Lab Workflow").strong());
                    ui.label("• Space: Play / Pause active dynamical population simulation");
                    ui.label("• Next Growth Phase button: Step automatically until the biological phase shifts (e.g., Lag -> Log -> Diauxic -> Stationary -> Senescence)");
                    ui.label("• F5 / F9: Instant SQLite Quicksave / Quickload of all 19 modules + custom sequence alignments");
                    ui.label("• Workset (.popworkset.json): Portable SHA-256 verified lab save-game files shareable across machines");
                    ui.separator();
                    ui.label(RichText::new("Empirical Sequence Header Convention").strong());
                    ui.label("Use '>SampleID|DemeName' or '>SampleID [Pop=DemeName]' in FASTA headers to automatically color Haplotype Network pie slices and compute pairwise Hudson F_ST / D_XY.");
                });
        }

        match self.screen {
            AppScreen::FrontMenu => {
                egui::CentralPanel::default().show(ctx, |ui| {
                    let act = render_front_menu(
                        ui,
                        &mut self.workset,
                        &mut self.filter_query,
                        &mut self.new_slot_name,
                        &self.saved_slots,
                    );
                    match act {
                        FrontMenuAction::None => {}
                        FrontMenuAction::LaunchModule(m) => {
                            self.workset.active_module = m;
                            self.screen = AppScreen::ModuleWorkspace;
                            self.active_tab = WorkspaceTab::VisualizationAndDistribution;
                            self.playing = !WorksetState::is_empirical_module(m);
                            self.steps_per_frame = 2;
                        }
                        FrontMenuAction::SaveNamedSlot(name) => self.save_slot(&name),
                        FrontMenuAction::LoadNamedSlot(name) => self.load_slot(&name),
                        FrontMenuAction::DeleteNamedSlot(name) => {
                            if let Some(db) = &self.db {
                                match db.delete_slot(&name) {
                                    Ok(()) => {
                                        self.status_message =
                                            format!("Deleted Workset checkpoint '{name}'.");
                                        self.refresh_slots();
                                    }
                                    Err(e) => self.status_message = e,
                                }
                            } else if let Some(err) = &self.db_error {
                                self.status_message =
                                    format!("Cannot delete slot '{name}': {err}");
                            }
                        }
                        FrontMenuAction::ExportWorksetFile => self.export_workset_dialog(),
                        FrontMenuAction::ImportWorksetFile => self.import_workset_dialog(),
                    }
                });
            }

            AppScreen::ModuleWorkspace => {
                egui::SidePanel::left("module_controls_panel")
                    .resizable(true)
                    .default_width(335.0)
                    .show(ctx, |ui| {
                        self.render_parameter_controls(ui);
                    });

                egui::CentralPanel::default().show(ctx, |ui| {
                    let mod_id = self.workset.active_module;
                    let meta = mod_id.metadata();

                    ui.horizontal_wrapped(|ui| {
                        ui.label(
                            RichText::new(format!("[{}] {}", meta.code, meta.title))
                                .size(17.5)
                                .strong(),
                        );
                        ui.label(
                            RichText::new(format!("({})", meta.organisms))
                                .italics()
                                .color(Color32::from_rgb(138, 201, 255)),
                        );
                    });
                    ui.separator();

                    // Empirical FASTA Suite Workspace Tabs (E1, E2, E3, E4)
                    if WorksetState::is_empirical_module(mod_id) {
                        ui.horizontal_wrapped(|ui| {
                            if ui
                                .selectable_label(
                                    mod_id == ModuleId::E1SequencePopGen,
                                    "1. Global PopGen, Neutrality & F_ST",
                                )
                                .clicked()
                            {
                                self.workset.active_module = ModuleId::E1SequencePopGen;
                            }
                            if ui
                                .selectable_label(
                                    mod_id == ModuleId::E2HaplotypeNetwork,
                                    "2. Haplotype Network Studio (MJ/MSN/TCS)",
                                )
                                .clicked()
                            {
                                self.workset.active_module = ModuleId::E2HaplotypeNetwork;
                            }
                            if ui
                                .selectable_label(
                                    mod_id == ModuleId::E3BayesianSkyline,
                                    "3. Bayesian Skyline & Demographic Expansion",
                                )
                                .clicked()
                            {
                                self.workset.active_module = ModuleId::E3BayesianSkyline;
                            }
                            if ui
                                .selectable_label(
                                    mod_id == ModuleId::E4SlidingWindowScan,
                                    "4. Sliding-Window Neutrality & Sweep Scan",
                                )
                                .clicked()
                            {
                                self.workset.active_module = ModuleId::E4SlidingWindowScan;
                            }
                        });
                        ui.separator();

                        match self.workset.active_module {
                            ModuleId::E1SequencePopGen => {
                                render_popgen_summary_view(ui, &mut self.workset.empirical);
                            }
                            ModuleId::E2HaplotypeNetwork => {
                                render_haplotype_network_view(
                                    ui,
                                    &mut self.workset.empirical,
                                    &mut self.dragged_hap_node,
                                );
                            }
                            ModuleId::E3BayesianSkyline => {
                                render_skyline_and_expansion_view(ui, &mut self.workset.empirical);
                            }
                            ModuleId::E4SlidingWindowScan => {
                                render_sliding_window_view(ui, &self.workset.empirical);
                            }
                            _ => {}
                        }
                        return;
                    }

                    // 4-Tab Bar for Dynamical & Theoretical Modules (V1..T4)
                    ui.horizontal_wrapped(|ui| {
                        ui.selectable_value(
                            &mut self.active_tab,
                            WorkspaceTab::VisualizationAndDistribution,
                            "1. Spatial Movement & Distribution",
                        );
                        ui.selectable_value(
                            &mut self.active_tab,
                            WorkspaceTab::TimeSeriesAndGrowthPhases,
                            "2. Population Graphs & Growth Phases",
                        );
                        ui.selectable_value(
                            &mut self.active_tab,
                            WorkspaceTab::PhasePortraitAndDiagnostics,
                            "3. Phase Portrait / Bifurcation / SFS",
                        );
                        ui.selectable_value(
                            &mut self.active_tab,
                            WorkspaceTab::EquationsAndExport,
                            "4. Governing Equations & CSV Export",
                        );
                    });
                    ui.separator();

                    match self.active_tab {
                        WorkspaceTab::VisualizationAndDistribution => {
                            render_spatial_canvas(ui, mod_id, &mut self.workset);
                        }
                        WorkspaceTab::TimeSeriesAndGrowthPhases => {
                            render_time_series_tab(ui, mod_id, &self.workset);
                        }
                        WorkspaceTab::PhasePortraitAndDiagnostics => {
                            render_phase_portrait_tab(ui, mod_id, &self.workset);
                        }
                        WorkspaceTab::EquationsAndExport => {
                            ui.label(RichText::new("Governing Mathematical Formulation").strong());
                            egui::Frame::group(ui.style())
                                .fill(Color32::from_rgb(18, 24, 34))
                                .inner_margin(Margin::same(10))
                                .show(ui, |ui| {
                                    ui.label(RichText::new(meta.equations_plain).monospace());
                                });
                            ui.add_space(8.0);
                            ui.label(RichText::new("Empirical Organism & Literature Provenance").strong());
                            ui.label(meta.summary);
                            ui.hyperlink_to(meta.citation_label, meta.citation_url);

                            ui.add_space(12.0);
                            if ui.button("Export Time-Series Trajectory to CSV...").clicked() {
                                if let Some(path) = rfd::FileDialog::new()
                                    .set_file_name(format!("{}_timeseries.csv", mod_id.slug()))
                                    .save_file()
                                {
                                    let csv = self.workset.export_module_csv(mod_id);
                                    if let Err(e) = fs::write(&path, csv) {
                                        self.status_message = format!("CSV write error: {e}");
                                    } else {
                                        self.status_message =
                                            format!("Exported CSV to {}", path.display());
                                    }
                                }
                            }
                        }
                    }
                });
            }
        }
    }
}

impl eframe::App for PopBioApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.ui_step(ctx);
    }
}

fn block_on_future<F: std::future::Future>(fut: F) -> F::Output {
    use std::sync::Arc;
    struct ThreadUnparker(std::thread::Thread);
    impl std::task::Wake for ThreadUnparker {
        fn wake(self: Arc<Self>) {
            self.0.unpark();
        }
    }
    let waker = std::task::Waker::from(Arc::new(ThreadUnparker(std::thread::current())));
    let mut cx = std::task::Context::from_waker(&waker);
    let mut fut = std::pin::pin!(fut);
    loop {
        match fut.as_mut().poll(&mut cx) {
            std::task::Poll::Ready(val) => return val,
            std::task::Poll::Pending => std::thread::park(),
        }
    }
}

/// Render an offscreen 1600x1000 RGBA PNG screenshot without opening any OS GUI window.
pub fn render_headless_screenshot(out_path: &std::path::Path, view: u8) -> Result<(), String> {
    use eframe::egui_wgpu::{wgpu, Renderer, ScreenDescriptor};

    const WIDTH: u32 = 1600;
    const HEIGHT: u32 = 1000;
    const BYTES_PER_ROW: u32 = WIDTH * 4;

    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
    let adapter = block_on_future(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        compatible_surface: None,
        force_fallback_adapter: false,
    }))
    .or_else(|| {
        block_on_future(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::LowPower,
            compatible_surface: None,
            force_fallback_adapter: true,
        }))
    })
    .ok_or_else(|| "No GPU adapter available for offscreen rendering".to_string())?;

    let (device, queue) =
        block_on_future(adapter.request_device(&wgpu::DeviceDescriptor::default(), None))
            .map_err(|e| format!("Failed to request offscreen GPU device: {e}"))?;

    let format = wgpu::TextureFormat::Rgba8Unorm;
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("popbio-offscreen-target"),
        size: wgpu::Extent3d {
            width: WIDTH,
            height: HEIGHT,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let texture_view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("popbio-offscreen-staging"),
        size: (BYTES_PER_ROW * HEIGHT) as u64,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let mut renderer = Renderer::new(&device, format, None, 1, false);
    let egui_ctx = egui::Context::default();
    egui_ctx.set_pixels_per_point(1.0);

    let mut ws = WorksetState::default();
    match view {
        2 => {
            ws.active_module = ModuleId::E3BayesianSkyline;
        }
        3 => {
            ws.active_module = ModuleId::E2HaplotypeNetwork;
        }
        4 => {
            ws.active_module = ModuleId::V1HareLynx;
        }
        _ => {
            ws.active_module = ModuleId::T5WrightLandscape;
            for _ in 0..35 {
                ws.step_module(ModuleId::T5WrightLandscape);
            }
        }
    }

    let mut app = PopBioApp::new_headless(&egui_ctx, ws, false);
    app.active_tab = WorkspaceTab::VisualizationAndDistribution;

    for frame_idx in 0..4_u32 {
        let raw_input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(WIDTH as f32, HEIGHT as f32),
            )),
            viewports: std::iter::once((
                egui::ViewportId::ROOT,
                egui::ViewportInfo {
                    native_pixels_per_point: Some(1.0),
                    ..Default::default()
                },
            ))
            .collect(),
            time: Some(frame_idx as f64 / 60.0),
            predicted_dt: 1.0 / 60.0,
            ..Default::default()
        };
        let full_output = egui_ctx.run(raw_input, |ctx| {
            app.ui_step(ctx);
        });
        for (id, image_delta) in &full_output.textures_delta.set {
            renderer.update_texture(&device, &queue, *id, image_delta);
        }
        let paint_jobs = egui_ctx.tessellate(full_output.shapes, full_output.pixels_per_point);
        let screen_descriptor = ScreenDescriptor {
            size_in_pixels: [WIDTH, HEIGHT],
            pixels_per_point: 1.0,
        };
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("popbio-offscreen-encoder"),
        });
        let user_cmd_bufs = renderer.update_buffers(
            &device,
            &queue,
            &mut encoder,
            &paint_jobs,
            &screen_descriptor,
        );
        {
            let render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("popbio-offscreen-pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &texture_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                occlusion_query_set: None,
                timestamp_writes: None,
            });
            let mut render_pass = render_pass.forget_lifetime();
            renderer.render(&mut render_pass, &paint_jobs, &screen_descriptor);
        }
        if frame_idx == 3 {
            encoder.copy_texture_to_buffer(
                wgpu::TexelCopyTextureInfo {
                    texture: &texture,
                    mip_level: 0,
                    origin: wgpu::Origin3d::ZERO,
                    aspect: wgpu::TextureAspect::All,
                },
                wgpu::TexelCopyBufferInfo {
                    buffer: &staging,
                    layout: wgpu::TexelCopyBufferLayout {
                        offset: 0,
                        bytes_per_row: Some(BYTES_PER_ROW),
                        rows_per_image: Some(HEIGHT),
                    },
                },
                wgpu::Extent3d {
                    width: WIDTH,
                    height: HEIGHT,
                    depth_or_array_layers: 1,
                },
            );
        }
        queue.submit(user_cmd_bufs.into_iter().chain(std::iter::once(encoder.finish())));
        device.poll(wgpu::Maintain::Wait);
        for id in &full_output.textures_delta.free {
            renderer.free_texture(id);
        }
    }

    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |res| {
        let _ = tx.send(res);
    });
    device.poll(wgpu::Maintain::Wait);
    rx.recv()
        .map_err(|e| format!("GPU readback channel error: {e}"))?
        .map_err(|e| format!("GPU buffer map failed: {e}"))?;
    let mapped = slice.get_mapped_range();

    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .map_err(|e| format!("Failed to create output directory: {e}"))?;
        }
    }
    let file = std::fs::File::create(out_path)
        .map_err(|e| format!("Failed to create PNG file {}: {e}", out_path.display()))?;
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), WIDTH, HEIGHT);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder
        .write_header()
        .map_err(|e| format!("PNG header error: {e}"))?;
    writer
        .write_image_data(&mapped)
        .map_err(|e| format!("PNG write error: {e}"))?;
    drop(mapped);
    staging.unmap();
    Ok(())
}
