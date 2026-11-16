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

use egui::{Align2, Color32, CornerRadius, FontId, Painter, Pos2, Rect, Sense, Stroke, Ui, Vec2};

use crate::model::{Agent2D, ModuleId};
use crate::workset::WorksetState;

fn draw_box_outline(painter: &Painter, rect: Rect, stroke: Stroke) {
    let lt = rect.left_top();
    let rt = rect.right_top();
    let rb = rect.right_bottom();
    let lb = rect.left_bottom();
    painter.line_segment([lt, rt], stroke);
    painter.line_segment([rt, rb], stroke);
    painter.line_segment([rb, lb], stroke);
    painter.line_segment([lb, lt], stroke);
}

fn draw_marginal_distribution_bars(painter: &Painter, arena: Rect, agents: &[Agent2D]) {
    if agents.is_empty() {
        return;
    }
    let bins = 16usize;
    let mut x_hist0 = [0usize; 16];
    let mut x_hist1 = [0usize; 16];
    for a in agents {
        let b = ((a.x.clamp(0.0, 0.9999) * (bins as f64)) as usize).min(bins - 1);
        if a.kind == 0 {
            x_hist0[b] += 1;
        } else {
            x_hist1[b] += 1;
        }
    }
    let max_bin = x_hist0
        .iter()
        .zip(x_hist1.iter())
        .map(|(a, b)| a + b)
        .max()
        .unwrap_or(1)
        .max(1) as f32;

    let strip_h = 26.0;
    let cell_w = arena.width() / (bins as f32);
    for i in 0..bins {
        let x0 = arena.left() + (i as f32) * cell_w;
        let h0 = ((x_hist0[i] as f32) / max_bin) * strip_h;
        let h1 = ((x_hist1[i] as f32) / max_bin) * strip_h;
        if h0 > 0.5 {
            let r0 = Rect::from_min_max(
                Pos2::new(x0 + 1.0, arena.bottom() - h0),
                Pos2::new(x0 + cell_w - 1.0, arena.bottom()),
            );
            painter.rect_filled(r0, CornerRadius::same(2), Color32::from_rgba_unmultiplied(72, 199, 142, 130));
        }
        if h1 > 0.5 {
            let r1 = Rect::from_min_max(
                Pos2::new(x0 + 1.0, arena.bottom() - h0 - h1),
                Pos2::new(x0 + cell_w - 1.0, arena.bottom() - h0),
            );
            painter.rect_filled(r1, CornerRadius::same(2), Color32::from_rgba_unmultiplied(242, 113, 78, 145));
        }
    }
}

/// Renders the interactive 2D spatial movement, patch distribution, microbial microcolony/plaque,
/// De Finetti HWE triangle, or Kingman coalescent genealogy view for the active simulation module.
pub fn render_spatial_canvas(ui: &mut Ui, mod_id: ModuleId, workset: &mut WorksetState) {
    let avail = ui.available_size_before_wrap();
    let canvas_size = Vec2::new(avail.x.max(360.0), (avail.y - 36.0).max(320.0));
    let (response, painter) = ui.allocate_painter(canvas_size, Sense::click_and_drag());
    let rect = response.rect;

    painter.rect_filled(rect, CornerRadius::same(6), Color32::from_rgb(13, 19, 28));
    draw_box_outline(&painter, rect, Stroke::new(1.2, Color32::from_rgb(56, 72, 96)));

    let pad = 18.0;
    let inner = rect.shrink(pad);

    match mod_id {
        ModuleId::V1HareLynx => {
            // Boreal forest grid + foraging hares & pursuing lynx + spatial X-density histogram
            for g in 1..8 {
                let gx = inner.left() + inner.width() * (g as f32 / 8.0);
                let gy = inner.top() + inner.height() * (g as f32 / 8.0);
                painter.line_segment(
                    [Pos2::new(gx, inner.top()), Pos2::new(gx, inner.bottom())],
                    Stroke::new(0.5, Color32::from_rgb(26, 38, 48)),
                );
                painter.line_segment(
                    [Pos2::new(inner.left(), gy), Pos2::new(inner.right(), gy)],
                    Stroke::new(0.5, Color32::from_rgb(26, 38, 48)),
                );
            }
            for a in &workset.v1_hare_lynx.agents {
                let px = inner.left() + (a.x as f32) * inner.width();
                let py = inner.top() + (a.y as f32) * (inner.height() - 30.0);
                if a.kind == 0 {
                    painter.circle_filled(Pos2::new(px, py), 4.0, Color32::from_rgb(138, 226, 174));
                } else {
                    painter.circle_filled(Pos2::new(px, py), 6.5, Color32::from_rgb(249, 115, 22));
                    painter.circle_stroke(
                        Pos2::new(px, py),
                        9.0,
                        Stroke::new(1.0, Color32::from_rgba_unmultiplied(249, 115, 22, 90)),
                    );
                }
            }
            draw_marginal_distribution_bars(&painter, inner, &workset.v1_hare_lynx.agents);
            ui.label(format!(
                "Legend: Green circles = Snowshoe Hares ({:.1}/km^2) | Orange halos = Canada Lynx ({:.1}/100 km^2) | Bottom bars = Spatial longitudinal density",
                workset.v1_hare_lynx.hares, workset.v1_hare_lynx.lynx
            ));
        }

        ModuleId::V2IsleRoyale => {
            // Lake Superior water backdrop + elongated Isle Royale island contour + pack territories
            painter.rect_filled(inner, CornerRadius::same(4), Color32::from_rgb(12, 34, 58));
            let island = Rect::from_min_max(
                Pos2::new(
                    inner.left() + 0.08 * inner.width(),
                    inner.top() + 0.20 * inner.height(),
                ),
                Pos2::new(
                    inner.left() + 0.92 * inner.width(),
                    inner.top() + 0.80 * inner.height(),
                ),
            );
            let fir_green = (45.0 + 75.0 * (workset.v2_isle_royale.fir_browse / 100.0)) as u8;
            painter.rect_filled(island, CornerRadius::same(28), Color32::from_rgb(26, fir_green, 42));

            // 3 Wolf pack territory rings (West, Central, East)
            for pack in 0..3 {
                let cx = inner.left() + (0.25 + 0.25 * pack as f32) * inner.width();
                let cy = inner.top() + (0.42 + 0.08 * pack as f32) * inner.height();
                painter.circle_stroke(
                    Pos2::new(cx, cy),
                    inner.width() * 0.11,
                    Stroke::new(1.0, Color32::from_rgba_unmultiplied(253, 224, 71, 80)),
                );
            }

            for a in &workset.v2_isle_royale.agents {
                let px = inner.left() + (a.x as f32) * inner.width();
                let py = inner.top() + (a.y as f32) * inner.height();
                if a.kind == 0 {
                    painter.circle_filled(Pos2::new(px, py), 3.6, Color32::from_rgb(214, 180, 132));
                } else {
                    painter.circle_filled(Pos2::new(px, py), 5.8, Color32::from_rgb(248, 113, 113));
                }
            }
            ui.label(format!(
                "Isle Royale Arena: Tan = Moose ({:.0}) | Red = Gray Wolves ({:.1}) across 3 pack territories | Island greenness = Balsam Fir Browse ({:.0}%)",
                workset.v2_isle_royale.moose, workset.v2_isle_royale.wolves, workset.v2_isle_royale.fir_browse
            ));
        }

        ModuleId::V3SoaySheep => {
            // Hirta island grazing pasture tint inversely related to sheep grazing pressure
            let graze_ratio =
                (workset.v3_soay_sheep.population / workset.v3_soay_sheep.carrying_scale_k.max(1.0))
                    .clamp(0.0, 2.0);
            let turf_g = (115.0 - 45.0 * graze_ratio) as u8;
            painter.rect_filled(inner, CornerRadius::same(16), Color32::from_rgb(38, turf_g, 44));
            for a in &workset.v3_soay_sheep.agents {
                let px = inner.left() + (a.x as f32) * inner.width();
                let py = inner.top() + (a.y as f32) * (inner.height() - 28.0);
                painter.circle_filled(Pos2::new(px, py), 4.2, Color32::from_rgb(241, 245, 249));
            }
            draw_marginal_distribution_bars(&painter, inner, &workset.v3_soay_sheep.agents);
            ui.label(format!(
                "Hirta (St. Kilda) Pasture: Census = {:.0} Soay sheep | Last Winter Storm Mortality = {:.1}% | Turf shade reflects Festuca forage standing crop",
                workset.v3_soay_sheep.population, workset.v3_soay_sheep.last_storm_mortality
            ));
        }

        ModuleId::I1TriboliumLpa => {
            // 3-horizon flour vial: Top = Surface Adults, Middle = Burrowing Feeding Larvae, Bottom = Quiescent Pupae
            let h3 = inner.height() / 3.0;
            let top_zone =
                Rect::from_min_max(inner.left_top(), Pos2::new(inner.right(), inner.top() + h3));
            let mid_zone = Rect::from_min_max(
                Pos2::new(inner.left(), inner.top() + h3),
                Pos2::new(inner.right(), inner.top() + 2.0 * h3),
            );
            let bot_zone = Rect::from_min_max(
                Pos2::new(inner.left(), inner.top() + 2.0 * h3),
                inner.right_bottom(),
            );
            painter.rect_filled(top_zone, CornerRadius::same(2), Color32::from_rgb(58, 48, 34));
            painter.rect_filled(mid_zone, CornerRadius::same(2), Color32::from_rgb(48, 39, 28));
            painter.rect_filled(bot_zone, CornerRadius::same(2), Color32::from_rgb(38, 30, 22));

            for (r_zone, lbl) in [
                (top_zone, "Upper Flour Layer: Adults A_t (cannibalize eggs & pupae)"),
                (mid_zone, "Middle Flour Layer: Feeding Larvae L_t (cannibalize eggs)"),
                (bot_zone, "Lower Flour Layer: Quiescent Pupae P_t"),
            ] {
                painter.text(
                    Pos2::new(r_zone.left() + 8.0, r_zone.top() + 10.0),
                    Align2::LEFT_TOP,
                    lbl,
                    FontId::proportional(12.0),
                    Color32::from_rgb(226, 212, 186),
                );
            }

            for a in &workset.i1_tribolium.agents {
                let px = inner.left() + (a.x as f32) * inner.width();
                let py = inner.top() + (a.y as f32) * inner.height();
                let (rad, col) = match a.kind {
                    0 => (3.5, Color32::from_rgb(253, 224, 71)), // Larvae
                    1 => (4.2, Color32::from_rgb(125, 211, 252)), // Pupae
                    _ => (5.0, Color32::from_rgb(217, 119, 6)),  // Adults
                };
                painter.circle_filled(Pos2::new(px, py), rad, col);
            }
            ui.label(format!(
                "Tribolium Stage Structure: Larvae = {:.1} (yellow) | Pupae = {:.1} (blue) | Adults = {:.1} (amber)",
                workset.i1_tribolium.larvae, workset.i1_tribolium.pupae, workset.i1_tribolium.adults
            ));
        }

        ModuleId::I2HuffakerMites => {
            // 6x4 Citrus tray grid showing each orange, Vaseline maze thickness, and local prey/predator densities
            let sim = &workset.i2_huffaker;
            let cw = inner.width() / (sim.cols as f32);
            let ch = inner.height() / (sim.rows as f32);
            let barrier_w = (sim.vaseline_barrier as f32) * 5.0;

            for r in 0..sim.rows {
                for c in 0..sim.cols {
                    let idx = r * sim.cols + c;
                    let cell = Rect::from_min_max(
                        Pos2::new(inner.left() + c as f32 * cw, inner.top() + r as f32 * ch),
                        Pos2::new(
                            inner.left() + (c + 1) as f32 * cw,
                            inner.top() + (r + 1) as f32 * ch,
                        ),
                    );
                    if barrier_w > 0.3 {
                        draw_box_outline(
                            &painter,
                            cell.shrink(1.0),
                            Stroke::new(barrier_w, Color32::from_rgb(148, 163, 184)),
                        );
                    }
                    let prey_frac = (sim.prey_patches[idx] / sim.k_patch.max(1.0)).clamp(0.0, 1.0);
                    let pred_frac = (sim.pred_patches[idx] / 25.0).clamp(0.0, 1.0);
                    let orange_col = Color32::from_rgb(
                        (70.0 + 175.0 * prey_frac) as u8,
                        (55.0 + 95.0 * prey_frac) as u8,
                        32,
                    );
                    let center = cell.center();
                    let rad = cw.min(ch) * 0.34;
                    painter.circle_filled(center, rad, orange_col);
                    if pred_frac > 0.02 {
                        painter.circle_stroke(
                            center,
                            rad + 3.5,
                            Stroke::new(1.5 + 3.5 * pred_frac as f32, Color32::from_rgb(239, 68, 68)),
                        );
                    }
                    painter.text(
                        center,
                        Align2::CENTER_CENTER,
                        format!("{:.0}/{:.0}", sim.prey_patches[idx], sim.pred_patches[idx]),
                        FontId::monospace(10.5),
                        Color32::WHITE,
                    );
                }
            }
            ui.label(format!(
                "Huffaker 6x4 Orange Tray: Bright orange = Eotetranychus prey density | Red halo = Typhlodromus predators | Gray border width = Vaseline barrier ({:.0}%)",
                sim.vaseline_barrier * 100.0
            ));
        }

        ModuleId::I3AphidParasitoid => {
            // 5 vertical host-plant stems with refuge zones at crown/axils
            for s_idx in 0..5 {
                let sx = inner.left() + (0.16 + 0.17 * s_idx as f32) * inner.width();
                painter.line_segment(
                    [
                        Pos2::new(sx, inner.top() + 16.0),
                        Pos2::new(sx, inner.bottom() - 16.0),
                    ],
                    Stroke::new(6.0, Color32::from_rgb(46, 125, 50)),
                );
            }
            for a in &workset.i3_aphid.agents {
                let px = inner.left() + (a.x as f32) * inner.width();
                let py = inner.top() + (a.y as f32) * inner.height();
                match a.kind {
                    0 => {
                        let col = if a.aux > 0.5 {
                            Color32::from_rgb(56, 189, 248) // Protected in stem refuge
                        } else {
                            Color32::from_rgb(163, 230, 53) // Exposed aphid
                        };
                        painter.circle_filled(Pos2::new(px, py), 3.8, col);
                    }
                    1 => {
                        painter.circle_filled(
                            Pos2::new(px, py),
                            4.8,
                            Color32::from_rgb(244, 63, 94),
                        );
                    }
                    _ => {
                        painter.circle_filled(
                            Pos2::new(px, py),
                            4.0,
                            Color32::from_rgb(180, 142, 92),
                        );
                    }
                }
            }
            ui.label(
                "Host Stems: Lime = Exposed Pea Aphids | Cyan = Refuge Aphids | Rose = Adult Aphidius Wasps | Tan = Parasitized Mummies",
            );
        }

        ModuleId::M1EcoliBatchChemostat => {
            // Microscope darkfield view of dividing E. coli rod cells + glucose/lactose turbidity tint
            let glc_tint = (workset.m1_ecoli.glucose_g_l * 6.0).clamp(0.0, 35.0) as u8;
            painter.rect_filled(inner, CornerRadius::same(8), Color32::from_rgb(14, 22 + glc_tint, 30));
            for a in &workset.m1_ecoli.agents {
                let cx = inner.left() + (a.x as f32) * inner.width();
                let cy = inner.top() + (a.y as f32) * inner.height();
                let ang = (a.aux as f32) * std::f32::consts::TAU;
                let dx = 5.5 * ang.cos();
                let dy = 5.5 * ang.sin();
                let col = if a.kind == 0 {
                    Color32::from_rgb(52, 211, 153)
                } else {
                    Color32::from_rgb(100, 116, 139)
                };
                painter.line_segment(
                    [Pos2::new(cx - dx, cy - dy), Pos2::new(cx + dx, cy + dy)],
                    Stroke::new(3.6, col),
                );
            }
            ui.label(format!(
                "Darkfield Microcolony View: Emerald rods = Viable E. coli ({:.2} x 10^7 CFU/mL) | Slate rods = Non-viable cells ({:.2}) | Glucose = {:.2} g/L | Enzyme Acclimation = {:.0}%",
                workset.m1_ecoli.viable_cells,
                workset.m1_ecoli.dead_cells,
                workset.m1_ecoli.glucose_g_l,
                workset.m1_ecoli.lag_enzyme_state * 100.0
            ));
        }

        ModuleId::M2GauseProtists => {
            let sim = &workset.m2_gause;
            if !sim.competition_mode && sim.sediment_refuge > 0.01 {
                let ref_top = inner.top() + inner.height() * 0.77;
                let ref_rect =
                    Rect::from_min_max(Pos2::new(inner.left(), ref_top), inner.right_bottom());
                painter.rect_filled(ref_rect, CornerRadius::same(2), Color32::from_rgb(48, 38, 24));
                painter.text(
                    Pos2::new(ref_rect.left() + 8.0, ref_rect.top() + 6.0),
                    Align2::LEFT_TOP,
                    "Oat-Medium Sediment Refuge (protects Paramecium caudatum from Didinium nasutum)",
                    FontId::proportional(12.0),
                    Color32::from_rgb(234, 210, 168),
                );
            }
            for a in &sim.agents {
                let px = inner.left() + (a.x as f32) * inner.width();
                let py = inner.top() + (a.y as f32) * inner.height();
                if a.kind == 0 {
                    painter.circle_filled(Pos2::new(px, py), 4.6, Color32::from_rgb(96, 165, 250));
                } else if sim.competition_mode {
                    painter.circle_filled(Pos2::new(px, py), 3.6, Color32::from_rgb(52, 211, 153));
                } else {
                    painter.circle_filled(Pos2::new(px, py), 5.8, Color32::from_rgb(244, 63, 94));
                }
            }
            if sim.competition_mode {
                ui.label(format!(
                    "Gause Competition Microcosm: Blue = P. caudatum ({:.1}) | Green = P. aurelia ({:.1})",
                    sim.p_caudatum, sim.p_aurelia
                ));
            } else {
                ui.label(format!(
                    "Gause Predation Microcosm: Blue = P. caudatum ({:.1}) | Crimson = Didinium nasutum ({:.1})",
                    sim.p_caudatum, sim.didinium
                ));
            }
        }

        ModuleId::M3PhagePlaque => {
            // 24x24 Petri-dish agar bacterial lawn showing radial T4 lytic plaques
            let sim = &workset.m3_phage;
            let side = sim.grid_side.max(1);
            let cw = inner.width() / (side as f32);
            let ch = inner.height() / (side as f32);
            for r in 0..side {
                for c in 0..side {
                    let idx = r * side + c;
                    let lawn = sim.agar_lawn[idx].clamp(0.0, 1.0);
                    let vir = (sim.agar_phage[idx] / 25.0).clamp(0.0, 1.0);
                    let cell_col = Color32::from_rgb(
                        (24.0 + 185.0 * vir) as u8,
                        (30.0 + 165.0 * lawn) as u8,
                        (42.0 + 60.0 * lawn) as u8,
                    );
                    let cell = Rect::from_min_max(
                        Pos2::new(inner.left() + c as f32 * cw, inner.top() + r as f32 * ch),
                        Pos2::new(
                            inner.left() + (c + 1) as f32 * cw - 0.8,
                            inner.top() + (r + 1) as f32 * ch - 0.8,
                        ),
                    );
                    painter.rect_filled(cell, CornerRadius::same(2), cell_col);
                }
            }
            ui.label(format!(
                "Agar Lawn Plaque Assay: Turbid green = E. coli lawn | Dark/amber lytic rings = Expanding T4 Bacteriophage plaques | Resistant R = {:.2}",
                sim.resistant_bact
            ));
        }

        ModuleId::T1LevinsMetapop => {
            let sim = &workset.t1_levins;
            for i in 0..sim.patches.len() {
                for j in (i + 1)..sim.patches.len() {
                    let dx = sim.patches[i].x - sim.patches[j].x;
                    let dy = sim.patches[i].y - sim.patches[j].y;
                    if dx * dx + dy * dy < 0.032 && !sim.patches[i].destroyed && !sim.patches[j].destroyed {
                        let p1 = Pos2::new(
                            inner.left() + sim.patches[i].x as f32 * inner.width(),
                            inner.top() + sim.patches[i].y as f32 * inner.height(),
                        );
                        let p2 = Pos2::new(
                            inner.left() + sim.patches[j].x as f32 * inner.width(),
                            inner.top() + sim.patches[j].y as f32 * inner.height(),
                        );
                        painter.line_segment(
                            [p1, p2],
                            Stroke::new(0.8, Color32::from_rgba_unmultiplied(100, 116, 139, 80)),
                        );
                    }
                }
            }
            for patch in &sim.patches {
                let center = Pos2::new(
                    inner.left() + patch.x as f32 * inner.width(),
                    inner.top() + patch.y as f32 * inner.height(),
                );
                if patch.destroyed {
                    painter.circle_stroke(
                        center,
                        7.0,
                        Stroke::new(1.5, Color32::from_rgb(127, 29, 29)),
                    );
                } else if patch.occupied {
                    painter.circle_filled(center, 8.0, Color32::from_rgb(52, 211, 153));
                } else {
                    painter.circle_stroke(
                        center,
                        7.5,
                        Stroke::new(1.5, Color32::from_rgb(148, 163, 184)),
                    );
                }
            }
            ui.label(format!(
                "64-Patch Landscape: Solid Green = Occupied | Hollow Gray = Empty Suitable | Dark Red = Destroyed Habitat (h = {:.1}%, critical h_crit = {:.1}%)",
                sim.destruction_h * 100.0,
                sim.critical_destruction_threshold() * 100.0
            ));
        }

        ModuleId::T2MayBifurcation => {
            // Live Cobweb Diagram (N_{t+1} vs N_t) side-by-side with density bar
            let sim = &workset.t2_may;
            let max_axis = (sim.carrying_k * if sim.ricker_map { 2.4 } else { 1.45 }).max(10.0);
            let to_screen = |nx: f64, ny: f64| -> Pos2 {
                let ux = (nx / max_axis).clamp(0.0, 1.0) as f32;
                let uy = (ny / max_axis).clamp(0.0, 1.0) as f32;
                Pos2::new(
                    inner.left() + ux * inner.width(),
                    inner.bottom() - uy * inner.height(),
                )
            };

            // Diagonal y = x line
            painter.line_segment(
                [to_screen(0.0, 0.0), to_screen(max_axis, max_axis)],
                Stroke::new(1.2, Color32::from_rgb(100, 116, 139)),
            );

            // Return map curve N_{t+1} = f(N_t)
            let mut prev_pt = to_screen(0.0, sim.map_eval(0.0, sim.r_param));
            for s in 1..=120 {
                let x_val = max_axis * (s as f64 / 120.0);
                let y_val = sim.map_eval(x_val, sim.r_param);
                let cur_pt = to_screen(x_val, y_val);
                painter.line_segment(
                    [prev_pt, cur_pt],
                    Stroke::new(2.2, Color32::from_rgb(56, 189, 248)),
                );
                prev_pt = cur_pt;
            }

            // Cobweb stair-steps
            for pair in sim.cobweb_points.windows(2) {
                let p1 = to_screen(pair[0].0, pair[0].1);
                let p2 = to_screen(pair[1].0, pair[1].1);
                painter.line_segment(
                    [p1, p2],
                    Stroke::new(1.3, Color32::from_rgba_unmultiplied(250, 204, 21, 185)),
                );
            }

            ui.label(format!(
                "Cobweb Map N_{{t+1}} vs N_t: Cyan = Return Curve | Gold = Cobweb Orbit | Current N_t = {:.2} | Lyapunov lambda = {:.4}",
                sim.n_current,
                sim.lyapunov_exponent()
            ));
        }

        ModuleId::T3WrightFisher => {
            // De Finetti Hardy-Weinberg Ternary Triangle showing all demes + HWE parabola
            let sim = &workset.t3_wright_fisher;
            let base_left = Pos2::new(inner.left() + inner.width() * 0.14, inner.bottom() - 24.0);
            let base_right = Pos2::new(inner.left() + inner.width() * 0.86, inner.bottom() - 24.0);
            let apex = Pos2::new(inner.center().x, inner.top() + 26.0);

            painter.line_segment([base_left, base_right], Stroke::new(1.8, Color32::GRAY));
            painter.line_segment([base_right, apex], Stroke::new(1.8, Color32::GRAY));
            painter.line_segment([apex, base_left], Stroke::new(1.8, Color32::GRAY));

            painter.text(base_left, Align2::RIGHT_TOP, "aa (p=0)", FontId::proportional(12.5), Color32::WHITE);
            painter.text(base_right, Align2::LEFT_TOP, "AA (p=1)", FontId::proportional(12.5), Color32::WHITE);
            painter.text(apex, Align2::CENTER_BOTTOM, "Aa Heterozygote (H=1)", FontId::proportional(12.5), Color32::WHITE);

            let ternary_pt = |p_freq: f64, het: f64| -> Pos2 {
                let x = base_left.x + (p_freq as f32) * (base_right.x - base_left.x);
                let y = base_left.y - (het as f32) * (base_left.y - apex.y);
                Pos2::new(x, y)
            };

            // HWE parabola H = 2p(1-p)
            let mut prev_hwe = ternary_pt(0.0, 0.0);
            for i in 1..=80 {
                let p = i as f64 / 80.0;
                let cur = ternary_pt(p, 2.0 * p * (1.0 - p));
                painter.line_segment(
                    [prev_hwe, cur],
                    Stroke::new(1.8, Color32::from_rgb(96, 165, 250)),
                );
                prev_hwe = cur;
            }

            let palette = [
                Color32::from_rgb(52, 211, 153),
                Color32::from_rgb(251, 146, 60),
                Color32::from_rgb(244, 114, 182),
                Color32::from_rgb(250, 204, 21),
                Color32::from_rgb(167, 139, 250),
                Color32::from_rgb(45, 212, 191),
            ];
            for (idx, &p) in sim.deme_freqs.iter().enumerate() {
                let het = 2.0 * p * (1.0 - p);
                let pt = ternary_pt(p, het);
                painter.circle_filled(pt, 6.5, palette[idx % palette.len()]);
            }

            let (p_bar, h_s, h_t, fst) = sim.summary_popgen();
            ui.label(format!(
                "De Finetti Ternary Diagram: Blue Curve = HWE 2p(1-p) | Colored Dots = {} Diploid Demes | p_bar = {:.3}, H_S = {:.3}, H_T = {:.3}, F_ST = {:.3}",
                sim.deme_freqs.len(), p_bar, h_s, h_t, fst
            ));
        }

        ModuleId::T4KingmanCoalescent => {
            // Render exact backward-in-time Kingman Coalescent Tree with infinite-sites mutations on branches
            let sim = &workset.t4_coalescent;
            let max_t = sim.t_mrca.max(1.0);
            let node_pos = |id: usize| -> Pos2 {
                let nd = &sim.nodes[id];
                let x = inner.left() + (nd.x_pos as f32) * inner.width();
                let y = inner.bottom() - 22.0 - ((nd.time_gen / max_t) as f32) * (inner.height() - 44.0);
                Pos2::new(x, y)
            };

            for br in &sim.branches {
                let pc = node_pos(br.child_id);
                let pp = node_pos(br.parent_id);
                let corner = Pos2::new(pc.x, pp.y);
                let stroke = Stroke::new(2.0, Color32::from_rgb(148, 163, 184));
                painter.line_segment([pc, corner], stroke);
                painter.line_segment([corner, pp], stroke);

                // Plot infinite-sites mutations along the vertical lineage branch
                for m in 0..br.mutations.min(12) {
                    let frac = (m as f32 + 1.0) / (br.mutations.min(12) as f32 + 1.0);
                    let my = pc.y + frac * (corner.y - pc.y);
                    painter.circle_filled(
                        Pos2::new(pc.x, my),
                        3.5,
                        Color32::from_rgb(250, 204, 21),
                    );
                }
            }
            for nd in &sim.nodes {
                let p = node_pos(nd.id);
                let col = if nd.left_child.is_none() {
                    Color32::from_rgb(56, 189, 248)
                } else {
                    Color32::from_rgb(244, 114, 182)
                };
                painter.circle_filled(p, 4.2, col);
            }
            ui.label(format!(
                "Kingman Coalescent Genealogy (n = {}): T_MRCA = {:.1} gen (Expected = {:.1}) | Gold ticks = Infinite-sites mutations",
                sim.sample_size_n,
                sim.t_mrca,
                sim.expected_tmrca()
            ));
        }

        ModuleId::T5WrightLandscape => {
            // Interactive mouse-drag 3D orbit camera
            let sim = &mut workset.t5_landscape;
            if response.dragged() {
                let delta = response.drag_delta();
                sim.camera_azimuth_deg = (sim.camera_azimuth_deg + (delta.x as f64) * 0.45).rem_euclid(360.0);
                sim.camera_elevation_deg = (sim.camera_elevation_deg + (delta.y as f64) * 0.35).clamp(12.0, 72.0);
            }

            let az = sim.camera_azimuth_deg.to_radians();
            let el = sim.camera_elevation_deg.to_radians();
            let (sin_az, cos_az) = az.sin_cos();
            let (sin_el, cos_el) = el.sin_cos();
            let cx_scr = inner.center().x;
            let cy_scr = inner.center().y + inner.height() * 0.14;
            let scale = inner.width().min(inner.height()) * 0.62;

            let project = |ux: f64, uy: f64, wz: f64| -> (Pos2, f64) {
                let x = ux - 0.5;
                let y = uy - 0.5;
                let z = (wz - 0.45) * 0.52;
                let rx = x * cos_az - y * sin_az;
                let ry = x * sin_az + y * cos_az;
                let depth = ry * cos_el + z * sin_el;
                let sy = ry * sin_el - z * cos_el;
                let persp = 1.0 / (1.0 + 0.22 * depth).max(0.4);
                (
                    Pos2::new(
                        cx_scr + (rx * persp) as f32 * scale,
                        cy_scr + (sy * persp) as f32 * scale,
                    ),
                    depth,
                )
            };

            // Build depth-sorted 22x22 surface quads (Painter's algorithm back-to-front)
            let g_n = 22usize;
            let mut quads: Vec<(f64, [Pos2; 4], Color32)> = Vec::with_capacity(g_n * g_n);
            for i in 0..g_n {
                let x0 = i as f64 / g_n as f64;
                let x1 = (i + 1) as f64 / g_n as f64;
                for j in 0..g_n {
                    let y0 = j as f64 / g_n as f64;
                    let y1 = (j + 1) as f64 / g_n as f64;
                    let w00 = sim.fitness_at(x0, y0);
                    let w10 = sim.fitness_at(x1, y0);
                    let w11 = sim.fitness_at(x1, y1);
                    let w01 = sim.fitness_at(x0, y1);
                    let (p00, d00) = project(x0, y0, w00);
                    let (p10, d10) = project(x1, y0, w10);
                    let (p11, d11) = project(x1, y1, w11);
                    let (p01, d01) = project(x0, y1, w01);
                    let mean_d = 0.25 * (d00 + d10 + d11 + d01);
                    let mean_w = 0.25 * (w00 + w10 + w11 + w01);
                    let t = ((mean_w - 0.45) / 1.15).clamp(0.0, 1.0);
                    let fill = Color32::from_rgb(
                        (22.0 + 215.0 * t.powf(1.3)) as u8,
                        (65.0 + 145.0 * (1.0 - (t - 0.55).abs() * 1.4).clamp(0.0, 1.0)) as u8,
                        (135.0 - 90.0 * t) as u8,
                    );
                    quads.push((mean_d, [p00, p10, p11, p01], fill));
                }
            }
            // Sort farthest depth first
            quads.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
            let grid_stroke = Stroke::new(0.6, Color32::from_rgba_unmultiplied(15, 23, 42, 130));
            for (_, pts, fill) in quads {
                painter.add(egui::Shape::convex_polygon(pts.to_vec(), fill, grid_stroke));
            }

            // Draw demes, vertical fitness drop-pins & 3D historical valley-crossing trails
            let palette = [
                Color32::from_rgb(248, 250, 252),
                Color32::from_rgb(56, 189, 248),
                Color32::from_rgb(244, 114, 182),
                Color32::from_rgb(163, 230, 53),
                Color32::from_rgb(251, 146, 60),
                Color32::from_rgb(192, 132, 252),
            ];
            for (idx, d) in sim.demes.iter().enumerate() {
                let col = palette[idx % palette.len()];
                for pair in d.trail.windows(2) {
                    let (pa, _) = project(pair[0][0], pair[0][1], pair[0][2] + 0.02);
                    let (pb, _) = project(pair[1][0], pair[1][1], pair[1][2] + 0.02);
                    painter.line_segment([pa, pb], Stroke::new(1.5, col));
                }
                let (p_base, _) = project(d.x, d.y, 0.45);
                let (p_top, _) = project(d.x, d.y, d.fitness + 0.05);
                painter.line_segment(
                    [p_base, p_top],
                    Stroke::new(1.1, Color32::from_rgba_unmultiplied(255, 255, 255, 110)),
                );
                painter.circle_filled(p_top, 5.8, col);
                painter.circle_stroke(p_top, 5.8, Stroke::new(1.4, Color32::BLACK));
            }

            let (p_pk1, _) = project(0.24, 0.26, sim.fitness_at(0.24, 0.26) + 0.12);
            let (p_pk2, _) = project(0.76, 0.74, sim.fitness_at(0.76, 0.74) + 0.12);
            painter.text(p_pk1, Align2::CENTER_BOTTOM, "Peak 1 (Local)", FontId::proportional(11.5), Color32::WHITE);
            painter.text(p_pk2, Align2::CENTER_BOTTOM, "Peak 2 (Global Optimum)", FontId::proportional(12.0), Color32::from_rgb(254, 240, 138));

            ui.label(format!(
                "Drag mouse on 3D surface to rotate (Azimuth = {:.0} deg, Elevation = {:.0} deg) | {} Demes (N_e = {:.0}) evolving via Drift + Steepest Ascent + Fitness-Weighted Dispersion",
                sim.camera_azimuth_deg,
                sim.camera_elevation_deg,
                sim.demes.len(),
                sim.effective_ne
            ));
        }

        ModuleId::T6LeslieLefkovitch => {
            // Left half: Current Stage Proportion vs. Stable Distribution w_i + Fisher Reproductive Value v_i
            // Right half: k x k Elasticity Matrix E_ij (%) & Transition Matrix A_ij Heatmap
            let sim = &workset.t6_matrix;
            let k = sim.stage_names.len().max(1);
            let mid_x = inner.left() + inner.width() * 0.48;

            // Left Panel: Stage Pyramid Bars
            let left_rect = Rect::from_min_max(inner.left_top(), Pos2::new(mid_x - 12.0, inner.bottom()));
            painter.text(
                Pos2::new(left_rect.left(), left_rect.top()),
                Align2::LEFT_TOP,
                "Current Cohort n_i/N (Cyan) vs. Stable w_i (Gold) & Reprod. Value v_i",
                FontId::proportional(12.0),
                Color32::WHITE,
            );
            let tot_n: f64 = sim.abundances.iter().sum::<f64>().max(1e-9);
            let row_h = ((left_rect.height() - 28.0) / (k as f32)).min(48.0);
            let max_v = sim
                .eigen
                .reproductive_value_v
                .iter()
                .cloned()
                .fold(1.0f64, f64::max);

            for i in 0..k {
                let y_top = left_rect.top() + 26.0 + (i as f32) * row_h;
                let cur_prop = (sim.abundances[i] / tot_n).clamp(0.0, 1.0) as f32;
                let w_prop = sim
                    .eigen
                    .stable_distribution_w
                    .get(i)
                    .copied()
                    .unwrap_or(0.0)
                    .clamp(0.0, 1.0) as f32;
                let v_val = sim
                    .eigen
                    .reproductive_value_v
                    .get(i)
                    .copied()
                    .unwrap_or(0.0);

                let bar_x0 = left_rect.left() + 135.0;
                let max_bw = (left_rect.right() - bar_x0 - 52.0).max(40.0);

                painter.text(
                    Pos2::new(left_rect.left(), y_top + row_h * 0.42),
                    Align2::LEFT_CENTER,
                    &sim.stage_names[i],
                    FontId::proportional(11.0),
                    Color32::from_rgb(203, 213, 225),
                );
                // Current proportion bar
                painter.rect_filled(
                    Rect::from_min_max(
                        Pos2::new(bar_x0, y_top + 3.0),
                        Pos2::new(bar_x0 + max_bw * cur_prop, y_top + row_h * 0.44),
                    ),
                    CornerRadius::same(2),
                    Color32::from_rgb(56, 189, 248),
                );
                // Analytical stable w_i bar
                painter.rect_filled(
                    Rect::from_min_max(
                        Pos2::new(bar_x0, y_top + row_h * 0.48),
                        Pos2::new(bar_x0 + max_bw * w_prop, y_top + row_h * 0.84),
                    ),
                    CornerRadius::same(2),
                    Color32::from_rgb(250, 204, 21),
                );
                painter.text(
                    Pos2::new(left_rect.right(), y_top + row_h * 0.42),
                    Align2::RIGHT_CENTER,
                    format!("v={:.1}", v_val.min(max_v)),
                    FontId::monospace(10.5),
                    Color32::from_rgb(244, 114, 182),
                );
            }

            // Right Panel: k x k Matrix & Elasticity Heatmap
            let right_rect = Rect::from_min_max(Pos2::new(mid_x + 8.0, inner.top()), inner.right_bottom());
            painter.text(
                Pos2::new(right_rect.left(), right_rect.top()),
                Align2::LEFT_TOP,
                format!(
                    "Projection Matrix A_ij & Elasticity E_ij (lambda_1 = {:.4}, r = {:+.4}, R_0 = {:.2})",
                    sim.eigen.dominant_lambda, sim.eigen.intrinsic_r, sim.eigen.net_reproductive_r0
                ),
                FontId::proportional(12.0),
                Color32::WHITE,
            );
            let mat = sim.active_matrix();
            let grid_top = right_rect.top() + 26.0;
            let cw = (right_rect.width() / (k as f32)).min(64.0);
            let ch = ((right_rect.height() - 30.0) / (k as f32)).min(52.0);
            let max_e = sim
                .eigen
                .elasticity
                .iter()
                .flatten()
                .cloned()
                .fold(0.01f64, f64::max);

            for r in 0..k {
                for c in 0..k {
                    let cell = Rect::from_min_max(
                        Pos2::new(right_rect.left() + (c as f32) * cw, grid_top + (r as f32) * ch),
                        Pos2::new(
                            right_rect.left() + ((c + 1) as f32) * cw - 2.0,
                            grid_top + ((r + 1) as f32) * ch - 2.0,
                        ),
                    );
                    let a_val = mat[r][c];
                    let e_val = sim.eigen.elasticity[r][c];
                    let intensity = (e_val / max_e).clamp(0.0, 1.0);
                    let bg = if a_val <= 0.0 {
                        Color32::from_rgb(18, 25, 38)
                    } else {
                        Color32::from_rgb(
                            (25.0 + 195.0 * intensity) as u8,
                            (55.0 + 115.0 * (1.0 - intensity * 0.4)) as u8,
                            (85.0 - 45.0 * intensity) as u8,
                        )
                    };
                    painter.rect_filled(cell, CornerRadius::same(3), bg);
                    if a_val > 0.0 {
                        painter.text(
                            Pos2::new(cell.center().x, cell.center().y - 5.0),
                            Align2::CENTER_CENTER,
                            if a_val >= 10.0 {
                                format!("{:.1}", a_val)
                            } else {
                                format!("{:.3}", a_val)
                            },
                            FontId::monospace(9.5),
                            Color32::WHITE,
                        );
                        painter.text(
                            Pos2::new(cell.center().x, cell.center().y + 7.0),
                            Align2::CENTER_CENTER,
                            format!("E:{:.0}%", e_val * 100.0),
                            FontId::monospace(9.0),
                            Color32::from_rgb(254, 240, 138),
                        );
                    }
                }
            }

            ui.label(format!(
                "Dominant Eigenvalue lambda_1 = {:.4} (Realized N_{{t+1}}/N_t = {:.4}, Damping Ratio rho = {:.3}) | Bright amber cells mark highest conservation elasticity E_ij",
                sim.eigen.dominant_lambda, sim.realized_lambda, sim.eigen.damping_ratio
            ));
        }

        ModuleId::E1SequencePopGen
        | ModuleId::E2HaplotypeNetwork
        | ModuleId::E3BayesianSkyline
        | ModuleId::E4SlidingWindowScan => {}
    }
}
