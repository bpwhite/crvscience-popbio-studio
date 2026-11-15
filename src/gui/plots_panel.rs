use egui::{
    Align2, Color32, CornerRadius, FontId, Painter, Pos2, Rect, RichText, Sense, Stroke, Ui, Vec2,
};

use crate::model::ModuleId;
use crate::workset::WorksetState;

pub struct PlotSeries {
    pub name: String,
    pub color: Color32,
    pub points: Vec<[f64; 2]>,
    pub width: f32,
    pub as_dots: bool,
    pub dot_radius: f32,
}

pub struct PlotBand {
    pub name: String,
    pub fill: Color32,
    /// Ordered `[x, y_low, y_high]` triples for credible intervals / envelopes
    pub bounds: Vec<[f64; 3]>,
}

pub struct PlotBar {
    pub x: f64,
    pub y: f64,
    pub width: f64,
    pub color: Color32,
}

pub struct PlotVMarker {
    pub x: f64,
    pub color: Color32,
}

pub struct PlotHMarker {
    pub y: f64,
    pub color: Color32,
}

pub fn draw_scientific_plot(
    ui: &mut Ui,
    id_salt: &str,
    height: f32,
    series: &[PlotSeries],
    bars: &[PlotBar],
    vmarkers: &[PlotVMarker],
) {
    draw_scientific_plot_ex(ui, id_salt, height, &[], series, bars, vmarkers, &[]);
}

/// Extended zero-dependency native `egui::Painter` scientific plot canvas supporting shaded
/// credible-interval polygon bands (`PlotBand`), horizontal significance thresholds (`PlotHMarker`),
/// vertical phase lines, bars, scatter/line layers, and live hover crosshair readouts.
pub fn draw_scientific_plot_ex(
    ui: &mut Ui,
    id_salt: &str,
    height: f32,
    bands: &[PlotBand],
    series: &[PlotSeries],
    bars: &[PlotBar],
    vmarkers: &[PlotVMarker],
    hmarkers: &[PlotHMarker],
) {
    ui.horizontal_wrapped(|ui| {
        for b in bands {
            if !b.name.is_empty() {
                ui.colored_label(b.fill, format!("▒ {}", b.name));
                ui.add_space(6.0);
            }
        }
        for s in series {
            if !s.name.is_empty() {
                ui.colored_label(s.color, format!("■ {}", s.name));
                ui.add_space(6.0);
            }
        }
    });

    let avail_w = ui.available_width().max(320.0);
    let (response, painter) =
        ui.allocate_painter(Vec2::new(avail_w, height.max(135.0)), Sense::hover());
    let outer = response.rect;
    painter.rect_filled(outer, CornerRadius::same(6), Color32::from_rgb(12, 18, 28));
    draw_outline(&painter, outer, Stroke::new(1.0, Color32::from_rgb(52, 68, 92)));

    let plot_rect = Rect::from_min_max(
        Pos2::new(outer.left() + 62.0, outer.top() + 12.0),
        Pos2::new(outer.right() - 14.0, outer.bottom() - 26.0),
    );
    if plot_rect.width() <= 10.0 || plot_rect.height() <= 10.0 {
        return;
    }

    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;

    for bnd in bands {
        for &[x, yl, yh] in &bnd.bounds {
            if x.is_finite() && yl.is_finite() && yh.is_finite() {
                x_min = x_min.min(x);
                x_max = x_max.max(x);
                y_min = y_min.min(yl);
                y_max = y_max.max(yh);
            }
        }
    }
    for s in series {
        for &[x, y] in &s.points {
            if x.is_finite() && y.is_finite() {
                x_min = x_min.min(x);
                x_max = x_max.max(x);
                y_min = y_min.min(y);
                y_max = y_max.max(y);
            }
        }
    }
    for b in bars {
        x_min = x_min.min(b.x - b.width);
        x_max = x_max.max(b.x + b.width);
        y_min = y_min.min(0.0);
        y_max = y_max.max(b.y);
    }
    for hm in hmarkers {
        y_min = y_min.min(hm.y);
        y_max = y_max.max(hm.y);
    }

    if !x_min.is_finite() || !y_min.is_finite() {
        let _ = id_salt;
        return;
    }
    if (x_max - x_min).abs() < 1e-9 {
        x_max = x_min + 1.0;
    }
    if (y_max - y_min).abs() < 1e-9 {
        y_max = y_min + 1.0;
    } else {
        let pad_y = 0.07 * (y_max - y_min);
        y_max += pad_y;
        if y_min >= 0.0 && y_min < 0.15 * y_max {
            y_min = 0.0;
        } else {
            y_min -= pad_y;
        }
    }

    let to_screen = |x: f64, y: f64| -> Pos2 {
        let ux = ((x - x_min) / (x_max - x_min)).clamp(0.0, 1.0) as f32;
        let uy = ((y - y_min) / (y_max - y_min)).clamp(0.0, 1.0) as f32;
        Pos2::new(
            plot_rect.left() + ux * plot_rect.width(),
            plot_rect.bottom() - uy * plot_rect.height(),
        )
    };

    for tick in 0..=4 {
        let frac = tick as f64 / 4.0;
        let xv = x_min + frac * (x_max - x_min);
        let yv = y_min + frac * (y_max - y_min);
        let sx = plot_rect.left() + (frac as f32) * plot_rect.width();
        let sy = plot_rect.bottom() - (frac as f32) * plot_rect.height();

        painter.line_segment(
            [Pos2::new(sx, plot_rect.top()), Pos2::new(sx, plot_rect.bottom())],
            Stroke::new(0.6, Color32::from_rgb(28, 38, 54)),
        );
        painter.line_segment(
            [Pos2::new(plot_rect.left(), sy), Pos2::new(plot_rect.right(), sy)],
            Stroke::new(0.6, Color32::from_rgb(28, 38, 54)),
        );

        painter.text(
            Pos2::new(sx, plot_rect.bottom() + 5.0),
            Align2::CENTER_TOP,
            format!("{:.1}", xv),
            FontId::monospace(10.0),
            Color32::from_rgb(148, 163, 184),
        );
        painter.text(
            Pos2::new(plot_rect.left() - 5.0, sy),
            Align2::RIGHT_CENTER,
            format!("{:.2}", yv),
            FontId::monospace(10.0),
            Color32::from_rgb(148, 163, 184),
        );
    }

    // Shaded credible interval ribbons (adjacent convex quads)
    for bnd in bands {
        for pair in bnd.bounds.windows(2) {
            let [x0, l0, h0] = pair[0];
            let [x1, l1, h1] = pair[1];
            let p_tl = to_screen(x0, h0.max(l0));
            let p_tr = to_screen(x1, h1.max(l1));
            let p_br = to_screen(x1, l1.min(h1));
            let p_bl = to_screen(x0, l0.min(h0));
            painter.add(egui::Shape::convex_polygon(
                vec![p_tl, p_tr, p_br, p_bl],
                bnd.fill,
                Stroke::NONE,
            ));
        }
    }

    // Horizontal threshold markers
    for hm in hmarkers {
        if hm.y >= y_min && hm.y <= y_max {
            let pl = to_screen(x_min, hm.y);
            let pr = to_screen(x_max, hm.y);
            painter.line_segment([pl, pr], Stroke::new(1.1, hm.color));
        }
    }

    // Vertical phase markers
    for vm in vmarkers {
        if vm.x >= x_min && vm.x <= x_max {
            let p_top = to_screen(vm.x, y_max);
            let p_bot = to_screen(vm.x, y_min);
            painter.line_segment([p_top, p_bot], Stroke::new(1.2, vm.color));
        }
    }

    // Histogram bars
    for b in bars {
        let p_tl = to_screen(b.x - 0.5 * b.width, b.y);
        let p_br = to_screen(b.x + 0.5 * b.width, 0.0);
        let r = Rect::from_min_max(
            Pos2::new(p_tl.x.min(p_br.x), p_tl.y.min(p_br.y)),
            Pos2::new(p_tl.x.max(p_br.x), p_tl.y.max(p_br.y)),
        );
        painter.rect_filled(r, CornerRadius::same(2), b.color);
    }

    // Line & scatter series
    for s in series {
        if s.points.is_empty() {
            continue;
        }
        if s.as_dots {
            let stride = (s.points.len() / 2200).max(1);
            for pt in s.points.iter().step_by(stride) {
                let p = to_screen(pt[0], pt[1]);
                painter.circle_filled(p, s.dot_radius, s.color);
            }
        } else {
            let stride = (s.points.len() / 700).max(1);
            let mut prev: Option<Pos2> = None;
            for pt in s.points.iter().step_by(stride) {
                let cur = to_screen(pt[0], pt[1]);
                if let Some(p0) = prev {
                    painter.line_segment([p0, cur], Stroke::new(s.width, s.color));
                }
                prev = Some(cur);
            }
            if let Some(&last_pt) = s.points.last() {
                let cur = to_screen(last_pt[0], last_pt[1]);
                if let Some(p0) = prev {
                    painter.line_segment([p0, cur], Stroke::new(s.width, s.color));
                }
            }
        }
    }

    if let Some(hp) = response.hover_pos() {
        if plot_rect.contains(hp) {
            painter.line_segment(
                [Pos2::new(hp.x, plot_rect.top()), Pos2::new(hp.x, plot_rect.bottom())],
                Stroke::new(0.8, Color32::from_rgba_unmultiplied(226, 232, 240, 110)),
            );
            painter.line_segment(
                [Pos2::new(plot_rect.left(), hp.y), Pos2::new(plot_rect.right(), hp.y)],
                Stroke::new(0.8, Color32::from_rgba_unmultiplied(226, 232, 240, 110)),
            );
            let hx = x_min + ((hp.x - plot_rect.left()) / plot_rect.width()) as f64 * (x_max - x_min);
            let hy = y_max - ((hp.y - plot_rect.top()) / plot_rect.height()) as f64 * (y_max - y_min);
            painter.text(
                Pos2::new(plot_rect.right() - 6.0, plot_rect.top() + 6.0),
                Align2::RIGHT_TOP,
                format!("x = {:.3}, y = {:.3}", hx, hy),
                FontId::monospace(11.0),
                Color32::from_rgb(248, 250, 252),
            );
        }
    }
}

fn draw_outline(painter: &Painter, rect: Rect, stroke: Stroke) {
    painter.line_segment([rect.left_top(), rect.right_top()], stroke);
    painter.line_segment([rect.right_top(), rect.right_bottom()], stroke);
    painter.line_segment([rect.right_bottom(), rect.left_bottom()], stroke);
    painter.line_segment([rect.left_bottom(), rect.left_top()], stroke);
}

pub fn render_time_series_tab(ui: &mut Ui, mod_id: ModuleId, workset: &WorksetState) {
    let (records, transitions, names) = workset.history_and_names(mod_id);
    if records.is_empty() {
        ui.label("No time-series points recorded yet.");
        return;
    }

    let colors = [
        Color32::from_rgb(74, 222, 128),
        Color32::from_rgb(251, 146, 60),
        Color32::from_rgb(56, 189, 248),
        Color32::from_rgb(244, 114, 182),
        Color32::from_rgb(250, 204, 21),
    ];

    ui.label(
        RichText::new(
            "Population Trajectory & Growth Phase Transitions (vertical markers indicate biological phase shifts)",
        )
        .strong(),
    );

    let mut series_list = Vec::new();
    for (s_idx, &name) in names.iter().enumerate() {
        if name.is_empty() {
            continue;
        }
        let pts: Vec<[f64; 2]> = records.iter().map(|r| [r.t, r.values[s_idx]]).collect();
        series_list.push(PlotSeries {
            name: name.to_string(),
            color: colors[s_idx % colors.len()],
            points: pts,
            width: 2.0,
            as_dots: false,
            dot_radius: 0.0,
        });
    }

    if mod_id == ModuleId::T3WrightFisher {
        for (d_idx, hist) in workset.t3_wright_fisher.deme_histories.iter().enumerate() {
            let pts: Vec<[f64; 2]> = hist
                .iter()
                .enumerate()
                .map(|(g, &p)| [g as f64, p])
                .collect();
            series_list.push(PlotSeries {
                name: format!("Deme {}", d_idx + 1),
                color: colors[(d_idx + 2) % colors.len()],
                points: pts,
                width: 1.1,
                as_dots: false,
                dot_radius: 0.0,
            });
        }
    }

    let vmarkers: Vec<PlotVMarker> = transitions
        .iter()
        .rev()
        .take(14)
        .map(|tr| {
            let [r, g, b] = tr.phase.badge_rgb();
            PlotVMarker {
                x: tr.t,
                color: Color32::from_rgba_unmultiplied(r, g, b, 130),
            }
        })
        .collect();

    let main_h = (ui.available_height() * 0.56).max(210.0);
    draw_scientific_plot(ui, "main_ts", main_h, &series_list, &[], &vmarkers);

    ui.add_space(6.0);
    ui.label(
        RichText::new("Per-Capita Growth Rate (1/N dN/dt) vs. Primary Abundance N (Density Dependence)")
            .strong(),
    );
    let dd_pts: Vec<[f64; 2]> = records
        .iter()
        .map(|r| [r.values[0], r.per_capita_growth])
        .collect();
    let dd_series = [PlotSeries {
        name: "Per-capita growth rate vs N".to_string(),
        color: Color32::from_rgb(167, 139, 250),
        points: dd_pts,
        width: 1.8,
        as_dots: false,
        dot_radius: 0.0,
    }];
    let sub_h = (ui.available_height() - 28.0).max(145.0);
    draw_scientific_plot(ui, "dd_plot", sub_h, &dd_series, &[], &[]);
}

pub fn render_phase_portrait_tab(ui: &mut Ui, mod_id: ModuleId, workset: &WorksetState) {
    let plot_h = (ui.available_height() - 48.0).max(280.0);

    match mod_id {
        ModuleId::T2MayBifurcation => {
            ui.label(
                RichText::new(
                    "Live Bifurcation Orbit Diagram (Asymptotic Attractors across Intrinsic Growth Rate r)",
                )
                .strong(),
            );
            let sim = &workset.t2_may;
            let mut orbit_pts = Vec::with_capacity(3200);
            let r_max = if sim.ricker_map { 4.2 } else { 3.0 };
            for step in 0..120 {
                let r = 0.5 + (r_max - 0.5) * (step as f64 / 119.0);
                let mut x = 0.25 * sim.carrying_k;
                for burn in 0..105 {
                    x = sim.map_eval(x, r);
                    if burn >= 80 {
                        orbit_pts.push([r, x]);
                    }
                }
            }
            let s = [PlotSeries {
                name: "Asymptotic Orbit N*".to_string(),
                color: Color32::from_rgb(56, 189, 248),
                points: orbit_pts,
                width: 0.0,
                as_dots: true,
                dot_radius: 1.6,
            }];
            let vm = [PlotVMarker {
                x: sim.r_param,
                color: Color32::from_rgb(250, 204, 21),
            }];
            draw_scientific_plot(ui, "may_bif", plot_h, &s, &[], &vm);
        }

        ModuleId::T5WrightLandscape => {
            ui.label(
                RichText::new(
                    "Topographic Phase Map: Deme Trajectories (X, Y) Climbing Local Peak (0.24, 0.26) vs. Global Peak (0.76, 0.74)",
                )
                .strong(),
            );
            let sim = &workset.t5_landscape;
            let colors = [
                Color32::from_rgb(74, 222, 128),
                Color32::from_rgb(56, 189, 248),
                Color32::from_rgb(251, 146, 60),
                Color32::from_rgb(244, 114, 182),
                Color32::from_rgb(250, 204, 21),
            ];
            let mut series_vec = Vec::new();
            for (idx, d) in sim.demes.iter().enumerate() {
                let pts: Vec<[f64; 2]> = d.trail.iter().map(|t| [t[0], t[1]]).collect();
                series_vec.push(PlotSeries {
                    name: format!("Deme {} (W={:.2})", idx + 1, d.fitness),
                    color: colors[idx % colors.len()],
                    points: pts,
                    width: 1.9,
                    as_dots: false,
                    dot_radius: 0.0,
                });
            }
            series_vec.push(PlotSeries {
                name: "Adaptive Peaks (Local & Global)".to_string(),
                color: Color32::from_rgb(254, 240, 138),
                points: vec![[0.24, 0.26], [0.76, 0.74]],
                width: 0.0,
                as_dots: true,
                dot_radius: 6.5,
            });
            draw_scientific_plot(ui, "t5_contour_map", plot_h, &series_vec, &[], &[]);
        }

        ModuleId::T6LeslieLefkovitch => {
            ui.label(
                RichText::new(
                    "Conservation Leverage Profile: Summed Stage Elasticity E_j (%) vs. Stable Stage Share w_j (%)",
                )
                .strong(),
            );
            let sim = &workset.t6_matrix;
            let mut bars = Vec::new();
            let mut w_pts = Vec::new();
            for (j, &el) in sim.eigen.stage_elasticity_totals.iter().enumerate() {
                let x = (j + 1) as f64;
                bars.push(PlotBar {
                    x,
                    y: el * 100.0,
                    width: 0.52,
                    color: Color32::from_rgb(251, 146, 60),
                });
                let wj = sim.eigen.stable_distribution_w.get(j).copied().unwrap_or(0.0);
                w_pts.push([x, wj * 100.0]);
            }
            let s = [PlotSeries {
                name: "Stable Stage Proportion w_j (%)".to_string(),
                color: Color32::from_rgb(56, 189, 248),
                points: w_pts,
                width: 2.3,
                as_dots: false,
                dot_radius: 0.0,
            }];
            draw_scientific_plot(ui, "t6_elasticity_bars", plot_h, &s, &bars, &[]);
        }

        ModuleId::T4KingmanCoalescent => {
            ui.label(
                RichText::new(
                    "Unfolded Site Frequency Spectrum (SFS): Simulated Branch Mutations vs. Neutral theta / i",
                )
                .strong(),
            );
            let sim = &workset.t4_coalescent;
            let reps = sim.replicate_index.max(1) as f64;
            let mut bars = Vec::new();
            let mut exp_pts = Vec::new();
            for (idx, &cum) in sim.cumulative_sfs.iter().enumerate() {
                let i_class = (idx + 1) as f64;
                bars.push(PlotBar {
                    x: i_class,
                    y: cum / reps,
                    width: 0.58,
                    color: Color32::from_rgb(56, 189, 248),
                });
                exp_pts.push([i_class, sim.theta_locus / i_class]);
            }
            let s = [PlotSeries {
                name: "Neutral Expectation theta / i".to_string(),
                color: Color32::from_rgb(250, 204, 21),
                points: exp_pts,
                width: 2.3,
                as_dots: false,
                dot_radius: 0.0,
            }];
            draw_scientific_plot(ui, "sfs_plot", plot_h, &s, &bars, &[]);
        }

        _ => {
            let (records, _, names) = workset.history_and_names(mod_id);
            let y_label = if names[1].is_empty() {
                "Per-Capita Growth"
            } else {
                names[1]
            };
            ui.label(
                RichText::new(format!(
                    "2D State-Space Phase Portrait: {} (X) vs. {} (Y)",
                    names[0], y_label
                ))
                .strong(),
            );
            let traj: Vec<[f64; 2]> = records
                .iter()
                .map(|r| {
                    let y = if names[1].is_empty() {
                        r.per_capita_growth
                    } else {
                        r.values[1]
                    };
                    [r.values[0], y]
                })
                .collect();

            let mut series_vec = vec![PlotSeries {
                name: "Phase Trajectory".to_string(),
                color: Color32::from_rgb(74, 222, 128),
                points: traj,
                width: 2.0,
                as_dots: false,
                dot_radius: 0.0,
            }];

            if let Some(last) = records.last() {
                let y = if names[1].is_empty() {
                    last.per_capita_growth
                } else {
                    last.values[1]
                };
                series_vec.push(PlotSeries {
                    name: "Current State".to_string(),
                    color: Color32::from_rgb(250, 204, 21),
                    points: vec![[last.values[0], y]],
                    width: 0.0,
                    as_dots: true,
                    dot_radius: 5.0,
                });
            }

            let mut vmarkers = Vec::new();
            if mod_id == ModuleId::V1HareLynx {
                let s = &workset.v1_hare_lynx;
                let mut prey_nc = Vec::new();
                for i in 1..=90 {
                    let h = s.k_hare * (i as f64 / 90.0);
                    let fh = s
                        .functional_response
                        .evaluate(h, s.attack_rate, s.handling_time);
                    if fh > 1e-6 {
                        let l_iso = (s.r_hare * h * (1.0 - h / s.k_hare)) / fh;
                        if l_iso >= 0.0 {
                            prey_nc.push([h, l_iso]);
                        }
                    }
                }
                series_vec.push(PlotSeries {
                    name: "Prey Nullcline (dH/dt = 0)".to_string(),
                    color: Color32::from_rgb(56, 189, 248),
                    points: prey_nc,
                    width: 1.6,
                    as_dots: false,
                    dot_radius: 0.0,
                });
                let denom =
                    s.attack_rate * (s.conversion_eff - s.lynx_mortality * s.handling_time);
                if denom > 1e-6 {
                    vmarkers.push(PlotVMarker {
                        x: s.lynx_mortality / denom,
                        color: Color32::from_rgb(248, 113, 113),
                    });
                }
            }

            draw_scientific_plot(ui, "phase_plane", plot_h, &series_vec, &[], &vmarkers);
        }
    }
}
