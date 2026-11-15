use std::fs;

use egui::{Align2, Color32, CornerRadius, FontId, Pos2, Rect, RichText, Sense, Stroke, Ui, Vec2};

use crate::gui::plots_panel::{
    draw_scientific_plot, draw_scientific_plot_ex, PlotBand, PlotBar, PlotHMarker, PlotSeries,
};
use crate::haplotype_network::NetworkAlgorithm;
use crate::popgen_alignment::{builtin_alignment_presets, GapHandling};
use crate::workset::EmpiricalWorkspaceState;

fn base_color(ch: char) -> Color32 {
    match ch {
        'A' => Color32::from_rgb(74, 222, 128),
        'C' => Color32::from_rgb(96, 165, 250),
        'G' => Color32::from_rgb(250, 204, 21),
        'T' => Color32::from_rgb(248, 113, 113),
        _ => Color32::from_rgb(148, 163, 184),
    }
}

/// Renders the Left Control Panel shared across Empirical Modules E1, E2, E3, and E4.
pub fn render_empirical_controls(
    ui: &mut Ui,
    emp: &mut EmpiricalWorkspaceState,
    status_msg: &mut String,
) {
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.label(RichText::new("Empirical FASTA / Alignment Controls").strong());
            ui.add_space(4.0);

            ui.label("Built-in Benchmark Alignments:");
            for (idx, (title, _)) in builtin_alignment_presets().iter().enumerate() {
                if ui.button(*title).clicked() {
                    emp.load_preset(idx);
                    *status_msg = format!("Loaded benchmark alignment: {title}");
                }
            }

            ui.add_space(6.0);
            if ui.button("Load FASTA / PHYLIP / NEXUS File...").clicked() {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter(
                        "Sequence Alignments",
                        &["fasta", "fa", "fas", "phy", "phylip", "nex", "nexus", "txt"],
                    )
                    .pick_file()
                {
                    match fs::read_to_string(&path) {
                        Ok(content) => {
                            emp.dataset_name = path
                                .file_name()
                                .and_then(|s| s.to_str())
                                .unwrap_or("Uploaded Alignment")
                                .to_string();
                            emp.raw_alignment_text = content;
                            emp.recompute();
                            *status_msg = format!("Loaded empirical alignment from {}", path.display());
                        }
                        Err(e) => {
                            *status_msg = format!("File read error: {e}");
                        }
                    }
                }
            }

            ui.separator();
            ui.label(RichText::new("Gap Handling & Network Options").strong());
            let mut changed = false;
            for mode in [GapHandling::PairwiseDeletion, GapHandling::CompleteDeletion] {
                if ui
                    .radio_value(&mut emp.gap_handling, mode, mode.label())
                    .changed()
                {
                    changed = true;
                }
            }

            for alg in [
                NetworkAlgorithm::BandeltMedianJoining,
                NetworkAlgorithm::MinimumSpanningNetwork,
                NetworkAlgorithm::TcsStatisticalParsimony,
            ] {
                if ui
                    .radio_value(&mut emp.network_algorithm, alg, alg.label())
                    .changed()
                {
                    changed = true;
                }
            }

            if emp.network_algorithm == NetworkAlgorithm::BandeltMedianJoining
                && ui
                    .add(egui::Slider::new(&mut emp.mj_epsilon, 0..=3).text("MJ epsilon"))
                    .changed()
            {
                changed = true;
            }
            ui.checkbox(&mut emp.show_mutation_labels, "Show edge site mutations");

            ui.separator();
            ui.label(RichText::new("Bayesian Skyline & Clock Calibration (E3)").strong());
            if ui
                .add(egui::Slider::new(&mut emp.bsp_groups, 1..=8).text("Skyline groups m"))
                .changed()
            {
                changed = true;
            }
            if ui
                .add(
                    egui::Slider::new(&mut emp.bsp_mcmc_steps, 600..=8000)
                        .logarithmic(true)
                        .text("MCMC iterations"),
                )
                .changed()
            {
                changed = true;
            }
            if ui
                .add(
                    egui::Slider::new(&mut emp.molecular_clock_rate, 1.0e-9..=1.0e-5)
                        .logarithmic(true)
                        .text("Clock mu (subs/site/gen)"),
                )
                .changed()
            {
                changed = true;
            }
            if ui
                .add(
                    egui::Slider::new(&mut emp.generation_time_years, 0.1..=25.0)
                        .logarithmic(true)
                        .text("Generation time (yr)"),
                )
                .changed()
            {
                changed = true;
            }

            ui.separator();
            ui.label(RichText::new("Sliding-Window Neutrality Scan (E4)").strong());
            let max_w = emp
                .sequences
                .first()
                .map(|s| s.sequence.len())
                .unwrap_or(180)
                .max(12);
            if ui
                .add(egui::Slider::new(&mut emp.sliding_window_bp, 8..=max_w).text("Window size (bp)"))
                .changed()
            {
                changed = true;
            }
            if ui
                .add(
                    egui::Slider::new(&mut emp.sliding_step_bp, 1..=max_w.min(60))
                        .text("Step size (bp)"),
                )
                .changed()
            {
                changed = true;
            }

            if changed || ui.button("Recompute All Empirical Modules").clicked() {
                emp.recompute();
            }

            if let Some(net) = &mut emp.network {
                if ui.button("Re-relax Haplotype Force Layout").clicked() {
                    net.relax_layout(180);
                }
            }

            ui.separator();
            ui.label(RichText::new("Export Empirical Artifacts").strong());
            ui.horizontal_wrapped(|ui| {
                if ui.button("SVG Network...").clicked() {
                    if let Some(net) = &emp.network {
                        if let Some(path) = rfd::FileDialog::new()
                            .set_file_name("haplotype_network.svg")
                            .save_file()
                        {
                            let _ = fs::write(&path, net.to_svg(1024, 720));
                            *status_msg = format!("Exported SVG network to {}", path.display());
                        }
                    }
                }
                if ui.button("GEXF 1.3...").clicked() {
                    if let Some(net) = &emp.network {
                        if let Some(path) = rfd::FileDialog::new()
                            .set_file_name("haplotype_network.gexf")
                            .save_file()
                        {
                            let _ = fs::write(&path, net.to_gexf());
                            *status_msg = format!("Exported GEXF graph to {}", path.display());
                        }
                    }
                }
            });
            ui.horizontal_wrapped(|ui| {
                if ui.button("Skyline CSV...").clicked() {
                    if let Some(sky) = &emp.skyline_report {
                        if let Some(path) = rfd::FileDialog::new()
                            .set_file_name("bayesian_skyline.csv")
                            .save_file()
                        {
                            let _ = fs::write(&path, sky.to_csv());
                            *status_msg = format!("Exported Bayesian Skyline CSV to {}", path.display());
                        }
                    }
                }
                if ui.button("Sliding-Window CSV...").clicked() {
                    if let Some(win) = &emp.sliding_report {
                        if let Some(path) = rfd::FileDialog::new()
                            .set_file_name("sliding_window_neutrality.csv")
                            .save_file()
                        {
                            let _ = fs::write(&path, win.to_csv());
                            *status_msg = format!("Exported Sliding-Window CSV to {}", path.display());
                        }
                    }
                }
            });

            if let Some(err) = &emp.last_error {
                ui.colored_label(Color32::from_rgb(248, 113, 113), format!("Error: {err}"));
            }
        });
}

/// Renders Module E1: Global alignment polymorphism strip, neutrality summary, LD, and Hudson F_ST.
pub fn render_popgen_summary_view(ui: &mut Ui, emp: &mut EmpiricalWorkspaceState) {
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.collapsing("Editable Multi-FASTA / PHYLIP / NEXUS Alignment Source", |ui| {
                let resp = ui.add(
                    egui::TextEdit::multiline(&mut emp.raw_alignment_text)
                        .font(egui::TextStyle::Monospace)
                        .desired_rows(7)
                        .desired_width(f32::INFINITY),
                );
                if resp.changed() {
                    emp.recompute();
                }
            });

            if !emp.sequences.is_empty() {
                ui.add_space(6.0);
                ui.label(
                    RichText::new(
                        "Aligned Nucleotide Matrix & Segregating Polymorphism Strip (first 80 bp)",
                    )
                    .strong(),
                );
                egui::ScrollArea::horizontal().show(ui, |ui| {
                    let ref_seq = emp.sequences[0].sequence.as_bytes();
                    let show_len = ref_seq.len().min(80);
                    egui::Grid::new("dna_matrix_preview")
                        .spacing([2.0, 2.0])
                        .show(ui, |ui| {
                            ui.label(RichText::new("Taxon | Deme").monospace().strong());
                            for c in 0..show_len {
                                if (c + 1) % 10 == 0 {
                                    ui.label(RichText::new("|").monospace().color(Color32::GRAY));
                                } else {
                                    ui.label(RichText::new(" ").monospace());
                                }
                            }
                            ui.end_row();

                            for (r_idx, seq) in emp.sequences.iter().take(18).enumerate() {
                                ui.label(
                                    RichText::new(format!("{} [{}]", seq.id, seq.population))
                                        .monospace()
                                        .size(11.5),
                                );
                                for (c, ch) in seq.sequence.chars().take(show_len).enumerate() {
                                    let is_diff = r_idx > 0 && (ch as u8) != ref_seq[c];
                                    let txt = if r_idx > 0 && !is_diff {
                                        RichText::new(".").monospace().color(Color32::DARK_GRAY)
                                    } else {
                                        RichText::new(ch.to_string())
                                            .monospace()
                                            .strong()
                                            .color(base_color(ch))
                                    };
                                    ui.label(txt);
                                }
                                ui.end_row();
                            }
                        });
                });
            }

            let Some(sum) = &emp.summary else {
                return;
            };

            ui.add_space(10.0);
            ui.label(
                RichText::new(format!(
                    "Population Genetics & Neutrality Summary — {}",
                    sum.dataset_name
                ))
                .size(16.0)
                .strong(),
            );

            egui::Grid::new("popgen_metrics_grid")
                .striped(true)
                .min_col_width(180.0)
                .show(ui, |ui| {
                    ui.label("Sample Size (n) / Length (L):");
                    ui.label(
                        RichText::new(format!(
                            "n = {} sequences, L = {} bp (L_eff = {} bp, GC = {:.1}%)",
                            sum.num_sequences,
                            sum.alignment_length,
                            sum.effective_sites,
                            sum.gc_content_pct
                        ))
                        .monospace(),
                    );
                    ui.end_row();

                    ui.label("Segregating Sites (S):");
                    ui.label(
                        RichText::new(format!(
                            "S = {}  (Parsimony-informative = {}, Singletons eta_1 = {})",
                            sum.segregating_sites_s,
                            sum.parsimony_informative_sites,
                            sum.singleton_sites_eta1
                        ))
                        .monospace(),
                    );
                    ui.end_row();

                    ui.label("Haplotypes (K) & Diversity (Hd):");
                    ui.label(
                        RichText::new(format!(
                            "K = {} distinct haplotypes,  Hd = {:.5}",
                            sum.num_haplotypes_k, sum.haplotype_diversity_hd
                        ))
                        .monospace(),
                    );
                    ui.end_row();

                    ui.label("Nucleotide Diversity (pi) & k_hat:");
                    ui.label(
                        RichText::new(format!(
                            "pi = {:.6} / site,  mean pairwise diffs k_hat = {:.4}",
                            sum.nucleotide_diversity_pi, sum.mean_pairwise_diffs_k
                        ))
                        .monospace(),
                    );
                    ui.end_row();

                    ui.label("Watterson's Estimator (theta_W):");
                    ui.label(
                        RichText::new(format!(
                            "theta_W = {:.4} / locus  ({:.6} / site)",
                            sum.watterson_theta_per_locus, sum.watterson_theta_per_site
                        ))
                        .monospace(),
                    );
                    ui.end_row();

                    ui.label("Tajima's D (1989):");
                    let td_str = match (
                        sum.tajimas_d,
                        sum.tajimas_d_pval_normal,
                        sum.tajimas_d_pval_beta,
                    ) {
                        (Some(d), Some(pn), Some(pb)) => {
                            format!("D = {:+.4}   (p_norm = {:.4}, p_beta = {:.4})", d, pn, pb)
                        }
                        (Some(d), Some(pn), None) => format!("D = {:+.4}   (p_norm = {:.4})", d, pn),
                        _ => "Insufficient polymorphism / sample size (requires n >= 4, S >= 1)"
                            .to_string(),
                    };
                    ui.label(RichText::new(td_str).monospace().strong());
                    ui.end_row();

                    ui.label("Fu & Li's D* and F* (1993):");
                    let fl_str = match (sum.fu_li_d_star, sum.fu_li_f_star) {
                        (Some(ds), Some(fs)) => format!("D* = {:+.4},   F* = {:+.4}", ds, fs),
                        _ => "N/A".to_string(),
                    };
                    ui.label(RichText::new(fl_str).monospace());
                    ui.end_row();

                    ui.label("Harpending Mismatch Raggedness (r):");
                    ui.label(
                        RichText::new(format!(
                            "r = {:.5}  (see Tab 3 for Bayesian Skyline & Rogers-Harpending Expansion fit)",
                            sum.harpending_raggedness
                        ))
                        .monospace(),
                    );
                    ui.end_row();
                });

            if !sum.pop_divergences.is_empty() {
                ui.add_space(10.0);
                ui.label(
                    RichText::new(
                        "Multi-Population Structure: Hudson (1992) F_ST, D_XY, D_a & Gene Flow Nm",
                    )
                    .strong(),
                );
                egui::Grid::new("pop_divergence_grid")
                    .striped(true)
                    .min_col_width(95.0)
                    .show(ui, |ui| {
                        for hdr in [
                            "Deme Pair",
                            "pi_A",
                            "pi_B",
                            "D_XY",
                            "D_a (net)",
                            "Hudson F_ST",
                            "Migrants Nm",
                        ] {
                            ui.label(RichText::new(hdr).strong());
                        }
                        ui.end_row();
                        for row in &sum.pop_divergences {
                            ui.label(format!("{} vs {}", row.pop_a, row.pop_b));
                            ui.label(format!("{:.5}", row.pi_within_a));
                            ui.label(format!("{:.5}", row.pi_within_b));
                            ui.label(format!("{:.5}", row.d_xy));
                            ui.label(format!("{:.5}", row.d_a));
                            ui.label(RichText::new(format!("{:.4}", row.hudson_fst)).strong());
                            ui.label(format!("{:.2}", row.effective_nm));
                            ui.end_row();
                        }
                    });
            }

            ui.add_space(10.0);
            ui.label(RichText::new("Pairwise Nucleotide Mismatch Distribution").strong());
            let bars: Vec<PlotBar> = sum
                .mismatch_histogram
                .iter()
                .map(|&(d, cnt)| PlotBar {
                    x: d as f64,
                    y: cnt as f64,
                    width: 0.65,
                    color: Color32::from_rgb(56, 189, 248),
                })
                .collect();
            draw_scientific_plot(ui, "mismatch_dist_plot", 165.0, &[], &bars, &[]);

            if !sum.ld_pairs.is_empty() {
                ui.add_space(8.0);
                ui.collapsing("Linkage Disequilibrium (r^2 & |D'| across Informative Sites)", |ui| {
                    egui::Grid::new("ld_grid").striped(true).show(ui, |ui| {
                        for hdr in ["Site A", "Site B", "Distance (bp)", "r^2", "|D'|"] {
                            ui.label(RichText::new(hdr).strong());
                        }
                        ui.end_row();
                        for ld in sum.ld_pairs.iter().take(20) {
                            ui.label(format!("bp {}", ld.site_a_1based));
                            ui.label(format!("bp {}", ld.site_b_1based));
                            ui.label(format!("{}", ld.distance_bp));
                            ui.label(format!("{:.4}", ld.r_squared));
                            ui.label(format!("{:.4}", ld.d_prime_abs));
                            ui.end_row();
                        }
                    });
                });
            }
        });
}

/// Renders Module E2: Interactive Haplotype Network Studio.
pub fn render_haplotype_network_view(
    ui: &mut Ui,
    emp: &mut EmpiricalWorkspaceState,
    dragged_node: &mut Option<usize>,
) {
    let Some(net) = &mut emp.network else {
        ui.label("No valid alignment loaded to construct a haplotype network.");
        return;
    };

    let palette = [
        Color32::from_rgb(56, 189, 248),
        Color32::from_rgb(251, 146, 60),
        Color32::from_rgb(52, 211, 153),
        Color32::from_rgb(167, 139, 250),
        Color32::from_rgb(250, 204, 21),
        Color32::from_rgb(244, 63, 94),
    ];

    ui.horizontal_wrapped(|ui| {
        ui.label(RichText::new(net.algorithm.label()).strong());
        ui.separator();
        ui.label(format!(
            "Nodes: {} ({} sampled haplotypes, {} median vectors) | Edges: {} | TCS 95% Limit: {} steps",
            net.nodes.len(),
            net.nodes.iter().filter(|n| !n.is_median_vector).count(),
            net.nodes.iter().filter(|n| n.is_median_vector).count(),
            net.edges.len(),
            net.tcs_parsimony_limit
        ));
    });

    ui.horizontal_wrapped(|ui| {
        ui.label("Demes:");
        for (idx, pop) in net.populations.iter().enumerate() {
            ui.colored_label(palette[idx % palette.len()], format!("● {pop}"));
        }
        ui.colored_label(Color32::from_rgb(239, 68, 68), "◆ Median Vector (mv)");
    });

    let avail = ui.available_size_before_wrap();
    let canvas_h = (avail.y - 125.0).max(290.0);
    let (response, painter) =
        ui.allocate_painter(Vec2::new(avail.x.max(360.0), canvas_h), Sense::click_and_drag());
    let rect = response.rect;
    painter.rect_filled(rect, CornerRadius::same(6), Color32::from_rgb(12, 18, 28));

    let to_pos = |nx: f64, ny: f64| -> Pos2 {
        Pos2::new(
            rect.left() + (nx as f32) * rect.width(),
            rect.top() + (ny as f32) * rect.height(),
        )
    };

    if let Some(pointer) = response.interact_pointer_pos() {
        if response.drag_started() {
            let mut closest = None;
            let mut best_d = 28.0f32;
            for node in &net.nodes {
                let p = to_pos(node.x, node.y);
                let d = p.distance(pointer);
                if d < best_d {
                    best_d = d;
                    closest = Some(node.id);
                }
            }
            *dragged_node = closest;
        }
        if response.dragged() {
            if let Some(nid) = *dragged_node {
                if let Some(node) = net.nodes.get_mut(nid) {
                    node.x = (((pointer.x - rect.left()) / rect.width()) as f64).clamp(0.06, 0.94);
                    node.y = (((pointer.y - rect.top()) / rect.height()) as f64).clamp(0.06, 0.94);
                }
            }
        }
    }
    if response.drag_stopped() {
        *dragged_node = None;
    }

    for edge in &net.edges {
        let nu = &net.nodes[edge.source];
        let nv = &net.nodes[edge.target];
        let p1 = to_pos(nu.x, nu.y);
        let p2 = to_pos(nv.x, nv.y);
        painter.line_segment([p1, p2], Stroke::new(2.2, Color32::from_rgb(148, 163, 184)));

        let dir = p2 - p1;
        let len = dir.length().max(1.0);
        let perp = Vec2::new(-dir.y / len, dir.x / len);
        let ticks = edge.mutational_steps.min(12);
        for k in 1..=ticks {
            let frac = 0.30 + 0.40 * ((k as f32 - 0.5) / (ticks as f32));
            let mid = p1 + dir * frac;
            painter.line_segment(
                [mid - perp * 5.5, mid + perp * 5.5],
                Stroke::new(2.0, Color32::from_rgb(250, 204, 21)),
            );
        }

        if emp.show_mutation_labels && !edge.mutations.is_empty() {
            let label_txt = if edge.mutations.len() <= 3 {
                edge.mutations.join(",")
            } else {
                format!("{} steps", edge.mutational_steps)
            };
            let mid = p1 + dir * 0.5 + perp * 11.0;
            painter.text(
                mid,
                Align2::CENTER_CENTER,
                label_txt,
                FontId::monospace(10.5),
                Color32::from_rgb(253, 224, 71),
            );
        }
    }

    for node in &net.nodes {
        let center = to_pos(node.x, node.y);
        if node.is_median_vector {
            painter.circle_filled(center, 6.5, Color32::from_rgb(239, 68, 68));
            painter.circle_stroke(center, 6.5, Stroke::new(1.5, Color32::WHITE));
            painter.text(
                center + Vec2::new(10.0, -8.0),
                Align2::LEFT_CENTER,
                &node.label,
                FontId::monospace(11.0),
                Color32::from_rgb(252, 165, 165),
            );
        } else {
            let radius = 11.0 + 5.5 * (node.sample_count as f32).sqrt();
            let active: Vec<(usize, usize)> = node
                .pop_counts
                .iter()
                .copied()
                .enumerate()
                .filter(|&(_, c)| c > 0)
                .collect();

            if active.len() <= 1 {
                let col = palette[active.first().map(|x| x.0).unwrap_or(0) % palette.len()];
                painter.circle_filled(center, radius, col);
            } else {
                let mut start_a = 0.0f32;
                for (p_idx, cnt) in active {
                    let sweep =
                        std::f32::consts::TAU * (cnt as f32) / (node.sample_count.max(1) as f32);
                    let steps = ((sweep * 8.0).ceil() as usize).max(3);
                    let col = palette[p_idx % palette.len()];
                    for s in 0..steps {
                        let a0 = start_a + sweep * (s as f32 / steps as f32);
                        let a1 = start_a + sweep * ((s + 1) as f32 / steps as f32);
                        let p_a = center + Vec2::new(radius * a0.cos(), radius * a0.sin());
                        let p_b = center + Vec2::new(radius * a1.cos(), radius * a1.sin());
                        painter.add(egui::Shape::convex_polygon(
                            vec![center, p_a, p_b],
                            col,
                            Stroke::NONE,
                        ));
                    }
                    start_a += sweep;
                }
            }
            painter.circle_stroke(center, radius, Stroke::new(1.8, Color32::WHITE));
            painter.text(
                center - Vec2::new(0.0, radius + 9.0),
                Align2::CENTER_CENTER,
                format!("{} (n={})", node.label, node.sample_count),
                FontId::proportional(12.0),
                Color32::WHITE,
            );
        }
    }

    let _ = Rect::NOTHING;
    egui::ScrollArea::vertical()
        .max_height(105.0)
        .show(ui, |ui| {
            egui::Grid::new("hap_membership_grid")
                .striped(true)
                .show(ui, |ui| {
                    ui.label(RichText::new("Haplotype").strong());
                    ui.label(RichText::new("n").strong());
                    ui.label(RichText::new("Member Sequence IDs").strong());
                    ui.end_row();
                    for node in &net.nodes {
                        if !node.is_median_vector {
                            ui.label(RichText::new(&node.label).monospace().strong());
                            ui.label(format!("{}", node.sample_count));
                            ui.label(node.member_ids.join(", "));
                            ui.end_row();
                        }
                    }
                });
        });
}

/// Renders Module E3: MCMC Bayesian Skyline Plot (Drummond et al. 2005), Strimmer-Pybus AICc
/// Generalized Skyline, and Rogers-Harpending (1992) Coalescent Sudden Demographic Expansion fit.
pub fn render_skyline_and_expansion_view(ui: &mut Ui, emp: &EmpiricalWorkspaceState) {
    let Some(sky) = &emp.skyline_report else {
        ui.label("Requires at least 4 aligned sequences to reconstruct a coalescent Skyline.");
        return;
    };

    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    RichText::new(format!("Inferred Demographic Trajectory: {}", sky.demographic_trend))
                        .size(15.5)
                        .strong()
                        .color(Color32::from_rgb(134, 239, 172)),
                );
                ui.separator();
                ui.label(format!(
                    "T_MRCA = {:.5} subs/site ({:.0} yr ago) | MCMC accept = {:.1}% | Generalized Skyline AICc = {:.2}",
                    sky.tmrca_subs_per_site,
                    sky.tmrca_scaled,
                    sky.mcmc_acceptance_rate * 100.0,
                    sky.generalized_aicc
                ));
            });

            ui.add_space(6.0);
            ui.label(
                RichText::new(
                    "Bayesian Skyline Plot: Effective Population Size N_e(t) Backward from Present (t = 0 -> T_MRCA)",
                )
                .strong(),
            );

            let cred_band = [PlotBand {
                name: "95% Bayesian Credible Interval".to_string(),
                fill: Color32::from_rgba_unmultiplied(56, 189, 248, 52),
                bounds: sky
                    .grid
                    .iter()
                    .map(|g| [g.scaled_time, g.ne_lower_95, g.ne_upper_95])
                    .collect(),
            }];

            let bsp_series = [
                PlotSeries {
                    name: "BSP Posterior Median N_e(t)".to_string(),
                    color: Color32::from_rgb(56, 189, 248),
                    points: sky
                        .grid
                        .iter()
                        .map(|g| [g.scaled_time, g.ne_median])
                        .collect(),
                    width: 2.5,
                    as_dots: false,
                    dot_radius: 0.0,
                },
                PlotSeries {
                    name: "Generalized Skyline (AICc Step)".to_string(),
                    color: Color32::from_rgb(250, 204, 21),
                    points: sky
                        .grid
                        .iter()
                        .map(|g| [g.scaled_time, g.generalized_skyline_ne])
                        .collect(),
                    width: 1.6,
                    as_dots: false,
                    dot_radius: 0.0,
                },
            ];
            draw_scientific_plot_ex(
                ui,
                "bsp_skyline_plot",
                235.0,
                &cred_band,
                &bsp_series,
                &[],
                &[],
                &[],
            );

            ui.add_space(10.0);
            let exp = &sky.expansion;
            ui.label(
                RichText::new(
                    "Rogers & Harpending (1992) Coalescent Sudden Demographic Expansion Wave Fit",
                )
                .size(15.0)
                .strong(),
            );

            egui::Grid::new("rogers_harpending_grid")
                .striped(true)
                .min_col_width(170.0)
                .show(ui, |ui| {
                    ui.label("Ancestral / Post-Expansion Theta:");
                    ui.label(
                        RichText::new(format!(
                            "theta_0 = {:.4}  ->  theta_1 = {:.4}  (2 * N * u_locus)",
                            exp.theta_0, exp.theta_1
                        ))
                        .monospace(),
                    );
                    ui.end_row();

                    ui.label("Expansion Horizon (tau = 2*u*t):");
                    ui.label(
                        RichText::new(format!(
                            "tau = {:.4} mutational units  (~{:.0} calibrated yr ago)",
                            exp.tau, exp.expansion_time_scaled
                        ))
                        .monospace(),
                    );
                    ui.end_row();

                    ui.label("Goodness-of-Fit SSD:");
                    ui.label(
                        RichText::new(format!(
                            "SSD_expansion = {:.5} (bootstrap p = {:.3})  vs  SSD_constant = {:.5}",
                            exp.ssd_expansion, exp.bootstrap_p_ssd, exp.ssd_constant
                        ))
                        .monospace()
                        .strong(),
                    );
                    ui.end_row();

                    ui.label("Harpending Raggedness (r):");
                    ui.label(
                        RichText::new(format!(
                            "r = {:.5} (coalescent bootstrap p = {:.3}; p >= 0.05 supports expansion wave)",
                            exp.observed_raggedness, exp.bootstrap_p_raggedness
                        ))
                        .monospace(),
                    );
                    ui.end_row();
                });

            let obs_bars: Vec<PlotBar> = exp
                .wave_comparison
                .iter()
                .map(|&(d, obs, _, _)| PlotBar {
                    x: d as f64,
                    y: obs,
                    width: 0.56,
                    color: Color32::from_rgba_unmultiplied(96, 165, 250, 170),
                })
                .collect();
            let wave_series = [
                PlotSeries {
                    name: "Fitted Sudden-Expansion Wave F_i(tau)".to_string(),
                    color: Color32::from_rgb(74, 222, 128),
                    points: exp
                        .wave_comparison
                        .iter()
                        .map(|&(d, _, e_exp, _)| [d as f64, e_exp])
                        .collect(),
                    width: 2.4,
                    as_dots: false,
                    dot_radius: 0.0,
                },
                PlotSeries {
                    name: "Constant-Population Equilibrium Geometric".to_string(),
                    color: Color32::from_rgb(248, 113, 113),
                    points: exp
                        .wave_comparison
                        .iter()
                        .map(|&(d, _, _, e_con)| [d as f64, e_con])
                        .collect(),
                    width: 1.7,
                    as_dots: false,
                    dot_radius: 0.0,
                },
            ];
            draw_scientific_plot(ui, "rh_wave_plot", 190.0, &wave_series, &obs_bars, &[]);

            ui.add_space(6.0);
            ui.collapsing("Ultrametric UPGMA Coalescent Intervals (k_i, Delta u_i, Classic Skyline theta_i)", |ui| {
                egui::Grid::new("coal_iv_grid").striped(true).show(ui, |ui| {
                    for h in ["Interval", "Lineages k_i", "t_start (subs)", "t_end (subs)", "Delta u_i", "Classic theta_i"] {
                        ui.label(RichText::new(h).strong());
                    }
                    ui.end_row();
                    for iv in &sky.intervals {
                        ui.label(format!("#{}", iv.index));
                        ui.label(format!("{}", iv.lineages_k));
                        ui.label(format!("{:.6}", iv.t_start_subs));
                        ui.label(format!("{:.6}", iv.t_end_subs));
                        ui.label(format!("{:.6}", iv.duration_subs));
                        ui.label(format!("{:.6}", iv.classic_theta));
                        ui.end_row();
                    }
                });
            });
        });
}

/// Renders Module E4: Sliding-Window Polymorphism (`pi`, `theta_W`), Neutrality (`Tajima's D`, `Fu & Li's D*`),
/// and Multi-Deme Divergence (`Hudson F_ST`, `D_XY`) Scanner along sequence coordinates.
pub fn render_sliding_window_view(ui: &mut Ui, emp: &EmpiricalWorkspaceState) {
    let Some(scan) = &emp.sliding_report else {
        ui.label("Requires at least 4 aligned sequences (length >= 10 bp) for sliding-window scanning.");
        return;
    };

    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    RichText::new(format!(
                        "Sliding-Window Scan: {} windows (Window = {} bp, Step = {} bp across {} bp)",
                        scan.windows.len(),
                        scan.window_bp,
                        scan.step_bp,
                        scan.alignment_length
                    ))
                    .size(15.5)
                    .strong(),
                );
            });

            ui.add_space(4.0);
            ui.label(
                RichText::new("Track 1 — Nucleotide Diversity pi(w), Watterson's theta_W(w) & Divergence D_XY(w)")
                    .strong(),
            );
            let div_series = [
                PlotSeries {
                    name: "pi / site".to_string(),
                    color: Color32::from_rgb(56, 189, 248),
                    points: scan
                        .windows
                        .iter()
                        .map(|w| [w.midpoint_bp, w.pi_per_site])
                        .collect(),
                    width: 2.2,
                    as_dots: false,
                    dot_radius: 0.0,
                },
                PlotSeries {
                    name: "theta_W / site".to_string(),
                    color: Color32::from_rgb(250, 204, 21),
                    points: scan
                        .windows
                        .iter()
                        .map(|w| [w.midpoint_bp, w.theta_w_per_site])
                        .collect(),
                    width: 1.9,
                    as_dots: false,
                    dot_radius: 0.0,
                },
                PlotSeries {
                    name: "Mean D_XY".to_string(),
                    color: Color32::from_rgb(244, 114, 182),
                    points: scan
                        .windows
                        .iter()
                        .map(|w| [w.midpoint_bp, w.mean_d_xy])
                        .collect(),
                    width: 1.6,
                    as_dots: false,
                    dot_radius: 0.0,
                },
            ];
            draw_scientific_plot(ui, "sw_diversity_track", 160.0, &div_series, &[], &[]);

            ui.add_space(6.0);
            ui.label(
                RichText::new(
                    "Track 2 — Neutrality Profile: Tajima's D(w) & Fu & Li's D*(w) (thresholds at 0 and +/- 1.96)",
                )
                .strong(),
            );
            let neut_series = [
                PlotSeries {
                    name: "Tajima's D(w)".to_string(),
                    color: Color32::from_rgb(74, 222, 128),
                    points: scan
                        .windows
                        .iter()
                        .map(|w| [w.midpoint_bp, w.tajimas_d])
                        .collect(),
                    width: 2.3,
                    as_dots: false,
                    dot_radius: 0.0,
                },
                PlotSeries {
                    name: "Fu & Li's D*(w)".to_string(),
                    color: Color32::from_rgb(251, 146, 60),
                    points: scan
                        .windows
                        .iter()
                        .map(|w| [w.midpoint_bp, w.fu_li_d_star])
                        .collect(),
                    width: 1.7,
                    as_dots: false,
                    dot_radius: 0.0,
                },
            ];
            let neut_hlines = [
                PlotHMarker {
                    y: 0.0,
                    color: Color32::from_rgba_unmultiplied(148, 163, 184, 150),
                },
                PlotHMarker {
                    y: 1.96,
                    color: Color32::from_rgba_unmultiplied(248, 113, 113, 120),
                },
                PlotHMarker {
                    y: -1.96,
                    color: Color32::from_rgba_unmultiplied(56, 189, 248, 120),
                },
            ];
            draw_scientific_plot_ex(
                ui,
                "sw_neutrality_track",
                175.0,
                &[],
                &neut_series,
                &[],
                &[],
                &neut_hlines,
            );

            ui.add_space(6.0);
            ui.label(
                RichText::new("Track 3 — Genomic Differentiation Profile: Mean Pairwise Hudson F_ST(w)")
                    .strong(),
            );
            let fst_series = [PlotSeries {
                name: "Hudson F_ST(w)".to_string(),
                color: Color32::from_rgb(167, 139, 250),
                points: scan
                    .windows
                    .iter()
                    .map(|w| [w.midpoint_bp, w.mean_hudson_fst])
                    .collect(),
                width: 2.2,
                as_dots: false,
                dot_radius: 0.0,
            }];
            draw_scientific_plot(ui, "sw_fst_track", 145.0, &fst_series, &[], &[]);

            ui.add_space(8.0);
            ui.label(RichText::new("Window-by-Window Selection & Divergence Table").strong());
            egui::Grid::new("sw_windows_table")
                .striped(true)
                .min_col_width(72.0)
                .show(ui, |ui| {
                    for hdr in [
                        "Coords (bp)",
                        "S",
                        "pi",
                        "theta_W",
                        "Tajima D",
                        "Fu-Li D*",
                        "F_ST",
                        "Regime Classification",
                    ] {
                        ui.label(RichText::new(hdr).strong());
                    }
                    ui.end_row();

                    for w in &scan.windows {
                        ui.label(format!("{}..{}", w.start_1based, w.end_1based));
                        ui.label(format!("{}", w.segregating_sites));
                        ui.label(format!("{:.4}", w.pi_per_site));
                        ui.label(format!("{:.4}", w.theta_w_per_site));
                        ui.label(RichText::new(format!("{:+.3}", w.tajimas_d)).monospace());
                        ui.label(RichText::new(format!("{:+.3}", w.fu_li_d_star)).monospace());
                        ui.label(format!("{:.3}", w.mean_hudson_fst));
                        let flag_col = if w.tajimas_d <= -1.65 {
                            Color32::from_rgb(56, 189, 248)
                        } else if w.tajimas_d >= 1.65 || w.mean_hudson_fst >= 0.55 {
                            Color32::from_rgb(250, 204, 21)
                        } else {
                            Color32::from_rgb(160, 174, 192)
                        };
                        ui.colored_label(flag_col, &w.regime_flag);
                        ui.end_row();
                    }
                });
        });
}
