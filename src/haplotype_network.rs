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

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::popgen_alignment::{is_valid_base, AlignedSequence};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NetworkAlgorithm {
    BandeltMedianJoining,
    MinimumSpanningNetwork,
    TcsStatisticalParsimony,
}

impl NetworkAlgorithm {
    pub fn label(&self) -> &'static str {
        match self {
            Self::BandeltMedianJoining => "Bandelt et al. (1999) Median-Joining (MJ)",
            Self::MinimumSpanningNetwork => "Kruskal Minimum Spanning Network (MSN with Ties)",
            Self::TcsStatisticalParsimony => "TCS 95% Statistical Parsimony Network",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HaplotypeNode {
    pub id: usize,
    pub label: String,
    pub is_median_vector: bool,
    pub sample_count: usize,
    pub member_ids: Vec<String>,
    /// Counts ordered parallel to `HaplotypeNetwork::populations`
    pub pop_counts: Vec<usize>,
    /// Condensed sequence across complete unambiguous sites
    pub sequence: Vec<u8>,
    /// Normalized 2D layout position in `[0.08, 0.92]^2`
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HaplotypeEdge {
    pub source: usize,
    pub target: usize,
    pub mutational_steps: usize,
    /// Alignment 1-based positions that mutated along this edge (e.g., `"A19G"`).
    pub mutations: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HaplotypeNetwork {
    pub algorithm: NetworkAlgorithm,
    pub epsilon: usize,
    pub tcs_parsimony_limit: usize,
    pub populations: Vec<String>,
    pub retained_site_1based: Vec<usize>,
    pub nodes: Vec<HaplotypeNode>,
    pub edges: Vec<HaplotypeEdge>,
}

fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

fn edge_mutations(a: &[u8], b: &[u8], orig_cols_1based: &[usize]) -> Vec<String> {
    let mut out = Vec::new();
    for (i, (&ba, &bb)) in a.iter().zip(b.iter()).enumerate() {
        if ba != bb {
            let pos = orig_cols_1based.get(i).copied().unwrap_or(i + 1);
            out.push(format!("{}{}{}", ba as char, pos, bb as char));
        }
    }
    out
}

struct DisjointSet {
    parent: Vec<usize>,
}

impl DisjointSet {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }
    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            self.parent[rb] = ra;
        }
    }
}

/// TCS (Templeton, Crandall & Sing 1992) 95% statistical parsimony step limit approximation.
pub fn compute_tcs_95_limit(alignment_len: usize) -> usize {
    if alignment_len < 15 {
        return 3;
    }
    let l = alignment_len as f64;
    let mut cumulative_prob = 1.0;
    let mut max_steps = 1usize;
    for j in 1..=35usize {
        let jf = j as f64;
        // Conservative homoplasy-free step connection probability approximation
        let p_j = 1.0 - (jf * (jf + 1.0)) / (2.5 * l).max(20.0);
        if p_j <= 0.0 {
            break;
        }
        cumulative_prob *= p_j;
        if cumulative_prob >= 0.95 {
            max_steps = j;
        } else {
            break;
        }
    }
    max_steps.max(2)
}

/// Kruskal Minimum Spanning Network retaining ALL equally parsimonious edges of identical weight
/// that bridge distinct connected components at each step threshold `d`.
fn build_msn_edges(seqs: &[Vec<u8>], max_step_cutoff: Option<usize>) -> Vec<(usize, usize, usize)> {
    let k = seqs.len();
    if k <= 1 {
        return Vec::new();
    }
    let mut dist_buckets: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
    for i in 0..k {
        for j in (i + 1)..k {
            let d = hamming_distance(&seqs[i], &seqs[j]);
            if let Some(lim) = max_step_cutoff {
                if d > lim {
                    continue;
                }
            }
            dist_buckets.entry(d).or_default().push((i, j));
        }
    }

    let mut dsu = DisjointSet::new(k);
    let mut chosen = Vec::new();

    for (d, pairs) in dist_buckets {
        // First collect all pairs at distance `d` whose endpoints lie in different components BEFORE merging
        let bridging: Vec<(usize, usize)> = pairs
            .into_iter()
            .filter(|&(u, v)| dsu.find(u) != dsu.find(v))
            .collect();
        for &(u, v) in &bridging {
            chosen.push((u, v, d));
        }
        for (u, v) in bridging {
            dsu.union(u, v);
        }
    }
    chosen
}

/// Consensus majority-rule median sequence for a triplet `(a, b, c)` per Bandelt et al. (1999).
/// Sites where all three differ are resolved conservatively toward `a` (lowest total cost tie-break).
fn triplet_consensus_median(a: &[u8], b: &[u8], c: &[u8]) -> Vec<u8> {
    let mut med = Vec::with_capacity(a.len());
    for i in 0..a.len() {
        let ba = a[i];
        let bb = b[i];
        let bc = c[i];
        if ba == bb || ba == bc {
            med.push(ba);
        } else if bb == bc {
            med.push(bb);
        } else {
            med.push(ba);
        }
    }
    med
}

/// Build an interactive Haplotype Network from empirical aligned DNA sequences.
pub fn build_haplotype_network(
    sequences: &[AlignedSequence],
    algorithm: NetworkAlgorithm,
    epsilon: usize,
) -> Result<HaplotypeNetwork, String> {
    if sequences.len() < 2 {
        return Err("At least 2 aligned sequences are required to build a haplotype network.".to_string());
    }
    let l_total = sequences[0].sequence.len();
    let bytes: Vec<&[u8]> = sequences.iter().map(|s| s.sequence.as_bytes()).collect();

    // Retain all columns where every sequence has a valid unambiguous A/C/G/T base
    let mut retained_cols = Vec::new();
    for c in 0..l_total {
        if bytes.iter().all(|row| is_valid_base(row[c])) {
            retained_cols.push(c);
        }
    }
    if retained_cols.is_empty() {
        return Err("No gap-free unambiguous ACGT columns found across all sequences.".to_string());
    }
    let retained_site_1based: Vec<usize> = retained_cols.iter().map(|&c| c + 1).collect();

    let mut pop_set: BTreeSet<String> = BTreeSet::new();
    for s in sequences {
        pop_set.insert(s.population.clone());
    }
    let populations: Vec<String> = pop_set.into_iter().collect();

    // Collapse identical condensed sequences into sampled haplotypes H1, H2, ...
    let mut hap_map: BTreeMap<Vec<u8>, (Vec<String>, Vec<usize>)> = BTreeMap::new();
    for s in sequences {
        let raw = s.sequence.as_bytes();
        let key: Vec<u8> = retained_cols.iter().map(|&c| raw[c]).collect();
        let pop_idx = populations.iter().position(|p| p == &s.population).unwrap_or(0);
        let entry = hap_map
            .entry(key)
            .or_insert_with(|| (Vec::new(), vec![0; populations.len()]));
        entry.0.push(s.id.clone());
        entry.1[pop_idx] += 1;
    }

    // Order sampled haplotypes by descending sample frequency
    let mut ordered_haps: Vec<(Vec<u8>, Vec<String>, Vec<usize>)> = hap_map
        .into_iter()
        .map(|(seq, (members, pcounts))| (seq, members, pcounts))
        .collect();
    ordered_haps.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(&b.0)));

    let mut nodes: Vec<HaplotypeNode> = Vec::new();
    for (idx, (seq, members, pcounts)) in ordered_haps.into_iter().enumerate() {
        nodes.push(HaplotypeNode {
            id: idx,
            label: format!("H{}", idx + 1),
            is_median_vector: false,
            sample_count: members.len(),
            member_ids: members,
            pop_counts: pcounts,
            sequence: seq,
            x: 0.5,
            y: 0.5,
        });
    }

    let tcs_limit = compute_tcs_95_limit(retained_cols.len());

    let raw_edges = match algorithm {
        NetworkAlgorithm::MinimumSpanningNetwork => {
            let seqs: Vec<Vec<u8>> = nodes.iter().map(|n| n.sequence.clone()).collect();
            build_msn_edges(&seqs, None)
        }
        NetworkAlgorithm::TcsStatisticalParsimony => {
            let seqs: Vec<Vec<u8>> = nodes.iter().map(|n| n.sequence.clone()).collect();
            build_msn_edges(&seqs, Some(tcs_limit))
        }
        NetworkAlgorithm::BandeltMedianJoining => {
            // Bandelt, Forster & Rohl (1999) Median-Joining (MJ) algorithm with homoplasy tolerance epsilon
            let eps = epsilon.min(3);
            let max_median_vectors = 18usize;
            let mut known_seqs: BTreeSet<Vec<u8>> =
                nodes.iter().map(|n| n.sequence.clone()).collect();
            let mut mv_counter = 0usize;

            for _iter in 0..4 {
                if mv_counter >= max_median_vectors || nodes.len() < 3 {
                    break;
                }
                let seqs: Vec<Vec<u8>> = nodes.iter().map(|n| n.sequence.clone()).collect();
                // Feasible epsilon-relaxed link graph (MSN + links within epsilon of component threshold)
                let base_links = build_msn_edges(&seqs, None);
                let max_msn_d = base_links.iter().map(|x| x.2).max().unwrap_or(0);
                let mut adj: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
                for &(u, v, _) in &base_links {
                    adj[u].push(v);
                    adj[v].push(u);
                }
                if eps > 0 {
                    for u in 0..nodes.len() {
                        for v in (u + 1)..nodes.len() {
                            let d = hamming_distance(&seqs[u], &seqs[v]);
                            if d <= max_msn_d + eps && !adj[u].contains(&v) {
                                adj[u].push(v);
                                adj[v].push(u);
                            }
                        }
                    }
                }

                // Inspect connected triplets (v - u - w) for parsimony-saving consensus median vectors
                let mut candidate_medians: Vec<(usize, Vec<u8>)> = Vec::new();
                for u in 0..nodes.len() {
                    let nbrs = &adj[u];
                    for a_i in 0..nbrs.len() {
                        for b_i in (a_i + 1)..nbrs.len() {
                            let v = nbrs[a_i];
                            let w = nbrs[b_i];
                            let med = triplet_consensus_median(&seqs[u], &seqs[v], &seqs[w]);
                            if !known_seqs.contains(&med) {
                                let d_uv = hamming_distance(&seqs[u], &seqs[v]);
                                let d_uw = hamming_distance(&seqs[u], &seqs[w]);
                                let d_vw = hamming_distance(&seqs[v], &seqs[w]);
                                let mst_cost = d_uv + d_uw + d_vw - d_uv.max(d_uw).max(d_vw);
                                let star_cost = hamming_distance(&med, &seqs[u])
                                    + hamming_distance(&med, &seqs[v])
                                    + hamming_distance(&med, &seqs[w]);
                                // Strict parsimony reduction: star tree through `med` saves mutational steps over the triplet MST
                                if star_cost < mst_cost {
                                    candidate_medians.push((star_cost, med));
                                }
                            }
                        }
                    }
                }
                if candidate_medians.is_empty() {
                    break;
                }
                candidate_medians.sort_by_key(|x| x.0);
                let min_cost = candidate_medians[0].0;
                let mut added_any = false;
                for (cost, med) in candidate_medians {
                    if cost <= min_cost + eps && !known_seqs.contains(&med) && mv_counter < max_median_vectors {
                        known_seqs.insert(med.clone());
                        mv_counter += 1;
                        let id = nodes.len();
                        nodes.push(HaplotypeNode {
                            id,
                            label: format!("mv{}", mv_counter),
                            is_median_vector: true,
                            sample_count: 0,
                            member_ids: Vec::new(),
                            pop_counts: vec![0; populations.len()],
                            sequence: med,
                            x: 0.5,
                            y: 0.5,
                        });
                        added_any = true;
                    }
                }
                if !added_any {
                    break;
                }
            }

            // Build MSN over sampled haplotypes + inferred median vectors, then prune obsolete degree <= 2 median vectors
            loop {
                let seqs: Vec<Vec<u8>> = nodes.iter().map(|n| n.sequence.clone()).collect();
                let mut links = build_msn_edges(&seqs, None);
                if eps > 0 {
                    let existing: BTreeSet<(usize, usize)> =
                        links.iter().map(|&(u, v, _)| (u.min(v), u.max(v))).collect();
                    let min_pos = links.iter().map(|x| x.2).min().unwrap_or(1);
                    for u in 0..nodes.len() {
                        for v in (u + 1)..nodes.len() {
                            if !existing.contains(&(u, v)) {
                                let d = hamming_distance(&seqs[u], &seqs[v]);
                                if d <= min_pos + eps {
                                    links.push((u, v, d));
                                }
                            }
                        }
                    }
                }
                let mut deg = vec![0usize; nodes.len()];
                for &(u, v, _) in &links {
                    deg[u] += 1;
                    deg[v] += 1;
                }
                if let Some(obsolete_idx) = nodes
                    .iter()
                    .position(|n| n.is_median_vector && deg[n.id] <= 2)
                {
                    nodes.remove(obsolete_idx);
                    for (idx, n) in nodes.iter_mut().enumerate() {
                        n.id = idx;
                    }
                } else {
                    break links;
                }
            }
        }
    };

    let edges: Vec<HaplotypeEdge> = raw_edges
        .into_iter()
        .map(|(u, v, steps)| HaplotypeEdge {
            source: u,
            target: v,
            mutational_steps: steps,
            mutations: edge_mutations(&nodes[u].sequence, &nodes[v].sequence, &retained_site_1based),
        })
        .collect();

    let mut net = HaplotypeNetwork {
        algorithm,
        epsilon,
        tcs_parsimony_limit: tcs_limit,
        populations,
        retained_site_1based,
        nodes,
        edges,
    };
    net.relax_layout(180);
    Ok(net)
}

/// Escape XML/SVG special text and attribute characters.
pub fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }
    out
}

impl HaplotypeNetwork {
    /// Deterministic circle initialization followed by Fruchterman-Reingold force-directed relaxation.
    pub fn relax_layout(&mut self, iterations: usize) {
        let n = self.nodes.len();
        if n == 0 {
            return;
        }
        if n == 1 {
            self.nodes[0].x = 0.5;
            self.nodes[0].y = 0.5;
            return;
        }
        for (i, node) in self.nodes.iter_mut().enumerate() {
            let angle = std::f64::consts::TAU * (i as f64) / (n as f64);
            node.x = 0.5 + 0.30 * angle.cos();
            node.y = 0.5 + 0.30 * angle.sin();
        }
        self.relax_layout_from_current(iterations);
    }

    /// Force-directed relaxation starting from current node coordinates, followed by
    /// centroid centering and aspect-preserving normalization into `[0.14, 0.86] x [0.16, 0.88]`.
    pub fn relax_layout_from_current(&mut self, iterations: usize) {
        let n = self.nodes.len();
        if n == 0 {
            return;
        }
        if n == 1 {
            self.nodes[0].x = 0.5;
            self.nodes[0].y = 0.5;
            return;
        }

        // Deterministically nudge exact coincident nodes so repulsion separates them cleanly.
        for i in 1..n {
            let coincident = (0..i).any(|k| {
                (self.nodes[i].x - self.nodes[k].x).abs() < 1e-9
                    && (self.nodes[i].y - self.nodes[k].y).abs() < 1e-9
            });
            if coincident {
                let angle = std::f64::consts::TAU * (i as f64) / (n as f64);
                self.nodes[i].x += 1e-4 * angle.cos();
                self.nodes[i].y += 1e-4 * angle.sin();
            }
        }

        let k_ideal = (0.45 / (n as f64).sqrt()).max(0.09);
        for iter in 0..iterations {
            let frac = iter as f64 / (iterations.max(1) as f64);
            let temp = 0.10 * (1.0 - frac).powi(2).max(0.002);
            let mut disp = vec![(0.0f64, 0.0f64); n];

            for i in 0..n {
                for j in (i + 1)..n {
                    let dx = self.nodes[i].x - self.nodes[j].x;
                    let dy = self.nodes[i].y - self.nodes[j].y;
                    let (ux, uy, dist) = if dx.abs() < 1e-9 && dy.abs() < 1e-9 {
                        let angle =
                            std::f64::consts::TAU * (((i * 7 + j * 13 + 1) % 360) as f64) / 360.0;
                        (angle.cos(), angle.sin(), 0.015)
                    } else {
                        let raw_dist = (dx * dx + dy * dy).sqrt();
                        (dx / raw_dist, dy / raw_dist, raw_dist.max(0.015))
                    };
                    let rep = (k_ideal * k_ideal) / dist;
                    let fx = ux * rep;
                    let fy = uy * rep;
                    disp[i].0 += fx;
                    disp[i].1 += fy;
                    disp[j].0 -= fx;
                    disp[j].1 -= fy;
                }
            }

            for edge in &self.edges {
                let u = edge.source;
                let v = edge.target;
                let dx = self.nodes[u].x - self.nodes[v].x;
                let dy = self.nodes[u].y - self.nodes[v].y;
                let dist = (dx * dx + dy * dy).sqrt().max(0.01);
                let target_len = k_ideal * (0.75 + 0.22 * (edge.mutational_steps.min(6) as f64));
                let attr = (dist - target_len) * 0.65;
                let fx = (dx / dist) * attr;
                let fy = (dy / dist) * attr;
                disp[u].0 -= fx;
                disp[u].1 -= fy;
                disp[v].0 += fx;
                disp[v].1 += fy;
            }

            for i in 0..n {
                // Balanced centering gravity (keeps disconnected TCS components together without wall pinning)
                disp[i].0 += 0.12 * (0.5 - self.nodes[i].x);
                disp[i].1 += 0.12 * (0.5 - self.nodes[i].y);
                let mag = (disp[i].0 * disp[i].0 + disp[i].1 * disp[i].1).sqrt().max(1e-9);
                let step = mag.min(temp);
                self.nodes[i].x = (self.nodes[i].x + (disp[i].0 / mag) * step).clamp(-4.0, 4.0);
                self.nodes[i].y = (self.nodes[i].y + (disp[i].1 / mag) * step).clamp(-4.0, 4.0);
            }
        }

        self.normalize_layout_bounds();
    }

    fn normalize_layout_bounds(&mut self) {
        let n = self.nodes.len();
        if n == 0 {
            return;
        }
        if n == 1 {
            self.nodes[0].x = 0.50;
            self.nodes[0].y = 0.52;
            return;
        }

        let mut min_x = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        for node in &self.nodes {
            if node.x.is_finite() && node.y.is_finite() {
                min_x = min_x.min(node.x);
                max_x = max_x.max(node.x);
                min_y = min_y.min(node.y);
                max_y = max_y.max(node.y);
            }
        }
        let span_x = max_x - min_x;
        let span_y = max_y - min_y;
        let max_span = span_x.max(span_y);

        if !max_span.is_finite() || max_span < 1e-4 {
            if n == 2 {
                self.nodes[0].x = 0.32;
                self.nodes[0].y = 0.52;
                self.nodes[1].x = 0.68;
                self.nodes[1].y = 0.52;
            } else {
                for (i, node) in self.nodes.iter_mut().enumerate() {
                    let angle = std::f64::consts::TAU * (i as f64) / (n as f64);
                    node.x = 0.50 + 0.26 * angle.cos();
                    node.y = 0.52 + 0.26 * angle.sin();
                }
            }
            return;
        }

        let cx = 0.5 * (min_x + max_x);
        let cy = 0.5 * (min_y + max_y);
        let target_span = if n == 2 { 0.44 } else { 0.70 };
        let scale = (target_span / max_span).min(3.5);
        for node in &mut self.nodes {
            node.x = (0.50 + (node.x - cx) * scale).clamp(0.14, 0.86);
            node.y = (0.52 + (node.y - cy) * scale).clamp(0.16, 0.88);
        }
    }

    /// Export publication-ready standalone SVG haplotype network with mutational hatch marks and population legend.
    pub fn to_svg(&self, width: u32, height: u32) -> String {
        let palette = [
            "#2E86AB", "#E07A5F", "#3D9970", "#9B5DE5", "#F2CC8F", "#D62828", "#4EA8DE",
        ];
        let wf = width as f64;
        let hf = height as f64;
        let mut svg = String::with_capacity(8192);
        svg.push_str(&format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\">\n",
            width, height, width, height
        ));
        svg.push_str("<rect width=\"100%\" height=\"100%\" fill=\"#0F172A\"/>\n");
        svg.push_str(&format!(
            "<text x=\"24\" y=\"32\" fill=\"#F8FAFC\" font-family=\"sans-serif\" font-size=\"16\" font-weight=\"bold\">{}</text>\n",
            xml_escape(self.algorithm.label())
        ));

        // Edges + mutational step hatch marks
        for e in &self.edges {
            let nu = &self.nodes[e.source];
            let nv = &self.nodes[e.target];
            let (x1, y1) = (nu.x * wf, nu.y * hf);
            let (x2, y2) = (nv.x * wf, nv.y * hf);
            svg.push_str(&format!(
                "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#94A3B8\" stroke-width=\"2.2\"/>\n",
                x1, y1, x2, y2
            ));
            let dx = x2 - x1;
            let dy = y2 - y1;
            let len = (dx * dx + dy * dy).sqrt().max(1.0);
            let nx = -dy / len;
            let ny = dx / len;
            let ticks = e.mutational_steps.min(12);
            for k in 1..=ticks {
                let frac = 0.32 + 0.36 * ((k as f64 - 0.5) / (ticks as f64));
                let tx = x1 + frac * dx;
                let ty = y1 + frac * dy;
                svg.push_str(&format!(
                    "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#FACC15\" stroke-width=\"2\"/>\n",
                    tx - 5.5 * nx,
                    ty - 5.5 * ny,
                    tx + 5.5 * nx,
                    ty + 5.5 * ny
                ));
            }
        }

        // Nodes
        for node in &self.nodes {
            let cx = node.x * wf;
            let cy = node.y * hf;
            if node.is_median_vector {
                svg.push_str(&format!(
                    "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"6.5\" fill=\"#EF4444\" stroke=\"#F8FAFC\" stroke-width=\"1.5\"/>\n",
                    cx, cy
                ));
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" fill=\"#FCA5A5\" font-family=\"monospace\" font-size=\"11\">{}</text>\n",
                    cx + 9.0,
                    cy + 4.0,
                    xml_escape(&node.label)
                ));
            } else {
                let r = 11.0 + 5.5 * (node.sample_count as f64).sqrt();
                let active_pops: Vec<(usize, usize)> = node
                    .pop_counts
                    .iter()
                    .copied()
                    .enumerate()
                    .filter(|&(_, c)| c > 0)
                    .collect();
                if active_pops.len() <= 1 {
                    let color = palette[active_pops.first().map(|x| x.0).unwrap_or(0) % palette.len()];
                    svg.push_str(&format!(
                        "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{:.1}\" fill=\"{}\" stroke=\"#F8FAFC\" stroke-width=\"2\"/>\n",
                        cx, cy, r, color
                    ));
                } else {
                    let mut start_angle = 0.0f64;
                    for (p_idx, cnt) in active_pops {
                        let sweep = std::f64::consts::TAU * (cnt as f64) / (node.sample_count.max(1) as f64);
                        let end_angle = start_angle + sweep;
                        let x_start = cx + r * start_angle.cos();
                        let y_start = cy + r * start_angle.sin();
                        let x_end = cx + r * end_angle.cos();
                        let y_end = cy + r * end_angle.sin();
                        let large_arc = if sweep > std::f64::consts::PI { 1 } else { 0 };
                        let color = palette[p_idx % palette.len()];
                        svg.push_str(&format!(
                            "<path d=\"M {:.1} {:.1} L {:.1} {:.1} A {:.1} {:.1} 0 {} 1 {:.1} {:.1} Z\" fill=\"{}\" stroke=\"#0F172A\" stroke-width=\"1\"/>\n",
                            cx, cy, x_start, y_start, r, r, large_arc, x_end, y_end, color
                        ));
                        start_angle = end_angle;
                    }
                    svg.push_str(&format!(
                        "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{:.1}\" fill=\"none\" stroke=\"#F8FAFC\" stroke-width=\"1.8\"/>\n",
                        cx, cy, r
                    ));
                }
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"#F8FAFC\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"bold\">{} (n={})</text>\n",
                    cx,
                    cy - r - 6.0,
                    xml_escape(&node.label),
                    node.sample_count
                ));
            }
        }

        // Legend
        for (i, pop) in self.populations.iter().enumerate() {
            let color = palette[i % palette.len()];
            let ly = 60.0 + 22.0 * (i as f64);
            svg.push_str(&format!(
                "<rect x=\"24\" y=\"{:.1}\" width=\"13\" height=\"13\" rx=\"3\" fill=\"{}\"/>\n",
                ly - 11.0,
                color
            ));
            svg.push_str(&format!(
                "<text x=\"44\" y=\"{:.1}\" fill=\"#E2E8F0\" font-family=\"sans-serif\" font-size=\"12\">{}</text>\n",
                ly,
                xml_escape(pop)
            ));
        }

        svg.push_str("</svg>\n");
        svg
    }

    /// Export Gephi GEXF 1.3 XML graph format.
    pub fn to_gexf(&self) -> String {
        let mut xml = String::from(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
             <gexf xmlns=\"http://gexf.net/1.3\" version=\"1.3\">\n\
             <graph mode=\"static\" defaultedgetype=\"undirected\">\n\
             <nodes>\n",
        );
        for n in &self.nodes {
            xml.push_str(&format!(
                "  <node id=\"{}\" label=\"{} (n={})\"/>\n",
                n.id,
                xml_escape(&n.label),
                n.sample_count
            ));
        }
        xml.push_str("</nodes>\n<edges>\n");
        for (eid, e) in self.edges.iter().enumerate() {
            xml.push_str(&format!(
                "  <edge id=\"{}\" source=\"{}\" target=\"{}\" weight=\"{}\"/>\n",
                eid, e.source, e.target, e.mutational_steps
            ));
        }
        xml.push_str("</edges>\n</graph>\n</gexf>\n");
        xml
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::popgen_alignment::{builtin_alignment_presets, parse_sequence_alignment};

    #[test]
    fn xml_escape_escapes_svg_and_gexf_outputs() {
        let raw_fasta = "\
            >Seq1|A&B <North> \"Prime\" '1'\nACGTACGT\n\
            >Seq2|C&D <South>\nACGTACGA\n\
            >Seq3|A&B <North> \"Prime\" '1'\nACGTACGC\n";
        let seqs = parse_sequence_alignment(raw_fasta).expect("parse special-char FASTA");
        let mut net =
            build_haplotype_network(&seqs, NetworkAlgorithm::BandeltMedianJoining, 0)
                .expect("build network");
        net.nodes[0].label = "H<1>&\"A\"'B'".to_string();

        let svg = net.to_svg(1024, 720);
        assert!(svg.contains("H&lt;1&gt;&amp;&quot;A&quot;&apos;B&apos;"));
        assert!(svg.contains("A&amp;B &lt;North&gt; &quot;Prime&quot; &apos;1&apos;"));
        assert!(!svg.contains("<North>"));

        let gexf = net.to_gexf();
        assert!(gexf.contains("label=\"H&lt;1&gt;&amp;&quot;A&quot;&apos;B&apos; (n="));
        assert!(!gexf.contains("label=\"H<1>"));
    }

    #[test]
    fn relax_layout_prevents_wall_pinning_across_presets_and_algorithms() {
        for (_title, raw_fasta) in builtin_alignment_presets() {
            let seqs = parse_sequence_alignment(raw_fasta).expect("parse preset");
            for alg in [
                NetworkAlgorithm::BandeltMedianJoining,
                NetworkAlgorithm::MinimumSpanningNetwork,
                NetworkAlgorithm::TcsStatisticalParsimony,
            ] {
                let net = build_haplotype_network(&seqs, alg, 0).expect("build net");
                assert!(net.nodes.len() >= 2);
                let mut min_x = f64::INFINITY;
                let mut max_x = f64::NEG_INFINITY;
                let mut min_y = f64::INFINITY;
                let mut max_y = f64::NEG_INFINITY;
                for node in &net.nodes {
                    assert!(node.x.is_finite() && node.y.is_finite());
                    assert!(
                        node.x >= 0.139 && node.x <= 0.861,
                        "x={} outside [0.14, 0.86]",
                        node.x
                    );
                    assert!(
                        node.y >= 0.159 && node.y <= 0.881,
                        "y={} outside [0.16, 0.88]",
                        node.y
                    );
                    assert!(
                        (node.x - 0.10).abs() > 0.01 && (node.x - 0.90).abs() > 0.01,
                        "x={} pinned near legacy wall",
                        node.x
                    );
                    assert!(
                        (node.y - 0.10).abs() > 0.01 && (node.y - 0.90).abs() > 0.01,
                        "y={} pinned near legacy wall",
                        node.y
                    );
                    min_x = min_x.min(node.x);
                    max_x = max_x.max(node.x);
                    min_y = min_y.min(node.y);
                    max_y = max_y.max(node.y);
                }
                if net.nodes.len() >= 3 {
                    assert!((max_x - min_x).max(max_y - min_y) > 0.25);
                }
            }
        }

        // Two-node edge case
        let two_fasta = ">A|P1\nAAAA\n>B|P2\nAAAT\n";
        let two_seqs = parse_sequence_alignment(two_fasta).expect("two seqs");
        let two_net =
            build_haplotype_network(&two_seqs, NetworkAlgorithm::BandeltMedianJoining, 0)
                .expect("two net");
        assert_eq!(two_net.nodes.len(), 2);
        let dx = (two_net.nodes[0].x - two_net.nodes[1].x).abs();
        let dy = (two_net.nodes[0].y - two_net.nodes[1].y).abs();
        assert!(dx.hypot(dy) > 0.20, "two nodes collapsed");
    }

    #[test]
    fn relax_layout_from_current_separates_coincident_nodes() {
        let fasta = "\
            >H1|PopA\nAAAAAA\n\
            >H2|PopB\nCAAAAA\n\
            >H3|PopC\nACAAAA\n\
            >H4|PopD\nAACAAA\n";
        let seqs = parse_sequence_alignment(fasta).expect("parse fasta");
        let mut net = build_haplotype_network(&seqs, NetworkAlgorithm::MinimumSpanningNetwork, 0)
            .expect("build net");
        assert!(net.nodes.len() >= 4);

        for node in &mut net.nodes {
            node.x = 0.5;
            node.y = 0.5;
        }
        net.relax_layout_from_current(180);

        for i in 0..net.nodes.len() {
            let ni = &net.nodes[i];
            assert!(ni.x.is_finite() && ni.y.is_finite());
            assert!(ni.x >= 0.139 && ni.x <= 0.861, "x={} out of bounds", ni.x);
            assert!(ni.y >= 0.159 && ni.y <= 0.881, "y={} out of bounds", ni.y);
            for j in (i + 1)..net.nodes.len() {
                let nj = &net.nodes[j];
                let d = (ni.x - nj.x).hypot(ni.y - nj.y);
                assert!(
                    d > 0.12,
                    "coincident nodes {i} and {j} failed to separate: d={d}"
                );
            }
        }
    }
}
