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

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// Single aligned DNA sequence record with optional population/deme assignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignedSequence {
    pub id: String,
    pub population: String,
    pub sequence: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GapHandling {
    CompleteDeletion,
    PairwiseDeletion,
}

impl GapHandling {
    pub fn label(&self) -> &'static str {
        match self {
            Self::CompleteDeletion => "Complete Deletion (exclude any column with '-' or 'N')",
            Self::PairwiseDeletion => "Pairwise Deletion (compare shared ACGT sites per pair)",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolymorphicSite {
    pub alignment_col_1based: usize,
    pub alleles: Vec<(char, usize)>,
    pub minor_allele_freq: f64,
    pub is_singleton: bool,
    pub is_parsimony_informative: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwisePopDivergence {
    pub pop_a: String,
    pub pop_b: String,
    pub pi_within_a: f64,
    pub pi_within_b: f64,
    pub d_xy: f64,
    pub d_a: f64,
    pub hudson_fst: f64,
    pub effective_nm: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LdPairStat {
    pub site_a_1based: usize,
    pub site_b_1based: usize,
    pub distance_bp: usize,
    pub r_squared: f64,
    pub d_prime_abs: f64,
}

/// Comprehensive population genetics summary computed from an empirical sequence alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PopGenSummary {
    pub dataset_name: String,
    pub num_sequences: usize,
    pub alignment_length: usize,
    pub effective_sites: usize,
    pub segregating_sites_s: usize,
    pub parsimony_informative_sites: usize,
    pub singleton_sites_eta1: usize,
    pub total_mutations_eta: usize,
    pub num_haplotypes_k: usize,
    pub haplotype_diversity_hd: f64,
    pub mean_pairwise_diffs_k: f64,
    pub nucleotide_diversity_pi: f64,
    pub watterson_theta_per_locus: f64,
    pub watterson_theta_per_site: f64,
    pub tajimas_d: Option<f64>,
    pub tajimas_d_pval_normal: Option<f64>,
    pub tajimas_d_pval_beta: Option<f64>,
    pub fu_li_d_star: Option<f64>,
    pub fu_li_f_star: Option<f64>,
    pub mismatch_histogram: Vec<(usize, usize)>,
    pub harpending_raggedness: f64,
    pub gc_content_pct: f64,
    pub polymorphic_sites: Vec<PolymorphicSite>,
    pub pop_divergences: Vec<PairwisePopDivergence>,
    pub ld_pairs: Vec<LdPairStat>,
}

/// Parse multi-FASTA, relaxed PHYLIP, or simple NEXUS DNA alignments automatically.
/// Population tags can be supplied in headers as `>Sample1|IsleRoyale`, `>Sample1 [Pop=North]`,
/// or inferred from prefix before `_`.
pub fn parse_sequence_alignment(raw_text: &str) -> Result<Vec<AlignedSequence>, String> {
    let trimmed = raw_text.trim().trim_start_matches('\u{feff}').trim();
    if trimmed.is_empty() {
        return Err("Alignment text is empty.".to_string());
    }

    let raw_records = if trimmed.starts_with('>') {
        parse_fasta(trimmed)?
    } else if trimmed
        .get(..6)
        .is_some_and(|p| p.eq_ignore_ascii_case("#NEXUS"))
    {
        parse_nexus(trimmed)?
    } else {
        parse_phylip(trimmed)?
    };

    if raw_records.len() < 2 {
        return Err("At least 2 aligned DNA sequences are required.".to_string());
    }

    let expected_len = raw_records[0].sequence.len();
    if expected_len == 0 {
        return Err("Sequences have length 0.".to_string());
    }
    for rec in &raw_records {
        if rec.sequence.len() != expected_len {
            return Err(format!(
                "Unaligned sequence length mismatch: '{}' has {} bp vs expected {} bp.",
                rec.id,
                rec.sequence.len(),
                expected_len
            ));
        }
    }
    Ok(raw_records)
}

fn extract_id_and_pop(header: &str) -> (String, String) {
    let h = header.trim();
    if let Some((left, right)) = h.split_once('|') {
        let pop = right
            .trim()
            .strip_prefix("Pop=")
            .or_else(|| right.trim().strip_prefix("pop="))
            .unwrap_or(right.trim());
        return (
            left.trim().to_string(),
            if pop.is_empty() {
                "Default".to_string()
            } else {
                pop.to_string()
            },
        );
    }
    if let Some((_, rest)) = h.split_once("Pop=") {
        let pop: String = rest
            .chars()
            .take_while(|c| c.is_alphanumeric() || *c == '_' || *c == '-')
            .collect();
        let id = h.split_whitespace().next().unwrap_or("Seq").to_string();
        if !pop.is_empty() {
            return (id, pop);
        }
    }
    let id = h.split_whitespace().next().unwrap_or("Seq").to_string();
    let pop = id
        .split_once('_')
        .map(|(prefix, _)| prefix.to_string())
        .unwrap_or_else(|| "Deme1".to_string());
    (id, pop)
}

fn sanitize_dna(raw: &str) -> String {
    raw.chars()
        .filter(|c| !c.is_whitespace())
        .map(|c| {
            let u = c.to_ascii_uppercase();
            match u {
                'U' => 'T',
                'A' | 'C' | 'G' | 'T' | '-' | 'N' => u,
                '.' | '?' => 'N',
                _ => 'N',
            }
        })
        .collect()
}

fn parse_fasta(text: &str) -> Result<Vec<AlignedSequence>, String> {
    let mut out = Vec::new();
    let mut cur_header: Option<String> = None;
    let mut cur_seq = String::new();

    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with(';') {
            continue;
        }
        if let Some(rest) = t.strip_prefix('>') {
            if let Some(h) = cur_header.take() {
                let (id, pop) = extract_id_and_pop(&h);
                out.push(AlignedSequence {
                    id,
                    population: pop,
                    sequence: sanitize_dna(&cur_seq),
                });
                cur_seq.clear();
            }
            cur_header = Some(rest.to_string());
        } else {
            cur_seq.push_str(t);
        }
    }
    if let Some(h) = cur_header {
        let (id, pop) = extract_id_and_pop(&h);
        out.push(AlignedSequence {
            id,
            population: pop,
            sequence: sanitize_dna(&cur_seq),
        });
    }
    Ok(out)
}

fn parse_phylip(text: &str) -> Result<Vec<AlignedSequence>, String> {
    let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());
    let header = lines
        .next()
        .ok_or_else(|| "Missing PHYLIP header line.".to_string())?;
    let mut parts = header.split_whitespace();
    let n_tax: usize = parts
        .next()
        .ok_or_else(|| "Invalid PHYLIP taxa count.".to_string())?
        .parse()
        .map_err(|_| "Invalid PHYLIP integer taxa count.".to_string())?;
    let n_sites: usize = parts
        .next()
        .ok_or_else(|| "Invalid PHYLIP site count.".to_string())?
        .parse()
        .map_err(|_| "Invalid PHYLIP integer site count.".to_string())?;

    if n_tax < 2 {
        return Err("PHYLIP taxa count must be at least 2.".to_string());
    }
    if n_sites == 0 {
        return Err("PHYLIP site count must be at least 1.".to_string());
    }

    let mut out: Vec<AlignedSequence> = Vec::with_capacity(n_tax.min(4096));
    for (idx, line) in lines.enumerate() {
        if idx < n_tax {
            let mut ws = line.split_whitespace();
            let name = ws
                .next()
                .ok_or_else(|| "Empty PHYLIP taxon row.".to_string())?;
            let seq_part: String = ws.collect();
            let (id, pop) = extract_id_and_pop(name);
            out.push(AlignedSequence {
                id,
                population: pop,
                sequence: sanitize_dna(&seq_part),
            });
        } else {
            // Interleaved continuation block
            let target = idx % n_tax;
            out[target].sequence.push_str(&sanitize_dna(line));
        }
    }
    Ok(out)
}

fn parse_nexus(text: &str) -> Result<Vec<AlignedSequence>, String> {
    let mut in_matrix = false;
    let mut concat: BTreeMap<String, (usize, String)> = BTreeMap::new();
    let mut order = 0usize;

    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('[') {
            continue;
        }
        let upper = t.to_ascii_uppercase();
        if upper.starts_with("MATRIX") {
            in_matrix = true;
            continue;
        }
        if in_matrix {
            if t.starts_with(';') || upper.starts_with("END;") || upper == "END" {
                break;
            }
            let clean = t.trim_end_matches(';');
            let mut ws = clean.split_whitespace();
            if let Some(name) = ws.next() {
                let seq_chunk: String = ws.collect();
                let entry = concat.entry(name.to_string()).or_insert_with(|| {
                    let o = order;
                    order += 1;
                    (o, String::new())
                });
                entry.1.push_str(&sanitize_dna(&seq_chunk));
            }
        }
    }
    let mut ordered: Vec<(usize, String, String)> = concat
        .into_iter()
        .map(|(k, (o, s))| (o, k, s))
        .collect();
    ordered.sort_by_key(|x| x.0);
    Ok(ordered
        .into_iter()
        .map(|(_, name, seq)| {
            let (id, pop) = extract_id_and_pop(&name);
            AlignedSequence {
                id,
                population: pop,
                sequence: seq,
            }
        })
        .collect())
}

#[inline]
pub fn is_valid_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Standard normal two-tailed p-value via Abramowitz & Stegun erfc approximation.
fn normal_two_tailed_pvalue(z: f64) -> f64 {
    let x = z.abs() / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let erfc_approx = poly * (-x * x).exp();
    erfc_approx.clamp(0.0, 1.0)
}

/// Tajima (1989) rescaled Beta distribution two-tailed p-value via numerical Simpson integration.
fn tajima_beta_two_tailed_pvalue(d_val: f64, n: usize) -> Option<f64> {
    if n < 4 || !d_val.is_finite() {
        return None;
    }
    let nf = n as f64;
    let a1: f64 = (1..n).map(|i| 1.0 / (i as f64)).sum();
    let a2: f64 = (1..n).map(|i| 1.0 / ((i * i) as f64)).sum();
    // Tajima (1989) Eqs. 41-47 theoretical D_min and D_max bounds
    let d_min = (2.0 / nf - 1.0 / a1)
        / (((nf + 1.0) / (3.0 * (nf - 1.0)) - 1.0 / a1) / a1).max(1e-9).sqrt();
    let d_max = ((nf + 1.0) / (2.0 * nf) - 1.0 / a1)
        / ((a2 / (a1 * a1) - 2.0 / (nf * a1) + (nf + 1.0) / (nf * (nf - 1.0) * a1)).max(1e-9).sqrt())
            .max(0.8);
    if d_max <= d_min {
        return None;
    }
    let alpha = -(1.0 + d_min * d_max) * d_max / (d_max - d_min);
    let beta = (1.0 + d_min * d_max) * d_min / (d_max - d_min);
    if alpha <= 0.5 || beta <= 0.5 {
        return Some(normal_two_tailed_pvalue(d_val));
    }
    let x_obs = ((d_val - d_min) / (d_max - d_min)).clamp(0.0001, 0.9999);
    // Numerical integration of x^(alpha-1) * (1 - x)^(beta-1) over 400 bins
    let steps = 400usize;
    let mut total_mass = 0.0;
    let mut cdf_obs = 0.0;
    for i in 0..steps {
        let u = (i as f64 + 0.5) / (steps as f64);
        let dens = ((alpha - 1.0) * u.ln() + (beta - 1.0) * (1.0 - u).ln()).exp();
        total_mass += dens;
        if u <= x_obs {
            cdf_obs += dens;
        }
    }
    if total_mass <= 0.0 {
        return Some(normal_two_tailed_pvalue(d_val));
    }
    let cdf = (cdf_obs / total_mass).clamp(0.0, 1.0);
    Some((2.0 * cdf.min(1.0 - cdf)).clamp(0.0, 1.0))
}

/// Compute comprehensive population genetic summary statistics on an aligned DNA dataset.
pub fn analyze_alignment(
    dataset_name: &str,
    sequences: &[AlignedSequence],
    gap_mode: GapHandling,
) -> Result<PopGenSummary, String> {
    let n = sequences.len();
    if n < 2 {
        return Err("Need at least 2 aligned sequences.".to_string());
    }
    let l_total = sequences[0].sequence.len();
    if l_total == 0 {
        return Err("Aligned sequence length is zero.".to_string());
    }
    let bytes: Vec<&[u8]> = sequences.iter().map(|s| s.sequence.as_bytes()).collect();

    // GC content across all unambiguous ACGT bases
    let mut gc_count = 0usize;
    let mut acgt_count = 0usize;
    for row in &bytes {
        for &b in row.iter() {
            if is_valid_base(b) {
                acgt_count += 1;
                if b == b'G' || b == b'C' {
                    gc_count += 1;
                }
            }
        }
    }
    let gc_content_pct = if acgt_count > 0 {
        100.0 * (gc_count as f64) / (acgt_count as f64)
    } else {
        0.0
    };

    // Determine retained site columns & polymorphic site metadata
    let mut complete_cols = Vec::with_capacity(l_total);
    let mut polymorphic_sites = Vec::new();
    let mut segregating_s = 0usize;
    let mut parsimony_inf = 0usize;
    let mut singletons_eta1 = 0usize;
    let mut total_mutations_eta = 0usize;

    for col in 0..l_total {
        let mut has_missing = false;
        let mut counts = [0usize; 4]; // A, C, G, T
        for row in &bytes {
            match row[col] {
                b'A' => counts[0] += 1,
                b'C' => counts[1] += 1,
                b'G' => counts[2] += 1,
                b'T' => counts[3] += 1,
                _ => has_missing = true,
            }
        }
        if !has_missing {
            complete_cols.push(col);
        }
        if gap_mode == GapHandling::CompleteDeletion && has_missing {
            continue;
        }
        let valid_n: usize = counts.iter().sum();
        if valid_n < 2 {
            continue;
        }
        let mut alleles: Vec<(char, usize)> = [('A', counts[0]), ('C', counts[1]), ('G', counts[2]), ('T', counts[3])]
            .into_iter()
            .filter(|&(_, c)| c > 0)
            .collect();
        alleles.sort_by(|a, b| b.1.cmp(&a.1));

        if alleles.len() >= 2 {
            segregating_s += 1;
            total_mutations_eta += alleles.len() - 1;
            let minor_count = alleles.last().map(|x| x.1).unwrap_or(0);
            let second_count = alleles.get(1).map(|x| x.1).unwrap_or(0);
            let is_sing = alleles.iter().any(|&(_, c)| c == 1);
            let is_pars = second_count >= 2;
            if is_sing {
                singletons_eta1 += 1;
            }
            if is_pars {
                parsimony_inf += 1;
            }
            polymorphic_sites.push(PolymorphicSite {
                alignment_col_1based: col + 1,
                alleles,
                minor_allele_freq: (minor_count as f64) / (valid_n as f64),
                is_singleton: is_sing,
                is_parsimony_informative: is_pars,
            });
        }
    }

    let effective_sites = match gap_mode {
        GapHandling::CompleteDeletion => complete_cols.len(),
        GapHandling::PairwiseDeletion => l_total,
    };

    // Collapse unique haplotypes across complete columns (or full alignment replacing gaps consistently)
    let mut hap_counts: BTreeMap<Vec<u8>, usize> = BTreeMap::new();
    for row in &bytes {
        let key: Vec<u8> = if !complete_cols.is_empty() {
            complete_cols.iter().map(|&c| row[c]).collect()
        } else {
            row.to_vec()
        };
        *hap_counts.entry(key).or_insert(0) += 1;
    }
    let num_haplotypes_k = hap_counts.len();
    let nf = n as f64;
    let sum_p2: f64 = hap_counts
        .values()
        .map(|&c| {
            let p = (c as f64) / nf;
            p * p
        })
        .sum();
    let haplotype_diversity_hd = (nf / (nf - 1.0)) * (1.0 - sum_p2);

    // Pairwise differences, nucleotide diversity pi, and mismatch distribution
    let n_pairs = n * (n - 1) / 2;
    let mut sum_diffs = 0.0;
    let mut sum_pi_pair = 0.0;
    let mut diff_hist_map: BTreeMap<usize, usize> = BTreeMap::new();

    for i in 0..n {
        for j in (i + 1)..n {
            let mut diffs = 0usize;
            let mut compared = 0usize;
            match gap_mode {
                GapHandling::CompleteDeletion => {
                    for &c in &complete_cols {
                        compared += 1;
                        if bytes[i][c] != bytes[j][c] {
                            diffs += 1;
                        }
                    }
                }
                GapHandling::PairwiseDeletion => {
                    for c in 0..l_total {
                        let bi = bytes[i][c];
                        let bj = bytes[j][c];
                        if is_valid_base(bi) && is_valid_base(bj) {
                            compared += 1;
                            if bi != bj {
                                diffs += 1;
                            }
                        }
                    }
                }
            }
            sum_diffs += diffs as f64;
            if compared > 0 {
                sum_pi_pair += (diffs as f64) / (compared as f64);
            }
            *diff_hist_map.entry(diffs).or_insert(0) += 1;
        }
    }

    let mean_pairwise_diffs_k = sum_diffs / (n_pairs as f64);
    let nucleotide_diversity_pi = sum_pi_pair / (n_pairs as f64);

    // Harpending's (1994) mismatch raggedness index r = sum_{i=1}^{d_max+1} (x_i - x_{i-1})^2
    let max_d = diff_hist_map.keys().copied().max().unwrap_or(0);
    let mut mismatch_histogram = Vec::with_capacity(max_d + 1);
    let mut prev_freq = 0.0;
    let mut raggedness = 0.0;
    for d in 0..=max_d {
        let count = diff_hist_map.get(&d).copied().unwrap_or(0);
        mismatch_histogram.push((d, count));
        let freq = (count as f64) / (n_pairs as f64);
        let delta = freq - prev_freq;
        raggedness += delta * delta;
        prev_freq = freq;
    }
    raggedness += prev_freq * prev_freq;

    // Watterson's theta & Tajima's D (1989)
    let a1: f64 = (1..n).map(|i| 1.0 / (i as f64)).sum();
    let a2: f64 = (1..n).map(|i| 1.0 / ((i * i) as f64)).sum();
    let s_f = segregating_s as f64;
    let theta_w_locus = if a1 > 0.0 { s_f / a1 } else { 0.0 };
    let theta_w_site = if effective_sites > 0 {
        theta_w_locus / (effective_sites as f64)
    } else {
        0.0
    };

    let (tajimas_d, tajimas_d_pval_normal, tajimas_d_pval_beta) = if n >= 4 && segregating_s > 0 {
        let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
        let b2 = (2.0 * (nf * nf + nf + 3.0)) / (9.0 * nf * (nf - 1.0));
        let c1 = b1 - 1.0 / a1;
        let c2 = b2 - (nf + 2.0) / (a1 * nf) + a2 / (a1 * a1);
        let e1 = c1 / a1;
        let e2 = c2 / (a1 * a1 + a2);
        let var_d = (e1 * s_f + e2 * s_f * (s_f - 1.0)).max(1e-12);
        let d_val = (mean_pairwise_diffs_k - theta_w_locus) / var_d.sqrt();
        (
            Some(d_val),
            Some(normal_two_tailed_pvalue(d_val)),
            tajima_beta_two_tailed_pvalue(d_val, n),
        )
    } else {
        (None, None, None)
    };

    // Fu & Li's (1993) D* and F* (unrooted intraspecific singleton neutrality tests)
    let (fu_li_d_star, fu_li_f_star) = if n >= 4 && segregating_s > 0 {
        let a_n_plus_1 = a1 + 1.0 / nf;
        let c_n = if n == 2 {
            1.0
        } else {
            2.0 * (nf * a1 - 2.0 * (nf - 1.0)) / ((nf - 1.0) * (nf - 2.0))
        };
        let d_n = c_n
            + (nf - 2.0) / ((nf - 1.0) * (nf - 1.0))
            + (2.0 / (nf - 1.0)) * (1.5 - (2.0 * a_n_plus_1 - 3.0) / (nf - 2.0) - 1.0 / nf);
        let v_d_star = (((nf / (nf - 1.0)).powi(2)) * a2
            + a1 * a1 * d_n
            - 2.0 * nf * a1 * (a1 + 1.0) / ((nf - 1.0) * (nf - 1.0)))
            / (a1 * a1 + a2);
        let u_d_star = (nf / (nf - 1.0)) * (a1 - nf / (nf - 1.0)) - v_d_star;

        let eta = s_f;
        let eta_s = singletons_eta1 as f64;
        let denom_d = (u_d_star * eta + v_d_star * eta * eta).max(1e-12).sqrt();
        let d_star = ((nf / (nf - 1.0)) * eta - a1 * eta_s) / denom_d;

        let v_f_star = (d_n
            + (2.0 * (nf * nf + nf + 3.0)) / (9.0 * nf * (nf - 1.0))
            - (2.0 / (nf - 1.0))
                * (4.0 * a2 - 6.0 + 8.0 / nf))
            / (a1 * a1 + a2);
        let u_f_star = (nf / (nf - 1.0)
            + (nf + 1.0) / (3.0 * (nf - 1.0))
            - 4.0 / (nf * (nf - 1.0))
            + ((2.0 * (nf + 1.0)) / ((nf - 1.0) * (nf - 1.0))) * (a_n_plus_1 - 2.0 * nf / (nf + 1.0)))
            / a1
            - v_f_star;
        let denom_f = (u_f_star.max(0.01) * eta + v_f_star.max(0.001) * eta * eta)
            .max(1e-12)
            .sqrt();
        let f_star = (mean_pairwise_diffs_k - ((nf - 1.0) / nf) * eta_s) / denom_f;
        (Some(d_star), Some(f_star))
    } else {
        (None, None)
    };

    // Multi-population Hudson (1992) F_ST, D_XY, D_a, and Nm
    let mut pop_groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, seq) in sequences.iter().enumerate() {
        pop_groups.entry(seq.population.clone()).or_default().push(idx);
    }
    let pop_names: Vec<String> = pop_groups.keys().cloned().collect();
    let mut pop_divergences = Vec::new();

    let pair_site_pi = |i: usize, j: usize| -> f64 {
        let mut diffs = 0usize;
        let mut comp = 0usize;
        for c in 0..l_total {
            let bi = bytes[i][c];
            let bj = bytes[j][c];
            if is_valid_base(bi) && is_valid_base(bj) {
                comp += 1;
                if bi != bj {
                    diffs += 1;
                }
            }
        }
        if comp > 0 {
            (diffs as f64) / (comp as f64)
        } else {
            0.0
        }
    };

    let within_pi = |idxs: &[usize]| -> f64 {
        if idxs.len() < 2 {
            return 0.0;
        }
        let mut s = 0.0;
        let mut cnt = 0usize;
        for a in 0..idxs.len() {
            for b in (a + 1)..idxs.len() {
                s += pair_site_pi(idxs[a], idxs[b]);
                cnt += 1;
            }
        }
        s / (cnt.max(1) as f64)
    };

    for a in 0..pop_names.len() {
        for b in (a + 1)..pop_names.len() {
            let idxs_a = &pop_groups[&pop_names[a]];
            let idxs_b = &pop_groups[&pop_names[b]];
            let pi_a = within_pi(idxs_a);
            let pi_b = within_pi(idxs_b);
            let mut sum_xy = 0.0;
            let mut cnt_xy = 0usize;
            for &ia in idxs_a {
                for &ib in idxs_b {
                    sum_xy += pair_site_pi(ia, ib);
                    cnt_xy += 1;
                }
            }
            let d_xy = sum_xy / (cnt_xy.max(1) as f64);
            let pi_within = 0.5 * (pi_a + pi_b);
            let d_a = (d_xy - pi_within).max(0.0);
            // Hudson, Slatkin & Maddison (1992) Eq. 3: F_ST = 1 - H_w / H_b = (D_XY - pi_within) / D_XY
            let fst = if d_xy > 1e-9 {
                ((d_xy - pi_within) / d_xy).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let nm = if fst > 1e-4 && fst < 0.9999 {
                (1.0 - fst) / (4.0 * fst)
            } else if fst <= 1e-4 {
                999.9
            } else {
                0.0
            };
            pop_divergences.push(PairwisePopDivergence {
                pop_a: pop_names[a].clone(),
                pop_b: pop_names[b].clone(),
                pi_within_a: pi_a,
                pi_within_b: pi_b,
                d_xy,
                d_a,
                hudson_fst: fst,
                effective_nm: nm,
            });
        }
    }

    // Linkage disequilibrium r^2 & |D'| across up to first 24 biallelic non-singleton sites
    let informative_cols: Vec<&PolymorphicSite> = polymorphic_sites
        .iter()
        .filter(|p| p.alleles.len() == 2 && p.is_parsimony_informative)
        .take(24)
        .collect();
    let mut ld_pairs = Vec::new();
    for i in 0..informative_cols.len() {
        for j in (i + 1)..informative_cols.len() {
            let ca = informative_cols[i].alignment_col_1based - 1;
            let cb = informative_cols[j].alignment_col_1based - 1;
            let maj_a = informative_cols[i].alleles[0].0 as u8;
            let maj_b = informative_cols[j].alleles[0].0 as u8;

            let mut valid = 0usize;
            let mut n_a = 0usize;
            let mut n_b = 0usize;
            let mut n_ab = 0usize;
            for row in &bytes {
                if is_valid_base(row[ca]) && is_valid_base(row[cb]) {
                    valid += 1;
                    let ia = row[ca] == maj_a;
                    let ib = row[cb] == maj_b;
                    if ia {
                        n_a += 1;
                    }
                    if ib {
                        n_b += 1;
                    }
                    if ia && ib {
                        n_ab += 1;
                    }
                }
            }
            if valid >= 4 {
                let vf = valid as f64;
                let pa = n_a as f64 / vf;
                let pb = n_b as f64 / vf;
                let pab = n_ab as f64 / vf;
                let d_ld = pab - pa * pb;
                let denom = pa * (1.0 - pa) * pb * (1.0 - pb);
                let r2 = if denom > 1e-9 {
                    ((d_ld * d_ld) / denom).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let d_max = if d_ld >= 0.0 {
                    (pa * (1.0 - pb)).min((1.0 - pa) * pb)
                } else {
                    (pa * pb).min((1.0 - pa) * (1.0 - pb))
                };
                let d_prime = if d_max > 1e-9 {
                    (d_ld.abs() / d_max).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                ld_pairs.push(LdPairStat {
                    site_a_1based: ca + 1,
                    site_b_1based: cb + 1,
                    distance_bp: cb.saturating_sub(ca),
                    r_squared: r2,
                    d_prime_abs: d_prime,
                });
            }
        }
    }

    Ok(PopGenSummary {
        dataset_name: dataset_name.to_string(),
        num_sequences: n,
        alignment_length: l_total,
        effective_sites,
        segregating_sites_s: segregating_s,
        parsimony_informative_sites: parsimony_inf,
        singleton_sites_eta1: singletons_eta1,
        total_mutations_eta,
        num_haplotypes_k,
        haplotype_diversity_hd,
        mean_pairwise_diffs_k,
        nucleotide_diversity_pi,
        watterson_theta_per_locus: theta_w_locus,
        watterson_theta_per_site: theta_w_site,
        tajimas_d,
        tajimas_d_pval_normal,
        tajimas_d_pval_beta,
        fu_li_d_star,
        fu_li_f_star,
        mismatch_histogram,
        harpending_raggedness: raggedness,
        gc_content_pct,
        polymorphic_sites,
        pop_divergences,
        ld_pairs,
    })
}

/// Built-in benchmark DNA alignments spanning Vertebrate mtDNA, Invertebrate COI island clades,
/// and Microbial rpoB experimental evolution so users can test PopGen & Haplotype Networks immediately.
pub fn builtin_alignment_presets() -> [(&'static str, &'static str); 4] {
    [
        (
            "Canis lupus mtDNA Control Region (Isle Royale vs. Mainland Superior)",
            ">IR01|IsleRoyale\nATGCCGTAACCTGATCGGATCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTT\n\
             >IR02|IsleRoyale\nATGCCGTAACCTGATCGGATCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTT\n\
             >IR03|IsleRoyale\nATGCCGTAACCTGATCGGATCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTT\n\
             >IR04|IsleRoyale\nATGCCGTAACCTGATCGGATCCGTAACCTAGGCTTAGCCATGACCGCTAGGCCTAACGTT\n\
             >IR05|IsleRoyale\nATGCCGTAACCTGATCGGATCCGTAACCTAGGCTTAGCCATGACCGCTAGGCCTAACGTT\n\
             >MN01|MainlandMN\nATGCCGTAACCTGATCGGATTCGTAACCTAGGCTTAGCCATGACCGCTAGGCCTAACGTT\n\
             >MN02|MainlandMN\nATGCCGTAACCTGATCGGATTCGTAACCTAGGCTTAGCCATGACCGCTAGGCCTAACGCT\n\
             >MN03|MainlandMN\nATGCCGTAACCTGATCGGATTCGTAACCTAGGCTTAGCCATGACCGCTAGGCCTAACGCT\n\
             >MN04|MainlandMN\nATGCCGTAACCTGATCGGATTCGTAACCTAGGCTCAGCCATGACCGCTAGGCCTAACGCT\n\
             >ON01|OntarioNorth\nATGCCGCAACCTGATCGGATTCGTAACCTAGGCTCAGCCATGACCGCTAGGCCTAACGCT\n\
             >ON02|OntarioNorth\nATGCCGCAACCTGATCGGATTCGTAACCTAGGCTCAGCCATGACCGCTAGGCCCAACGCT\n\
             >ON03|OntarioNorth\nATGCCGCAACCTGATCGGATTCGTAACCTAGGCTCAGCCATGACCGCTAGGCCCAACGCT\n\
             >ON04|OntarioNorth\nATGCCGCAACCTGATCGGATTCGTAACCTAGGCTCAGCCATGACCGCTAGGCCCAACGTT\n\
             >ON05|OntarioNorth\nATGCCGCAACCTGATCGGATCCGTAACCTAGGCTCAGCCATGACCGTTAGGCCCAACGTT",
        ),
        (
            "Tribolium castaneum COI Archipelago Structured Demes (3 Islands)",
            ">North01|NorthIsland\nGGAACCTTATATTTTATTTTTGGAGCATGAGCAGGAATAGTAGGAACTTCTTTAAGACTT\n\
             >North02|NorthIsland\nGGAACCTTATATTTTATTTTTGGAGCATGAGCAGGAATAGTAGGAACTTCTTTAAGACTT\n\
             >North03|NorthIsland\nGGAACCTTATATTTTATTTTTGGAGCATGAGCAGGAATAGTAGGAACTTCTTTAAGACTC\n\
             >North04|NorthIsland\nGGAACCTTATATTTTATTTTTGGAGCATGGGCAGGAATAGTAGGAACTTCTTTAAGACTC\n\
             >Central01|CentralAtoll\nGGAACCTTATATTTTATTTTTGGAGCATGGGCAGGAATAGTAGGAACTTCTTTAAGACTC\n\
             >Central02|CentralAtoll\nGGAACCTTATATTTTATTTTTGGAGCATGGGCAGGAATAGTGGGAACTTCTTTAAGACTC\n\
             >Central03|CentralAtoll\nGGAACCTTATATTTTATTTTTGGAGCATGGGCAGGAATAGTGGGAACTTCTTTAAGACTC\n\
             >Central04|CentralAtoll\nGGAACCTTATATTTTATTTTTGGAGCATGGGCAGGAATAGTGGGAACCTCTTTAAGACTC\n\
             >South01|SouthReef\nGGAACCTTATACTTTATTTTTGGAGCATGGGCAGGAATAGTGGGAACCTCTTTAAGACTC\n\
             >South02|SouthReef\nGGAACCTTATACTTTATTTTTGGAGCATGGGCAGGAATAGTGGGAACCTCTTTAAGTCTC\n\
             >South03|SouthReef\nGGAACCTTATACTTTATTTTTGGAGCATGGGCAGGAATAGTGGGAACCTCTTTAAGTCTC\n\
             >South04|SouthReef\nGGAACCTTATACTTTATTTTTGGAGCATGGGCAGGAATAGTGGGAACCTCTTTAAGTCTC\n\
             >South05|SouthReef\nGGAACCTTATACTTTATTTTTGGAGCATGAGCAGGAATAGTGGGAACCTCTTTAAGTCTT",
        ),
        (
            "Escherichia coli rpoB Star-Phylogeny Selective Sweep Sample",
            ">Clone01|ChemostatA\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n\
             >Clone02|ChemostatA\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n\
             >Clone03|ChemostatA\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n\
             >Clone04|ChemostatA\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n\
             >Clone05|ChemostatA\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n\
             >Clone06|ChemostatA\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTG\n\
             >Clone07|ChemostatB\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n\
             >Clone08|ChemostatB\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGTTA\n\
             >Clone09|ChemostatB\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTGGCTA\n\
             >Clone10|ChemostatB\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTCGCTAGCTA\n\
             >Clone11|ChemostatB\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTTGCTAGCTAGCTA\n\
             >Clone12|ChemostatB\nCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAAGCTAGCTAGCTA",
        ),
        (
            "Chromosomal 180 bp Locus: Selective Sweep (bp 1-60) & Balancing F_ST Peak (bp 121-180)",
            ">Sample01|ValleyNorth\nATGAGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAACCTTATATTTTATTTTTGGAGCATGAGCAGGAATAGTAGGAACTTCTTTAAGACTT\n\
             >Sample02|ValleyNorth\nATGCGTACCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAACCTTATATTTTATTTTTGGAGCATGAGCAGGAATAGTAGGAACTTCTTTAAGACTT\n\
             >Sample03|ValleyNorth\nATGCGTACGTAGCCAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTGAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAACCTTATATTTTATTTTTGGAGCATGAGCAGGAATAGTAGGAACTTCTTTAAGACTT\n\
             >Sample04|ValleyNorth\nATGCGTACGTAGCTAGCTCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTGAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAACCTTATATTTTATTTTTGGAGCATGAGCAGGAATAGTAGGAACTTCTTTAAGACTT\n\
             >Sample05|ValleyNorth\nATGCGTACGTAGCTAGCTAGCTACCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTGAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAACCTTATATTTTATTTTTGGAGCATGAGCAGGAATAGTAGGAACTTCTTTAAGACTT\n\
             >Sample06|ValleyNorth\nATGCGTACGTAGCTAGCTAGCTAGCTAGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAACCTTATATTTTATTTTTGGAGCATGAGCAGGAATAGTAGGAACTTCTTTAAGACTT\n\
             >Sample07|RidgeSouth\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCCAGCTAGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAAGCTTATCTTTTACTTTTGCAGCATCAGCAGCAATAGCAGGAAGTTCTTCAAGAGTT\n\
             >Sample08|RidgeSouth\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTCGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGATATGCCGTAACCTGATCGGATGGAAGCTTATCTTTTACTTTTGCAGCATCAGCAGCAATAGCAGGAAGTTCTTCAAGAGTT\n\
             >Sample09|RidgeSouth\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTACCTAGCTAGCTAGCTAGCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGATATGCCGTAACCTGATCGGATGGAAGCTTATCTTTTACTTTTGCAGCATCAGCAGCAATAGCAGGAAGTTCTTCAAGAGTT\n\
             >Sample10|RidgeSouth\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGATAGCTAGCTAGCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAAGCTTATCTTTTACTTTTGCAGCATCAGCAGCAATAGCAGGAAGTTCTTCAAGAGTT\n\
             >Sample11|RidgeSouth\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAAGCTTATCTTTTACTTTTGCAGCATCAGCAGCAATAGCAGGAAGTTCTTCAAGAGTT\n\
             >Sample12|RidgeSouth\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCGTAACCTAGGCTTAGCCATGACCGTTAGGCCTAACGTTATGCCGTAACCTGATCGGATGGAAGCTTATCTTTTACTTTTGCAGCATCAGCAGCAATAGCAGGAAGTTCTTCAAGAGTT",
        ),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_sequence_alignment_strips_utf8_bom_across_formats() {
        let bom_fasta = "\u{feff}>S1|PopA\nACGTACGT\n>S2|PopB\nACGTACGA\n";
        let fasta_seqs = parse_sequence_alignment(bom_fasta).expect("BOM-prefixed FASTA");
        assert_eq!(fasta_seqs.len(), 2);
        assert_eq!(fasta_seqs[0].id, "S1");
        assert_eq!(fasta_seqs[0].population, "PopA");
        let summary =
            analyze_alignment("bom_fasta", &fasta_seqs, GapHandling::CompleteDeletion)
                .expect("analyze BOM FASTA");
        assert_eq!(summary.segregating_sites_s, 1);

        let bom_phylip = "\u{feff}2 8\nS1_PopA   ACGTACGT\nS2_PopB   ACGTACGA\n";
        let phylip_seqs = parse_sequence_alignment(bom_phylip).expect("BOM-prefixed PHYLIP");
        assert_eq!(phylip_seqs.len(), 2);
        assert_eq!(phylip_seqs[1].sequence, "ACGTACGA");
    }

    #[test]
    fn parse_sequence_alignment_rejects_multibyte_utf8_and_degenerate_phylip_without_panicking() {
        assert!(parse_sequence_alignment("aaa\u{1f9ec}AAAA\n").is_err());
        assert!(parse_sequence_alignment("\u{1f9ec}\u{1f9ec}\nACGT\n").is_err());

        let zero_tax = parse_sequence_alignment("0 4\nS1 ACGT\n")
            .expect_err("0-taxa PHYLIP should return Err");
        assert!(zero_tax.contains("at least 2"), "unexpected err: {zero_tax}");

        let one_tax = parse_sequence_alignment("1 4\nS1 ACGT\n")
            .expect_err("1-taxon PHYLIP should return Err");
        assert!(one_tax.contains("at least 2"), "unexpected err: {one_tax}");

        let zero_sites = parse_sequence_alignment("2 0\nS1 ACGT\nS2 ACGT\n")
            .expect_err("0-site PHYLIP should return Err");
        assert!(
            zero_sites.contains("at least 1"),
            "unexpected err: {zero_sites}"
        );
    }
}
