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

use crate::popgen_alignment::{analyze_alignment, is_valid_base, AlignedSequence, GapHandling};
use crate::rng::SimRng;

/// One coalescent interval `i` (`1..=n-1`, ordered tips-to-root) extracted from an ultrametric genealogy.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoalescentInterval {
    /// 1-based interval index from present (`0.0`) toward `T_MRCA`
    pub index: usize,
    /// Number of active ancestral lineages `k_i = n - i + 1` during this interval
    pub lineages_k: usize,
    /// Start time (substitutions/site backward from present)
    pub t_start_subs: f64,
    /// End time (substitutions/site backward from present)
    pub t_end_subs: f64,
    /// Waiting duration `Delta u_i` (substitutions/site)
    pub duration_subs: f64,
    /// Classic Skyline point estimate `theta_i = k_i * (k_i - 1) * Delta u_i / 2` (`= 2 * N_e * mu`)
    pub classic_theta: f64,
}

/// Time-grid point along the Bayesian Skyline Plot (`BSP`) reconstruction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SkylineGridPoint {
    /// Evolutionary distance backward from present (substitutions/site)
    pub time_subs_per_site: f64,
    /// Scaled calendar / generation time (`time_subs_per_site / clock_rate * gen_time`)
    pub scaled_time: f64,
    /// Posterior median `N_e` (or `theta / (2 * clock_rate)`)
    pub ne_median: f64,
    /// Lower 2.5% credible bound
    pub ne_lower_95: f64,
    /// Upper 97.5% credible bound
    pub ne_upper_95: f64,
    /// Generalized Skyline (Strimmer & Pybus 2001) step-function estimate at this time
    pub generalized_skyline_ne: f64,
}

/// Rogers & Harpending (1992) coalescent sudden demographic expansion fit on the mismatch wave.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DemographicExpansionFit {
    /// Pre-expansion scaled population size `theta_0 = 2 * N_0 * u` (locus-wide)
    pub theta_0: f64,
    /// Post-expansion scaled population size `theta_1 = 2 * N_1 * u` (locus-wide)
    pub theta_1: f64,
    /// Mutational expansion time `tau = 2 * u * t` (in locus mutational units)
    pub tau: f64,
    /// Estimated calendar/generation time since sudden expansion `t = tau / (2 * u_locus)`
    pub expansion_time_scaled: f64,
    /// Sum of Squared Deviations between observed and expected sudden-expansion mismatch wave
    pub ssd_expansion: f64,
    /// Sum of Squared Deviations under constant-population equilibrium geometric distribution
    pub ssd_constant: f64,
    /// Parametric coalescent bootstrap p-value for `SSD` (`p >= 0.05` means sudden expansion is NOT rejected)
    pub bootstrap_p_ssd: f64,
    /// Observed Harpending raggedness index `r`
    pub observed_raggedness: f64,
    /// Parametric coalescent bootstrap p-value for raggedness `r`
    pub bootstrap_p_raggedness: f64,
    /// `(diffs_d, observed_freq, expected_expansion_freq, expected_constant_freq)`
    pub wave_comparison: Vec<(usize, f64, f64, f64)>,
}

/// Combined Bayesian Skyline + Generalized Skyline + Rogers-Harpending Demographic Expansion result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SkylineAndExpansionReport {
    pub num_tips: usize,
    pub tmrca_subs_per_site: f64,
    pub tmrca_scaled: f64,
    pub bsp_groups_m: usize,
    pub mcmc_steps: usize,
    pub mcmc_acceptance_rate: f64,
    pub generalized_aicc: f64,
    pub inferred_fold_change: f64,
    pub demographic_trend: String,
    pub intervals: Vec<CoalescentInterval>,
    pub grid: Vec<SkylineGridPoint>,
    pub expansion: DemographicExpansionFit,
}

/// Single window record from the Sliding-Window Polymorphism, Neutrality & Divergence Scan.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlidingWindowRecord {
    pub window_index: usize,
    pub start_1based: usize,
    pub end_1based: usize,
    pub midpoint_bp: f64,
    pub segregating_sites: usize,
    pub pi_per_site: f64,
    pub theta_w_per_site: f64,
    pub tajimas_d: f64,
    pub fu_li_d_star: f64,
    pub mean_hudson_fst: f64,
    pub mean_d_xy: f64,
    pub regime_flag: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlidingWindowScanReport {
    pub window_bp: usize,
    pub step_bp: usize,
    pub alignment_length: usize,
    pub windows: Vec<SlidingWindowRecord>,
}

/// Pairwise Jukes-Cantor (1969) distance per site (`d_JC = -0.75 * ln(1 - 4/3 * p)`).
fn jukes_cantor_distance(a: &[u8], b: &[u8]) -> f64 {
    let mut diffs = 0usize;
    let mut comp = 0usize;
    for (&ba, &bb) in a.iter().zip(b.iter()) {
        if is_valid_base(ba) && is_valid_base(bb) {
            comp += 1;
            if ba != bb {
                diffs += 1;
            }
        }
    }
    if comp == 0 {
        return 0.0;
    }
    let p = ((diffs as f64) / (comp as f64)).min(0.74);
    let inside = (1.0 - (4.0 / 3.0) * p).max(1e-6);
    (-0.75 * inside.ln()).max(0.0)
}

/// Reconstruct an ultrametric UPGMA coalescent genealogy and return the `n - 1` coalescent intervals.
pub fn extract_ultrametric_coalescent_intervals(
    sequences: &[AlignedSequence],
    l_eff_hint: usize,
) -> Vec<CoalescentInterval> {
    let n = sequences.len();
    if n < 2 {
        return Vec::new();
    }
    let bytes: Vec<&[u8]> = sequences.iter().map(|s| s.sequence.as_bytes()).collect();
    let min_res = 1.0 / ((l_eff_hint.max(sequences[0].sequence.len()).max(50) as f64) * 8.0);

    // Active UPGMA clusters: (cluster_size, height_subs, member_tip_indices)
    let mut clusters: Vec<(usize, f64, Vec<usize>)> =
        (0..n).map(|i| (1usize, 0.0f64, vec![i])).collect();

    // Precompute tip-to-tip JC69 distance matrix
    let mut tip_d = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = jukes_cantor_distance(bytes[i], bytes[j]);
            tip_d[i][j] = d;
            tip_d[j][i] = d;
        }
    }

    let mut coal_heights = Vec::with_capacity(n - 1);
    let mut last_height = 0.0f64;

    while clusters.len() > 1 {
        let m = clusters.len();
        let mut best_pair = (0usize, 1usize);
        let mut best_dist = f64::INFINITY;
        for a in 0..m {
            for b in (a + 1)..m {
                let mut sum = 0.0;
                for &ia in &clusters[a].2 {
                    for &ib in &clusters[b].2 {
                        sum += tip_d[ia][ib];
                    }
                }
                let avg = sum / ((clusters[a].0 * clusters[b].0) as f64);
                if avg < best_dist {
                    best_dist = avg;
                    best_pair = (a, b);
                }
            }
        }
        let (a, b) = best_pair;
        // Ultrametric node height u = d_AB / 2, regularized strictly non-decreasing with minimal resolution
        let raw_h = 0.5 * best_dist;
        let node_h = raw_h.max(last_height + min_res);
        last_height = node_h;
        coal_heights.push(node_h);

        let cb = clusters.remove(b);
        let ca = &mut clusters[a];
        ca.0 += cb.0;
        ca.1 = node_h;
        ca.2.extend(cb.2);
    }

    let mut intervals = Vec::with_capacity(n - 1);
    let mut t_prev = 0.0;
    for (idx, &t_end) in coal_heights.iter().enumerate() {
        let k = n - idx;
        let dur = (t_end - t_prev).max(min_res * 0.5);
        let binom = ((k * (k - 1)) as f64) * 0.5;
        let classic_theta = (binom * dur).max(1e-8);
        intervals.push(CoalescentInterval {
            index: idx + 1,
            lineages_k: k,
            t_start_subs: t_prev,
            t_end_subs: t_end,
            duration_subs: dur,
            classic_theta,
        });
        t_prev = t_end;
    }
    intervals
}

/// Strimmer & Pybus (2001) Generalized Skyline grouping via small-sample AICc optimization.
fn compute_generalized_skyline(intervals: &[CoalescentInterval]) -> (Vec<f64>, f64) {
    let m_total = intervals.len();
    if m_total == 0 {
        return (Vec::new(), 0.0);
    }
    let max_dur = intervals
        .iter()
        .map(|iv| iv.duration_subs)
        .fold(0.0f64, f64::max);

    let mut best_thetas = intervals.iter().map(|iv| iv.classic_theta).collect::<Vec<_>>();
    let mut best_aicc = f64::INFINITY;

    for step in 0..=16 {
        let eps = max_dur * (step as f64 / 24.0);
        // Pool consecutive intervals until cumulative duration >= eps
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut cur_grp = Vec::new();
        let mut cur_dur = 0.0;
        for (i, iv) in intervals.iter().enumerate() {
            cur_grp.push(i);
            cur_dur += iv.duration_subs;
            if cur_dur >= eps || i + 1 == m_total {
                groups.push(cur_grp);
                cur_grp = Vec::new();
                cur_dur = 0.0;
            }
        }
        if groups.len() > 1 && groups.last().map(|g| g.len()).unwrap_or(0) == 1 {
            // Merge trailing tiny group if its duration < eps
            let last_idx = groups.last().unwrap()[0];
            if intervals[last_idx].duration_subs < eps {
                let popped = groups.pop().unwrap();
                groups.last_mut().unwrap().extend(popped);
            }
        }

        let k_params = groups.len();
        let mut assigned_theta = vec![0.0; m_total];
        let mut log_lik = 0.0;
        for grp in &groups {
            let c_events = grp.len() as f64;
            let weighted_wait: f64 = grp
                .iter()
                .map(|&idx| {
                    let k = intervals[idx].lineages_k as f64;
                    0.5 * k * (k - 1.0) * intervals[idx].duration_subs
                })
                .sum();
            let theta_hat = (weighted_wait / c_events).max(1e-8);
            for &idx in grp {
                assigned_theta[idx] = theta_hat;
                let k = intervals[idx].lineages_k as f64;
                let binom = 0.5 * k * (k - 1.0);
                log_lik +=
                    (binom / theta_hat).ln() - (binom * intervals[idx].duration_subs) / theta_hat;
            }
        }

        let kf = k_params as f64;
        let nf = m_total as f64;
        let penalty = if nf - kf - 1.0 > 0.5 {
            2.0 * kf + (2.0 * kf * (kf + 1.0)) / (nf - kf - 1.0)
        } else {
            2.0 * kf + 100.0
        };
        let aicc = -2.0 * log_lik + penalty;
        if aicc < best_aicc {
            best_aicc = aicc;
            best_thetas = assigned_theta;
        }
    }
    (best_thetas, best_aicc)
}

/// Drummond, Rambaut, Shapiro & Pybus (2005) piecewise-constant Bayesian Skyline log-posterior:
/// exact coalescent log-likelihood + exponential Markov smoothing prior across `m` grouped intervals.
fn bsp_log_posterior(
    intervals: &[CoalescentInterval],
    group_sizes: &[usize],
    group_thetas: &[f64],
) -> f64 {
    let mut ll = 0.0;
    let mut iv_idx = 0usize;
    for (g_idx, &g_count) in group_sizes.iter().enumerate() {
        let theta = group_thetas[g_idx].max(1e-10);
        for _ in 0..g_count {
            let iv = &intervals[iv_idx];
            let k = iv.lineages_k as f64;
            let binom = 0.5 * k * (k - 1.0);
            ll += (binom / theta).ln() - (binom * iv.duration_subs) / theta;
            iv_idx += 1;
        }
    }
    // Jeffreys scale prior `-ln(theta_0)` on present-day group + exponential Markov chain smoothing prior
    // `theta_j ~ Exp(mean = theta_{j-1})` -> `log p(theta_j | theta_{j-1}) = -ln(theta_{j-1}) - theta_j / theta_{j-1}`
    let mut lp = -group_thetas[0].max(1e-10).ln();
    for j in 1..group_thetas.len() {
        let prev = group_thetas[j - 1].max(1e-10);
        let curr = group_thetas[j].max(1e-10);
        lp += -prev.ln() - (curr / prev);
    }
    ll + lp
}

/// Rogers & Harpending (1992) theoretical mismatch probability distribution `F_0..F_{d_max}`
/// under a sudden demographic expansion from `theta_0` to `theta_1` at mutational time `tau = 2*u*t`.
pub fn rogers_harpending_wave(
    theta_0: f64,
    theta_1: f64,
    tau: f64,
    max_d: usize,
) -> Vec<f64> {
    let t0 = theta_0.max(1e-4);
    let t1 = theta_1.max(t0 + 1e-4);
    let tau_pos = tau.max(0.0);

    // Pre-expansion equilibrium geometric distribution F_i(0)
    let q0 = t0 / (1.0 + t0);
    let mut f_init = vec![0.0; max_d + 1];
    let mut pow_q = 1.0;
    for f_val in f_init.iter_mut() {
        *f_val = pow_q / (1.0 + t0);
        pow_q *= q0;
    }

    // Exact Rogers & Harpending (1992) Eq. 4 closed-form post-expansion mismatch wave:
    // F_i(tau) = F_i_hat(theta_1) + exp(-tau * (1 + 1/theta_1)) * sum_{j=0}^i (tau^j / j!) * [F_{i-j}(0) - F_{i-j}_hat(theta_1)]
    // where F_k_hat(theta_1) = (theta_1 / (1 + theta_1))^k / (1 + theta_1) is the new equilibrium toward which the wave relaxes!
    let q1 = t1 / (1.0 + t1);
    let mut f_eq1 = vec![0.0; max_d + 1];
    let mut pow_q1 = 1.0;
    for f_val in f_eq1.iter_mut() {
        *f_val = pow_q1 / (1.0 + t1);
        pow_q1 *= q1;
    }

    let rate = 1.0 + 1.0 / t1;
    let damp = (-tau_pos * rate).exp();
    // Note: in mutational time tau = 2*u*t with coalescence hazard 1/theta_1 per 2u,
    // the exact Rogers & Harpending (1992) Eq. 7 wave argument uses `tau * (1 + 1/theta_1)`
    let mut pois_term = vec![0.0; max_d + 1];
    pois_term[0] = 1.0;
    for j in 1..=max_d {
        pois_term[j] = pois_term[j - 1] * tau_pos / (j as f64);
    }

    let mut out = vec![0.0; max_d + 1];
    for i in 0..=max_d {
        let mut transient = 0.0;
        for j in 0..=i {
            transient += pois_term[j] * (f_init[i - j] - f_eq1[i - j]);
        }
        out[i] = (f_eq1[i] + damp * transient).max(0.0);
    }
    let norm: f64 = out.iter().sum();
    if norm > 1e-12 {
        for v in &mut out {
            *v /= norm;
        }
    }
    out
}

fn raggedness_from_freqs(freqs: &[f64]) -> f64 {
    let mut prev = 0.0;
    let mut rag = 0.0;
    for &f in freqs {
        let d = f - prev;
        rag += d * d;
        prev = f;
    }
    rag + prev * prev
}

/// Fit Rogers & Harpending (1992) sudden demographic expansion parameters `(theta_0, theta_1, tau)`
/// and evaluate parametric coalescent bootstrap p-values (`B = 100`).
pub fn fit_demographic_expansion(
    sequences: &[AlignedSequence],
    mismatch_histogram: &[(usize, usize)],
    clock_rate_per_site: f64,
    generation_time: f64,
) -> DemographicExpansionFit {
    let n = sequences.len();
    let l_bp = sequences.first().map(|s| s.sequence.len()).unwrap_or(1).max(1) as f64;
    let max_d = mismatch_histogram
        .iter()
        .map(|&(d, _)| d)
        .max()
        .unwrap_or(0)
        .max(8);

    let total_pairs: usize = mismatch_histogram.iter().map(|&(_, c)| c).sum::<usize>().max(1);
    let mut obs_freq = vec![0.0; max_d + 1];
    let mut mean_m = 0.0;
    for &(d, cnt) in mismatch_histogram {
        if d <= max_d {
            let f = (cnt as f64) / (total_pairs as f64);
            obs_freq[d] = f;
            mean_m += (d as f64) * f;
        }
    }
    let mut var_m = 0.0;
    for (d, &f) in obs_freq.iter().enumerate() {
        let diff = (d as f64) - mean_m;
        var_m += diff * diff * f;
    }

    // Rogers (1995) method-of-moments initial seed: theta_0_init = sqrt(max(v - m, 0)), tau_init = m - theta_0
    let mom_theta0 = (var_m - mean_m).max(0.0).sqrt().clamp(0.05, mean_m.max(0.5));
    let mom_tau = (mean_m - mom_theta0).max(0.25);

    // Constant-size geometric baseline with theta = mean_m
    let theta_const = mean_m.max(0.05);
    let q_c = theta_const / (1.0 + theta_const);
    let mut const_freq = vec![0.0; max_d + 1];
    let mut pq = 1.0;
    for val in const_freq.iter_mut() {
        *val = pq / (1.0 + theta_const);
        pq *= q_c;
    }
    let norm_c: f64 = const_freq.iter().sum::<f64>().max(1e-12);
    for val in &mut const_freq {
        *val /= norm_c;
    }
    let ssd_constant: f64 = obs_freq
        .iter()
        .zip(const_freq.iter())
        .map(|(o, e)| (o - e).powi(2))
        .sum();

    // Multi-resolution least-squares refinement of (theta_0, theta_1, tau)
    let mut best_params = (mom_theta0, (mean_m * 12.0).max(5.0), mom_tau);
    let mut best_ssd = f64::INFINITY;

    for &t0_mult in &[0.1, 0.4, 0.8, 1.0, 1.5, 2.5] {
        let cand_t0 = (mom_theta0 * t0_mult).clamp(0.01, mean_m.max(1.0) * 1.5);
        for &t1 in &[2.0, 8.0, 25.0, 80.0, 250.0] {
            if t1 <= cand_t0 {
                continue;
            }
            for step_tau in 1..=24 {
                let cand_tau = (mean_m.max(1.0) * 2.2) * (step_tau as f64 / 24.0);
                let wave = rogers_harpending_wave(cand_t0, t1, cand_tau, max_d);
                let ssd: f64 = obs_freq
                    .iter()
                    .zip(wave.iter())
                    .take(max_d + 1)
                    .map(|(o, e)| (o - e).powi(2))
                    .sum();
                if ssd < best_ssd {
                    best_ssd = ssd;
                    best_params = (cand_t0, t1, cand_tau);
                }
            }
        }
    }

    // Fine coordinate polish around best_params
    let mut step_scale = 0.25;
    for _ in 0..10 {
        let (c0, c1, ctau) = best_params;
        for &(d0, d1, dtau) in &[
            (-1.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, -1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, -1.0),
            (0.0, 0.0, 1.0),
        ] {
            let n0 = (c0 * (1.0 + step_scale * d0)).clamp(0.005, 100.0);
            let n1 = (c1 * (1.0 + step_scale * d1)).clamp(n0 + 0.05, 1000.0);
            let ntau = (ctau * (1.0 + step_scale * dtau)).clamp(0.05, (max_d as f64) * 2.5);
            let wave = rogers_harpending_wave(n0, n1, ntau, max_d);
            let ssd: f64 = obs_freq
                .iter()
                .zip(wave.iter())
                .map(|(o, e)| (o - e).powi(2))
                .sum();
            if ssd < best_ssd {
                best_ssd = ssd;
                best_params = (n0, n1, ntau);
            }
        }
        step_scale *= 0.65;
    }

    let (theta_0, theta_1, tau) = best_params;
    let fitted_wave = rogers_harpending_wave(theta_0, theta_1, tau, max_d);
    let obs_rag = raggedness_from_freqs(&obs_freq);

    // Parametric coalescent bootstrap (B = 80 replicates of an n-tip sudden-expansion coalescent tree)
    let mut rng = SimRng::new(1992_0552);
    let b_reps = 80usize;
    let mut ge_ssd = 0usize;
    let mut ge_rag = 0usize;
    let n_sim = n.clamp(4, 30);
    let sim_pairs = (n_sim * (n_sim - 1) / 2).max(1);

    for _ in 0..b_reps {
        // Simulate backward coalescent tree where theta(u) = theta_1 for 2*u*t < tau, and theta_0 before tau
        let mut active: Vec<Vec<usize>> = (0..n_sim).map(|i| vec![i]).collect();
        let mut pair_muts = vec![vec![0usize; n_sim]; n_sim];
        let mut u_mut = 0.0; // measured in tau units (2*u*t)

        while active.len() > 1 {
            let k = active.len() as f64;
            let theta_cur = if u_mut < tau { theta_1 } else { theta_0 };
            // Coalescence hazard per unit of tau (= 2*u*t): rate = k*(k - 1) / (2 * theta_cur)
            let rate = (k * (k - 1.0)) / (2.0 * theta_cur.max(0.01));
            let mut dt_tau = rng.exponential(rate);
            if u_mut < tau && u_mut + dt_tau > tau {
                // Cross expansion horizon seamlessly into ancestral smaller population size theta_0
                let rem_before = tau - u_mut;
                // Add Poisson mutations along all k branches over rem_before (branch mutation rate = 0.5 per unit tau!)
                for cluster in &active {
                    let mut m_cnt = 0usize;
                    let lim = (-0.5 * rem_before).exp();
                    let mut p = rng.next_f64();
                    while p > lim && m_cnt < 40 {
                        m_cnt += 1;
                        p *= rng.next_f64();
                    }
                    if m_cnt > 0 {
                        for &tip_in in cluster {
                            for tip_out in 0..n_sim {
                                if !cluster.contains(&tip_out) {
                                    pair_muts[tip_in][tip_out] += m_cnt;
                                }
                            }
                        }
                    }
                }
                u_mut = tau;
                let rate_anc = (k * (k - 1.0)) / (2.0 * theta_0.max(0.01));
                dt_tau = rng.exponential(rate_anc);
            }

            // Accumulate branch mutations over dt_tau
            for cluster in &active {
                let mut m_cnt = 0usize;
                let lim = (-0.5 * dt_tau.min(40.0)).exp();
                let mut p = rng.next_f64();
                while p > lim && m_cnt < 50 {
                    m_cnt += 1;
                    p *= rng.next_f64();
                }
                if m_cnt > 0 {
                    for &tip_in in cluster {
                        for tip_out in 0..n_sim {
                            if !cluster.contains(&tip_out) {
                                pair_muts[tip_in][tip_out] += m_cnt;
                            }
                        }
                    }
                }
            }
            u_mut += dt_tau;
            let ia = rng.gen_index(active.len());
            let ca = active.remove(ia);
            let ib = rng.gen_index(active.len());
            active[ib].extend(ca);
        }

        let mut sim_freq = vec![0.0; max_d + 1];
        for i in 0..n_sim {
            for j in (i + 1)..n_sim {
                let d = pair_muts[i][j].min(max_d);
                sim_freq[d] += 1.0 / (sim_pairs as f64);
            }
        }
        let sim_ssd: f64 = sim_freq
            .iter()
            .zip(fitted_wave.iter())
            .map(|(s, e)| (s - e).powi(2))
            .sum();
        let sim_rag = raggedness_from_freqs(&sim_freq);
        if sim_ssd >= best_ssd {
            ge_ssd += 1;
        }
        if sim_rag >= obs_rag {
            ge_rag += 1;
        }
    }

    let u_locus = (clock_rate_per_site.max(1e-12)) * l_bp;
    let expansion_time_scaled = (tau / (2.0 * u_locus)) * generation_time.max(1e-6);

    let wave_comparison = (0..=max_d)
        .map(|d| (d, obs_freq[d], fitted_wave[d], const_freq[d]))
        .collect();

    DemographicExpansionFit {
        theta_0,
        theta_1,
        tau,
        expansion_time_scaled,
        ssd_expansion: best_ssd,
        ssd_constant,
        bootstrap_p_ssd: (ge_ssd as f64) / (b_reps as f64),
        observed_raggedness: obs_rag,
        bootstrap_p_raggedness: (ge_rag as f64) / (b_reps as f64),
        wave_comparison,
    }
}

/// Run full Bayesian Skyline MCMC (`Drummond et al. 2005`), Generalized Skyline (`Strimmer & Pybus 2001`),
/// and Rogers & Harpending (`1992`) sudden expansion analysis on an empirical sequence alignment.
pub fn run_skyline_and_expansion(
    sequences: &[AlignedSequence],
    mismatch_histogram: &[(usize, usize)],
    bsp_groups: usize,
    mcmc_steps: usize,
    clock_rate_per_site: f64,
    generation_time: f64,
) -> Result<SkylineAndExpansionReport, String> {
    if sequences.len() < 4 {
        return Err("Bayesian Skyline inference requires at least 4 aligned sequences.".to_string());
    }
    let l_bp = sequences[0].sequence.len();
    let intervals = extract_ultrametric_coalescent_intervals(sequences, l_bp);
    let m_iv = intervals.len();
    let m_groups = bsp_groups.clamp(1, m_iv.min(10));
    let tmrca_subs = intervals.last().map(|iv| iv.t_end_subs).unwrap_or(0.0);
    let mu = clock_rate_per_site.max(1e-12);
    let gen_t = generation_time.max(1e-6);

    let (gen_thetas, gen_aicc) = compute_generalized_skyline(&intervals);

    // Initialize balanced group sizes summing to `m_iv`
    let mut group_sizes = vec![m_iv / m_groups; m_groups];
    for g in 0..(m_iv % m_groups) {
        group_sizes[g] += 1;
    }

    // Initialize group_thetas from group means of Generalized Skyline
    let mut group_thetas = vec![0.0; m_groups];
    let mut cursor = 0usize;
    for g in 0..m_groups {
        let cnt = group_sizes[g].max(1);
        let mean_th: f64 = gen_thetas[cursor..cursor + cnt].iter().sum::<f64>() / (cnt as f64);
        group_thetas[g] = mean_th.max(1e-6);
        cursor += cnt;
    }

    let mut cur_lp = bsp_log_posterior(&intervals, &group_sizes, &group_thetas);
    let mut rng = SimRng::new(2005_1185);
    let total_iters = mcmc_steps.clamp(400, 15_000);
    let burn_in = total_iters / 4;
    let sample_stride = ((total_iters - burn_in) / 200).max(1);

    let grid_n = 56usize;
    let grid_times: Vec<f64> = (0..grid_n)
        .map(|i| tmrca_subs * (i as f64 / ((grid_n - 1) as f64)))
        .collect();
    let mut posterior_samples: Vec<Vec<f64>> = vec![Vec::new(); grid_n];

    let mut proposals = 0usize;
    let mut accepted = 0usize;

    for iter in 0..total_iters {
        if m_groups == 1 || rng.next_f64() < 0.78 {
            // Log-normal scale multiplier proposal on a random group's theta_g
            let g = rng.gen_index(m_groups);
            let delta = rng.standard_normal() * 0.32;
            let u = rng.next_f64();
            proposals += 1;
            if bsp_try_scale_move(&intervals, &group_sizes, &mut group_thetas, g, delta, u, &mut cur_lp) {
                accepted += 1;
            }
        } else {
            // Shift a group boundary by +/- 1 coalescent interval between adjacent groups `g` and `g+1`
            let g = rng.gen_index(m_groups - 1);
            let move_right = rng.next_f64() < 0.5;
            let u = rng.next_f64();
            match bsp_try_boundary_move(&intervals, &mut group_sizes, &group_thetas, g, move_right, u, &mut cur_lp) {
                Some(true) => {
                    proposals += 1;
                    accepted += 1;
                }
                Some(false) => {
                    proposals += 1;
                }
                None => {} // Skipped invalid boundary move: do not count in `proposals`
            }
        }

        if iter >= burn_in && (iter - burn_in) % sample_stride == 0 {
            // Map each interval to its active group theta, then convert theta = 2 * N_e * mu -> N_e = theta / (2 * mu)
            let mut iv_theta = Vec::with_capacity(m_iv);
            for (g_idx, &cnt) in group_sizes.iter().enumerate() {
                for _ in 0..cnt {
                    iv_theta.push(group_thetas[g_idx]);
                }
            }
            for (gi, &t_q) in grid_times.iter().enumerate() {
                let iv_pos = intervals
                    .iter()
                    .position(|iv| t_q <= iv.t_end_subs)
                    .unwrap_or(m_iv - 1);
                let ne_est = iv_theta[iv_pos] / (2.0 * mu);
                posterior_samples[gi].push(ne_est);
            }
        }
    }

    let mut grid = Vec::with_capacity(grid_n);
    for (gi, &t_q) in grid_times.iter().enumerate() {
        let samples = &mut posterior_samples[gi];
        samples.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let sn = samples.len().max(1);
        let med = samples[sn / 2];
        let lo = samples[(sn * 25 / 1000).min(sn - 1)];
        let hi = samples[(sn * 975 / 1000).min(sn - 1)];
        let iv_pos = intervals
            .iter()
            .position(|iv| t_q <= iv.t_end_subs)
            .unwrap_or(m_iv - 1);
        let gen_ne = gen_thetas[iv_pos] / (2.0 * mu);

        grid.push(SkylineGridPoint {
            time_subs_per_site: t_q,
            scaled_time: (t_q / mu) * gen_t,
            ne_median: med,
            ne_lower_95: lo,
            ne_upper_95: hi,
            generalized_skyline_ne: gen_ne,
        });
    }

    let ne_present = grid.first().map(|g| g.ne_median).unwrap_or(1.0).max(1e-6);
    let ne_ancestral = grid.last().map(|g| g.ne_median).unwrap_or(1.0).max(1e-6);
    let fold_change = ne_present / ne_ancestral;
    let trend = if fold_change > 1.65 {
        format!("Demographic Expansion ({:.2}x growth from T_MRCA to present)", fold_change)
    } else if fold_change < 0.60 {
        format!("Demographic Contraction / Bottleneck ({:.2}x of ancestral N_e)", fold_change)
    } else {
        format!("Near-Constant Effective Size ({:.2}x ratio)", fold_change)
    };

    let expansion = fit_demographic_expansion(sequences, mismatch_histogram, mu, gen_t);

    Ok(SkylineAndExpansionReport {
        num_tips: sequences.len(),
        tmrca_subs_per_site: tmrca_subs,
        tmrca_scaled: (tmrca_subs / mu) * gen_t,
        bsp_groups_m: m_groups,
        mcmc_steps: total_iters,
        mcmc_acceptance_rate: (accepted as f64) / (proposals.max(1) as f64),
        generalized_aicc: gen_aicc,
        inferred_fold_change: fold_change,
        demographic_trend: trend,
        intervals,
        grid,
        expansion,
    })
}

/// Sliding-window scan across an empirical DNA alignment computing `pi`, `theta_W`,
/// `Tajima's D`, `Fu & Li's D*`, and multi-deme `Hudson F_ST` / `D_XY` along sequence coordinates.
pub fn run_sliding_window_scan(
    sequences: &[AlignedSequence],
    window_bp: usize,
    step_bp: usize,
    gap_mode: GapHandling,
) -> Result<SlidingWindowScanReport, String> {
    if sequences.len() < 4 {
        return Err("Sliding-window neutrality scan requires at least 4 sequences.".to_string());
    }
    let l_total = sequences[0].sequence.len();
    if l_total < 10 {
        return Err("Alignment must be at least 10 bp long for sliding-window scanning.".to_string());
    }
    let win = window_bp.clamp(8, l_total);
    let step = step_bp.clamp(1, win);

    let mut windows = Vec::new();
    let mut start = 0usize;
    let mut w_idx = 0usize;

    while start < l_total {
        let end = (start + win).min(l_total);
        if end - start < 6 && !windows.is_empty() {
            break;
        }
        let sliced: Vec<AlignedSequence> = sequences
            .iter()
            .map(|s| AlignedSequence {
                id: s.id.clone(),
                population: s.population.clone(),
                sequence: s.sequence[start..end].to_string(),
            })
            .collect();

        if let Ok(sum) = analyze_alignment("window", &sliced, gap_mode) {
            let td = sum.tajimas_d.unwrap_or(0.0);
            let fld = sum.fu_li_d_star.unwrap_or(0.0);
            let (mean_fst, mean_dxy) = if sum.pop_divergences.is_empty() {
                (0.0, 0.0)
            } else {
                let k = sum.pop_divergences.len() as f64;
                let sf: f64 = sum.pop_divergences.iter().map(|d| d.hudson_fst).sum();
                let sd: f64 = sum.pop_divergences.iter().map(|d| d.d_xy).sum();
                (sf / k, sd / k)
            };

            let flag = if td <= -1.65 {
                "Putative Selective Sweep / Expansion Excess Singletons".to_string()
            } else if td >= 1.65 {
                "Putative Balancing Selection / Admixture Peak".to_string()
            } else if mean_fst >= 0.55 {
                "Divergence Island (Elevated F_ST)".to_string()
            } else {
                "Consistent with Neutrality".to_string()
            };

            windows.push(SlidingWindowRecord {
                window_index: w_idx + 1,
                start_1based: start + 1,
                end_1based: end,
                midpoint_bp: 0.5 * ((start + 1 + end) as f64),
                segregating_sites: sum.segregating_sites_s,
                pi_per_site: sum.nucleotide_diversity_pi,
                theta_w_per_site: sum.watterson_theta_per_site,
                tajimas_d: td,
                fu_li_d_star: fld,
                mean_hudson_fst: mean_fst,
                mean_d_xy: mean_dxy,
                regime_flag: flag,
            });
            w_idx += 1;
        }

        if end == l_total {
            break;
        }
        start += step;
    }

    let _ = BTreeMap::<u8, u8>::new();
    Ok(SlidingWindowScanReport {
        window_bp: win,
        step_bp: step,
        alignment_length: l_total,
        windows,
    })
}

impl SkylineAndExpansionReport {
    pub fn to_csv(&self) -> String {
        let mut csv = String::from(
            "time_subs_per_site,scaled_time,bsp_ne_median,bsp_ne_lower_95,bsp_ne_upper_95,generalized_skyline_ne\n",
        );
        for pt in &self.grid {
            csv.push_str(&format!(
                "{:.8},{:.4},{:.4},{:.4},{:.4},{:.4}\n",
                pt.time_subs_per_site,
                pt.scaled_time,
                pt.ne_median,
                pt.ne_lower_95,
                pt.ne_upper_95,
                pt.generalized_skyline_ne
            ));
        }
        csv
    }
}

/// Counts valid single-interval boundary shifts across adjacent skyline groups.
pub fn bsp_count_valid_boundary_moves(sizes: &[usize]) -> usize {
    let mut cnt = 0usize;
    for idx in 0..sizes.len().saturating_sub(1) {
        if sizes[idx] > 1 {
            cnt += 1;
        }
        if sizes[idx + 1] > 1 {
            cnt += 1;
        }
    }
    cnt
}

/// Log-normal scale proposal on `group_thetas[g]` rejecting out-of-bounds draws (`[1e-8, 50.0]`)
/// instead of clamping, preserving detailed balance with exact Jacobian `+delta`.
pub fn bsp_try_scale_move(
    intervals: &[CoalescentInterval],
    group_sizes: &[usize],
    group_thetas: &mut [f64],
    g: usize,
    delta: f64,
    u: f64,
    cur_lp: &mut f64,
) -> bool {
    let old_val = group_thetas[g];
    let proposed = old_val * delta.exp();
    if !proposed.is_finite() || !(1e-8..=50.0).contains(&proposed) {
        return false;
    }
    group_thetas[g] = proposed;
    let new_lp = bsp_log_posterior(intervals, group_sizes, group_thetas);
    let log_alpha = (new_lp - *cur_lp) + delta;
    if log_alpha >= 0.0 || u.max(1e-300).ln() < log_alpha {
        *cur_lp = new_lp;
        true
    } else {
        group_thetas[g] = old_val;
        false
    }
}

/// Adjacent-group boundary shift with exact Hastings ratio `ln(N_valid(cur) / N_valid(next))`.
/// Returns `None` when the selected donor group has size 1 (skipped move, excluded from `proposals`).
pub fn bsp_try_boundary_move(
    intervals: &[CoalescentInterval],
    group_sizes: &mut [usize],
    group_thetas: &[f64],
    g: usize,
    move_right: bool,
    u: f64,
    cur_lp: &mut f64,
) -> Option<bool> {
    let (donor, recv) = if move_right { (g, g + 1) } else { (g + 1, g) };
    if group_sizes[donor] <= 1 {
        return None;
    }
    let n_valid_cur = bsp_count_valid_boundary_moves(group_sizes);
    group_sizes[donor] -= 1;
    group_sizes[recv] += 1;
    let n_valid_next = bsp_count_valid_boundary_moves(group_sizes);
    let new_lp = bsp_log_posterior(intervals, group_sizes, group_thetas);
    let log_hastings = (n_valid_cur as f64).ln() - (n_valid_next as f64).ln();
    let log_alpha = (new_lp - *cur_lp) + log_hastings;
    if log_alpha >= 0.0 || u.max(1e-300).ln() < log_alpha {
        *cur_lp = new_lp;
        Some(true)
    } else {
        group_sizes[donor] += 1;
        group_sizes[recv] -= 1;
        Some(false)
    }
}

impl SlidingWindowScanReport {
    pub fn to_csv(&self) -> String {
        let mut csv = String::from(
            "window_index,start_bp,end_bp,midpoint_bp,segregating_S,pi_per_site,theta_w_per_site,tajimas_d,fu_li_d_star,mean_hudson_fst,mean_d_xy,regime_flag\n",
        );
        for w in &self.windows {
            csv.push_str(&format!(
                "{},{},{},{:.1},{},{:.6},{:.6},{:.4},{:.4},{:.4},{:.6},\"{}\"\n",
                w.window_index,
                w.start_1based,
                w.end_1based,
                w.midpoint_bp,
                w.segregating_sites,
                w.pi_per_site,
                w.theta_w_per_site,
                w.tajimas_d,
                w.fu_li_d_star,
                w.mean_hudson_fst,
                w.mean_d_xy,
                w.regime_flag
            ));
        }
        csv
    }
}
