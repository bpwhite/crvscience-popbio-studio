# PopBio Studio — Scientific, Evolutionary & Coalescent Methods (`v1.2.0`)

PopBio Studio (`GPL-3.0-or-later`) couples deterministic/stochastic ecological systems with population-genetic and coalescent inference across 19 modules.

---

## 1. Ecological & Stage-Structured Dynamical Systems (`V1`–`V3`, `I1`–`I3`, `M1`–`M3`)

1. **Rosenzweig-MacArthur Predator-Prey (`V1`)**:
   Supports Holling Type I (`a*H`), Type II (`a*H / (1 + a*T_h*H)`), and Type III (`a*H^2 / (1 + a*T_h*H^2)`) functional responses with logistic prey growth `r_H * H * (1 - H / K)` and optional Wiener demographic noise.
2. **Overcompensatory & Stage-Structured Maps (`V3`, `I1`)**:
   - **Soay Sheep (`V3`)**: Maynard Smith-Slatkin map `N_{t+1} = lambda * N_t / (1 + (a * N_t)^beta)` where `beta > 2` enters oscillatory/chaotic overcompensation alongside Atlantic storm pulses.
   - **Tribolium LPA (`I1`)**: Discrete 2-week Larvae-Pupae-Adults cannibalism system (`L_{t+1} = b*A_t*exp(-c_el*L_t - c_ea*A_t)`, `P_{t+1} = (1 - mu_l)*L_t`, `A_{t+1} = P_t*exp(-c_pa*A_t) + (1 - mu_a)*A_t`) capturing Costantino et al. (1997) bifurcations.
3. **Leslie / Lefkovitch Matrices (`T6`)**:
   Computes dominant eigenvalue `lambda_1`, intrinsic rate `r = ln(lambda_1)`, net reproductive rate `R_0`, damping ratio `rho = lambda_1 / |lambda_2|`, right/left eigenvectors `w, v`, sensitivity `S_ij = v_i * w_j / <v, w>`, and elasticity `E_ij = (A_ij / lambda_1) * S_ij` (`sum_{i,j} E_ij = 1`).

---

## 2. Evolutionary Genetics & Adaptive Landscapes (`T1`–`T5`)

- **Wright-Fisher Multi-Deme Island Model (`T3`)**: Combines binomial genetic drift (`2 * N_e` draws), genotype viability fitnesses `(w_AA, w_Aa, w_aa)`, bidirectional mutation `(mu, nu)`, and island migration `m`, reporting subpopulation `H_S`, total `H_T`, and `F_ST = (H_T - H_S) / H_T`.
- **Sewall Wright 3D Shifting Balance (`T5`)**: Demes `(x_i, y_i)` explore a smooth multi-peak surface `W_bar(x, y)` via Phase I finite-size drift variance `sigma_drift ~ 1 / sqrt(2 * N_e)`, Phase II gradient selection `delta_q ~ s * grad W_bar`, and Phase III fitness-weighted propagule migration `m`.

---

## 3. Empirical Sequence Polymorphism & Neutrality (`E1`, `E4`)

Let `n` be aligned sequences (`n >= 4`), `a_1 = sum_{k=1}^{n-1} 1/k`, `a_2 = sum_{k=1}^{n-1} 1/k^2`, `k_bar` mean pairwise differences, `S` segregating sites, `eta = sum_col (k_alleles - 1)` total mutations, and `eta_1` external singleton mutation count (`c == 1` alleles bounded by `k_alleles - 1` per valid column when `valid_n > 2`, e.g. tri-allelic `{A:6, C:1, G:1}` yields `eta = 2, eta_1 = 2`).

1. **Harpending's Mismatch Raggedness (`r`)**:
   With empirical pairwise difference frequencies `x_d` for `d = 0..=d_max` and `x_{d_max + 1} = 0`:

       r = sum_{i=1}^{d_max + 1} (x_i - x_{i-1})^2

   yielding exact `r = 1.0` when all `n` sequences are identical (`d_max = 0, x_0 = 1.0`).

2. **Fu & Li's `D*` and `F*` (Fu & Li 1993; Simonsen, Churchill & Aquadro 1995)**:
   Using `a_{n+1} = a_1 + 1/n`, `c_n = 2*(n*a_1 - 2*(n - 1)) / ((n - 1)*(n - 2))`, and `d_n = c_n + (n - 2)/(n - 1)^2 + (2/(n - 1))*(3/2 - (2*a_{n+1} - 3)/(n - 2) - 1/n)`:

       v_D* = [ (n/(n - 1))^2 * a_2 + a_1^2 * d_n - 2*n*a_1*(a_1 + 1) / (n - 1)^2 ] / (a_1^2 + a_2)
       u_D* = (n/(n - 1)) * (a_1 - n/(n - 1)) - v_D*
       D*   = [ (n/(n - 1))*eta - a_1*eta_1 ] / sqrt(u_D* * eta + v_D* * eta^2)

       A_theta(n) = (n + 1) / (3*(n - 1)) + (n - 5) / n
       B_theta^2(n) = ((n - 1)/n)^2 * d_n + 2*(n^2 + n + 3) / (9*n*(n - 1)) - 2*(4*a_1 + 2 - n) / n^2
       v_F* = B_theta^2(n) / (a_1^2 + a_2)
       u_F* = (A_theta(n) / a_1) - v_F*
       F*   = [ k_bar - ((n - 1)/n)*eta_1 ] / sqrt(u_F* * eta + v_F* * eta^2)

---

## 4. Bayesian Skyline Coalescent MCMC & Detailed Balance (`E3`)

Given ultrametric UPGMA Jukes-Cantor coalescent waiting durations `Delta u_i` and active lineages `k_i`, intervals are partitioned into `m` piecewise-constant demographic groups `(group_sizes, group_thetas)` with `theta_g = 2 * N_e * mu` under an exponential Markov smoothing prior (`Drummond et al. 2005`):

1. **Scale Move on `theta_g`**: Draw `delta ~ Normal(0, 0.32)` and propose `theta' = theta_g * exp(delta)`. Whenever `theta'` lies outside `[1e-8, 50.0]`, the move is rejected immediately without clamping, preserving detailed balance with exact Jacobian term `+delta`:

       log(alpha) = (log_posterior(theta') - log_posterior(theta_g)) + delta

2. **Group-Boundary Shift**: Randomly pick boundary `g in 0..(m - 2)` and direction; skip without incrementing `proposals` if the donor group already has size `1`. Otherwise compute valid-move counts `N_valid(cur)` and `N_valid(next)` and include the asymmetric boundary Hastings ratio:

       log(alpha) = (log_posterior(next) - log_posterior(cur)) + ln(N_valid(cur) / N_valid(next))
