# PopBio Studio (`popbio-studio` v1.2.0)

Native safe-Rust (`#![forbid(unsafe_code)]`) educational population biology, theoretical evolutionary genetics, Bayesian Skyline demographic inference, sliding-window neutrality, and empirical haplotype network desktop workstation built with `eframe`/`egui` `0.31`, a built-in hover-crosshair scientific plotting engine, and bundled SQLite WAL (`rusqlite`).

## Quick Start

```bash
# Launch the interactive Front Menu & GUI Workstation (19 Modules across 5 Groups)
cargo run --release

# Jump directly into Bayesian Skyline Demographic Expansion (E3) or Sliding-Window Scan (E4)
cargo run --release -- --module e3-bayesian-skyline
cargo run --release -- --module e4-sliding-window

# Restore a portable SHA-256 verified Workset ("Save Game")
cargo run --release -- --load-workset session.popworkset.json

# Headless CLI batch simulation -> CSV export
cargo run --release -- simulate --module v1-hare-lynx --steps 600 --out hare_lynx.csv

# Headless CLI empirical FASTA PopGen + Bayesian Skyline CSV + Sliding-Window CSV + SVG/GEXF Network
cargo run --release -- analyze-fasta examples/chromosomal_sweep_balancing_180bp.fasta \
  --bsp-groups 4 --bsp-mcmc-steps 3000 --window-bp 48 --step-bp 12 \
  --out-stats popgen.json --out-skyline-csv skyline.csv \
  --out-windows-csv windows.csv --out-svg mj_network.svg

# Run unit & integration test suite
cargo test
```

## Module Catalog (19 Modules Organized into 5 Front-Menu Groups)

### I. Vertebrate Predator-Prey & Herbivore Dynamics (`src/sim_vertebrates.rs`)
* **[V1] Snowshoe Hare (*Lepus americanus*) & Canada Lynx (*Lynx canadensis*)**: Boreal Rosenzweig-MacArthur predator-prey system with selectable Holling Type I, II, and III functional responses, demographic stochasticity, nullcline phase portraits, and 2D pursuit/evasion movement ([Stenseth et al. 1997 *PNAS*](https://pubmed.ncbi.nlm.nih.gov/9144205/); [Holling 1959 *Can. Entomol.*](https://www.cambridge.org/core/journals/canadian-entomologist/article/some-characteristics-of-simple-types-of-predation-and-parasitism1/9E1E7D2CCC314766A424680444F4EA9F)).
* **[V2] Isle Royale Moose (*Alces alces*) & Gray Wolf (*Canis lupus*)**: Ratio-dependent pack predation, balsam fir browse carrying capacity feedback, winter NAO severity slider, pack inbreeding/Allee threshold, and Lake Superior ice-bridge genetic rescue toggle ([Isle Royale Wolf-Moose Project](https://www.isleroyalewolf.org/population-dynamics)).
* **[V3] St. Kilda Soay Sheep (*Ovis aries*)**: Discrete-time Maynard Smith-Slatkin overcompensatory density regulation (`beta` above 2 produces crash-recovery cycles) coupled with stochastic North Atlantic winter gales ([Coulson et al. 2001 *Science*](https://pubmed.ncbi.nlm.nih.gov/11375487/)).

### II. Invertebrate Stage-Structure & Spatial Dynamics (`src/sim_invertebrates.rs`)
* **[I1] Flour Beetles (*Tribolium castaneum* & *Tribolium confusum*)**: Stage-structured Larvae-Pupae-Adults (`LPA`) egg/pupa cannibalism bifurcations (`mu_a`, `c_pa`) with 3-layer flour column visualization and Thomas Park two-species interspecific competition toggle ([Costantino et al. 1997 *Science*](https://pmc.ncbi.nlm.nih.gov/articles/PMC314198/)).
* **[I2] Huffaker Citrus-Tray Mites (*Eotetranychus sexmaculatus* & *Typhlodromus occidentalis*)**: 24-patch (`6x4`) orange tray metapopulation demonstrating how partial Vaseline dispersal barriers generate asynchronous spatial predator-prey waves ([Mühlbauer et al. 2020 *Ecol. Evol.*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7713957/)).
* **[I3] Pea Aphid (*Acyrthosiphon pisum*) & Parasitoid Wasp (*Aphidius ervi*)**: Nicholson-Bailey host-parasitoid interaction stabilized by Negative-Binomial spatial clumping (`k` below 1) and physical legume-stem refuges ([Holling 1959 *Can. Entomol.*](https://www.cambridge.org/core/journals/canadian-entomologist/article/some-characteristics-of-simple-types-of-predation-and-parasitism1/9E1E7D2CCC314766A424680444F4EA9F)).

### III. Microbiology Growth Phases, Chemostat & Phage (`src/sim_microbiology.rs`)
* **[M1] *Escherichia coli* K-12 Growth Phases & Monod Chemostat**: Tracks intracellular lag acclimation, Exponential (Log) binary fission, Glucose-Lactose Diauxic Shift, Stationary substrate depletion, and Senescence/Death phases in batch culture or continuous chemostat dilution `D`, paired with a darkfield dividing rod-cell microcolony canvas ([Monod 1949 *Annu. Rev. Microbiol.*](https://www.annualreviews.org/content/journals/10.1146/annurev.mi.03.100149.002103)).
* **[M2] Gause Protist Microcosms (*Paramecium caudatum*, *P. aurelia* & *Didinium nasutum*)**: Toggle between competitive exclusion and protist predator-prey dynamics with oat-medium sediment refuges and periodic 3-day immigration ([Mühlbauer et al. 2020 *Ecol. Evol.*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7713957/)).
* **[M3] Bacteriophage T4 & *Escherichia coli* Lytic Plaque Wave**: Couples Campbell susceptible/eclipse-infected/resistant/virion ODEs with a `24x24` radial agar-lawn plaque assay ([Campbell 1961 *Evolution*](https://www.scilit.com/publications/3e46b60b6e9b6b089247f59bdcbe9dff)).

### IV. Theoretical Ecology & Evolutionary Genetics (`src/sim_theoretical.rs`)
* **[T1] Levins Metapopulation & Habitat Destruction Network**: Compares deterministic patch occupancy `dp/dt = c*p*(1 - h - p) - e*p` against a stochastic 64-patch landscape graph around the critical habitat destruction threshold `h_crit = 1 - e/c` ([Levins 1969 *Bull. Entomol. Soc. Am.*](https://academic.oup.com/ae/article-abstract/15/3/237/255899)).
* **[T2] May Discrete Logistic & Ricker Bifurcation Lab**: Interactive return-map cobweb plot, period-doubling cascade orbit diagram, and running Lyapunov exponent `lambda` ([May 1976 *Nature*](https://pubmed.ncbi.nlm.nih.gov/934280/)).
* **[T3] Wright-Fisher Drift, Viability Selection, Mutation & Island `F_ST`**: Multi-deme diploid (`2N_e`) allele trajectories, relative fitnesses (`w_AA`, `w_Aa`, `w_aa`), bidirectional mutation, migration `m`, real-time `H_S`, `H_T`, `F_ST`, and a De Finetti Hardy-Weinberg ternary diagram ([Wright 1931 *Genetics*](https://pmc.ncbi.nlm.nih.gov/articles/PMC1201091/)).
* **[T4] Kingman Coalescent Genealogy & Site Frequency Spectrum**: Backward-in-time `n`-sample coalescent tree simulator with branch mutations, `T_MRCA` tracking, and observed vs. neutral `theta / i` unfolded SFS bar charts ([Kingman 1982 *Stoch. Process. Appl.*](https://ideas.repec.org/a/eee/spapps/v13y1982i3p235-248.html)).
* **[T5] Sewall Wright 3D Adaptive Fitness Landscape (`v1.2.0`)**: Interactive mouse-rotatable 3D multi-peak adaptive surface `W_bar(x, y)` modeling Wright's three-phase Shifting Balance theory—exploratory genetic drift (`1 / sqrt(2*N_e)`) across maladaptive epistatic valleys, steepest-ascent intrademic selection, and fitness-weighted interdeme dispersion `m` ([Wright 1932 Adaptive Landscape Analysis](https://www.cambridge.org/core/journals/philosophy-of-science/article/heuristic-role-of-sewall-wrights-1932-adaptive-landscape-diagram/4FB15417FA20A1FE32301B147AF28CB9)).
* **[T6] Age- & Stage-Structured Leslie / Lefkovitch Matrix Projection (`v1.2.0`)**: Projects `n_{t+1} = A * n_t` across loggerhead sea turtle (*Caretta caretta*, 7-stage Lefkovitch), songbird (4-age strict Leslie), and woodland perennial (5-stage Lefkovitch) presets, solving dominant eigenvalue `lambda_1`, intrinsic `r`, net reproductive rate `R_0`, damping ratio `rho`, right eigenvector stable stage distribution `w`, left eigenvector Fisher reproductive value `v`, sensitivity `S_ij`, and conservation elasticity `E_ij` ([Leslie 1945 *Biometrika*](https://pubmed.ncbi.nlm.nih.gov/21006835/); [Crouse, Crowder & Caswell 1987 *Ecology*](https://web.abo.fi/fak/mnf/mate/kurser/matriser/Turtlearticle.pdf)).

### V. Empirical Sequence Alignments, Demographic Skylines & Haplotype Networks (`src/popgen_alignment.rs`, `src/popgen_skyline.rs`, `src/haplotype_network.rs`)
* **[E1] Sequence Alignment Polymorphism, Neutrality & `F_ST` Suite**: Parses multi-FASTA, relaxed PHYLIP, and NEXUS DNA alignments; renders a color-coded polymorphism strip; and computes segregating sites `S`, haplotype diversity `Hd`, nucleotide diversity `pi`, Watterson's `theta_W` ([Watterson 1975 *Theor. Popul. Biol.*](https://pubmed.ncbi.nlm.nih.gov/1145509/)), Tajima's *D* with normal and Beta p-values ([Tajima 1989 *Genetics*](https://pubmed.ncbi.nlm.nih.gov/2513255/)), Fu & Li's `D*` and `F*` ([Fu & Li 1993 *Genetics*](https://pmc.ncbi.nlm.nih.gov/articles/PMC1206737/)), Harpending's mismatch raggedness `r`, pairwise linkage disequilibrium (`r^2`, `|D'|`), and multi-deme Hudson `F_ST`, `D_XY`, `D_a`, and `Nm` ([Hudson et al. 1992 *Genetics*](https://pubmed.ncbi.nlm.nih.gov/1427045/)).
* **[E2] Interactive Haplotype Network Studio**: Constructs Bandelt Median-Joining (`MJ`) networks with adjustable `epsilon` (`0..=3`) and inferred median vectors (`mv1`, `mv2`, ...) ([Bandelt, Forster & Röhl 1999 *Mol. Biol. Evol.*](https://pubmed.ncbi.nlm.nih.gov/10331250/)), Kruskal Minimum Spanning Networks (`MSN`), and TCS 95% Statistical Parsimony networks, with proportional population pie wedges, mutational hatch marks, interactive node dragging, and SVG / GEXF 1.3 export.
* **[E3] Bayesian Skyline & Coalescent Demographic Expansion Lab (`v1.2.0`)**: Reconstructs ultrametric Jukes-Cantor UPGMA coalescent intervals (`k_i`, `Delta u_i`), samples piecewise-constant `N_e(t)` trajectories and 95% credible ribbons via Metropolis-Hastings MCMC under an exponential Markov smoothing prior ([Drummond et al. 2005 *Mol. Biol. Evol.*](https://pubmed.ncbi.nlm.nih.gov/15703244/)), overlays small-sample `AICc` Generalized Skylines ([Strimmer & Pybus 2001 *Mol. Biol. Evol.*](https://pubmed.ncbi.nlm.nih.gov/11719579/)), and fits sudden demographic expansion waves (`theta_0`, `theta_1`, `tau = 2*u*t`) with parametric coalescent bootstrap `SSD` and raggedness p-values ([Rogers & Harpending 1992 *Mol. Biol. Evol.*](https://pubmed.ncbi.nlm.nih.gov/1316531/)).
* **[E4] Sliding-Window Neutrality, Sweep & Divergence Scanner (`v1.2.0`)**: Scans configurable window widths (`W` bp) and step sizes (`S` bp) across aligned sequences to plot synchronized tracks of `pi(w)`, `theta_W(w)`, `D_XY(w)`, `Tajima's D(w)`, `Fu & Li's D*(w)`, and multi-deme `Hudson F_ST(w)`, flagging putative selective sweep troughs (`D` below `-1.65`), balancing selection peaks (`D` above `+1.65`), and divergence islands ([Librado & Rozas 2009 *Bioinformatics*](https://academic.oup.com/bioinformatics/article/25/11/1451/332507)).

## Workset ("Save Game") Architecture (`src/workset.rs`)
* **SQLite WAL Checkpoints (`popbio_worksets.sqlite3`)**: Press `F5` (`Quicksave`), `F9` (`Quickload`), or create named checkpoint slots directly on the Front Menu to persist all 19 modules, RNG states, custom parameters, time-series histories, and empirical Skyline/sliding-window configurations.
* **Portable `.popworkset.json` Envelopes**: Export or import self-contained JSON worksets verified via SHA-256 checksums for classroom lab assignments or reproducible study states.
