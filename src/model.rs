use serde::{Deserialize, Serialize};

/// Top-level taxonomic and methodological groupings displayed on the Front Menu.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ModuleGroup {
    Vertebrates,
    Invertebrates,
    Microbiology,
    Theoretical,
    EmpiricalPopGen,
}

impl ModuleGroup {
    pub const ALL: [ModuleGroup; 5] = [
        ModuleGroup::Vertebrates,
        ModuleGroup::Invertebrates,
        ModuleGroup::Microbiology,
        ModuleGroup::Theoretical,
        ModuleGroup::EmpiricalPopGen,
    ];

    pub fn title(&self) -> &'static str {
        match self {
            ModuleGroup::Vertebrates => "I. Vertebrate Predator-Prey & Herbivore Dynamics",
            ModuleGroup::Invertebrates => "II. Invertebrate Stage-Structure & Spatial Dynamics",
            ModuleGroup::Microbiology => "III. Microbiology Growth Phases, Chemostat & Phage",
            ModuleGroup::Theoretical => "IV. Theoretical Ecology & Evolutionary Genetics",
            ModuleGroup::EmpiricalPopGen => {
                "V. Empirical Sequence Alignments, Demographic Skylines & Haplotype Networks"
            }
        }
    }

    pub fn short_badge(&self) -> &'static str {
        match self {
            ModuleGroup::Vertebrates => "VERTEBRATES",
            ModuleGroup::Invertebrates => "INVERTEBRATES",
            ModuleGroup::Microbiology => "MICROBIOLOGY",
            ModuleGroup::Theoretical => "THEO POPBIO",
            ModuleGroup::EmpiricalPopGen => "EMPIRICAL POPGEN",
        }
    }

    pub fn subtitle(&self) -> &'static str {
        match self {
            ModuleGroup::Vertebrates => {
                "Boreal hare-lynx limit cycles, Isle Royale wolf-moose pack dynamics, and St. Kilda Soay sheep crashes"
            }
            ModuleGroup::Invertebrates => {
                "Tribolium LPA cannibalism bifurcations, Huffaker citrus-mite trays, and host-parasitoid refuges"
            }
            ModuleGroup::Microbiology => {
                "Escherichia coli 4-phase batch & Monod chemostat cultures, Gause Paramecium-Didinium, and T4 phage plaques"
            }
            ModuleGroup::Theoretical => {
                "Levins metapopulations, May chaos, Wright-Fisher drift, Kingman coalescents, Sewall Wright 3D fitness landscapes, and Leslie/Lefkovitch matrices"
            }
            ModuleGroup::EmpiricalPopGen => {
                "FASTA/PHYLIP/NEXUS polymorphism, Bandelt Median-Joining networks, Bayesian Skyline / Rogers-Harpending expansion, and sliding-window neutrality scans"
            }
        }
    }

    pub fn accent_rgb(&self) -> [u8; 3] {
        match self {
            ModuleGroup::Vertebrates => [56, 142, 60],
            ModuleGroup::Invertebrates => [216, 122, 28],
            ModuleGroup::Microbiology => [0, 151, 167],
            ModuleGroup::Theoretical => [123, 82, 171],
            ModuleGroup::EmpiricalPopGen => [194, 59, 90],
        }
    }
}

/// Unique identifier for all 19 interactive modules in PopBio Studio v1.2.0.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ModuleId {
    V1HareLynx,
    V2IsleRoyale,
    V3SoaySheep,
    I1TriboliumLpa,
    I2HuffakerMites,
    I3AphidParasitoid,
    M1EcoliBatchChemostat,
    M2GauseProtists,
    M3PhagePlaque,
    T1LevinsMetapop,
    T2MayBifurcation,
    T3WrightFisher,
    T4KingmanCoalescent,
    T5WrightLandscape,
    T6LeslieLefkovitch,
    E1SequencePopGen,
    E2HaplotypeNetwork,
    E3BayesianSkyline,
    E4SlidingWindowScan,
}

impl ModuleId {
    pub const ALL: [ModuleId; 19] = [
        ModuleId::V1HareLynx,
        ModuleId::V2IsleRoyale,
        ModuleId::V3SoaySheep,
        ModuleId::I1TriboliumLpa,
        ModuleId::I2HuffakerMites,
        ModuleId::I3AphidParasitoid,
        ModuleId::M1EcoliBatchChemostat,
        ModuleId::M2GauseProtists,
        ModuleId::M3PhagePlaque,
        ModuleId::T1LevinsMetapop,
        ModuleId::T2MayBifurcation,
        ModuleId::T3WrightFisher,
        ModuleId::T4KingmanCoalescent,
        ModuleId::T5WrightLandscape,
        ModuleId::T6LeslieLefkovitch,
        ModuleId::E1SequencePopGen,
        ModuleId::E2HaplotypeNetwork,
        ModuleId::E3BayesianSkyline,
        ModuleId::E4SlidingWindowScan,
    ];

    pub fn from_slug(slug: &str) -> Option<Self> {
        let s = slug.trim().to_ascii_lowercase();
        match s.as_str() {
            "v1" | "v1-hare-lynx" | "hare-lynx" => Some(Self::V1HareLynx),
            "v2" | "v2-isle-royale" | "isle-royale" => Some(Self::V2IsleRoyale),
            "v3" | "v3-soay-sheep" | "soay-sheep" => Some(Self::V3SoaySheep),
            "i1" | "i1-tribolium" | "tribolium" => Some(Self::I1TriboliumLpa),
            "i2" | "i2-huffaker-mites" | "mites" => Some(Self::I2HuffakerMites),
            "i3" | "i3-aphid-parasitoid" | "aphid" => Some(Self::I3AphidParasitoid),
            "m1" | "m1-ecoli" | "ecoli" => Some(Self::M1EcoliBatchChemostat),
            "m2" | "m2-gause-protists" | "gause" => Some(Self::M2GauseProtists),
            "m3" | "m3-phage-plaque" | "phage" => Some(Self::M3PhagePlaque),
            "t1" | "t1-levins" | "levins" => Some(Self::T1LevinsMetapop),
            "t2" | "t2-may-chaos" | "bifurcation" => Some(Self::T2MayBifurcation),
            "t3" | "t3-wright-fisher" | "wright-fisher" => Some(Self::T3WrightFisher),
            "t4" | "t4-coalescent" | "coalescent" => Some(Self::T4KingmanCoalescent),
            "t5" | "t5-wright-landscape" | "landscape" => Some(Self::T5WrightLandscape),
            "t6" | "t6-leslie-matrix" | "leslie" | "lefkovitch" => Some(Self::T6LeslieLefkovitch),
            "e1" | "e1-sequence-popgen" | "popgen" => Some(Self::E1SequencePopGen),
            "e2" | "e2-haplotype-network" | "network" => Some(Self::E2HaplotypeNetwork),
            "e3" | "e3-bayesian-skyline" | "skyline" => Some(Self::E3BayesianSkyline),
            "e4" | "e4-sliding-window" | "sliding-window" => Some(Self::E4SlidingWindowScan),
            _ => None,
        }
    }

    pub fn slug(&self) -> &'static str {
        match self {
            Self::V1HareLynx => "v1-hare-lynx",
            Self::V2IsleRoyale => "v2-isle-royale",
            Self::V3SoaySheep => "v3-soay-sheep",
            Self::I1TriboliumLpa => "i1-tribolium",
            Self::I2HuffakerMites => "i2-huffaker-mites",
            Self::I3AphidParasitoid => "i3-aphid-parasitoid",
            Self::M1EcoliBatchChemostat => "m1-ecoli",
            Self::M2GauseProtists => "m2-gause-protists",
            Self::M3PhagePlaque => "m3-phage-plaque",
            Self::T1LevinsMetapop => "t1-levins",
            Self::T2MayBifurcation => "t2-may-chaos",
            Self::T3WrightFisher => "t3-wright-fisher",
            Self::T4KingmanCoalescent => "t4-coalescent",
            Self::T5WrightLandscape => "t5-wright-landscape",
            Self::T6LeslieLefkovitch => "t6-leslie-matrix",
            Self::E1SequencePopGen => "e1-sequence-popgen",
            Self::E2HaplotypeNetwork => "e2-haplotype-network",
            Self::E3BayesianSkyline => "e3-bayesian-skyline",
            Self::E4SlidingWindowScan => "e4-sliding-window",
        }
    }

    pub fn group(&self) -> ModuleGroup {
        match self {
            Self::V1HareLynx | Self::V2IsleRoyale | Self::V3SoaySheep => ModuleGroup::Vertebrates,
            Self::I1TriboliumLpa | Self::I2HuffakerMites | Self::I3AphidParasitoid => {
                ModuleGroup::Invertebrates
            }
            Self::M1EcoliBatchChemostat | Self::M2GauseProtists | Self::M3PhagePlaque => {
                ModuleGroup::Microbiology
            }
            Self::T1LevinsMetapop
            | Self::T2MayBifurcation
            | Self::T3WrightFisher
            | Self::T4KingmanCoalescent
            | Self::T5WrightLandscape
            | Self::T6LeslieLefkovitch => ModuleGroup::Theoretical,
            Self::E1SequencePopGen
            | Self::E2HaplotypeNetwork
            | Self::E3BayesianSkyline
            | Self::E4SlidingWindowScan => ModuleGroup::EmpiricalPopGen,
        }
    }

    pub fn metadata(&self) -> ModuleMetadata {
        match self {
            Self::V1HareLynx => ModuleMetadata {
                code: "V1",
                title: "Snowshoe Hare & Canada Lynx Boreal Cycle",
                organisms: "Lepus americanus & Lynx canadensis",
                time_unit: "Years",
                summary: "Rosenzweig-MacArthur predator-prey dynamics with Holling Type I/II/III functional response and boreal foraging movement.",
                equations_plain: "dH/dt = r_H * H * (1 - H / K) - f(H) * L\ndL/dt = e * f(H) * L - m_L * L\nType II f(H) = (a * H) / (1 + a * T_h * H)",
                citation_label: "Stenseth et al. (1997) PNAS 94:5147-5152; Holling (1959) Can. Entomol. 91:385-398",
                citation_url: "https://pubmed.ncbi.nlm.nih.gov/9144205/",
            },
            Self::V2IsleRoyale => ModuleMetadata {
                code: "V2",
                title: "Isle Royale Moose & Gray Wolf Pack Dynamics",
                organisms: "Alces alces & Canis lupus",
                time_unit: "Years",
                summary: "Ratio-dependent pack predation, winter severity forcing, and small-pack inbreeding / Allee rescue on Lake Superior's Isle Royale.",
                equations_plain: "dM/dt = r_M * (1 - 0.35 * W_sev) * M * (1 - M / K_fir) - (a * M * P) / (P + a * T_h * M)\ndP/dt = e_W * g(M, P) * P * A_allee(P) - m_W * P",
                citation_label: "Vucetich & Peterson — Isle Royale Wolf-Moose Project (1958-present)",
                citation_url: "https://www.isleroyalewolf.org/population-dynamics",
            },
            Self::V3SoaySheep => ModuleMetadata {
                code: "V3",
                title: "St. Kilda Soay Sheep Overcompensatory Crashes",
                organisms: "Ovis aries (Hirta, St. Kilda)",
                time_unit: "Years",
                summary: "Discrete-time Maynard Smith-Slatkin overcompensatory herbivore regulation coupled with stochastic North Atlantic winter gales.",
                equations_plain: "N_{t+1} = (N_t * lambda * exp(-w_storm)) / (1 + (N_t / K)^beta)\nbeta > 2 induces overcompensatory 2-point cycles and crash-recovery phases",
                citation_label: "Coulson, Clutton-Brock et al. (2001) Science 292:1528-1531",
                citation_url: "https://pubmed.ncbi.nlm.nih.gov/11375487/",
            },
            Self::I1TriboliumLpa => ModuleMetadata {
                code: "I1",
                title: "Flour Beetle Stage-Structured LPA & Park Competition",
                organisms: "Tribolium castaneum & Tribolium confusum",
                time_unit: "2-Week Censuses",
                summary: "Nonlinear Larvae-Pupae-Adults (LPA) egg/pupa cannibalism bifurcations and two-species flour microcosm competition.",
                equations_plain: "L_{t+1} = b * A_t * exp(-c_el * L_t - c_ea * A_t)\nP_{t+1} = (1 - mu_l) * L_t\nA_{t+1} = P_t * exp(-c_pa * A_t) + (1 - mu_a) * A_t",
                citation_label: "Costantino, Desharnais, Cushing & Dennis (1997) Science 275:389-391",
                citation_url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC314198/",
            },
            Self::I2HuffakerMites => ModuleMetadata {
                code: "I2",
                title: "Huffaker Spatial Citrus-Tray Mite Metapopulation",
                organisms: "Eotetranychus sexmaculatus & Typhlodromus occidentalis",
                time_unit: "Weeks",
                summary: "24-patch (6x4) orange tray lattice showing how partial dispersal barriers rescue locally unstable predator-prey interactions via spatial waves.",
                equations_plain: "dN_i/dt = r * N_i * (1 - N_i / K_patch) - a * N_i * P_i / (1 + a * T_h * N_i) + D_N * (1 - B_vaseline) * Lap(N)_i\ndP_i/dt = e * f(N_i) * P_i - m * P_i + D_P * (1 - B_vaseline) * Lap(P)_i",
                citation_label: "Huffaker (1958) Hilgardia; Muhlbauer et al. (2020) Ecol. Evol. (PMC7713957)",
                citation_url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC7713957/",
            },
            Self::I3AphidParasitoid => ModuleMetadata {
                code: "I3",
                title: "Pea Aphid & Parasitoid Wasp Host-Refuge Dynamics",
                organisms: "Acyrthosiphon pisum & Aphidius ervi",
                time_unit: "Generations",
                summary: "Nicholson-Bailey host-parasitoid interaction stabilized by Negative-Binomial spatial aggregation (clumping k) and physical plant-stem refuges.",
                equations_plain: "H_{t+1} = lambda * H_t * exp(-r_intra * H_t / K) * [refuge + (1 - refuge) * (1 + a * W_t / k)^(-k)]\nW_{t+1} = c * H_t * (1 - escape_fraction)",
                citation_label: "Holling (1959) Can. Entomol. 91:385-398 (Predation & Parasitism)",
                citation_url: "https://www.cambridge.org/core/journals/canadian-entomologist/article/some-characteristics-of-simple-types-of-predation-and-parasitism1/9E1E7D2CCC314766A424680444F4EA9F",
            },
            Self::M1EcoliBatchChemostat => ModuleMetadata {
                code: "M1",
                title: "Escherichia coli Growth Phases & Monod Chemostat",
                organisms: "Escherichia coli K-12",
                time_unit: "Hours",
                summary: "Explicit Lag, Exponential (Log), Diauxic Shift, Stationary, and Death/Senescence phases in batch culture or continuous chemostat dilution.",
                equations_plain: "dA_enz/dt = k_lag * (1 - A_enz)   [Lag physiological adaptation]\nmu(S1, S2) = A_enz * [mu1 * S1 / (Ks1 + S1) + mu2 * S2 / (Ks2 + S2) * (Ki / (Ki + S1))] / (1 + I_tox)\ndX_live/dt = (mu - k_death * I_tox - D) * X_live",
                citation_label: "Monod (1949) Annual Review of Microbiology 3:371-394",
                citation_url: "https://www.annualreviews.org/content/journals/10.1146/annurev.mi.03.100149.002103",
            },
            Self::M2GauseProtists => ModuleMetadata {
                code: "M2",
                title: "Gause Microcosms: Paramecium Competition & Didinium Predation",
                organisms: "Paramecium caudatum, Paramecium aurelia & Didinium nasutum",
                time_unit: "Days",
                summary: "Toggle between Gause's competitive exclusion (P. aurelia vs. P. caudatum) and P. caudatum + Didinium nasutum predation with sediment refuge.",
                equations_plain: "Competition: dN1/dt = r1*N1*(1 - (N1 + alpha12*N2)/K1), dN2/dt = r2*N2*(1 - (N2 + alpha21*N1)/K2)\nPredation: dN_free/dt = r*N*(1 - N/K) - a*(1 - ref)*N*D / (1 + a*Th*N) + I_immig",
                citation_label: "Gause (1934) Struggle for Existence; Muhlbauer et al. (2020) PMC7713957",
                citation_url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC7713957/",
            },
            Self::M3PhagePlaque => ModuleMetadata {
                code: "M3",
                title: "T4 Bacteriophage & E. coli Lytic Plaque Wave",
                organisms: "Escherichia virus T4 & Escherichia coli B",
                time_unit: "Hours",
                summary: "Campbell phage-bacteria chemostat/lawn dynamics tracking susceptible cells, infected eclipse cells, free virions, resistant mutants, and 2D agar plaques.",
                equations_plain: "dS_cell/dt = mu(R)*S_cell - phi*S_cell*V - u_mut*S_cell\ndI_cell/dt = phi*S_cell*V - k_lysis*I_cell\ndV/dt = beta*k_lysis*I_cell - phi*S_cell*V - d_V*V",
                citation_label: "Campbell (1961) Evolution 15:153-165 (Conditions for Existence of Bacteriophage)",
                citation_url: "https://www.scilit.com/publications/3e46b60b6e9b6b089247f59bdcbe9dff",
            },
            Self::T1LevinsMetapop => ModuleMetadata {
                code: "T1",
                title: "Levins Metapopulation & Habitat Fragmentation Threshold",
                organisms: "Patch-Network Metapopulation (N = 64 patches)",
                time_unit: "Seasons",
                summary: "Couples deterministic Levins occupancy ODE with a stochastic 64-patch spatial network subject to habitat destruction fraction h.",
                equations_plain: "dp/dt = c * p * (1 - h - p) - e * p\nEquilibrium p* = (1 - h) - (e / c); Persistence requires h < 1 - (e / c)",
                citation_label: "Levins (1969) Bull. Entomol. Soc. Am. 15:237-240",
                citation_url: "https://academic.oup.com/ae/article-abstract/15/3/237/255899",
            },
            Self::T2MayBifurcation => ModuleMetadata {
                code: "T2",
                title: "May Discrete Logistic & Ricker Bifurcation / Chaos Lab",
                organisms: "Discrete Non-Overlapping Generation Population",
                time_unit: "Generations",
                summary: "Explore monotonic approach, damped oscillations, 2^k period-doubling cascades, chaotic attractors, cobweb maps, and Lyapunov exponents.",
                equations_plain: "May Logistic: N_{t+1} = N_t * [1 + r * (1 - N_t / K)]\nRicker Map:   N_{t+1} = N_t * exp(r * (1 - N_t / K))\nLyapunov lambda = mean(ln |f'(N_t)|)",
                citation_label: "May (1976) Nature 261:459-467",
                citation_url: "https://pubmed.ncbi.nlm.nih.gov/934280/",
            },
            Self::T3WrightFisher => ModuleMetadata {
                code: "T3",
                title: "Wright-Fisher Genetic Drift, Selection & Island F_ST",
                organisms: "Multi-Deme Diploid Mendelian Populations (2N_e)",
                time_unit: "Generations",
                summary: "Simultaneous diploid demes evolving under genetic drift, viability selection (w_AA, w_Aa, w_aa), mutation, and island migration with De Finetti ternary plot.",
                equations_plain: "p' = (p^2 * w_AA + p * q * w_Aa) / w_bar;  w_bar = p^2*w_AA + 2*p*q*w_Aa + q^2*w_aa\np'' = (1 - m)*p_mut + m*p_mean;  k_{t+1} ~ Binomial(2 * N_e, p'')\nH_t = H_0 * (1 - 1 / (2 * N_e))^t",
                citation_label: "Wright (1931) Genetics 16:97-159",
                citation_url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC1201091/",
            },
            Self::T4KingmanCoalescent => ModuleMetadata {
                code: "T4",
                title: "Kingman Coalescent Genealogy & Site Frequency Spectrum",
                organisms: "Neutral Haploid Sample (n = 4..40 chromosomes)",
                time_unit: "Generations (Backward)",
                summary: "Simulates backward-in-time coalescent trees, waiting times E[T_k] = 4*N_e / (k*(k-1)), infinite-sites mutations, and empirical vs. neutral SFS.",
                equations_plain: "Rate lambda_k = k * (k - 1) / (4 * N_e)\nE[T_MRCA] = 4 * N_e * (1 - 1 / n);  E[xi_i] = theta / i",
                citation_label: "Kingman (1982) Stochastic Processes and their Applications 13:235-248",
                citation_url: "https://ideas.repec.org/a/eee/spapps/v13y1982i3p235-248.html",
            },
            Self::T5WrightLandscape => ModuleMetadata {
                code: "T5",
                title: "Sewall Wright 3D Adaptive Fitness Landscape",
                organisms: "Multi-Deme Shifting Balance Metapopulation",
                time_unit: "Generations",
                summary: "Interactive 3D rotatable multi-peak fitness surface W_bar(x, y) illustrating Phase I drift across epistatic valleys, Phase II hill-climbing, and Phase III interdeme dispersion.",
                equations_plain: "dz_i = (s / W_bar) * grad(W_bar) + m * (z_weighted - z_i) + N(0, 1 / sqrt(2 * N_e))\nW_bar(x, y) = W_0 + h_1 * exp(-d_1^2 / sigma_1) + h_2 * exp(-d_2^2 / sigma_2) - v_epistasis",
                citation_label: "Wright (1932) Proc. 6th Int. Congr. Genet. 1:356-366 (Adaptive Landscape)",
                citation_url: "https://www.cambridge.org/core/journals/philosophy-of-science/article/heuristic-role-of-sewall-wrights-1932-adaptive-landscape-diagram/4FB15417FA20A1FE32301B147AF28CB9",
            },
            Self::T6LeslieLefkovitch => ModuleMetadata {
                code: "T6",
                title: "Leslie & Lefkovitch Age/Stage Matrix Projection",
                organisms: "Caretta caretta (7-Stage), Passerine (4-Age) & Perennial Plant (5-Stage)",
                time_unit: "Projection Steps",
                summary: "Projects age- and stage-structured matrices n_{t+1} = A * n_t with dominant eigenvalue lambda_1, stable stage distribution w, Fisher reproductive value v, and elasticity E_ij.",
                equations_plain: "n_{t+1} = A * n_t;  A * w = lambda_1 * w;  A^T * v = lambda_1 * v\nSensitivity S_ij = (v_i * w_j) / <v, w>;  Elasticity E_ij = (a_ij / lambda_1) * S_ij",
                citation_label: "Leslie (1945) Biometrika 33:183-212; Crouse, Crowder & Caswell (1987) Ecology 68:1412-1423",
                citation_url: "https://pubmed.ncbi.nlm.nih.gov/21006835/",
            },
            Self::E1SequencePopGen => ModuleMetadata {
                code: "E1",
                title: "Sequence Alignment Polymorphism, Neutrality & F_ST Suite",
                organisms: "Empirical Multi-FASTA / PHYLIP / NEXUS DNA Alignments",
                time_unit: "Base Pairs (bp)",
                summary: "Computes S, Hd, pi, Watterson's theta_W, Tajima's D, Fu & Li's D*/F*, mismatch raggedness, LD (r^2, |D'|), and multi-population Hudson F_ST / D_XY.",
                equations_plain: "theta_W = S / a_1;  pi = mean_ij(d_ij / L_eff)\nTajima's D = (k_hat - S / a_1) / sqrt(e_1 * S + e_2 * S * (S - 1))\nHudson F_ST = 1 - (pi_within / D_XY)",
                citation_label: "Tajima (1989) Genetics 123:585-595; Watterson (1975); Hudson et al. (1992) Genetics 132:583-589",
                citation_url: "https://pubmed.ncbi.nlm.nih.gov/2513255/",
            },
            Self::E2HaplotypeNetwork => ModuleMetadata {
                code: "E2",
                title: "Bandelt Median-Joining & TCS Haplotype Network Studio",
                organisms: "Intraspecific Mitochondrial / Chloroplast / Nuclear Haplotypes",
                time_unit: "Mutational Steps",
                summary: "Constructs Minimum Spanning, TCS 95% Parsimony, and Bandelt Median-Joining haplotype networks with median vectors, deme pie wedges, and SVG/GEXF export.",
                equations_plain: "Median-Joining (MJ): Union of d <= d_min + epsilon links -> consensus triplet median vectors -> final MST pruning",
                citation_label: "Bandelt, Forster & Rohl (1999) Molecular Biology and Evolution 16:37-48",
                citation_url: "https://pubmed.ncbi.nlm.nih.gov/10331250/",
            },
            Self::E3BayesianSkyline => ModuleMetadata {
                code: "E3",
                title: "Bayesian Skyline & Coalescent Demographic Expansion Lab",
                organisms: "Ultrametric Coalescent Genealogies & Mismatch Waves",
                time_unit: "Substitutions / Generations Ago",
                summary: "MCMC Bayesian Skyline N_e(t) with 95% credible ribbon, Strimmer-Pybus AICc Generalized Skyline, and Rogers-Harpending sudden expansion wave + bootstrap p-values.",
                equations_plain: "log L(g | Theta, A) = sum_i [ ln(k_i*(k_i-1)/(2*theta_{h(i)})) - k_i*(k_i-1)*du_i/(2*theta_{h(i)}) ]\nSudden Expansion Wave: F_i(tau) from (theta_0, theta_1, tau = 2*u*t) with coalescent bootstrap SSD & raggedness",
                citation_label: "Drummond et al. (2005) MBE 22:1185-1192; Strimmer & Pybus (2001); Rogers & Harpending (1992) MBE 9:552-569",
                citation_url: "https://pubmed.ncbi.nlm.nih.gov/15703244/",
            },
            Self::E4SlidingWindowScan => ModuleMetadata {
                code: "E4",
                title: "Sliding-Window Neutrality, Sweep & Divergence Scanner",
                organisms: "Chromosomal / Multi-Gene Aligned Loci",
                time_unit: "Alignment Coordinate (bp)",
                summary: "Scans configurable bp windows and steps across FASTA alignments to localize selective sweeps (Tajima's D < 0), balancing selection peaks (D > 0), and F_ST islands.",
                equations_plain: "Window [x, x + W): pi(w), theta_W(w), Tajima's D(w), Fu & Li's D*(w), Hudson F_ST(w), D_XY(w)",
                citation_label: "Librado & Rozas (2009) Bioinformatics 25:1451-1452 (DnaSP Sliding-Window Polymorphism)",
                citation_url: "https://academic.oup.com/bioinformatics/article/25/11/1451/332507",
            },
        }
    }
}

#[derive(Debug, Clone)]
pub struct ModuleMetadata {
    pub code: &'static str,
    pub title: &'static str,
    pub organisms: &'static str,
    pub time_unit: &'static str,
    pub summary: &'static str,
    pub equations_plain: &'static str,
    pub citation_label: &'static str,
    pub citation_url: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GrowthPhase {
    LagPhase,
    ExponentialLog,
    DiauxicShift,
    DeceleratingTransition,
    StationaryCarryingCapacity,
    PredatorOvershoot,
    PopulationCrashBottleneck,
    DeathSenescence,
    LimitCycleOscillation,
    ChaoticRegime,
    DriftSelectionTrajectory,
    FixationOrLoss,
    EmpiricalAnalysis,
}

impl GrowthPhase {
    pub fn label(&self) -> &'static str {
        match self {
            Self::LagPhase => "Lag Phase (Metabolic / Demographic Acclimation)",
            Self::ExponentialLog => "Exponential / Log Growth Phase",
            Self::DiauxicShift => "Diauxic Metabolic Shift (Secondary Substrate)",
            Self::DeceleratingTransition => "Density-Dependent Deceleration Phase",
            Self::StationaryCarryingCapacity => "Stationary / Carrying Capacity Equilibrium",
            Self::PredatorOvershoot => "Predator / Consumer Peak & Prey Suppression",
            Self::PopulationCrashBottleneck => "Population Crash / Trough Bottleneck",
            Self::DeathSenescence => "Senescence / Exponential Death Phase",
            Self::LimitCycleOscillation => "Sustained Limit-Cycle Oscillation",
            Self::ChaoticRegime => "Aperiodic / Chaotic Fluctuation Regime",
            Self::DriftSelectionTrajectory => "Segregating Polymorphism (Drift / Selection)",
            Self::FixationOrLoss => "Allelic Fixation / Extinction Boundary",
            Self::EmpiricalAnalysis => "Empirical Alignment, Skyline & Neutrality Ready",
        }
    }

    pub fn badge_rgb(&self) -> [u8; 3] {
        match self {
            Self::LagPhase => [120, 144, 156],
            Self::ExponentialLog => [46, 160, 67],
            Self::DiauxicShift => [0, 150, 136],
            Self::DeceleratingTransition => [210, 153, 34],
            Self::StationaryCarryingCapacity => [56, 139, 253],
            Self::PredatorOvershoot => [219, 109, 40],
            Self::PopulationCrashBottleneck => [218, 54, 51],
            Self::DeathSenescence => [140, 45, 65],
            Self::LimitCycleOscillation => [137, 87, 229],
            Self::ChaoticRegime => [209, 59, 138],
            Self::DriftSelectionTrajectory => [63, 136, 197],
            Self::FixationOrLoss => [110, 118, 129],
            Self::EmpiricalAnalysis => [194, 59, 90],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FunctionalResponse {
    TypeILinear,
    TypeIIHyperbolic,
    TypeIIISigmoidal,
}

impl FunctionalResponse {
    pub fn label(&self) -> &'static str {
        match self {
            Self::TypeILinear => "Holling Type I (Linear)",
            Self::TypeIIHyperbolic => "Holling Type II (Saturating)",
            Self::TypeIIISigmoidal => "Holling Type III (Prey-Switching Sigmoidal)",
        }
    }

    pub fn evaluate(&self, prey: f64, a: f64, th: f64) -> f64 {
        let n = prey.max(0.0);
        match self {
            Self::TypeILinear => a * n,
            Self::TypeIIHyperbolic => (a * n) / (1.0 + a * th * n).max(1e-9),
            Self::TypeIIISigmoidal => {
                let n2 = n * n;
                let half_sat_sq = (1.0 / (a * th.max(0.01)).max(1e-6)).max(1.0);
                (a * n2) / (half_sat_sq + a * th * n2).max(1e-9)
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Agent2D {
    pub x: f64,
    pub y: f64,
    pub vx: f64,
    pub vy: f64,
    pub kind: u8,
    pub aux: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeRecord {
    pub t: f64,
    pub values: [f64; 5],
    pub per_capita_growth: f64,
    pub phase: GrowthPhase,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhaseTransitionEvent {
    pub t: f64,
    pub phase: GrowthPhase,
}
