#![forbid(unsafe_code)]
//! PopBio Studio (`popbio-studio` v1.2.0)
//! Native safe-Rust educational population biology, Sewall Wright 3D adaptive landscapes,
//! Leslie/Lefkovitch stage matrices, Bayesian Skyline demographic inference, and haplotype network workstation.

pub mod gui;
pub mod haplotype_network;
pub mod model;
pub mod popgen_alignment;
pub mod popgen_skyline;
pub mod rng;
pub mod sim_invertebrates;
pub mod sim_landscapes_matrices;
pub mod sim_microbiology;
pub mod sim_theoretical;
pub mod sim_vertebrates;
pub mod workset;
