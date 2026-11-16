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
