//! 3D CFD solvers and algorithms
//!
//! This crate provides various 3D computational fluid dynamics solvers including:
//! - Finite Element Method (FEM)
//! - Spectral methods
//! - Mesh integration with CSGrs via cfd-mesh
//! - Immersed Boundary Method (IBM)
//! - Level Set Method
//! - Volume of Fluid (VOF)

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod fem;
pub mod spectral; // Refactored into submodules for proper separation
pub mod ibm;
pub mod level_set;
pub mod vof;

// Re-export commonly used types
pub use fem::{FemSolver, FemConfig, FluidProperties, FluidElement};
pub use spectral::{SpectralSolver, SpectralConfig, SpectralBasis};
pub use ibm::{IbmSolver, IbmConfig, LagrangianPoint};
pub use level_set::{LevelSetSolver, LevelSetConfig};
pub use vof::{VofSolver, VofConfig};

// Use CSGrs integration from cfd-mesh
#[cfg(feature = "csg")]
pub use cfd_mesh::csg::CsgMeshAdapter;

// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface