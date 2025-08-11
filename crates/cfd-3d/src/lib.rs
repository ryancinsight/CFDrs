//! 3D CFD solvers and algorithms
//!
//! This crate provides various 3D computational fluid dynamics solvers including:
//! - Finite Element Method (FEM)
//! - Spectral methods
//! - Mesh integration with CSGrs
//! - Immersed Boundary Method (IBM)
//! - Level Set Method
//! - Volume of Fluid (VOF)

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod fem;
pub mod spectral;
pub mod mesh_integration;
pub mod ibm;
pub mod level_set;
pub mod vof;

// Re-export commonly used types
pub use fem::{FemSolver, FemConfig};
pub use spectral::{SpectralSolver, SpectralConfig};
pub use mesh_integration::{MeshAdapter, MeshQualityReport, StlAdapter, CsgMeshAdapter};
pub use ibm::{IbmSolver, IbmConfig, LagrangianPoint};
pub use level_set::{LevelSetSolver, LevelSetConfig};
pub use vof::{VofSolver, VofConfig};


// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface