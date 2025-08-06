//! 3D CFD solvers with CSGrs mesh integration.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod fem;
pub mod spectral;
pub mod ibm;
pub mod level_set;
pub mod vof;
pub mod mesh_integration;

pub use fem::{FiniteElementMethod, FemSolver, Element};
pub use spectral::{SpectralMethod, SpectralSolver};
pub use ibm::{ImmersedBoundaryMethod, IbmSolver};
pub use level_set::{LevelSetMethod, LevelSet};
pub use vof::{VolumeOfFluid, VofSolver};
pub use mesh_integration::{MeshAdapter, CsgMeshAdapter};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        fem::FemSolver,
        spectral::SpectralSolver,
        ibm::IbmSolver,
        level_set::LevelSetMethod,
        vof::VofSolver,
        mesh_integration::MeshAdapter,
    };
}