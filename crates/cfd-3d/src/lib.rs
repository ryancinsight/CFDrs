//! 3D CFD simulations with CSGrs integration.

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

// TODO: Implement these exports
// pub use fem::{FiniteElementMethod, FemSolver, Element};
// pub use spectral::{SpectralMethod, SpectralSolver};
// pub use ibm::{ImmersedBoundaryMethod, IbmSolver};
// pub use level_set::{LevelSetMethod, LevelSet};
// pub use vof::{VolumeOfFluid, VofSolver};
// pub use mesh_integration::{MeshAdapter, CsgMeshAdapter};

/// Common 3D CFD types and traits
pub mod prelude {
    // TODO: Add exports when implemented
    // pub use crate::{
    //     fem::FemSolver,
    //     spectral::SpectralSolver,
    //     ibm::IbmSolver,
    //     level_set::LevelSetMethod,
    //     vof::VofSolver,
    //     mesh_integration::MeshAdapter,
    // };
}