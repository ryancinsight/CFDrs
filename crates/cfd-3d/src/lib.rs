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

// Export implemented functionality
pub use fem::{
    FemSolver, FemConfig, Element, ElementType, MaterialProperties, Tetrahedron4
};
pub use mesh_integration::{
    MeshAdapter, MeshQualityReport, StlAdapter, CsgMeshAdapter
};
pub use spectral::{
    SpectralSolver, SpectralConfig, SpectralBasis, SpectralSolution
};

// TODO: Implement these exports when modules are completed
// pub use ibm::{ImmersedBoundaryMethod, IbmSolver};
// pub use level_set::{LevelSetMethod, LevelSet};
// pub use vof::{VolumeOfFluid, VofSolver};

/// Common 3D CFD types and traits
pub mod prelude {
    // Export implemented functionality
    pub use crate::{
        fem::{FemSolver, FemConfig, Element, ElementType, MaterialProperties, Tetrahedron4},
        mesh_integration::{MeshAdapter, MeshQualityReport, StlAdapter, CsgMeshAdapter},
        spectral::{SpectralSolver, SpectralConfig, SpectralBasis, SpectralSolution},
    };

    // TODO: Add exports when implemented
    // pub use crate::{
    //     ibm::IbmSolver,
    //     level_set::LevelSetMethod,
    //     vof::VofSolver,
    // };
}