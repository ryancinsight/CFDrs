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



/// 3D CFD domain-specific prelude for advanced mesh operations
///
/// This prelude exports 3D-specific functionality not available in the main prelude.
/// Use this when working extensively with 3D meshes, CSG operations, or advanced
/// numerical methods like spectral or immersed boundary methods.
///
/// For basic 3D functionality, prefer `cfd_suite::prelude::*`.
pub mod prelude {
    // === Advanced Mesh Operations ===
    // Detailed mesh functionality beyond basic adapters
    pub use crate::mesh_integration::{
        MeshQualityReport, StlAdapter, CsgMeshAdapter
    };

    // === Specialized Numerical Methods ===
    // Advanced 3D-specific solvers and techniques
    pub use crate::{
        fem::{Element, ElementType, MaterialProperties, Tetrahedron4},
        spectral::{SpectralSolution, SpectralBasis},
    };
}