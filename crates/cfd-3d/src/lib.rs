//! 3D CFD simulations with CSGrs integration.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod fem;
pub mod spectral;
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



// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface