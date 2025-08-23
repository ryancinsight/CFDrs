//! 3D CFD simulations.

#![warn(missing_docs)]

pub mod fem;
pub mod ibm;
pub mod level_set;
pub mod spectral;
pub mod vof;

// Export FEM functionality
pub use fem::{FemSolver, FemConfig, StokesFlowProblem};

// Export spectral functionality
pub use spectral::{
    SpectralSolver, SpectralConfig, SpectralSolution, SpectralBasis,
    BasisFunction, FourierTransform, ChebyshevPolynomial,
};

// Export IBM functionality
pub use ibm::{IbmSolver, IbmConfig, LagrangianPoint};

// Export level set functionality
pub use level_set::{LevelSetSolver, LevelSetConfig};

// Export VOF functionality
pub use vof::{VofSolver, VofConfig};

// Use CSGrs integration from cfd-mesh
#[cfg(feature = "csg")]
pub use cfd_mesh::csg::CsgMeshAdapter;