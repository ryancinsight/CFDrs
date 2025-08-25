//! 3D CFD simulations.

#![warn(missing_docs)]

pub mod fem;
pub mod ibm;
pub mod level_set;
pub mod spectral;
pub mod vof;

// Export FEM functionality
pub use fem::{FemConfig, FemSolver, StokesFlowProblem};

// Export spectral functionality
pub use spectral::{
    BasisFunction, ChebyshevPolynomial, FourierTransform, SpectralBasis, SpectralConfig,
    SpectralSolution, SpectralSolver,
};

// Export IBM functionality
pub use ibm::{IbmConfig, IbmSolver, LagrangianPoint};

// Export level set functionality
pub use level_set::{LevelSetConfig, LevelSetSolver};

// Export VOF functionality
pub use vof::{VofConfig, VofSolver};

// Use CSGrs integration from cfd-mesh
#[cfg(feature = "csg")]
pub use cfd_mesh::csg::CsgMeshAdapter;
