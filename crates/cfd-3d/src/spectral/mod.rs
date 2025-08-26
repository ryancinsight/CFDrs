//! Spectral methods for 3D CFD simulations
//!
//! This module provides high-order spectral methods following
//! proper separation of concerns.

pub mod basis;
pub mod chebyshev;
pub mod fourier;
pub mod poisson;
pub mod solver;
pub use basis::{BasisFunction, SpectralBasis};
pub use chebyshev::{ChebyshevDifferentiation, ChebyshevPolynomial};
pub use fourier::{FourierTransform, SpectralDerivative};
pub use poisson::{PoissonBoundaryCondition, PoissonSolver};
pub use solver::{SpectralConfig, SpectralSolution, SpectralSolver};
