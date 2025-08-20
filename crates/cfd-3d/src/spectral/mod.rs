//! Spectral methods for 3D CFD simulations
//!
//! This module provides high-order spectral methods following
//! proper separation of concerns.

pub mod basis;
pub mod chebyshev;
pub mod fourier;
pub mod poisson;
pub mod solver;

pub use basis::{SpectralBasis, BasisFunction};
pub use chebyshev::{ChebyshevPolynomial, ChebyshevDifferentiation};
pub use fourier::{FourierTransform, SpectralDerivative};
pub use poisson::{PoissonSolver, PoissonBoundaryCondition};
pub use solver::{SpectralSolver, SpectralConfig, SpectralSolution};