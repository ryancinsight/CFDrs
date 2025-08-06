//! Mathematical utilities and numerical methods for CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod interpolation;
pub mod linear_solver;
pub mod sparse;
pub mod differentiation;
pub mod integration;

pub use interpolation::{Interpolation, LinearInterpolation, CubicSplineInterpolation};
pub use linear_solver::{LinearSolver, ConjugateGradient, GMRES, BiCGSTAB};
pub use sparse::{SparseMatrix, SparseMatrixBuilder};
pub use differentiation::{FiniteDifference, Gradient};
pub use integration::{Quadrature, GaussQuadrature};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        interpolation::{Interpolation, LinearInterpolation},
        linear_solver::{LinearSolver, ConjugateGradient},
        sparse::{SparseMatrix, SparseMatrixBuilder},
        differentiation::FiniteDifference,
        integration::Quadrature,
    };
}