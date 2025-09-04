//! Mathematical utilities and numerical methods for CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod differentiation;
pub mod error;
pub mod integration;
pub mod interpolation;
pub mod iterators;
pub mod linear_solver;
pub mod simd;
pub mod sparse;
pub mod vector_ops;
pub mod vectorization;

// --- Curated Top-Level API ---
// Only expose a very small number of absolutely fundamental traits or structs.
// Users should interact with the module hierarchy for most types.

pub use self::interpolation::Interpolation;
pub use self::sparse::SparseMatrix;

// The primary API is through the public modules themselves.
// This creates a hierarchical, self-documenting structure.
// Example usage:
//   use cfd_math::linear_solver::BiCGSTAB;
//   use cfd_math::interpolation::CubicSplineInterpolation;
//   use cfd_math::integration::GaussQuadrature;

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        differentiation::FiniteDifference,
        integration::Quadrature,
        interpolation::{Interpolation, LinearInterpolation},
        linear_solver::ConjugateGradient,
        sparse::{SparseMatrix, SparseMatrixBuilder},
    };
}
