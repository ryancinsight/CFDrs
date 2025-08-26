#![allow(dead_code)]
//! Mathematical utilities and numerical methods for CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
pub mod differentiation;
pub mod integration;
pub mod interpolation;
pub mod iterators;
pub mod linear_solver;
pub mod sparse;
pub mod vector_ops;
pub mod vectorization;
pub use interpolation::{CubicSplineInterpolation, Interpolation, LinearInterpolation};
pub use linear_solver::{
    BiCGSTAB, ConjugateGradient, IdentityPreconditioner, JacobiPreconditioner, LinearSolver,
    Preconditioner, SORPreconditioner,
};
// LinearSolverConfig is re-exported from cfd_core in linear_solver module
pub use differentiation::{FiniteDifference, Gradient};
pub use integration::{GaussQuadrature, Quadrature};
pub use iterators::{
    MathIteratorExt, NormIteratorExt, ParallelIteratorExt, StatisticsIteratorExt, StencilIterator,
    StencilPattern, StridedWindowIterator, WindowIterator,
pub use sparse::{SparseMatrix, SparseMatrixBuilder, SparseMatrixExt, SparsePatterns};
pub use vector_ops::{sparse_matvec, SimdVectorOps};
pub use vectorization::{StencilOps, VectorizedOps};
/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        differentiation::FiniteDifference,
        integration::Quadrature,
        interpolation::{Interpolation, LinearInterpolation},
        linear_solver::{ConjugateGradient, LinearSolver},
        sparse::{SparseMatrix, SparseMatrixBuilder},
    };


}
}
