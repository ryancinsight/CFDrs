#![allow(dead_code)]
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
pub mod iterators;
pub mod vector_ops;
pub mod vectorization;

pub use interpolation::{Interpolation, LinearInterpolation, CubicSplineInterpolation};
pub use linear_solver::{
    LinearSolver, LinearSolverConfig, ConjugateGradient, BiCGSTAB, 
    Preconditioner, IdentityPreconditioner, JacobiPreconditioner, 
    SORPreconditioner
};
pub use sparse::{SparseMatrix, SparseMatrixBuilder};
pub use differentiation::{FiniteDifference, Gradient};
pub use integration::{Quadrature, GaussQuadrature};
pub use iterators::{MathIteratorExt, NormIteratorExt, StatisticsIteratorExt, 
                    WindowIterator, StridedWindowIterator, StencilIterator, StencilPattern,
                    ParallelIteratorExt};
pub use vectorization::{VectorizedOps, StencilOps};
pub use vector_ops::{SimdVectorOps, sparse_matvec};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        interpolation::{Interpolation, LinearInterpolation},
        linear_solver::{LinearSolver, LinearSolverConfig, ConjugateGradient},
        sparse::{SparseMatrix, SparseMatrixBuilder},
        differentiation::FiniteDifference,
        integration::Quadrature,
    };
}