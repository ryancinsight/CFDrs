//! Mathematical utilities and numerical methods for CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod constants;
pub mod interpolation;
pub mod linear_solver;
pub mod sparse;
pub mod differentiation;
pub mod integration;
pub mod iterators;
pub mod vectorization;

pub use constants::{factors, fem, channel, physical};
pub use interpolation::{Interpolation, LinearInterpolation, CubicSplineInterpolation};
pub use linear_solver::{
    LinearSolver, LinearSolverConfig, ConjugateGradient, GMRES, BiCGSTAB, 
    Preconditioner, IdentityPreconditioner, JacobiPreconditioner, 
    SORPreconditioner, ILU0Preconditioner
};
pub use sparse::{SparseMatrix, SparseMatrixBuilder};
pub use differentiation::{FiniteDifference, Gradient};
pub use integration::{Quadrature, GaussQuadrature};
pub use iterators::{MathIteratorExt, VectorOps, SliceOps, CfdIteratorChain};
pub use vectorization::{VectorizedOps, StencilOps};

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