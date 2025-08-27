//! Numerical methods domain - Mathematical algorithms and discretization schemes.
//!
//! This module encapsulates numerical method knowledge following DDD principles.
//! It provides abstractions for discretization, time integration, and linear system solving.

pub mod discretization;
pub mod linear_solvers;
pub mod service;
pub mod time_integration;
pub mod traits;

// Re-export main types
pub use discretization::finite_difference;
pub use linear_solvers::{ConjugateGradient, DirectSolver, GaussSeidel, Jacobi};
pub use service::NumericalMethodsService;
pub use time_integration::time_schemes;
pub use traits::{DiscretizationScheme, LinearSystemSolver};
