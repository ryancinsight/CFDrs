//! Trifurcation solvers and validation for 3-way branching flows.
//!
//! # Theorem — Generalized Murray Law for Three Daughters
//!
//! Optimal branching under minimum work satisfies
//! $D_0^3 = D_1^3 + D_2^3 + D_3^3$, used as a physiological consistency
//! constraint in trifurcation geometry/validation paths.

pub mod geometry;
pub mod solver;
pub mod validation;

pub use geometry::TrifurcationGeometry3D;
pub use solver::{TrifurcationConfig3D, TrifurcationSolver3D, TrifurcationSolution3D};
pub use validation::{TrifurcationValidator3D, TrifurcationValidationResult3D};

#[cfg(test)]
mod tests;
