//! Generalized Minimal Residual (GMRES) method — compatibility wrapper module.
//!
//! The solver implementation now delegates to [`leto_ops::GMRES`] (SSOT).

mod solver;

pub use solver::GMRES;
