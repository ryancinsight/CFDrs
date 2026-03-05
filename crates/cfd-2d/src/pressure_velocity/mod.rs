//! Pressure-velocity coupling using Semi-Implicit Method for Pressure-Linked Equations
//!
//! Modularized implementation following SLAP and SOC principles.
//! Reference: Patankar, S.V. (1980). "Numerical Heat Transfer and Fluid Flow"
//!
//! # Theorem (Pressure Correction M-Matrix)
//!
//! The pressure correction equation $\nabla \cdot (\mathbf{d}\,\nabla p') = \nabla \cdot \mathbf{u}^*$
//! produces a symmetric negative-definite coefficient matrix (M-matrix) when
//! $d_f = V_f / A_P > 0$, guaranteeing unique solvability by CG or GMRES.
//!
//! **Proof sketch**:
//! Each face coefficient $a_f = d_f A_f / \delta_f > 0$, and $a_P = \sum_f a_f$,
//! yielding strict diagonal dominance with non-positive off-diagonals.
//! guarantees stability and physical realism.

pub mod coefficients;
pub mod config;
mod correction;
pub mod pressure;
pub mod rhie_chow;
pub mod solver;

pub use coefficients::CellCoefficients;
pub use config::{PressureLinearSolver, PressureVelocityConfig};
pub use pressure::PressureCorrectionSolver;
pub use rhie_chow::RhieChowInterpolation;
pub use solver::PressureVelocitySolver;
