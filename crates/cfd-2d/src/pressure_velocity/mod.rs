//! Pressure-velocity coupling using Semi-Implicit Method for Pressure-Linked Equations
//!
//! Modularized implementation following SLAP and SOC principles.
//! Reference: Patankar, S.V. (1980). "Numerical Heat Transfer and Fluid Flow"
//!
//! # Theorem (Pressure Correction M-Matrix)
//!
//! On an orthogonal collocated grid with one pressure anchor and positive
//! momentum diagonals, the pressure correction equation
//! $\nabla \cdot (\mathbf{d}\,\nabla p') = \nabla \cdot \mathbf{u}^*$
//! assembles into a symmetric M-matrix, and the anchored reduced system is
//! symmetric positive definite.
//!
//! **Proof sketch**:
//! Each face coefficient $a_f = d_f A_f / \delta_f$ is positive when
//! $d_f = V_f / A_P > 0$. The discrete Laplacian therefore has non-positive
//! off-diagonals and positive row sums, while anchoring one pressure degree of
//! freedom removes the constant null-space. The reduced system is then suitable
//! for CG; GMRES/BiCGSTAB remain safe fallback solvers for non-ideal variants.

pub mod coefficients;
pub(crate) mod boundary;
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
