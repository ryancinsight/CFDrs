//! `NavierStokesSolver2D` — SIMPLE pressure-velocity coupling solver.
//!
//! Implements the Semi-Implicit Method for Pressure-Linked Equations (SIMPLE;
//! Patankar, 1980) for 2D incompressible flow with non-Newtonian blood rheology
//! on a Cartesian staggered grid.
//!
//! ## Algorithm
//! 1. Solve u-momentum (Gauss-Seidel + under-relaxation)
//! 2. Solve v-momentum (Gauss-Seidel + under-relaxation)
//! 3. Solve pressure-correction equation from continuity residual
//! 4. Correct velocities + pressure
//! 5. Apply Rhie-Chow face-velocity interpolation (via `cfd-core`)
//! 6. Update viscosity from shear rate (non-Newtonian)
//! 7. Test convergence
//!
//! ## References
//! - Patankar (1980): §6.3–6.7, §7.4
//! - Rhie & Chow (1983): Pressure-velocity interpolation
//! - Versteeg & Malalasekera (2007): §11.5
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

mod momentum;
mod monitoring;
mod pressure;
mod solve;
mod types;
mod velocity_interpolation;
pub use types::NavierStokesSolver2D;
