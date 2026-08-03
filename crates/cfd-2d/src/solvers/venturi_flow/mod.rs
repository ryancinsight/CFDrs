//! 2D Venturi throat flow solver with pressure recovery validation
//!
//! This module implements a validated solver for Venturi throat flows, which are
//! critical for understanding pressure drops and flow rate measurement in microfluidics.
//!
//! # Physics Background
//!
//! ## Venturi Effect
//!
//! A Venturi is a tapered channel that creates:
//! 1. **Acceleration**: Flow speeds up as area decreases
//! 2. **Pressure drop**: Pressure decreases in the throat (Bernoulli principle)
//! 3. **Recovery**: Pressure recovers in the diffuser as flow decelerates
//!
//! ## Bernoulli Equation (Frictionless)
//!
//! ```text
//! P₁ + (1/2)ρu₁² + ρgh₁ = P₂ + (1/2)ρu₂² + ρgh₂
//! ```
//!
//! At the throat (z₂ = z₁), mass conservation (A₁u₁ = A₂u₂) gives:
//! ```text
//! P₂ = P₁ + (1/2)ρ(u₁² - u₂²)
//!    = P₁ + (1/2)ρu₁²(1 - (A₁/A₂)²)
//! ```
//!
//! ## Pressure Recovery Coefficient
//!
//! The pressure recovery coefficient Cp quantifies the pressure drop:
//! ```text
//! Cp = (P₂ - P₁) / ((1/2)ρu₁²)
//! ```
//!
//! For ideal (frictionless) Venturi:
//! - Minimum Cp in throat (most negative)
//! - Cp → 0 at outlet (full recovery)
//! - Cp_ideal = 1 - (A₁/A₂)² (from Bernoulli)
//!
//! For real (viscous) Venturi:
//! - Cp_actual = (1 - ε) · Cp_ideal, where ε ≈ 0.1-0.3 (recovery loss)
//!
//! # Validation Strategy
//!
//! 1. **Analytical**: Compare against Bernoulli predictions
//! 2. **Convergence**: Grid refinement study with Richardson extrapolation
//! 3. **Literature**: Benchmark against ISO 5167 standards for Venturi meters
//! 4. **Energy**: Verify energy conservation with viscous dissipation
//!
//! # References
//!
//! - Shapiro, A.H. (1953). "The Dynamics and Thermodynamics of Compressible Fluid Flow"
//! - ISO 5167-1:2003. "Measurement of fluid flow by means of pressure differential devices"
//! - White, F.M. (2011). "Fluid Mechanics" (7th ed.)
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

mod analytical;
mod geometry;
mod solution;
mod solver;

pub use analytical::{BernoulliVenturi, ViscousVenturi};
pub use geometry::VenturiGeometry;
pub use solution::VenturiFlowSolution;
pub use solver::{VenturiSolver2D, VenturiValidationResult, VenturiValidator};
