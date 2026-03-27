//! Energy equation solver for temperature transport.
//!
//! Solves the energy equation for incompressible flows:
//! ∂T/∂t + (u·∇)T = α∇²T + Q/(ρCp)
//! where α = k/(ρCp) is thermal diffusivity
//!
//! # Theorem (Energy Conservation — First Law of Thermodynamics)
//!
//! The discrete energy equation preserves total thermal energy: for an adiabatic
//! system with no internal heat generation (Q = 0), the total enthalpy
//! H = ∫_Ω ρ C_p T dV is conserved to machine precision.
//!
//! **Proof sketch**:
//! The FVM discretisation of ∂T/∂t + ∇·(uT) = α∇²T
//! yields telescoping face fluxes (as in the momentum FVM). With no-flux BCs,
//! all boundary contributions vanish, so Σ_i V_i (T_i^{n+1} - T_i^n) / Δt = 0.
//! The explicit time integration is stable under the thermal CFL condition
//! αΔt/Δx² ≤ 1/4 (2D), ensuring positivity of the temperature
//! update stencil (all coefficients non-negative).

/// Primary solver algorithms.
pub mod solver;
/// Viscous dissipation definitions.
pub mod viscous_dissipation;

#[cfg(test)]
mod tests;

pub use solver::*;
pub use viscous_dissipation::*;

/// Constants for energy equation
pub mod constants {
    /// Default Prandtl number for air
    pub const DEFAULT_PRANDTL: f64 = 0.71;
    /// Stefan-Boltzmann constant [W/(m²·K⁴)]
    pub const STEFAN_BOLTZMANN: f64 = 5.67e-8;
    /// Denominator coefficient for central difference: (u_{i+1} - u_{i-1}) / (CENTRAL_DIFF_COEFF · Δx)
    pub const CENTRAL_DIFF_COEFF: f64 = 2.0;
}
