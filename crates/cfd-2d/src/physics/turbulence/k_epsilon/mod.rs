//! # k-ε Turbulence Model
//!
//! Implements the standard k-ε (Launder & Spalding 1974) and Realizable k-ε
//! (Shih, Zhu & Lumley 1995) turbulence models, with optional Kato-Launder
//! (1993) vorticity-strain production.
//!
//! ## Governing Equations
//!
//! **Turbulent Kinetic Energy (k):**
//! ```text
//! ∂k/∂t + U_j ∂k/∂x_j = ∂/∂x_j[(ν + ν_t/σ_k) ∂k/∂x_j] + P_k − ε
//! ```
//!
//! **Dissipation Rate (ε):**
//! ```text
//! ∂ε/∂t + U_j ∂ε/∂x_j = ∂/∂x_j[(ν + ν_t/σ_ε) ∂ε/∂x_j] + C_ε1 ε/k P_k − C_ε2 ε²/k
//! ```
//!
//! **Eddy Viscosity (Boussinesq, 1877):**
//! ```text
//! ν_t = C_μ k² / ε
//! ```
//!
//! ## Standard Constants (Launder & Spalding 1974)
//!
//! | Constant | Value | Source |
//! |----------|-------|--------|
//! | C_μ      | 0.09  | Launder & Spalding (1974), Table 1 |
//! | C_ε1     | 1.44  | '' |
//! | C_ε2     | 1.92  | '' |
//! | σ_k      | 1.0   | '' |
//! | σ_ε      | 1.3   | '' |
//!
//! ## Realizability Theorem
//!
//! The Reynolds stress tensor τ_ij = −ρ u'_i u'_j must be positive
//! semi-definite for physical consistency. This requires k ≥ 0 and
//! ε ≥ 0. The implementation enforces these bounds after every update step
//! via clipping: k = max(k, k_min), ε = max(ε, ε_min).
//!
//! ## Literature References
//!
//! - Launder & Spalding (1974). *Comp. Meth. Appl. Mech. Eng.*, 3(2), 269–289.
//! - Jones & Launder (1972). *Int. J. Heat Mass Transfer*, 15(2), 301–314.
//! - Shih, Zhu & Lumley (1995). *Computers & Fluids*, 24(3), 227–238.
//! - Kato & Launder (1993). *Proc. 9th Symp. Turb. Shear Flows*, Kyoto.

/// Core k-ε model struct and `TurbulenceModel` trait implementation.
pub mod model;

/// Realizable k-ε variant (Shih, Zhu & Lumley 1995): strain-dependent C_μ.
pub mod realizable;

/// Kato-Launder (1993) vorticity-strain production modification.
pub mod kato_launder;

pub use model::KEpsilonModel;

#[cfg(test)]
mod tests;
