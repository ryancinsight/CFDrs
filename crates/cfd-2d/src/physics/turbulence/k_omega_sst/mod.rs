//! # k-ω SST Turbulence Model
//!
//! ## Mathematical Foundation
//!
//! The k-ω SST (Shear Stress Transport) model was developed by Menter (1994).
//! It combines the robust near-wall performance of the k-ω model with
//! the freestream independence of the k-ε model through blending functions.
//!
//! ## Module Structure
//!
//! | Module | Responsibility |
//! |--------|---------------|
//! | [`blending`] | F1/F2 blending functions, cross-diffusion, coefficient blending |
//! | [`limiter`] | Menter (2003) SST production limiter |
//! | [`model`] | `KOmegaSSTModel` struct and `TurbulenceModel` trait impl |
//!
//! ## Governing Equations
//!
//! **k equation**: `∂k/∂t + Uⱼ ∂k/∂xⱼ = ∂/∂xⱼ[(ν + νₜ/σₖ) ∂k/∂xⱼ] + Pₖ − β* k ω`
//!
//! **ω equation**: `∂ω/∂t + Uⱼ ∂ω/∂xⱼ = ∂/∂xⱼ[(ν + νₜ/σω) ∂ω/∂xⱼ] + γ Pₖ/νₜ − β ω² + CD_kω`
//!
//! **Eddy viscosity**: `νₜ = a₁ k / max(a₁ ω, S F₂)` (Bradshaw assumption)
//!
//! ## Realizability
//!
//! The model enforces `k ≥ k_min` and `ω ≥ ω_min` at every time step,
//! ensuring the Reynolds stress tensor remains positive semi-definite.
//!
//! ## References
//! - Menter, F.R. (1994). AIAA Journal, 32(8), 1598-1605.
//! - Menter, F.R., Kuntz, M. & Langtry, R. (2003). Turbulence, Heat and
//!   Mass Transfer 4:625-632.
//! - Wilcox, D.C. (2008). *Turbulence Modeling for CFD*. DCW Industries.

pub mod blending;
pub mod limiter;
mod model;

pub use limiter::limit_production;
pub use model::KOmegaSSTModel;

#[cfg(test)]
mod tests;
