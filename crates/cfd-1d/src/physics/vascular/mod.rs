//! Vascular flow models for hemodynamic simulations
//!
//! This module provides 1D blood flow models for arterial and venous networks,
//! including pulsatile flow solutions and bifurcation geometry optimization.
//!
//! # Mathematical Foundations
//!
//! ## Womersley Flow
//! Analytical solution for pulsatile flow in rigid cylindrical vessels:
//! ```text
//! u(r,t) = Re{ (P̂/iρω) · [1 - J₀(α·r/R·i^(3/2)) / J₀(α·i^(3/2))] · e^(iωt) }
//! ```
//! Where α = R√(ωρ/μ) is the Womersley number.
//!
//! ## Murray's Law
//! Optimal bifurcation geometry minimizing metabolic cost:
//! ```text
//! D₀³ = D₁³ + D₂³  (parent = sum of cubed daughter diameters)
//! ```
//!
//! ## 1D Wave Propagation
//! Characteristic equations for arterial pressure waves.
//!
//! # References
//! - Womersley, J.R. (1955) "Method for the calculation of velocity, rate of flow
//!   and viscous drag in arteries when the pressure gradient is known"
//! - Murray, C.D. (1926) "The Physiological Principle of Minimum Work"
//! - Olufsen, M.S. (1999) "Structured tree outflow condition for blood flow in
//!   larger systemic arteries"

pub mod bessel;
pub mod bifurcation;
pub mod murrays_law;
pub mod womersley;

pub use bifurcation::{Bifurcation, BifurcationNetwork, JunctionType};
pub use murrays_law::{MurraysLaw, OptimalBifurcation};
pub use womersley::{WomersleyFlow, WomersleyNumber, WomersleyProfile};
