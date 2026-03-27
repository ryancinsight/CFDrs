//! Channel flow solver implementations.
//!
//! This module decomposes channel flow computation into three sub-modules
//! following the Single Responsibility Principle:
//!
//! | Module | Responsibility |
//! |--------|---------------|
//! | [`geometry_impls`] | `ChannelGeometry` cross-section methods (`area`, `D_h`, `P`) |
//! | [`shape_factors`] | Poiseuille number `Po = f·Re` for each cross-section shape |
//! | [`flow_resistance`] | `Channel` resistance computation (Stokes/laminar/turbulent/slip) |
//!
//! # Governing Physics
//!
//! The hydraulic resistance `R = ΔP / Q` for each flow regime is derived from
//! the Darcy-Weisbach equation and the Poiseuille number:
//!
//! ```text
//! R = Po · μ · L / (2 · A · D_h²)    (laminar / Stokes)
//! R = f · ρ · L · V / (2 · D_h · A)  (turbulent, f from Haaland 1983)
//! R_slip = R_lam / (1 + 4·Kn)        (Beskok-Karniadakis 1999)
//! ```

mod flow_resistance;
mod geometry_impls;
pub mod shape_factors;

#[cfg(test)]
mod tests;
