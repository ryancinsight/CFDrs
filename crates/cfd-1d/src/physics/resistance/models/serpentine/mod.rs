//! Serpentine channel resistance model with Dean flow corrections.
//!
//! ## Dean Flow in Curved Channels
//!
//! **Dean Number**: When fluid flows through a curved channel, centrifugal forces
//! create secondary flow vortices (Dean vortices). The Dean number quantifies
//! the ratio of centrifugal to viscous forces:
//!
//! De = Re √(`D_h` / `2R_c`)
//!
//! ### Friction Factor Enhancement
//!
//! The friction factor in curved channels is higher than in straight channels
//! due to the secondary flow. The enhancement factor depends on the Dean number:
//!
//! **Laminar regime** (De < 17): Negligible curvature effect
//! **Low Dean number** (17 ≤ De < 370): Ito (1959) or Dean (1928) perturbation
//! **High Dean number** (De ≥ 370): Ito (1959) boundary layer limit
//!
//! ### Bend Loss Coefficient
//!
//! Each 180° bend introduces an additional minor loss:
//! `K_bend` = C₁ + C₂/Re (Idelchik 2007, §6.2)
//!
//! ## Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`model`] | `SerpentineModel<T>` struct, `ResistanceModel<T>` impl |
//! | [`analysis`] | `SerpentineAnalysis<T>`, `bayat_rezai_enhancement()` |
//!
//! ## References
//!
//! - Dean, W. R. (1928). *Proc. R. Soc. Lond. A*, 121(787), 402-420.
//! - Ito, H. (1959). *ASME J. Basic Eng.*, 81(2), 123-134.
//! - Shah & London (1978). *Laminar Flow Forced Convection in Ducts*.
//! - Bayat, P. & Rezai, P. (2017). *Sci. Rep.* 7:13655.
//! - Idelchik, I. E. (2007). *Handbook of Hydraulic Resistance* (4th ed.), §6.1-6.4.

pub mod analysis;
pub mod model;
#[cfg(test)]
mod tests;
mod types;

use super::traits;

// ── Re-exports ──────────────────────────────────────────────────────────────

pub use analysis::{bayat_rezai_enhancement, SerpentineAnalysis};
pub use model::SerpentineModel;
pub use types::{BendType, SerpentineCrossSection};
