//! Cavitation regime classification: stable vs. inertial cavitation.
//!
//! ## Cavitation Regimes
//!
//! ### Stable Cavitation
//! Bubbles oscillate about equilibrium radius without collapse.
//! Less damaging to materials and biological cells.
//!
//! ### Inertial Cavitation
//! Bubbles grow rapidly and collapse violently. Produces shock waves,
//! microjets, and high local temperatures.
//!
//! ## Mathematical Criteria
//!
//! ### Blake Threshold
//! ```math
//! P_Blake = P_v + 4σ/(3R_c)
//! ```
//!
//! ### Inertial Cavitation Threshold (Apfel & Holland 1991)
//! ```math
//! P_threshold = P_v + √(8σ/(3R_0)) · (P_∞ + 2σ/R_0)^(1/2)
//! ```

/// Cavitation regime analysis results and reporting.
mod analysis;
/// Cavitation regime classifier.
mod classifier;
/// Cavitation regime types.
mod types;

pub use analysis::CavitationRegimeAnalysis;
pub use classifier::CavitationRegimeClassifier;
pub use types::CavitationRegime;

#[cfg(test)]
mod tests;
