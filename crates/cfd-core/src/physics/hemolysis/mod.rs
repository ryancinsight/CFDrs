//! Hemolysis and blood damage models for microfluidic and millifluidic applications.
//!
//! Provides mathematical models for predicting red blood cell (RBC) damage
//! and hemolysis in flow systems, particularly relevant for:
//! - Cardiovascular devices (pumps, valves, oxygenators)
//! - Microfluidic blood processing
//! - Millifluidic diagnostic devices
//! - Venturi cavitation systems
//!
//! ## Mathematical Foundation
//!
//! ### Power Law Model (Giersiepen et al. 1990)
//! ```math
//! D = C · τ^α · t^β
//! ```
//!
//! ### Normalized Index of Hemolysis (NIH)
//! ```math
//! NIH = (100 − Hct) / Hct · ΔHb / Hb₀ · 100%
//! ```
//!
//! ## References
//! - Giersiepen, M. et al. (1990). Estimation of shear stress-related blood damage.
//! - Zhang, T. et al. (2011). Study of flow-induced hemolysis.

/// Hemolysis calculator with blood properties and clinical indices.
mod calculator;
/// Hemolysis model definitions and damage index calculations.
mod models;
/// Blood trauma assessment, platelet activation, and severity classification.
mod trauma;

pub use calculator::HemolysisCalculator;
pub use models::{
    CAVITATION_HI_SLOPE, GIERSIEPEN_MILLIFLUIDIC_C, GIERSIEPEN_MILLIFLUIDIC_STRESS,
    GIERSIEPEN_MILLIFLUIDIC_TIME, HemolysisModel,
};
pub use trauma::{BloodTrauma, BloodTraumaSeverity, PlateletActivation};

#[cfg(test)]
mod tests;
