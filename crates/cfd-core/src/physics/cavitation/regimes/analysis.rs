//! Cavitation regime analysis results and reporting.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

use crate::error::Result;

use super::classifier::CavitationRegimeClassifier;
use super::types::CavitationRegime;

/// Cavitation regime analysis results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CavitationRegimeAnalysis<T: RealField + Copy> {
    /// Identified regime
    pub regime: CavitationRegime,
    /// Blake threshold pressure (Pa)
    pub blake_threshold: T,
    /// Inertial threshold pressure amplitude (Pa)
    pub inertial_threshold: T,
    /// Cavitation number
    pub cavitation_number: T,
    /// Mechanical index (for acoustic cavitation)
    pub mechanical_index: T,
    /// Estimated maximum bubble radius (m)
    pub max_bubble_radius: T,
    /// Sonoluminescence probability (0-1)
    pub sonoluminescence_probability: T,
    /// Material damage potential (0-1)
    pub damage_potential: T,
    /// Hemolysis risk (0-1)
    pub hemolysis_risk: T,
}

impl<T: RealField + Copy + FromPrimitive> CavitationRegimeClassifier<T> {
    /// Perform comprehensive regime analysis
    pub fn analyze(&self) -> Result<CavitationRegimeAnalysis<T>> {
        let regime = self.classify_regime();
        let blake_threshold = self.blake_threshold();
        let inertial_threshold = self.inertial_threshold();

        let current_pressure = if let Some(p_ac) = self.acoustic_pressure {
            self.ambient_pressure - p_ac
        } else {
            self.ambient_pressure
        };

        let two = T::from_f64(2.0).unwrap_or_else(|| T::one());
        let v_sq =
            (self.ambient_pressure - current_pressure) * two / self.bubble_model.liquid_density;
        let current_velocity = if v_sq > T::zero() {
            v_sq.sqrt()
        } else {
            T::from_f64(1e-6).unwrap_or_else(|| T::one())
        };

        let cavitation_number = self.cavitation_number(current_pressure, current_velocity);
        let mechanical_index = self.mechanical_index()?;
        let max_bubble_radius = self.estimate_max_radius(regime);
        let sonoluminescence_probability = self.sonoluminescence_probability();
        let damage_potential = self.damage_potential();
        let hemolysis_risk = self.hemolysis_risk();

        Ok(CavitationRegimeAnalysis {
            regime,
            blake_threshold,
            inertial_threshold,
            cavitation_number,
            mechanical_index,
            max_bubble_radius,
            sonoluminescence_probability,
            damage_potential,
            hemolysis_risk,
        })
    }
}

impl<T: RealField + Copy + std::fmt::Display> std::fmt::Display for CavitationRegimeAnalysis<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Cavitation Regime Analysis:\n\
             Regime: {:?}\n\
             Blake Threshold: {} Pa\n\
             Inertial Threshold: {} Pa\n\
             Cavitation Number: {}\n\
             Mechanical Index: {}\n\
             Max Bubble Radius: {} μm\n\
             Sonoluminescence Probability: {:.1}%\n\
             Damage Potential: {:.1}%\n\
             Hemolysis Risk: {:.1}%",
            self.regime,
            self.blake_threshold,
            self.inertial_threshold,
            self.cavitation_number,
            self.mechanical_index,
            self.max_bubble_radius,
            self.sonoluminescence_probability,
            self.damage_potential,
            self.hemolysis_risk,
        )
    }
}
