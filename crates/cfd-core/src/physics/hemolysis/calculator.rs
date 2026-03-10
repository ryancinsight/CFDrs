//! Hemolysis calculator with blood properties and clinical indices.
//!
//! # Theorem — Normalized Index of Hemolysis (ASTM F1841-97)
//!
//! The Normalized Index of Hemolysis (NIH) converts the measured increase in
//! plasma-free haemoglobin $\Delta Hb$ into a standardised damage metric:
//!
//! ```text
//! NIH = (1 − Hct) / Hct · ΔHb / Hb₀ · 100   [g / 100 L]
//! ```
//!
//! where $Hct$ is the haematocrit (volume fraction of packed red cells) and
//! $Hb_0$ is the initial total haemoglobin concentration.
//!
//! **Proof sketch.** The free haemoglobin resides in the plasma fraction
//! $(1 - Hct)$; dividing by $Hct$ normalises to the erythrocyte volume that
//! could release haemoglobin. The further division by $Hb_0$ gives a
//! dimensionless damage fraction, and the factor 100 re-scales to the
//! conventional g / 100 L reporting unit.
//!
//! # Theorem — Modified Index of Hemolysis (MIH)
//!
//! The MIH accounts for flow rate $Q$ and exposure time $t_e$:
//!
//! ```text
//! MIH = V / (Q · t_e) · ΔHb   [g / 100 L]
//! ```
//!
//! enabling direct comparison across devices with different priming volumes.
//!
//! ## References
//!
//! - ASTM F1841-97 (2017). "Standard Practice for Assessment of Hemolysis
//!   in Continuous Flow Blood Pumps."
//! - Taskin, M. E. et al. (2012). "Evaluation of Eulerian and Lagrangian
//!   models for hemolysis estimation." *ASAIO J.* 58:363–372.

use super::models::HemolysisModel;
use crate::error::{Error, Result};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Hemolysis calculator with blood properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HemolysisCalculator<T: RealField + Copy> {
    /// Hemolysis model
    model: HemolysisModel,
    /// Hematocrit (volume fraction, 0-1)
    hematocrit: T,
    /// Initial hemoglobin concentration (g/dL)
    hemoglobin_initial: T,
    /// Blood volume (m³)
    blood_volume: T,
    /// Flow rate (m³/s)
    flow_rate: T,
}

impl<T: RealField + Copy + FromPrimitive> HemolysisCalculator<T> {
    /// Create new hemolysis calculator
    pub fn new(
        model: HemolysisModel,
        hematocrit: T,
        hemoglobin_initial: T,
        blood_volume: T,
        flow_rate: T,
    ) -> Result<Self> {
        if hematocrit < T::zero() || hematocrit > T::one() {
            return Err(Error::InvalidConfiguration(
                "Hematocrit must be between 0 and 1".to_string(),
            ));
        }
        if hemoglobin_initial <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Initial hemoglobin must be positive".to_string(),
            ));
        }
        if blood_volume <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Blood volume must be positive".to_string(),
            ));
        }
        if flow_rate <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Flow rate must be positive".to_string(),
            ));
        }

        Ok(Self {
            model,
            hematocrit,
            hemoglobin_initial,
            blood_volume,
            flow_rate,
        })
    }

    /// Calculate Normalized Index of Hemolysis (NIH)
    pub fn normalized_index(&self, delta_hemoglobin: T) -> T {
        let hundred = T::from_f64(100.0).unwrap_or_else(|| T::one());
        let one_minus_hct = (hundred - self.hematocrit * hundred) / hundred;

        if self.hemoglobin_initial > T::zero() {
            one_minus_hct / self.hematocrit * delta_hemoglobin / self.hemoglobin_initial * hundred
        } else {
            T::zero()
        }
    }

    /// Calculate Modified Index of Hemolysis (MIH) in g/(100L)
    pub fn modified_index(&self, delta_hemoglobin: T, exposure_time: T) -> T {
        let hundred = T::from_f64(100.0).unwrap_or_else(|| T::one());

        if self.flow_rate > T::zero() && exposure_time > T::zero() {
            self.blood_volume / (self.flow_rate * exposure_time)
                * delta_hemoglobin
                * (hundred - self.hematocrit * hundred)
                / hundred
        } else {
            T::zero()
        }
    }

    /// Calculate hemoglobin release from damage index
    pub fn hemoglobin_release(&self, damage_index: T) -> T {
        let scaling = self.hematocrit / (T::one() - self.hematocrit);
        self.hemoglobin_initial * damage_index * scaling
    }

    /// Estimate exposure time from geometry
    pub fn estimate_exposure_time(&self, characteristic_length: T) -> T {
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(|| T::one());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::one());

        let area = pi * characteristic_length * characteristic_length / four;

        if area > T::zero() {
            let velocity = self.flow_rate / area;
            if velocity > T::zero() {
                characteristic_length / velocity
            } else {
                T::zero()
            }
        } else {
            T::zero()
        }
    }

    /// Calculate cumulative damage over path with varying shear stress
    pub fn cumulative_damage(&self, stress_time_history: &[(f64, f64)]) -> Result<f64> {
        let mut total_damage = 0.0;

        for &(stress, time) in stress_time_history {
            let damage = self.model.damage_index(stress, time)?;
            total_damage += damage;
        }

        Ok(total_damage)
    }

    /// Get a reference to the hemolysis model
    pub fn model(&self) -> &HemolysisModel {
        &self.model
    }

    /// Calculate critical shear stress for onset of hemolysis
    pub fn critical_shear_stress(&self, exposure_time: f64) -> Result<f64> {
        let threshold_damage = 0.01;

        match self.model {
            HemolysisModel::PowerLaw {
                coefficient,
                stress_exponent,
                time_exponent,
            } => {
                let time_term = exposure_time.powf(time_exponent);
                if time_term > 0.0 && coefficient > 0.0 {
                    Ok((threshold_damage / (coefficient * time_term)).powf(1.0 / stress_exponent))
                } else {
                    Ok(f64::INFINITY)
                }
            }
            HemolysisModel::HeuserOpitz { threshold, .. } => Ok(threshold),
            HemolysisModel::Zhang { .. } => Ok(1000.0 * 0.0035),
        }
    }
}
