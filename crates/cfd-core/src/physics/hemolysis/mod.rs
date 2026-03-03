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
    HemolysisModel,
    CAVITATION_HI_SLOPE,
    GIERSIEPEN_MILLIFLUIDIC_C,
    GIERSIEPEN_MILLIFLUIDIC_STRESS,
    GIERSIEPEN_MILLIFLUIDIC_TIME,
};
pub use trauma::{BloodTrauma, BloodTraumaSeverity, PlateletActivation};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_giersiepen_damage_increases_with_stress() {
        let model = HemolysisModel::giersiepen_standard();
        let t = 1.0;

        let d1 = model.damage_index(50.0, t).unwrap();
        let d2 = model.damage_index(100.0, t).unwrap();
        let d3 = model.damage_index(200.0, t).unwrap();

        assert!(d2 > d1);
        assert!(d3 > d2);
    }

    #[test]
    fn test_damage_increases_with_time() {
        let model = HemolysisModel::giersiepen_standard();
        let tau = 100.0;

        let d1 = model.damage_index(tau, 0.1).unwrap();
        let d2 = model.damage_index(tau, 1.0).unwrap();
        let d3 = model.damage_index(tau, 10.0).unwrap();

        assert!(d2 > d1);
        assert!(d3 > d2);
    }

    #[test]
    fn test_heuser_opitz_threshold() {
        let model = HemolysisModel::heuser_opitz();

        let d_below = model.damage_index(100.0, 1.0).unwrap();
        let d_above = model.damage_index(200.0, 1.0).unwrap();

        assert_eq!(d_below, 0.0);
        assert!(d_above > 0.0);
    }

    #[test]
    fn test_nih_calculation() {
        let calc = HemolysisCalculator::new(
            HemolysisModel::default(),
            0.45,
            15.0,
            5e-3,
            1e-4,
        )
        .unwrap();

        let delta_hb = 0.1;
        let nih = calc.normalized_index(delta_hb);

        assert!(nih > 0.0);
        assert!(nih < 100.0);
    }

    #[test]
    fn test_critical_shear_stress() {
        let calc = HemolysisCalculator::new(
            HemolysisModel::giersiepen_standard(),
            0.45,
            15.0,
            5e-3,
            1e-4,
        )
        .unwrap();

        let tau_crit = calc.critical_shear_stress(1.0).unwrap();

        assert!(tau_crit > 10.0);
        assert!(tau_crit < 1000.0);
    }

    #[test]
    fn test_blood_trauma_severity() {
        let trauma_minimal = BloodTrauma {
            hemolysis_level: 5.0,
            platelet_activation: 2.0,
            thrombosis_risk: 0.01,
            max_shear_stress: 80.0,
            avg_exposure_time: 0.1,
        };

        assert_eq!(trauma_minimal.severity(), BloodTraumaSeverity::Minimal);
        assert!(trauma_minimal.meets_fda_guidance());

        let trauma_severe = BloodTrauma {
            hemolysis_level: 100.0,
            platelet_activation: 30.0,
            thrombosis_risk: 0.5,
            max_shear_stress: 400.0,
            avg_exposure_time: 2.0,
        };

        assert_eq!(trauma_severe.severity(), BloodTraumaSeverity::Severe);
        assert!(!trauma_severe.meets_fda_guidance());
    }
}
