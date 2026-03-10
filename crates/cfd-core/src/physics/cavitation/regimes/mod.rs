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
//! P_Blake = P_v + 2σ/R_c · (1 + 2σ/(3R_c(P_∞ − P_v)))
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
mod tests {
    use super::*;
    use crate::physics::cavitation::rayleigh_plesset::RayleighPlesset;

    fn create_test_bubble() -> RayleighPlesset<f64> {
        RayleighPlesset {
            initial_radius: 1e-6,
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        }
    }

    #[test]
    fn test_regime_classification_none() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(bubble, 1e5, None, None);

        let regime = classifier.classify_regime();
        assert_eq!(regime, CavitationRegime::None);
    }

    #[test]
    fn test_regime_classification_stable() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(bubble, 3e4, Some(2e4), Some(20e3));

        let regime = classifier.classify_regime();
        assert_eq!(regime, CavitationRegime::Stable);
    }

    #[test]
    fn test_regime_classification_inertial() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(bubble, 1e5, Some(1e6), Some(20e3));

        let regime = classifier.classify_regime();
        assert_eq!(regime, CavitationRegime::Inertial);
    }

    #[test]
    fn test_mechanical_index() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(bubble, 1e5, Some(1e6), Some(1e6));

        let mi = classifier.mechanical_index().unwrap();
        assert!(mi > 0.0);
        assert!((mi - 1e6).abs() < 1.0);
    }

    #[test]
    fn test_damage_potential_ordering() {
        let bubble = create_test_bubble();

        let classifier_none = CavitationRegimeClassifier::new(bubble, 1e5, None, None);
        let classifier_stable = CavitationRegimeClassifier::new(bubble, 3e4, Some(2e4), Some(20e3));
        let classifier_inertial =
            CavitationRegimeClassifier::new(bubble, 1e5, Some(1e6), Some(20e3));

        let damage_none = classifier_none.damage_potential();
        let damage_stable = classifier_stable.damage_potential();
        let damage_inertial = classifier_inertial.damage_potential();

        assert!(damage_none < damage_stable);
        assert!(damage_stable < damage_inertial);
    }

    #[test]
    fn test_hemolysis_risk_inertial_highest() {
        let bubble = create_test_bubble();

        let classifier_stable = CavitationRegimeClassifier::new(bubble, 1e5, Some(1e4), Some(20e3));
        let classifier_inertial =
            CavitationRegimeClassifier::new(bubble, 1e5, Some(1e6), Some(20e3));

        let risk_stable = classifier_stable.hemolysis_risk();
        let risk_inertial = classifier_inertial.hemolysis_risk();

        assert!(risk_inertial > risk_stable);
        assert!(risk_inertial > 0.5);
    }

    #[test]
    fn test_full_analysis() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(bubble, 1e5, Some(5e5), Some(20e3));

        let analysis = classifier.analyze().unwrap();

        assert_eq!(analysis.regime, CavitationRegime::Inertial);
        assert!(analysis.blake_threshold > 0.0);
        assert!(analysis.inertial_threshold > 0.0);
        assert!(analysis.max_bubble_radius > bubble.initial_radius);
        assert!(analysis.sonoluminescence_probability > 0.5);
    }
}
