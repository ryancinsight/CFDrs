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
mod tests {
    use super::*;
    use crate::physics::cavitation::rayleigh_plesset::RayleighPlesset;
    use aequitas::systems::si::quantities::{
        DynamicViscosity, Frequency, Length, MassDensity, Pressure, SurfaceTension,
    };

    fn pressure(value: f64) -> Pressure<f64> {
        Pressure::from_base(value)
    }

    fn frequency(value: f64) -> Frequency<f64> {
        Frequency::from_base(value)
    }

    fn create_test_bubble() -> RayleighPlesset<f64> {
        RayleighPlesset {
            initial_radius: Length::from_base(1e-6),
            liquid_density: MassDensity::from_base(998.0),
            liquid_viscosity: DynamicViscosity::from_base(1.002e-3),
            surface_tension: SurfaceTension::from_base(0.0728),
            vapor_pressure: Pressure::from_base(2339.0),
            polytropic_index: 1.4,
        }
    }

    #[test]
    fn test_regime_classification_none() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(bubble, pressure(1e5), None, None);

        let regime = classifier.classify_regime();
        assert_eq!(regime, CavitationRegime::None);
    }

    #[test]
    fn test_regime_classification_stable() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(
            bubble,
            pressure(3e4),
            Some(pressure(2e4)),
            Some(frequency(20e3)),
        );

        let regime = classifier.classify_regime();
        assert_eq!(regime, CavitationRegime::Stable);
    }

    #[test]
    fn test_regime_classification_inertial() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(
            bubble,
            pressure(1e5),
            Some(pressure(1e6)),
            Some(frequency(20e3)),
        );

        let regime = classifier.classify_regime();
        assert_eq!(regime, CavitationRegime::Inertial);
    }

    #[test]
    fn test_mechanical_index() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(
            bubble,
            pressure(1e5),
            Some(pressure(1e6)),
            Some(frequency(1e6)),
        );

        let mi = classifier.mechanical_index().expect("expected value");
        assert!(mi > 0.0);
        assert!((mi - 1e6).abs() < 1.0);
    }

    #[test]
    fn test_blake_threshold_matches_critical_radius_contract() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(
            bubble,
            pressure(1e5),
            Some(pressure(1e6)),
            Some(frequency(1e6)),
        );

        let critical_radius = bubble.blake_critical_radius(1e5);
        let expected = bubble.vapor_pressure.into_base()
            + (4.0 / 3.0) * bubble.surface_tension.into_base() / critical_radius;
        let actual = classifier.blake_threshold();

        assert!(critical_radius > 0.0);
        assert!((actual.into_base() - expected).abs() < 1e-12);
    }

    #[test]
    fn test_damage_potential_ordering() {
        let bubble = create_test_bubble();

        let classifier_none = CavitationRegimeClassifier::new(bubble, pressure(1e5), None, None);
        let classifier_stable = CavitationRegimeClassifier::new(
            bubble,
            pressure(3e4),
            Some(pressure(2e4)),
            Some(frequency(20e3)),
        );
        let classifier_inertial = CavitationRegimeClassifier::new(
            bubble,
            pressure(1e5),
            Some(pressure(1e6)),
            Some(frequency(20e3)),
        );

        let damage_none = classifier_none.damage_potential();
        let damage_stable = classifier_stable.damage_potential();
        let damage_inertial = classifier_inertial.damage_potential();

        assert!(damage_none < damage_stable);
        assert!(damage_stable < damage_inertial);
    }

    #[test]
    fn test_hemolysis_risk_inertial_highest() {
        let bubble = create_test_bubble();

        let classifier_stable = CavitationRegimeClassifier::new(
            bubble,
            pressure(1e5),
            Some(pressure(1e4)),
            Some(frequency(20e3)),
        );
        let classifier_inertial = CavitationRegimeClassifier::new(
            bubble,
            pressure(1e5),
            Some(pressure(1e6)),
            Some(frequency(20e3)),
        );

        let risk_stable = classifier_stable.hemolysis_risk();
        let risk_inertial = classifier_inertial.hemolysis_risk();

        assert!(risk_inertial > risk_stable);
        assert!(risk_inertial > 0.5);
    }

    #[test]
    fn test_full_analysis() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(
            bubble,
            pressure(1e5),
            Some(pressure(5e5)),
            Some(frequency(20e3)),
        );

        let analysis = classifier.analyze().expect("expected value");

        assert_eq!(analysis.regime, CavitationRegime::Inertial);
        assert!(analysis.blake_threshold.into_base() > 0.0);
        assert!(analysis.inertial_threshold.into_base() > 0.0);
        assert!(analysis.max_bubble_radius > bubble.initial_radius);
        assert!(analysis.sonoluminescence_probability > 0.5);
    }
}
