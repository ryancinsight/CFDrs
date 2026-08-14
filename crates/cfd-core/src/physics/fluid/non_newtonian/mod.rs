//! Non-Newtonian fluid models.
//!
//! References:
//! - Bird, R.B., Armstrong, R.C., Hassager, O. (1987) "Dynamics of Polymeric Liquids"
//! - Chhabra, R.P., Richardson, J.F. (2008) "Non-Newtonian Flow and Applied Rheology"

#![cfg_attr(
    test,
    expect(
        clippy::unwrap_used,
        reason = "ratchet CFDRS-UNWRAP-1: pre-existing debt"
    )
)]

/// Bingham plastic fluid model.
mod bingham;
/// Carreau–Yasuda fluid model for blood and polymer solutions.
mod carreau_yasuda;
/// Casson fluid model for blood.
mod casson;
/// Herschel–Bulkley generalised Bingham plastic model.
mod herschel_bulkley;
/// Power-law (Ostwald–de Waele) fluid model.
mod power_law;

pub use bingham::BinghamPlastic;
pub use carreau_yasuda::CarreauYasuda;
pub use casson::Casson;
pub use herschel_bulkley::HerschelBulkley;
pub use power_law::PowerLawFluid;

#[cfg(test)]
mod tests {
    use super::*;
    use aequitas::systems::si::quantities::{
        Dimensionless, DynamicViscosity, MassDensity, MolarEnergy, Pressure, ReciprocalTime,
        SpecificHeatCapacity, ThermalConductivity, ThermodynamicTemperature, Time, Velocity,
    };
    use eunomia::assert_relative_eq;

    #[test]
    fn test_herschel_bulkley_constant() {
        let fluid = HerschelBulkley::<f64>::new(
            "Test Fluid".to_string(),
            MassDensity::from_base(1000.0),
            Pressure::from_base(10.0),
            5.0,
            Dimensionless::from_base(0.5),
            SpecificHeatCapacity::from_base(4000.0),
            ThermalConductivity::from_base(0.6),
            Velocity::from_base(1500.0),
            ReciprocalTime::from_base(10.0),
        );

        use crate::physics::fluid::traits::Fluid;
        let props = fluid.properties_at(300.0, 101325.0).unwrap();
        assert_relative_eq!(props.dynamic_viscosity.into_base(), 2.581_138_830_084_19);
        assert!(!fluid.is_temperature_dependent());
    }

    #[test]
    fn test_herschel_bulkley_temperature_dependent() {
        let t_ref = 300.0;
        let ea_k = 5000.0;

        let fluid = HerschelBulkley::<f64>::new(
            "Test Fluid".to_string(),
            MassDensity::from_base(1000.0),
            Pressure::from_base(10.0),
            5.0,
            Dimensionless::from_base(0.5),
            SpecificHeatCapacity::from_base(4000.0),
            ThermalConductivity::from_base(0.6),
            Velocity::from_base(1500.0),
            ReciprocalTime::from_base(10.0),
        )
        .with_temperature_dependence(
            ThermodynamicTemperature::from_base(t_ref),
            Some(MolarEnergy::from_base(ea_k)),
            None,
        );

        use crate::physics::fluid::traits::Fluid;
        assert!(fluid.is_temperature_dependent());
        assert_eq!(
            fluid
                .reference_temperature()
                .map(aequitas::Quantity::into_base),
            Some(t_ref)
        );

        let props_ref = fluid.properties_at(t_ref, 101325.0).unwrap();
        assert_relative_eq!(
            props_ref.dynamic_viscosity.into_base(),
            2.581_138_830_084_19
        );

        let t_high = 350.0;
        let props_high = fluid.properties_at(t_high, 101325.0).unwrap();

        let r = 8.314462618;
        let arg = (ea_k / r) * (1.0 / t_high - 1.0 / t_ref);
        let k_high = 5.0 * arg.exp();

        let shear_rate = 10.0_f64;
        let power_law_term = k_high * shear_rate.powf(0.5 - 1.0);
        let yield_term = 10.0 / shear_rate;
        let expected_viscosity = yield_term + power_law_term;

        assert_relative_eq!(props_high.dynamic_viscosity.into_base(), expected_viscosity);
        assert!(props_high.dynamic_viscosity.into_base() < props_ref.dynamic_viscosity.into_base());
    }

    #[test]
    fn test_power_law_constant() {
        let fluid = PowerLawFluid::<f64>::new(
            "Test Fluid".to_string(),
            MassDensity::from_base(1000.0),
            5.0,
            Dimensionless::from_base(0.5),
            SpecificHeatCapacity::from_base(4000.0),
            ThermalConductivity::from_base(0.6),
            Velocity::from_base(1500.0),
            ReciprocalTime::from_base(10.0),
        );

        use crate::physics::fluid::traits::Fluid;
        let props = fluid.properties_at(300.0, 101325.0).unwrap();
        assert_relative_eq!(props.dynamic_viscosity.into_base(), 1.5811388300841898);
        assert!(!fluid.is_temperature_dependent());
    }

    #[test]
    fn test_power_law_temperature_dependent() {
        let t_ref = 300.0;
        let ea_k = 5000.0;

        let fluid = PowerLawFluid::<f64>::new(
            "Test Fluid".to_string(),
            MassDensity::from_base(1000.0),
            5.0,
            Dimensionless::from_base(0.5),
            SpecificHeatCapacity::from_base(4000.0),
            ThermalConductivity::from_base(0.6),
            Velocity::from_base(1500.0),
            ReciprocalTime::from_base(10.0),
        )
        .with_temperature_dependence(
            ThermodynamicTemperature::from_base(t_ref),
            Some(MolarEnergy::from_base(ea_k)),
        );

        use crate::physics::fluid::traits::Fluid;
        assert!(fluid.is_temperature_dependent());
        assert_eq!(
            fluid
                .reference_temperature()
                .map(aequitas::Quantity::into_base),
            Some(t_ref)
        );

        let props_ref = fluid.properties_at(t_ref, 101325.0).unwrap();
        assert_relative_eq!(props_ref.dynamic_viscosity.into_base(), 1.5811388300841898);

        let t_high = 350.0;
        let props_high = fluid.properties_at(t_high, 101325.0).unwrap();

        let r = 8.314462618;
        let arg = (ea_k / r) * (1.0 / t_high - 1.0 / t_ref);
        let k_high = 5.0 * arg.exp();

        let shear_rate = 10.0_f64;
        let expected_viscosity = k_high * shear_rate.powf(0.5 - 1.0);

        assert_relative_eq!(props_high.dynamic_viscosity.into_base(), expected_viscosity);
        assert!(props_high.dynamic_viscosity.into_base() < props_ref.dynamic_viscosity.into_base());
    }

    #[test]
    fn fixed_dimension_model_parameters_preserve_value_semantics() {
        let bingham = BinghamPlastic::new(
            "Bingham".to_string(),
            MassDensity::from_base(1000.0),
            Pressure::from_base(10.0),
            DynamicViscosity::from_base(5.0),
            SpecificHeatCapacity::from_base(4000.0),
            ThermalConductivity::from_base(0.6),
            Velocity::from_base(1500.0),
            ReciprocalTime::from_base(10.0),
        );
        assert_relative_eq!(bingham.apparent_viscosity(10.0), 6.0);
        assert!(bingham.is_yielded(10.1));

        let casson = Casson::new(
            "Casson".to_string(),
            MassDensity::from_base(1000.0),
            Pressure::from_base(10.0),
            DynamicViscosity::from_base(5.0),
            SpecificHeatCapacity::from_base(4000.0),
            ThermalConductivity::from_base(0.6),
            Velocity::from_base(1500.0),
        );
        assert_relative_eq!(
            casson.apparent_viscosity(10.0),
            (5.0_f64.sqrt() + 1.0).powi(2)
        );

        let carreau = CarreauYasuda::new(
            "Carreau".to_string(),
            MassDensity::from_base(1060.0),
            DynamicViscosity::from_base(0.056),
            DynamicViscosity::from_base(0.0035),
            Time::from_base(3.313),
            Dimensionless::from_base(0.3568),
            Dimensionless::from_base(2.0),
            SpecificHeatCapacity::from_base(3600.0),
            ThermalConductivity::from_base(0.5),
            Velocity::from_base(1540.0),
        );
        let apparent_viscosity = carreau.apparent_viscosity(100.0);
        assert!(apparent_viscosity > 0.0035);
        assert!(apparent_viscosity < 0.056);
    }
}
