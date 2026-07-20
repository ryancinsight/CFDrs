use aequitas::systems::si::quantities::{MassDensity, SpecificHeatCapacity, ThermalConductivity};
use eunomia::RealField;
use proteus::ThermophysicalProperties;

use crate::error::Error;

pub(super) fn thermal_diffusivity<T: RealField>(
    density: T,
    specific_heat: T,
    thermal_conductivity: T,
) -> Result<T, Error> {
    ThermophysicalProperties::try_from_quantities(
        MassDensity::from_base(density),
        SpecificHeatCapacity::from_base(specific_heat),
        ThermalConductivity::from_base(thermal_conductivity),
    )
    .map(|properties| properties.thermal_diffusivity().into_base())
    .map_err(|error| Error::InvalidInput(error.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delegates_the_dimensional_diffusivity_law_to_proteus() {
        let actual =
            thermal_diffusivity(1_000.0_f64, 4_000.0, 0.6).expect("positive continuum properties");
        let expected = 0.6_f64 / (1_000.0 * 4_000.0);
        assert_eq!(actual.to_bits(), expected.to_bits());
    }

    #[test]
    fn reports_the_proteus_property_that_violated_its_contract() {
        let error = thermal_diffusivity(-1.0_f64, 4_000.0, 0.6)
            .expect_err("negative density is outside the material domain");
        match error {
            Error::InvalidInput(message) => {
                assert!(message.contains("MassDensity"));
                assert!(message.contains("FiniteNonNegative"));
            }
            other => panic!("unexpected error variant: {other}"),
        }
    }
}
