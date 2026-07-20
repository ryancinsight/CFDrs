use aequitas::systems::si::quantities::{
    MassDensity, ReciprocalTemperature, SpecificHeatCapacity, ThermalConductivity,
    ThermodynamicTemperature,
};
use eunomia::RealField;
use proteus::{
    ConstantResponse, ConstitutiveLaw, LinearResponse, ResponseSet, TemperatureLaw,
    ThermophysicalProperties,
};

use crate::error::Error;

/// Validate and compute thermal diffusivity via the Proteus thermophysical contract.
///
/// Returns `Ok(α)` where `α = k / (ρ·cₚ)` \[m²/s], or an error if any
/// quantity violates the `FiniteNonNegative` constraint.
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

/// Validate the thermophysical subset (density, specific heat, thermal conductivity)
/// via the Proteus `FiniteNonNegative` property contract, without computing any
/// derived quantity.
///
/// Returns `Ok(())` when all three quantities satisfy the Proteus constraint, or a
/// descriptive `Error::InvalidInput` naming the violated property.
pub(super) fn validate_thermophysical_subset<T: RealField>(
    density: T,
    specific_heat: T,
    thermal_conductivity: T,
) -> Result<(), Error> {
    ThermophysicalProperties::try_from_quantities(
        MassDensity::from_base(density),
        SpecificHeatCapacity::from_base(specific_heat),
        ThermalConductivity::from_base(thermal_conductivity),
    )
    .map(|_| ())
    .map_err(|error| Error::InvalidInput(error.to_string()))
}

/// Evaluate a linearly temperature-dependent density through Proteus.
///
/// The CFD thermal-expansion convention
/// `rho(T) = rho_ref * (1 - beta * (T - T_ref))` maps to a Proteus
/// [`LinearResponse`] coefficient of `-beta`. Heat capacity and conductivity
/// remain invariant ZST responses.
pub(super) fn linear_density_at<T: RealField>(
    density_ref: T,
    thermal_expansion: T,
    reference_temperature: T,
    specific_heat: T,
    thermal_conductivity: T,
    temperature: T,
) -> Result<T, Error> {
    let reference_properties = ThermophysicalProperties::try_from_quantities(
        MassDensity::from_base(density_ref),
        SpecificHeatCapacity::from_base(specific_heat),
        ThermalConductivity::from_base(thermal_conductivity),
    )
    .map_err(invalid_input)?;
    let density_response =
        LinearResponse::new(ReciprocalTemperature::from_base(-thermal_expansion))
            .map_err(invalid_input)?;
    let law = TemperatureLaw::new(
        reference_properties,
        ThermodynamicTemperature::from_base(reference_temperature),
        ResponseSet::new(density_response, ConstantResponse, ConstantResponse),
    )
    .map_err(invalid_input)?;
    let evaluation_temperature = ThermodynamicTemperature::from_base(temperature);

    law.properties(&evaluation_temperature)
        .map(|properties| *properties.density().quantity().as_base())
        .map_err(invalid_input)
}

fn invalid_input(error: impl ToString) -> Error {
    Error::InvalidInput(error.to_string())
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

    #[test]
    fn linear_density_matches_the_closed_form_oracle() {
        let actual = linear_density_at(1_000.0_f64, 2.0e-4, 300.0, 4_180.0, 0.6, 310.0)
            .expect("finite positive thermophysical state");
        let expected = 1_000.0_f64 * (1.0 - 2.0e-4 * (310.0 - 300.0));
        // The Proteus response uses one FMA while the independent oracle uses
        // separate product and subtraction. Four native roundings bound both
        // evaluation paths and the final density scaling.
        let rounding = 4.0 * f64::EPSILON * expected.abs();

        assert!((actual - expected).abs() <= rounding);
    }

    #[test]
    fn linear_density_rejects_a_nonphysical_response() {
        let error = linear_density_at(1_000.0_f64, 0.2, 300.0, 4_180.0, 0.6, 310.0)
            .expect_err("the response produces a negative mass density");

        match error {
            Error::InvalidInput(message) => {
                assert!(message.contains("MassDensity"));
                assert!(message.contains("FiniteNonNegative"));
            }
            other => panic!("unexpected error variant: {other}"),
        }
    }
}
