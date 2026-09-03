//! Solid material implementations

use super::traits::SolidProperties;
use crate::error::{Error, Result};
use aequitas::systems::si::quantities::{
    Dimensionless, MassDensity, Pressure, ReciprocalTemperature, SpecificHeatCapacity,
    ThermalConductivity,
};
use eunomia::RealField;
use proteus::{IsotropicSolid, NamedIsotropicSolid};

/// Elastic solid with constant isotropic properties.
///
/// Elastic state is held as a [`proteus::IsotropicSolid`], so the
/// `(E, nu) <-> (lambda, mu) <-> (c_p, c_s)` conversion contract and the named
/// material catalog have one owner. CFDrs adds only the thermal properties its
/// flow and conjugate-heat paths consume, and stores no material constant.
#[derive(Debug, Clone)]
pub struct ElasticSolid<T> {
    /// Validated isotropic elastic state and density.
    elastic: IsotropicSolid<T>,
    /// Thermal conductivity [W/(m·K)]
    thermal_conductivity: ThermalConductivity<T>,
    /// Specific heat capacity [J/(kg·K)]
    specific_heat: SpecificHeatCapacity<T>,
    /// Thermal expansion coefficient [1/K]
    thermal_expansion: ReciprocalTemperature<T>,
}

impl<T: RealField> ElasticSolid<T> {
    /// Build a solid from a Proteus catalog entry.
    ///
    /// Elastic constants, density, and thermal properties all come from the
    /// provider.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the scalar type cannot represent
    /// the entry's published constants.
    pub fn from_catalog(entry: NamedIsotropicSolid) -> Result<Self> {
        let elastic = entry.solid::<T>().map_err(|invalid| {
            Error::InvalidInput(format!("{} elastic state: {invalid}", entry.name()))
        })?;
        let thermophysical = entry.thermophysical::<T>().map_err(|invalid| {
            Error::InvalidInput(format!("{} thermal state: {invalid}", entry.name()))
        })?;

        Ok(Self {
            elastic,
            thermal_conductivity: *thermophysical.thermal_conductivity().quantity(),
            specific_heat: *thermophysical.specific_heat_capacity().quantity(),
            thermal_expansion: entry.thermal_expansion::<T>(),
        })
    }

    /// Compose a solid from an explicit elastic state and thermal properties.
    #[must_use]
    pub const fn new(
        elastic: IsotropicSolid<T>,
        thermal_conductivity: ThermalConductivity<T>,
        specific_heat: SpecificHeatCapacity<T>,
        thermal_expansion: ReciprocalTemperature<T>,
    ) -> Self {
        Self {
            elastic,
            thermal_conductivity,
            specific_heat,
            thermal_expansion,
        }
    }

    /// Borrow the validated elastic state.
    #[must_use]
    pub const fn elastic(&self) -> &IsotropicSolid<T> {
        &self.elastic
    }
}

impl<T: RealField> SolidProperties<T> for ElasticSolid<T> {
    fn density(&self) -> MassDensity<T> {
        *self.elastic.density().quantity()
    }

    fn youngs_modulus(&self) -> Pressure<T> {
        self.elastic.moduli().youngs_modulus()
    }

    fn poissons_ratio(&self) -> Dimensionless<T> {
        self.elastic.moduli().poissons_ratio()
    }

    fn shear_modulus(&self) -> Pressure<T> {
        *self.elastic.moduli().shear_modulus()
    }

    fn thermal_conductivity(&self) -> ThermalConductivity<T> {
        self.thermal_conductivity
    }

    fn specific_heat(&self) -> SpecificHeatCapacity<T> {
        self.specific_heat
    }

    fn thermal_expansion(&self) -> ReciprocalTemperature<T> {
        self.thermal_expansion
    }
}

#[cfg(test)]
mod tests {
    use super::{ElasticSolid, NamedIsotropicSolid, SolidProperties};

    fn tolerance(expected: f64) -> f64 {
        expected.abs() * f64::EPSILON * 16.0
    }

    #[test]
    fn shear_modulus_matches_the_isotropic_definition() {
        // An independent oracle: mu = E / (2 (1 + nu)) derived from whatever
        // (E, nu) the provider publishes, never from literals restating them.
        let steel = ElasticSolid::<f64>::from_catalog(NamedIsotropicSolid::CarbonSteel)
            .expect("catalog entry is representable at f64");
        let young = steel.youngs_modulus().into_base();
        let poisson = steel.poissons_ratio().into_base();
        let expected = young / (2.0 * (1.0 + poisson));

        assert!((steel.shear_modulus().into_base() - expected).abs() <= tolerance(expected));
    }

    #[test]
    fn every_accessor_relays_the_provider_faithfully() {
        // The contract under test is relay fidelity, not the provider's
        // values. Asserting literals here would duplicate the provider's own
        // catalog test and break whenever a published constant is corrected -
        // which is exactly what happened when 6061-T6 moved to its published
        // 68.9 GPa.
        for entry in [
            NamedIsotropicSolid::CarbonSteel,
            NamedIsotropicSolid::StainlessSteel316L,
            NamedIsotropicSolid::Aluminium6061,
            NamedIsotropicSolid::TitaniumGrade5,
        ] {
            let solid = ElasticSolid::<f64>::from_catalog(entry).expect("representable at f64");
            let provider = entry.solid::<f64>().expect("representable at f64");
            let thermophysical = entry.thermophysical::<f64>().expect("representable at f64");

            for (actual, expected, label) in [
                (
                    solid.density().into_base(),
                    provider.density().quantity().into_base(),
                    "density",
                ),
                (
                    solid.youngs_modulus().into_base(),
                    provider.moduli().youngs_modulus().into_base(),
                    "E",
                ),
                (
                    solid.poissons_ratio().into_base(),
                    provider.moduli().poissons_ratio().into_base(),
                    "nu",
                ),
                (
                    solid.shear_modulus().into_base(),
                    *provider.moduli().shear_modulus().as_base(),
                    "mu",
                ),
                (
                    solid.thermal_conductivity().into_base(),
                    thermophysical.thermal_conductivity().quantity().into_base(),
                    "k",
                ),
                (
                    solid.specific_heat().into_base(),
                    thermophysical
                        .specific_heat_capacity()
                        .quantity()
                        .into_base(),
                    "c_p",
                ),
                (
                    solid.thermal_expansion().into_base(),
                    entry.thermal_expansion::<f64>().into_base(),
                    "alpha",
                ),
            ] {
                assert!(
                    (actual - expected).abs() <= tolerance(expected),
                    "{} {label} does not relay the provider value",
                    entry.name()
                );
            }
        }
    }

    #[test]
    fn distinct_steel_grades_are_not_conflated() {
        // The catalog this replaced said "steel" and meant carbon steel;
        // kwavers said "steel" and meant 316L. Keying by grade keeps that
        // substitution unrepresentable.
        let carbon = ElasticSolid::<f64>::from_catalog(NamedIsotropicSolid::CarbonSteel)
            .expect("representable");
        let austenitic = ElasticSolid::<f64>::from_catalog(NamedIsotropicSolid::StainlessSteel316L)
            .expect("representable");

        assert!((carbon.density().into_base() - austenitic.density().into_base()).abs() > 1.0);
    }
}
