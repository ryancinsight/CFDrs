//! Solid material implementations

use super::traits::SolidProperties;
use aequitas::systems::si::quantities::{
    Dimensionless, MassDensity, Pressure, ReciprocalTemperature, SpecificHeatCapacity,
    ThermalConductivity,
};
use eunomia::FloatElement;
use serde::{Deserialize, Serialize};

/// Elastic solid with constant isotropic properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElasticSolid<T> {
    /// Solid density [kg/m³]
    pub density: MassDensity<T>,
    /// Young's modulus \[Pa]
    pub youngs_modulus: Pressure<T>,
    /// Poisson's ratio [-]
    pub poissons_ratio: Dimensionless<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal expansion coefficient [1/K]
    pub thermal_expansion: ReciprocalTemperature<T>,
}

impl<T: FloatElement + Copy> SolidProperties<T> for ElasticSolid<T> {
    fn density(&self) -> MassDensity<T> {
        self.density
    }

    fn youngs_modulus(&self) -> Pressure<T> {
        self.youngs_modulus
    }

    fn poissons_ratio(&self) -> Dimensionless<T> {
        self.poissons_ratio
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

impl<T: FloatElement + Copy> ElasticSolid<T> {
    /// Create steel properties (typical carbon steel)
    #[must_use]
    pub fn steel() -> Self {
        Self {
            density: MassDensity::from_base(<T as FloatElement>::from_f64(7850.0)),
            youngs_modulus: Pressure::from_base(<T as FloatElement>::from_f64(200e9)),
            poissons_ratio: Dimensionless::from_base(<T as FloatElement>::from_f64(0.3)),
            thermal_conductivity: ThermalConductivity::from_base(<T as FloatElement>::from_f64(
                50.0,
            )),
            specific_heat: SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(460.0)),
            thermal_expansion: ReciprocalTemperature::from_base(<T as FloatElement>::from_f64(
                12e-6,
            )),
        }
    }

    /// Create aluminum properties (typical 6061 alloy)
    #[must_use]
    pub fn aluminum() -> Self {
        Self {
            density: MassDensity::from_base(<T as FloatElement>::from_f64(2700.0)),
            youngs_modulus: Pressure::from_base(<T as FloatElement>::from_f64(70e9)),
            poissons_ratio: Dimensionless::from_base(<T as FloatElement>::from_f64(0.33)),
            thermal_conductivity: ThermalConductivity::from_base(<T as FloatElement>::from_f64(
                237.0,
            )),
            specific_heat: SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(900.0)),
            thermal_expansion: ReciprocalTemperature::from_base(<T as FloatElement>::from_f64(
                23e-6,
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{ElasticSolid, SolidProperties};

    #[test]
    fn steel_shear_modulus_matches_isotropic_definition() {
        let steel = ElasticSolid::<f64>::steel();
        let expected: f64 = 200e9 / (2.0 * (1.0 + 0.3));
        let tolerance = expected.abs() * f64::EPSILON * 4.0;

        assert!((steel.shear_modulus().into_base() - expected).abs() <= tolerance);
    }

    #[test]
    fn aluminum_constructor_preserves_material_constants() {
        let aluminum = ElasticSolid::<f64>::aluminum();

        assert_eq!(aluminum.density.into_base(), 2700.0);
        assert_eq!(aluminum.youngs_modulus.into_base(), 70e9);
        assert_eq!(aluminum.poissons_ratio.into_base(), 0.33);
        assert_eq!(aluminum.thermal_conductivity.into_base(), 237.0);
        assert_eq!(aluminum.specific_heat.into_base(), 900.0);
        assert_eq!(aluminum.thermal_expansion.into_base(), 23e-6);
    }
}
