//! Solid material implementations

use super::traits::SolidProperties;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Elastic solid with constant isotropic properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElasticSolid<T: RealField + Copy> {
    /// Solid density [kg/m³]
    pub density: T,
    /// Young's modulus [Pa]
    pub youngs_modulus: T,
    /// Poisson's ratio [-]
    pub poissons_ratio: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal expansion coefficient [1/K]
    pub thermal_expansion: T,
}

impl<T: RealField + Copy> SolidProperties<T> for ElasticSolid<T> {
    fn density(&self) -> T {
        self.density
    }

    fn youngs_modulus(&self) -> T {
        self.youngs_modulus
    }

    fn poissons_ratio(&self) -> T {
        self.poissons_ratio
    }

    fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity
    }

    fn specific_heat(&self) -> T {
        self.specific_heat
    }

    fn thermal_expansion(&self) -> T {
        self.thermal_expansion
    }
}

impl<T: RealField + Copy> ElasticSolid<T> {
    /// Create steel properties (typical carbon steel)
    #[must_use]
    pub fn steel() -> Self
    where
        T: nalgebra::RealField + Copy + num_traits::FromPrimitive,
    {
        Self {
            density: T::from_f64(7850.0).unwrap(),
            youngs_modulus: T::from_f64(200e9).unwrap(),
            poissons_ratio: T::from_f64(0.3).unwrap(),
            thermal_conductivity: T::from_f64(50.0).unwrap(),
            specific_heat: T::from_f64(460.0).unwrap(),
            thermal_expansion: T::from_f64(12e-6).unwrap(),
        }
    }

    /// Create aluminum properties (typical 6061 alloy)
    #[must_use]
    pub fn aluminum() -> Self
    where
        T: nalgebra::RealField + Copy + num_traits::FromPrimitive,
    {
        Self {
            density: T::from_f64(2700.0).unwrap(),
            youngs_modulus: T::from_f64(70e9).unwrap(),
            poissons_ratio: T::from_f64(0.33).unwrap(),
            thermal_conductivity: T::from_f64(237.0).unwrap(),
            specific_heat: T::from_f64(900.0).unwrap(),
            thermal_expansion: T::from_f64(23e-6).unwrap(),
        }
    }
}
