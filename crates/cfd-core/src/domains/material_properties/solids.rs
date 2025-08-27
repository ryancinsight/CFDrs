//! Solid material implementations

use super::traits::SolidProperties;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Elastic solid implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElasticSolid<T: RealField + Copy> {
    /// Solid density
    pub density: T,
    /// Young's modulus
    pub youngs_modulus: T,
    /// Poisson's ratio
    pub poissons_ratio: T,
    /// Thermal conductivity
    pub thermal_conductivity: T,
    /// Specific heat capacity
    pub specific_heat: T,
    /// Thermal expansion coefficient
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
    /// Create steel properties
    pub fn steel() -> Self
    where
        T: From<f64>,
    {
        Self {
            density: T::from(7850.0),
            youngs_modulus: T::from(200e9),
            poissons_ratio: T::from(0.3),
            thermal_conductivity: T::from(50.0),
            specific_heat: T::from(460.0),
            thermal_expansion: T::from(12e-6),
        }
    }

    /// Create aluminum properties
    pub fn aluminum() -> Self
    where
        T: From<f64>,
    {
        Self {
            density: T::from(2700.0),
            youngs_modulus: T::from(70e9),
            poissons_ratio: T::from(0.33),
            thermal_conductivity: T::from(237.0),
            specific_heat: T::from(900.0),
            thermal_expansion: T::from(23e-6),
        }
    }

    /// Create concrete properties
    pub fn concrete() -> Self
    where
        T: From<f64>,
    {
        Self {
            density: T::from(2400.0),
            youngs_modulus: T::from(30e9),
            poissons_ratio: T::from(0.2),
            thermal_conductivity: T::from(1.7),
            specific_heat: T::from(880.0),
            thermal_expansion: T::from(10e-6),
        }
    }
}
