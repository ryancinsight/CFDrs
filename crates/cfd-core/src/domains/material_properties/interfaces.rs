//! Interface properties between materials

use super::traits::InterfaceProperties;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Wetting properties for fluid-solid interfaces
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WettingProperties<T: RealField + Copy> {
    /// Contact angle in radians
    pub contact_angle: T,
    /// Advancing contact angle
    pub advancing_angle: T,
    /// Receding contact angle
    pub receding_angle: T,
}

/// Fluid-solid interface implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FluidSolidInterface<T: RealField + Copy> {
    /// Surface tension
    pub surface_tension: T,
    /// Contact angle
    pub contact_angle: T,
    /// Wetting properties
    pub wetting: WettingProperties<T>,
}

impl<T: RealField + Copy> InterfaceProperties<T> for FluidSolidInterface<T> {
    fn surface_tension(&self) -> T {
        self.surface_tension
    }

    fn contact_angle(&self) -> T {
        self.contact_angle
    }
}

impl<T: RealField + Copy> FluidSolidInterface<T> {
    /// Create water-air interface at 20°C
    pub fn water_air() -> Self
    where
        T: From<f64>,
    {
        let contact_angle = T::from(std::f64::consts::PI / 2.0); // 90 degrees
        Self {
            surface_tension: T::from(0.0728), // N/m
            contact_angle,
            wetting: WettingProperties {
                contact_angle,
                advancing_angle: T::from(std::f64::consts::PI * 100.0 / 180.0),
                receding_angle: T::from(std::f64::consts::PI * 80.0 / 180.0),
            },
        }
    }

    /// Create oil-water interface
    pub fn oil_water() -> Self
    where
        T: From<f64>,
    {
        let contact_angle = T::from(std::f64::consts::PI * 140.0 / 180.0);
        Self {
            surface_tension: T::from(0.030), // N/m
            contact_angle,
            wetting: WettingProperties {
                contact_angle,
                advancing_angle: T::from(std::f64::consts::PI * 150.0 / 180.0),
                receding_angle: T::from(std::f64::consts::PI * 130.0 / 180.0),
            },
        }
    }
}
