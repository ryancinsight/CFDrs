//! Interface properties between materials

use super::traits::InterfaceProperties;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Wetting properties for fluid-solid interfaces
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WettingProperties<T: RealField + Copy> {
    /// Static contact angle [rad]
    pub contact_angle: T,
    /// Advancing contact angle [rad]
    pub advancing_angle: T,
    /// Receding contact angle [rad]
    pub receding_angle: T,
}

/// Fluid-solid interface implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FluidSolidInterface<T: RealField + Copy> {
    /// Surface tension [N/m]
    pub surface_tension: T,
    /// Static contact angle [rad]
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
    #[must_use]
    pub fn water_air() -> Self
    where
        T: nalgebra::RealField + Copy + num_traits::FromPrimitive,
    {
        let pi = T::from_f64(std::f64::consts::PI).unwrap();
        let contact_angle = pi / T::from_f64(2.0).unwrap(); // 90 degrees
        Self {
            surface_tension: T::from_f64(0.0728).unwrap(),
            contact_angle,
            wetting: WettingProperties {
                contact_angle,
                advancing_angle: pi * T::from_f64(100.0 / 180.0).unwrap(),
                receding_angle: pi * T::from_f64(80.0 / 180.0).unwrap(),
            },
        }
    }
}
