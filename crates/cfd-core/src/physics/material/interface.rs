//! Interface properties between materials

use super::traits::InterfaceProperties;
use aequitas::systems::si::quantities::{Angle, SurfaceTension};
use eunomia::FloatElement;
use serde::{Deserialize, Serialize};

/// Wetting properties for fluid-solid interfaces
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WettingProperties<T: FloatElement + Copy> {
    /// Static contact angle \[rad]
    pub contact_angle: Angle<T>,
    /// Advancing contact angle \[rad]
    pub advancing_angle: Angle<T>,
    /// Receding contact angle \[rad]
    pub receding_angle: Angle<T>,
}

/// Fluid-solid interface implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FluidSolidInterface<T: FloatElement + Copy> {
    /// Surface tension [N/m]
    pub surface_tension: SurfaceTension<T>,
    /// Static contact angle \[rad]
    pub contact_angle: Angle<T>,
    /// Wetting properties
    pub wetting: WettingProperties<T>,
}

impl<T: FloatElement + Copy> InterfaceProperties<T> for FluidSolidInterface<T> {
    fn surface_tension(&self) -> SurfaceTension<T> {
        self.surface_tension
    }

    fn contact_angle(&self) -> Angle<T> {
        self.contact_angle
    }
}

impl<T: FloatElement + Copy> FluidSolidInterface<T> {
    /// Create water-air interface at 20°C
    #[must_use]
    pub fn water_air() -> Self {
        let pi = <T as FloatElement>::from_f64(std::f64::consts::PI);
        let contact_angle = Angle::from_base(pi / <T as FloatElement>::from_f64(2.0));
        Self {
            surface_tension: SurfaceTension::from_base(<T as FloatElement>::from_f64(
                crate::physics::cavitation::constants::SURFACE_TENSION_WATER_VALUE,
            )),
            contact_angle,
            wetting: WettingProperties {
                contact_angle,
                advancing_angle: Angle::from_base(
                    pi * <T as FloatElement>::from_f64(100.0 / 180.0),
                ),
                receding_angle: Angle::from_base(pi * <T as FloatElement>::from_f64(80.0 / 180.0)),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{FluidSolidInterface, InterfaceProperties};

    #[test]
    fn water_air_interface_preserves_contact_angles() {
        let interface = FluidSolidInterface::<f64>::water_air();

        assert_eq!(interface.surface_tension.into_base(), 0.0728);
        assert_eq!(
            interface.contact_angle.into_base(),
            std::f64::consts::PI / 2.0
        );
        assert_eq!(
            interface.wetting.advancing_angle.into_base(),
            std::f64::consts::PI * (100.0 / 180.0)
        );
        assert_eq!(
            interface.wetting.receding_angle.into_base(),
            std::f64::consts::PI * (80.0 / 180.0)
        );
    }

    #[test]
    fn adhesion_energy_uses_eunomia_trigonometry() {
        let interface = FluidSolidInterface::<f64>::water_air();
        let tolerance = interface.surface_tension.into_base() * f64::EPSILON * 16.0;

        assert!(interface.adhesion_energy().into_base().abs() <= tolerance);
    }
}
