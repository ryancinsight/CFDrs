//! Surface properties and wettability characteristics

use crate::scalar::Cfd1dScalar;
use aequitas::systems::si::quantities::{Angle, EnergyPerArea, Length};
use serde::{Deserialize, Serialize};

/// Surface properties affecting flow
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurfaceProperties<T: Cfd1dScalar + Copy> {
    /// Surface roughness \[m]
    pub roughness: Length<T>,
    /// Contact angle \[radians]
    pub contact_angle: Option<Angle<T>>,
    /// Surface energy [J/m²]
    pub surface_energy: Option<EnergyPerArea<T>>,
    /// Hydrophobic/hydrophilic nature
    pub wettability: Wettability,
}

/// Surface wettability characteristics
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Wettability {
    /// Hydrophilic surface
    Hydrophilic,
    /// Hydrophobic surface
    Hydrophobic,
    /// Superhydrophilic surface
    Superhydrophilic,
    /// Superhydrophobic surface
    Superhydrophobic,
}

#[cfg(test)]
mod tests {
    use super::{SurfaceProperties, Wettability};
    use aequitas::systems::si::quantities::{Angle, EnergyPerArea, Length};

    #[test]
    fn physical_surface_metrics_preserve_si_values() {
        let surface = SurfaceProperties::<f64> {
            roughness: Length::from_base(2.5e-7),
            contact_angle: Some(Angle::from_base(1.0)),
            surface_energy: Some(EnergyPerArea::from_base(0.12)),
            wettability: Wettability::Hydrophilic,
        };

        assert_eq!(surface.roughness.into_base(), 2.5e-7);
        assert_eq!(
            surface
                .contact_angle
                .expect("contact angle is present")
                .into_base(),
            1.0
        );
        assert_eq!(
            surface
                .surface_energy
                .expect("surface energy is present")
                .into_base(),
            0.12
        );
    }
}
