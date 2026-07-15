//! Surface properties and wettability characteristics

use crate::scalar::Cfd1dScalar;
use serde::{Deserialize, Serialize};

/// Surface properties affecting flow
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurfaceProperties<T: Cfd1dScalar + Copy> {
    /// Surface roughness \[m]
    pub roughness: T,
    /// Contact angle \[radians]
    pub contact_angle: Option<T>,
    /// Surface energy [J/m²]
    pub surface_energy: Option<T>,
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
