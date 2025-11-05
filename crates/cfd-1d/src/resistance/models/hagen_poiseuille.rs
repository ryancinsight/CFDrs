//! Hagen-Poiseuille resistance model for laminar flow in circular pipes.
//!
//! ## Hagen-Poiseuille Theorem
//!
//! **Theorem**: For steady, fully developed, laminar flow of an incompressible Newtonian fluid
//! in a straight circular pipe of constant cross-section, the pressure drop ΔP over a length L
//! is given by:
//!
//! ΔP = (32 μ L V) / D²
//!
//! where:
//! - μ is the dynamic viscosity [Pa·s]
//! - L is the pipe length [m]
//! - V is the average velocity [m/s]
//! - D is the pipe diameter [m]
//!
//! **Hydraulic Resistance Form**: R = ΔP / Q = (128 μ L) / (π D⁴)
//!
//! where Q is the volumetric flow rate [m³/s].
//!
//! ### Derivation
//!
//! Starting from the Navier-Stokes equations for axisymmetric flow:
//!
//! (1/ρ) ∇P = ν ∇²u + body forces
//!
//! For fully developed laminar flow (∂u/∂z = 0), the axial momentum equation reduces to:
//!
//! dP/dz = μ d²u/dr² + (μ/r) du/dr
//!
//! Integrating twice with boundary conditions u(r=R) = 0 and du/dr(r=0) = 0:
//!
//! u(r) = (1/(4μ)) (-dP/dz) (R² - r²)
//!
//! The Hagen-Poiseuille velocity profile is parabolic:
//!
//! u(r) = u_max (1 - (r/R)²)
//!
//! where u_max = (R²/(4μ)) (-dP/dz)
//!
//! Average velocity: V = (1/2) u_max = (R²/(8μ)) (-dP/dz)
//!
//! Pressure drop: ΔP = - (8μ L V) / R² = (32 μ L V) / D²
//!
//! ### Validity Conditions
//!
//! 1. **Laminar Flow**: Re < 2300 (typically Re < 2000 for safety)
//! 2. **Fully Developed Flow**: L/D > 10 (entrance effects negligible)
//! 3. **Newtonian Fluid**: Constant viscosity, no shear-thinning/thickening
//! 4. **Incompressible**: ρ constant, Mach < 0.3
//! 5. **No Body Forces**: Gravity and centrifugal forces negligible
//! 6. **Constant Properties**: Temperature constant, no property variations
//! 7. **Straight Pipe**: No curvature or bends
//! 8. **Smooth Walls**: No surface roughness effects
//!
//! ### References
//!
//! - Hagen, G. (1839). "Über die Bewegung der Flüssigkeiten in cylindrischen Röhren."
//!   *Poggendorff's Annalen der Physik und Chemie*, 46, 423-442.
//! - Poiseuille, J. L. M. (1840). "Recherches expérimentales sur le mouvement des liquides
//!   dans les tubes de très petits diamètres." *Mémoires présentés par divers savants à
//!   l'Académie Royale des Sciences de l'Institut de France*, 9, 433-544.
//! - White, F. M. (2006). *Viscous Fluid Flow* (3rd ed.). McGraw-Hill. Eq. 3-52.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::Result;
use cfd_core::fluid::Fluid;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants
const HAGEN_POISEUILLE_COEFFICIENT: f64 = 128.0;
const POWER_FOUR: f64 = 4.0;

/// Hagen-Poiseuille resistance model for circular channels
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HagenPoiseuilleModel<T: RealField + Copy> {
    /// Channel diameter [m]
    pub diameter: T,
    /// Channel length [m]
    pub length: T,
}

impl<T: RealField + Copy> HagenPoiseuilleModel<T> {
    /// Create a new Hagen-Poiseuille model
    pub fn new(diameter: T, length: T) -> Self {
        Self { diameter, length }
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> ResistanceModel<T>
    for HagenPoiseuilleModel<T>
{
    fn calculate_resistance(&self, fluid: &Fluid<T>, _conditions: &FlowConditions<T>) -> Result<T> {
        let viscosity = fluid.viscosity;
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(|| T::zero());

        let coefficient = T::from_f64(HAGEN_POISEUILLE_COEFFICIENT).unwrap_or_else(|| T::zero());

        // R = (128 * μ * L) / (π * D^4)
        let d4 = num_traits::Float::powf(
            self.diameter,
            T::from_f64(POWER_FOUR).unwrap_or_else(|| T::zero()),
        );
        let resistance = coefficient * viscosity * self.length / (pi * d4);

        Ok(resistance)
    }

    fn model_name(&self) -> &'static str {
        "Hagen-Poiseuille"
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::zero(),
            T::from_f64(cfd_core::constants::dimensionless::reynolds::PIPE_CRITICAL_LOWER)
                .unwrap_or_else(|| T::zero()),
        )
    }
}
