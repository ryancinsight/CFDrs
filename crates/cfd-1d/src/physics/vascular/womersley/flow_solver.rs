//! Complete Womersley flow solver for arterial segments.
//!
//! Provides time-varying flow solutions for vessel segments with
//! given inlet conditions and geometry.

use super::profile::WomersleyProfile;
use super::WomersleyNumber;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Complete Womersley flow solver for arterial segments
///
/// Provides time-varying flow solutions for vessel segments with
/// given inlet conditions and geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WomersleyFlow<T: RealField + Copy> {
    /// Vessel radius [m]
    pub radius: T,
    /// Vessel length [m]
    pub length: T,
    /// Fluid density [kg/m³]
    pub density: T,
    /// Dynamic viscosity [Pa·s]
    pub viscosity: T,
    /// Angular frequency [rad/s]
    pub omega: T,
    /// Inlet pressure amplitude [Pa]
    pub inlet_pressure_amplitude: T,
    /// Mean pressure gradient [Pa/m]
    pub mean_pressure_gradient: T,
}

impl<T: RealField + FromPrimitive + Copy> WomersleyFlow<T> {
    /// Create new Womersley flow solver
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        radius: T,
        length: T,
        density: T,
        viscosity: T,
        omega: T,
        inlet_pressure_amplitude: T,
        mean_pressure_gradient: T,
    ) -> Self {
        Self {
            radius,
            length,
            density,
            viscosity,
            omega,
            inlet_pressure_amplitude,
            mean_pressure_gradient,
        }
    }

    /// Get Womersley number for this flow configuration
    pub fn womersley_number(&self) -> WomersleyNumber<T> {
        WomersleyNumber::new(self.radius, self.omega, self.density, self.viscosity)
    }

    /// Get velocity profile calculator
    pub fn profile(&self) -> WomersleyProfile<T> {
        let pressure_gradient_amplitude = self.inlet_pressure_amplitude / self.length;
        WomersleyProfile::new(self.womersley_number(), pressure_gradient_amplitude)
    }

    /// Calculate total velocity (mean + pulsatile) at position and time
    pub fn velocity(&self, xi: T, t: T) -> T {
        // Mean (steady) component - Poiseuille
        let r = self.radius;
        let mu = self.viscosity;
        let four = T::one() + T::one() + T::one() + T::one();
        let u_mean = -self.mean_pressure_gradient * r * r / (four * mu) * (T::one() - xi * xi);

        // Pulsatile component
        let u_pulsatile = self.profile().velocity(xi, t);

        u_mean + u_pulsatile
    }

    /// Calculate impedance magnitude |Z| for this segment
    ///
    /// Z = ΔP / Q (complex impedance)
    pub fn impedance_magnitude(&self) -> T {
        let alpha = self.womersley_number().value();
        let r = self.radius;
        let mu = self.viscosity;
        let rho = self.density;
        let omega = self.omega;
        let pi = T::pi();
        let eight = T::from_f64(8.0).expect("Mathematical constant conversion compromised");

        if alpha < T::one() {
            // Low α: Z ≈ 8μL/(πR⁴) (Poiseuille resistance dominates)
            eight * mu * self.length / (pi * r.powi(4))
        } else {
            // High α: Z ≈ ρωL/(πR²) (inertance dominates)
            rho * omega * self.length / (pi * r * r)
        }
    }
}
