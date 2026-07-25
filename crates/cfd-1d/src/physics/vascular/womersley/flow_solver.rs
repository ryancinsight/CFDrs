//! Complete Womersley flow solver for arterial segments.
//!
//! Provides time-varying flow solutions for vessel segments with
//! given inlet conditions and geometry.

use super::profile::WomersleyProfile;
use super::WomersleyNumber;
use crate::scalar::Cfd1dScalar;
use aequitas::systems::si::quantities::{
    DynamicViscosity, HydraulicResistance, Length, MassDensity, Pressure, PressureGradient,
    ReciprocalTime, Time, Velocity,
};
use eunomia::FloatElement;
use serde::{Deserialize, Serialize};

/// Complete Womersley flow solver for arterial segments
///
/// Provides time-varying flow solutions for vessel segments with
/// given inlet conditions and geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WomersleyFlow<T: Cfd1dScalar + Copy> {
    /// Vessel radius \[m]
    pub radius: Length<T>,
    /// Vessel length \[m]
    pub length: Length<T>,
    /// Fluid density [kg/m³]
    pub density: MassDensity<T>,
    /// Dynamic viscosity [Pa·s]
    pub viscosity: DynamicViscosity<T>,
    /// Angular frequency \[rad/s]
    pub omega: ReciprocalTime<T>,
    /// Inlet pressure amplitude \[Pa]
    pub inlet_pressure_amplitude: Pressure<T>,
    /// Mean pressure gradient [Pa/m]
    pub mean_pressure_gradient: PressureGradient<T>,
}

impl<T: Cfd1dScalar + FloatElement + Copy> WomersleyFlow<T> {
    /// Create new Womersley flow solver
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        radius: Length<T>,
        length: Length<T>,
        density: MassDensity<T>,
        viscosity: DynamicViscosity<T>,
        omega: ReciprocalTime<T>,
        inlet_pressure_amplitude: Pressure<T>,
        mean_pressure_gradient: PressureGradient<T>,
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
        let pressure_gradient_amplitude = PressureGradient::from_base(
            self.inlet_pressure_amplitude.into_base() / self.length.into_base(),
        );
        WomersleyProfile::new(self.womersley_number(), pressure_gradient_amplitude)
    }

    /// Calculate total velocity (mean + pulsatile) at position and time
    pub fn velocity(&self, xi: T, t: Time<T>) -> Velocity<T> {
        // Mean (steady) component - Poiseuille
        let r = self.radius.into_base();
        let mu = self.viscosity.into_base();
        let four = T::one() + T::one() + T::one() + T::one();
        let u_mean =
            -self.mean_pressure_gradient.into_base() * r * r / (four * mu) * (T::one() - xi * xi);

        // Pulsatile component
        let u_pulsatile = self.profile().velocity(xi, t).into_base();

        Velocity::from_base(u_mean + u_pulsatile)
    }

    /// Calculate impedance magnitude |Z| for this segment
    ///
    /// Z = ΔP / Q (complex impedance)
    pub fn impedance_magnitude(&self) -> HydraulicResistance<T> {
        let alpha = self.womersley_number().value().into_base();
        let r = self.radius.into_base();
        let mu = self.viscosity.into_base();
        let rho = self.density.into_base();
        let omega = self.omega.into_base();
        let pi = T::pi();
        let eight = <T as FloatElement>::from_f64(8.0);

        let resistance = if alpha < T::one() {
            // Low α: Z ≈ 8μL/(πR⁴) (Poiseuille resistance dominates)
            eight * mu * self.length.into_base() / (pi * <T as FloatElement>::powi(r, 4))
        } else {
            // High α: Z ≈ ρωL/(πR²) (inertance dominates)
            rho * omega * self.length.into_base() / (pi * r * r)
        };
        HydraulicResistance::from_base(resistance)
    }
}
