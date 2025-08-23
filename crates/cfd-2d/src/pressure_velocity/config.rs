//! Configuration for STANDARD algorithm

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use cfd_core::constants;

/// Pressure-velocity coupling configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PressureVelocityConfig<T: RealField + Copy> {
    /// Base solver configuration
    pub base: cfd_core::solver::SolverConfig<T>,
    /// Time step (for unsteady problems)
    pub dt: T,
    /// Velocity under-relaxation factor (0 < α < 1)
    pub alpha_u: T,
    /// Pressure under-relaxation factor (0 < α < 1)
    pub alpha_p: T,
    /// Use Rhie-Chow interpolation for colocated grid
    pub use_rhie_chow: bool,
    /// Convection scheme
    pub convection_scheme: crate::schemes::SpatialScheme,
    /// Use implicit momentum solver
    pub implicit_momentum: bool,
}

impl<T: RealField + Copy + FromPrimitive + Copy> PressureVelocityConfig<T> {
    /// Create new configuration with validation
    pub fn new() -> cfd_core::error::Result<Self> {
        Ok(Self {
            base: cfd_core::solver::SolverConfig::builder()
                .max_iterations(100)
                .tolerance(T::from_f64(cfd_core::constants::numerical::solver::CONVERGENCE_TOLERANCE)
                    .ok_or_else(|| cfd_core::error::Error::InvalidConfiguration(
                        "Cannot convert tolerance".into()
                    ))?)
                .build(),
            dt: T::from_f64(cfd_core::constants::numerical::discretization::CFL_EXPLICIT * cfd_core::constants::numerical::time::SAFETY_FACTOR)
                .ok_or_else(|| cfd_core::error::Error::InvalidConfiguration(
                    "Cannot convert time step".into()
                ))?,
            alpha_u: T::from_f64(cfd_core::constants::numerical::relaxation::VELOCITY)
                .ok_or_else(|| cfd_core::error::Error::InvalidConfiguration(
                    "Cannot convert velocity relaxation".into()
                ))?,
            alpha_p: T::from_f64(cfd_core::constants::numerical::relaxation::PRESSURE)
                .ok_or_else(|| cfd_core::error::Error::InvalidConfiguration(
                    "Cannot convert pressure relaxation".into()
                ))?,
            use_rhie_chow: true,
            convection_scheme: crate::schemes::SpatialScheme::SecondOrderUpwind,
            implicit_momentum: true,
        })
    }
    
    /// Validate configuration parameters
    pub fn validate(&self) -> cfd_core::error::Result<()> {
        if self.alpha_u <= T::zero() || self.alpha_u >= T::one() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Velocity relaxation factor must be in (0, 1)".into()
            ));
        }
        if self.alpha_p <= T::zero() || self.alpha_p >= T::one() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Pressure relaxation factor must be in (0, 1)".into()
            ));
        }
        if self.dt <= T::zero() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Time step must be positive".into()
            ));
        }
        Ok(())
    }
}