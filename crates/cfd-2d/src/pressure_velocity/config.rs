//! Configuration for STANDARD algorithm

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Linear solver choice for pressure Poisson equation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PressureLinearSolver {
    /// Conjugate Gradient (symmetric matrices only)
    ConjugateGradient,
    /// BiCGSTAB (general non-symmetric matrices)
    BiCGSTAB,
    /// GMRES(m) with restart (industry standard for SIMPLE/PISO)
    /// Reference: Saad & Schultz (1986), Saad (2003) §6.5
    GMRES { 
        /// Restart dimension (typically 20-50 for CFD, default 30)
        restart_dim: usize 
    },
}

impl Default for PressureLinearSolver {
    fn default() -> Self {
        // GMRES is the industry standard for pressure correction equations
        // restart_dim=30 is standard for CFD applications
        Self::GMRES { restart_dim: 30 }
    }
}

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
    /// Linear solver for pressure Poisson equation
    #[serde(default)]
    pub pressure_linear_solver: PressureLinearSolver,
}

impl<T: RealField + Copy + FromPrimitive + Copy> PressureVelocityConfig<T> {
    /// Create new configuration with validation
    pub fn new() -> cfd_core::error::Result<Self> {
        Ok(Self {
            base: cfd_core::solver::SolverConfig::builder()
                .max_iterations(crate::constants::solver::LOG_INTERVAL)
                .tolerance(
                    T::from_f64(cfd_core::constants::numerical::solver::CONVERGENCE_TOLERANCE)
                        .ok_or_else(|| {
                            cfd_core::error::Error::InvalidConfiguration(
                                "Cannot convert tolerance".into(),
                            )
                        })?,
                )
                .build(),
            dt: T::from_f64(
                cfd_core::constants::numerical::time::DEFAULT_CFL
                    * cfd_core::constants::numerical::time::TIME_STEP_SAFETY,
            )
            .ok_or_else(|| {
                cfd_core::error::Error::InvalidConfiguration("Cannot convert time step".into())
            })?,
            alpha_u: T::from_f64(cfd_core::constants::numerical::relaxation::VELOCITY_RELAXATION)
                .ok_or_else(|| {
                cfd_core::error::Error::InvalidConfiguration(
                    "Cannot convert velocity relaxation".into(),
                )
            })?,
            alpha_p: T::from_f64(cfd_core::constants::numerical::relaxation::PRESSURE_RELAXATION)
                .ok_or_else(|| {
                cfd_core::error::Error::InvalidConfiguration(
                    "Cannot convert pressure relaxation".into(),
                )
            })?,
            use_rhie_chow: true,
            convection_scheme: crate::schemes::SpatialScheme::SecondOrderUpwind,
            implicit_momentum: true,
            pressure_linear_solver: PressureLinearSolver::default(),
        })
    }

    /// Validate configuration parameters
    pub fn validate(&self) -> cfd_core::error::Result<()> {
        if self.alpha_u <= T::zero() || self.alpha_u >= T::one() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Velocity relaxation factor must be in (0, 1)".into(),
            ));
        }
        if self.alpha_p <= T::zero() || self.alpha_p >= T::one() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Pressure relaxation factor must be in (0, 1)".into(),
            ));
        }
        if self.dt <= T::zero() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Time step must be positive".into(),
            ));
        }
        Ok(())
    }
}
