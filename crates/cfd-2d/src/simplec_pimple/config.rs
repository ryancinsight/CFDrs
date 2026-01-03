//! Configuration for SIMPLEC and PIMPLE algorithms

use nalgebra::RealField;
use num_traits::FromPrimitive;
use crate::schemes::SpatialScheme;
use crate::pressure_velocity::PressureLinearSolver;

/// Algorithm selection
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlgorithmType {
    /// SIMPLEC: Consistent pressure correction
    Simplec,
    /// PIMPLE: Merged PISO-SIMPLE algorithm
    Pimple,
}

/// Configuration for SIMPLEC/PIMPLE algorithms
#[derive(Debug, Clone)]
pub struct SimplecPimpleConfig<T: RealField + Copy> {
    /// Algorithm type
    pub algorithm: AlgorithmType,
    /// Time step size
    pub dt: T,
    /// Velocity under-relaxation factor (0 < α_u ≤ 1)
    pub alpha_u: T,
    /// Pressure under-relaxation factor (0 < α_p ≤ 1)
    pub alpha_p: T,
    /// Number of outer PIMPLE iterations (only used for PIMPLE)
    pub n_outer_correctors: usize,
    /// Number of inner SIMPLE corrections per outer iteration (only used for PIMPLE)
    pub n_inner_correctors: usize,
    /// Convergence tolerance for inner iterations
    pub tolerance: T,
    /// Maximum iterations for inner loops
    pub max_inner_iterations: usize,
    /// Use Rhie-Chow interpolation for momentum interpolation
    pub use_rhie_chow: bool,
    /// Convection scheme for momentum equations
    pub convection_scheme: SpatialScheme,
    /// Linear solver for pressure Poisson equation
    pub pressure_linear_solver: PressureLinearSolver,
}

impl<T: RealField + Copy + FromPrimitive> Default for SimplecPimpleConfig<T> {
    fn default() -> Self {
        Self {
            algorithm: AlgorithmType::Simplec,
            dt: T::from_f64(0.01).unwrap_or_else(T::one),
            alpha_u: T::from_f64(0.7).unwrap_or_else(|| T::from_f64(0.7).unwrap()),
            alpha_p: T::from_f64(0.3).unwrap_or_else(|| T::from_f64(0.3).unwrap()),
            n_outer_correctors: 2,
            n_inner_correctors: 1,
            tolerance: T::from_f64(1e-6).unwrap_or_else(|| T::from_f64(1e-6).unwrap()),
            max_inner_iterations: 50,
            use_rhie_chow: true,
            convection_scheme: SpatialScheme::SecondOrderUpwind,
            pressure_linear_solver: PressureLinearSolver::default(),
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> SimplecPimpleConfig<T> {
    /// Create configuration for SIMPLEC algorithm
    pub fn simplec() -> Self {
        Self {
            algorithm: AlgorithmType::Simplec,
            ..Default::default()
        }
    }

    /// Create configuration for PIMPLE algorithm
    pub fn pimple() -> Self {
        Self {
            algorithm: AlgorithmType::Pimple,
            ..Default::default()
        }
    }

    /// Validate configuration parameters
    pub fn validate(&self) -> cfd_core::error::Result<()> {
        if self.alpha_u <= T::zero() || self.alpha_u > T::one() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "alpha_u must be in (0, 1]".to_string(),
            ));
        }
        if self.alpha_p <= T::zero() || self.alpha_p > T::one() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "alpha_p must be in (0, 1]".to_string(),
            ));
        }
        if self.n_outer_correctors == 0 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "n_outer_correctors must be > 0".to_string(),
            ));
        }
        if self.n_inner_correctors == 0 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "n_inner_correctors must be > 0".to_string(),
            ));
        }
        Ok(())
    }
}
