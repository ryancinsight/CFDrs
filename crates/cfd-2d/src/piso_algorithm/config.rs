//! PISO solver configuration

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use cfd_core::solver::SolverConfig;

/// PISO solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PisoConfig<T: RealField + Copy> {
    /// Base solver configuration
    pub base: SolverConfig<T>,
    /// Time step for transient simulation
    pub time_step: T,
    /// Number of pressure corrector steps (typically 2, max 3)
    pub num_correctors: usize,
    /// Number of non-orthogonal corrector iterations
    pub non_orthogonal_correctors: usize,
    /// Velocity under-relaxation factor
    pub velocity_relaxation: T,
    /// Pressure under-relaxation factor  
    pub pressure_relaxation: T,
    /// Maximum iterations for pressure equation
    pub pressure_max_iterations: usize,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for PisoConfig<T> {
    fn default() -> Self {
        Self {
            base: SolverConfig::default(),
            time_step: T::from_f64(0.01).unwrap_or_else(|| T::from_f64(0.01).unwrap()),
            num_correctors: 2,  // Standard PISO uses 2 correctors
            non_orthogonal_correctors: 1,
            velocity_relaxation: T::from_f64(0.7).unwrap_or_else(|| T::from_f64(0.7).unwrap()),
            pressure_relaxation: T::from_f64(0.3).unwrap_or_else(|| T::from_f64(0.3).unwrap()),
            pressure_max_iterations: 100,
        }
    }
}