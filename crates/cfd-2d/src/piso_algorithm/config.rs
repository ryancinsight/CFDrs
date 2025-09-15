//! PISO solver configuration

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Configuration for PISO algorithm
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PisoConfig<T: RealField + Copy> {
    /// Number of pressure corrector steps
    pub n_correctors: usize,

    /// Number of non-orthogonal corrector loops
    pub n_non_orthogonal_correctors: usize,

    /// Time step size
    pub time_step: T,

    /// Under-relaxation factors
    pub velocity_relaxation: T,
    /// Under-relaxation factor for pressure correction
    pub pressure_relaxation: T,

    /// Logging frequency (None disables logging, Some(n) logs every n steps)
    pub log_frequency: Option<usize>,
}

impl<T: RealField + Copy + FromPrimitive> Default for PisoConfig<T> {
    fn default() -> Self {
        Self {
            n_correctors: 2,
            n_non_orthogonal_correctors: 1,
            time_step: T::from_f64(0.01)
                .expect("Failed to represent default time step (0.01) in numeric type T"),
            velocity_relaxation: T::from_f64(0.7)
                .expect("Failed to represent velocity relaxation factor (0.7) in numeric type T"),
            pressure_relaxation: T::from_f64(0.3)
                .expect("Failed to represent pressure relaxation factor (0.3) in numeric type T"),
            log_frequency: Some(100), // Default to logging every 100 steps
        }
    }
}
