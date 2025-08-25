//! PISO solver configuration

use cfd_core::solver::SolverConfig;
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
    pub pressure_relaxation: T,
}

impl<T: RealField + Copy + FromPrimitive> Default for PisoConfig<T> {
    fn default() -> Self {
        Self {
            n_correctors: 2,
            n_non_orthogonal_correctors: 1,
            time_step: T::from_f64(0.01).unwrap_or_else(|| {
                T::from_f64(1e-2).unwrap_or_else(|| {
                    T::from_usize(1).unwrap_or_else(T::one)
                        / T::from_usize(100).unwrap_or_else(T::one)
                })
            }),
            velocity_relaxation: T::from_f64(0.7).unwrap_or_else(|| {
                T::from_usize(7).unwrap_or_else(T::one) / T::from_usize(10).unwrap_or_else(T::one)
            }),
            pressure_relaxation: T::from_f64(0.3).unwrap_or_else(|| {
                T::from_usize(3).unwrap_or_else(T::one) / T::from_usize(10).unwrap_or_else(T::one)
            }),
        }
    }
}
