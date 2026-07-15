//! PISO solver configuration
//!
//! # Invariant
//!
//! Default parameter ranges are chosen to satisfy the PISO stability constraints:
//! $n_{\text{correctors}} \ge 2$ (Issa 1986), CFL $\le 1$ for explicit time
//! advancement, and convergence tolerance $> 0$.

use crate::scalar;
use crate::scalar::Cfd2dScalar;
use eunomia::FloatElement;
use serde::{Deserialize, Serialize};

/// Configuration for PISO algorithm
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PisoConfig<T: Cfd2dScalar + Copy> {
    /// Number of pressure corrector steps
    pub n_correctors: usize,

    /// Number of non-orthogonal corrector loops
    pub n_non_orthogonal_correctors: usize,

    /// Time step size
    pub time_step: T,

    /// Under-relaxation factors
    pub velocity_relaxation: T,
    /// Under-relaxation factor for pressure equation (typically 0.3-0.8)
    pub pressure_relaxation: T,

    /// Logging frequency (None disables logging, Some(n) logs every n steps)
    pub log_frequency: Option<usize>,
}

impl<T: Cfd2dScalar + Copy + FloatElement> Default for PisoConfig<T> {
    fn default() -> Self {
        Self {
            n_correctors: 2,
            n_non_orthogonal_correctors: 1,
            time_step: scalar::from_f64(0.01),
            velocity_relaxation: scalar::from_f64(0.7),
            pressure_relaxation: scalar::from_f64(0.3),
            log_frequency: Some(100), // Default to logging every 100 steps
        }
    }
}
