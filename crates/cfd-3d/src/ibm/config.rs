//! IBM configuration and constants

use serde::{Deserialize, Serialize};

// Named constants for IBM
const DEFAULT_SMOOTHING_WIDTH: f64 = 1.5; // Delta function width in grid cells
const DEFAULT_FORCE_SCALE: f64 = 1.0;
const DEFAULT_MAX_ITERATIONS: usize = 100;
const DEFAULT_TOLERANCE: f64 = 1e-6;

/// IBM configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IbmConfig {
    /// Width of smoothing kernel (in grid cells)
    pub smoothing_width: f64,
    /// Force scaling factor
    pub force_scale: f64,
    /// Maximum iterations for force calculation
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: f64,
    /// Use direct forcing method
    pub use_direct_forcing: bool,
}

impl Default for IbmConfig {
    fn default() -> Self {
        Self {
            smoothing_width: DEFAULT_SMOOTHING_WIDTH,
            force_scale: DEFAULT_FORCE_SCALE,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: DEFAULT_TOLERANCE,
            use_direct_forcing: true,
        }
    }
}
