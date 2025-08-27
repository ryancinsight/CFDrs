//! Level Set configuration

use serde::{Deserialize, Serialize};

// Named constants for Level Set
const DEFAULT_REINITIALIZATION_INTERVAL: usize = 5;
const DEFAULT_BAND_WIDTH: f64 = 5.0; // Width of narrow band in grid cells
const DEFAULT_CFL_NUMBER: f64 = 0.5;
const DEFAULT_TOLERANCE: f64 = 1e-6;
const DEFAULT_MAX_ITERATIONS: usize = 100;

/// Level Set configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LevelSetConfig {
    /// Reinitialization interval (timesteps)
    pub reinitialization_interval: usize,
    /// Narrow band width
    pub band_width: f64,
    /// CFL number for time stepping
    pub cfl_number: f64,
    /// Convergence tolerance
    pub tolerance: f64,
    /// Maximum iterations for reinitialization
    pub max_iterations: usize,
    /// Use narrow band method
    pub use_narrow_band: bool,
    /// Use WENO scheme for advection
    pub use_weno: bool,
}

impl Default for LevelSetConfig {
    fn default() -> Self {
        Self {
            reinitialization_interval: DEFAULT_REINITIALIZATION_INTERVAL,
            band_width: DEFAULT_BAND_WIDTH,
            cfl_number: DEFAULT_CFL_NUMBER,
            tolerance: DEFAULT_TOLERANCE,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            use_narrow_band: false,
            use_weno: true,
        }
    }
}
