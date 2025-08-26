//! Configuration for Level Set method

use serde::{Deserialize, Serialize};

// Named constants for Level Set - NO MAGIC NUMBERS!
/// Default reinitialization interval in timesteps
pub const DEFAULT_REINITIALIZATION_INTERVAL: usize = 5;
/// Width of narrow band in grid cells
pub const DEFAULT_BAND_WIDTH: f64 = 5.0;
/// CFL number for time stepping stability
pub const DEFAULT_CFL_NUMBER: f64 = 0.5;
/// Convergence tolerance for iterative methods
pub const DEFAULT_TOLERANCE: f64 = 1e-6;
/// Maximum iterations for reinitialization
pub const DEFAULT_MAX_ITERATIONS: usize = 100;
/// Interface thickness for smoothing (in grid cells)
pub const EPSILON_SMOOTHING: f64 = 1.5;
/// Exponent for power calculations
pub const POWER_EXPONENT_TWO: f64 = 2.0;
/// Half value constant
pub const HALF_VALUE: f64 = 0.5;
/// Small epsilon for avoiding division by zero
pub const EPSILON_DIVISION: f64 = 1e-10;

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