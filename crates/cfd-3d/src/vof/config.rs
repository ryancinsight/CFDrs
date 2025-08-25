//! VOF configuration and constants

use serde::{Deserialize, Serialize};

// Named constants for VOF - SSOT principle
const DEFAULT_MAX_ITERATIONS: usize = 100;
const DEFAULT_TOLERANCE: f64 = 1e-6;
const DEFAULT_CFL_NUMBER: f64 = 0.3;
pub const VOF_EPSILON: f64 = 1e-10;  // Small value to avoid division by zero
pub const INTERFACE_THICKNESS: f64 = 1.5;  // Interface thickness in cells
pub const VOF_INTERFACE_LOWER: f64 = 0.01;  // Lower bound for interface cells
pub const VOF_INTERFACE_UPPER: f64 = 0.99;  // Upper bound for interface cells

/// VOF solver configuration constants
pub mod constants {
    /// Maximum iterations for PLIC reconstruction
    pub const PLIC_MAX_ITERATIONS: usize = 10;
    /// Tolerance for PLIC convergence
    pub const PLIC_TOLERANCE: f64 = 1e-6;
}

/// VOF configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VofConfig {
    /// Maximum iterations for interface reconstruction
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: f64,
    /// CFL number for time stepping
    pub cfl_number: f64,
    /// Use PLIC (Piecewise Linear Interface Calculation)
    pub use_plic: bool,
    /// Use geometric advection
    pub use_geometric_advection: bool,
    /// Enable compression to sharpen interface
    pub enable_compression: bool,
}

impl Default for VofConfig {
    fn default() -> Self {
        Self {
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: DEFAULT_TOLERANCE,
            cfl_number: DEFAULT_CFL_NUMBER,
            use_plic: true,
            use_geometric_advection: true,
            enable_compression: false,
        }
    }
}