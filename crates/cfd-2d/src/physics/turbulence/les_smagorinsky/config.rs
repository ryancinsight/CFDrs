//! Smagorinsky LES model configuration
//!
//! Defines the configuration parameters and defaults for the Smagorinsky
//! Large Eddy Simulation model.

/// Smagorinsky LES model configuration
#[derive(Debug, Clone)]
pub struct SmagorinskyConfig {
    /// Smagorinsky constant (C_S)
    pub smagorinsky_constant: f64,
    /// Use dynamic procedure to compute C_S locally
    pub dynamic_procedure: bool,
    /// Wall-damping coefficient
    pub wall_damping: bool,
    /// Van Driest damping constant
    pub van_driest_constant: f64,
    /// Minimum SGS viscosity (numerical stability)
    pub min_sgs_viscosity: f64,
    /// Use GPU acceleration if available
    pub use_gpu: bool,
}

impl Default for SmagorinskyConfig {
    fn default() -> Self {
        Self {
            smagorinsky_constant: 0.1,      // Standard value for homogeneous turbulence
            dynamic_procedure: false,       // Use fixed constant by default
            wall_damping: true,             // Enable wall damping
            van_driest_constant: 0.4,       // Standard van Driest constant
            min_sgs_viscosity: 1e-8,        // Prevent division by zero
            use_gpu: cfg!(feature = "gpu"), // Use GPU if feature is enabled
        }
    }
}

impl SmagorinskyConfig {
    /// Create a new configuration with custom Smagorinsky constant
    pub const fn with_constant(smagorinsky_constant: f64) -> Self {
        Self {
            smagorinsky_constant,
            dynamic_procedure: false,
            wall_damping: true,
            van_driest_constant: 0.4,
            min_sgs_viscosity: 1e-8,
            use_gpu: cfg!(feature = "gpu"),
        }
    }

    /// Enable dynamic procedure
    pub const fn with_dynamic(mut self) -> Self {
        self.dynamic_procedure = true;
        self
    }

    /// Disable wall damping
    pub const fn without_wall_damping(mut self) -> Self {
        self.wall_damping = false;
        self
    }
}

