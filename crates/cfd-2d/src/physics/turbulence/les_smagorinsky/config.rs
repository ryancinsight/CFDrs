//! Smagorinsky LES model configuration
//!
//! Defines the configuration parameters and defaults for the Smagorinsky
//! Large Eddy Simulation model.
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \overline{u_i^\prime u_j^\prime}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\overline{u_i^\prime u_i^\prime} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

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
    /// Optional lower bound on SGS viscosity.
    ///
    /// The default is zero so laminar or strain-free resolved states remain
    /// exactly laminar. Nonzero values are explicit user-selected numerical
    /// regularization and are not part of the Smagorinsky closure.
    pub min_sgs_viscosity: f64,
    /// Use GPU acceleration if available
    pub use_gpu: bool,
}

impl Default for SmagorinskyConfig {
    fn default() -> Self {
        Self {
            smagorinsky_constant: 0.1, // Standard value for homogeneous turbulence
            dynamic_procedure: false,  // Use fixed constant by default
            wall_damping: true,        // Enable wall damping
            van_driest_constant: 0.4,  // Standard van Driest constant
            min_sgs_viscosity: 0.0,
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
            min_sgs_viscosity: 0.0,
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
