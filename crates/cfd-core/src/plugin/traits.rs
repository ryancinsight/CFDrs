//! Core plugin traits for CFD solvers.

use crate::error::Result;

/// Core plugin trait - the single abstraction for all plugin functionality
pub trait Plugin: Send + Sync {
    /// Plugin name (must be unique)
    fn name(&self) -> &str;

    /// Plugin version
    fn version(&self) -> &str;

    /// Plugin description
    fn description(&self) -> &str;

    /// Plugin author (optional)
    fn author(&self) -> Option<&str> {
        None
    }

    /// Plugin license (optional)
    fn license(&self) -> Option<&str> {
        None
    }

    /// Plugin dependencies (other plugin names)
    fn dependencies(&self) -> Vec<&str> {
        vec![]
    }

    /// Plugin capabilities
    fn capabilities(&self) -> Vec<&str> {
        vec![]
    }

    /// Check if plugin supports parallel execution
    fn supports_parallel(&self) -> bool {
        false
    }

    /// Check if plugin supports adaptive time stepping
    fn supports_adaptive_timestep(&self) -> bool {
        false
    }

    /// Check if plugin supports multi-phase flows
    fn supports_multiphase(&self) -> bool {
        false
    }

    /// Check if plugin supports turbulence modeling
    fn supports_turbulence(&self) -> bool {
        false
    }
}

/// Extended trait for simulation plugins
pub trait SimulationPlugin: Plugin {
    /// Initialize the plugin with given configuration
    ///
    /// # Errors
    /// Returns error if configuration is invalid or initialization fails
    fn initialize(&mut self, config: &str) -> Result<()>;

    /// Execute the plugin's main functionality
    ///
    /// # Errors
    /// Returns error if execution fails or encounters numerical issues
    fn execute(&mut self, timestep: f64) -> Result<()>;

    /// Clean up resources
    ///
    /// # Errors
    /// Returns error if cleanup fails or resources cannot be released
    fn cleanup(&mut self) -> Result<()>;

    /// Get current state for checkpointing
    ///
    /// # Errors
    /// Returns error if state cannot be serialized or accessed
    fn get_state(&self) -> Result<String>;

    /// Restore from checkpoint
    ///
    /// # Errors
    /// Returns error if state cannot be deserialized or restoration fails
    fn set_state(&mut self, state: &str) -> Result<()>;

    /// Validate plugin configuration
    ///
    /// # Errors
    /// Returns error if configuration is invalid or incompatible
    fn validate_config(&self, config: &str) -> Result<()> {
        // Default implementation: accept any config
        let _ = config;
        Ok(())
    }

    /// Get performance metrics
    fn get_metrics(&self) -> Result<String> {
        Ok(String::new())
    }
}
