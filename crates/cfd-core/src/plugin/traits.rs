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
    fn initialize(&mut self, config: &str) -> Result<()>;

    /// Execute the plugin's main functionality
    fn execute(&mut self, timestep: f64) -> Result<()>;

    /// Clean up resources
    fn cleanup(&mut self) -> Result<()>;

    /// Get current state for checkpointing
    fn get_state(&self) -> Result<String>;

    /// Restore from checkpoint
    fn set_state(&mut self, state: &str) -> Result<()>;

    /// Validate plugin configuration
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
