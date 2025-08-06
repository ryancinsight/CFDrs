//! Plugin system for extensible CFD solvers.

use crate::{Error, Result};
use std::any::Any;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

/// Trait for simulation plugins
pub trait SimulationPlugin: Send + Sync + 'static {
    /// Configuration type for this plugin
    type Config: Send + Sync + 'static;
    /// State type maintained by this plugin
    type State: Send + Sync + 'static;
    /// Output type produced by this plugin
    type Output: Send + Sync + 'static;

    /// Initialize the plugin with given configuration
    fn initialize(&self, config: Self::Config) -> Result<Self::State>;

    /// Perform one simulation step
    fn step(&self, state: &mut Self::State, dt: f64) -> Result<()>;

    /// Generate output from current state
    fn output(&self, state: &Self::State) -> Self::Output;

    /// Get plugin metadata
    fn metadata(&self) -> PluginMetadata {
        PluginMetadata::default()
    }
}

/// Generic plugin trait for type erasure
pub trait Plugin: Any + Send + Sync {
    /// Get plugin name
    fn name(&self) -> &str;

    /// Get plugin version
    fn version(&self) -> &str;

    /// Get plugin description
    fn description(&self) -> &str;

    /// As any for downcasting
    fn as_any(&self) -> &dyn Any;
}

/// Plugin metadata
#[derive(Debug, Clone, Default)]
pub struct PluginMetadata {
    /// Plugin name
    pub name: String,
    /// Plugin version
    pub version: String,
    /// Plugin description
    pub description: String,
    /// Plugin author
    pub author: Option<String>,
    /// Plugin license
    pub license: Option<String>,
}

/// Plugin registry for managing available plugins
#[derive(Clone)]
pub struct PluginRegistry {
    plugins: Arc<RwLock<HashMap<String, Arc<dyn Plugin>>>>,
}

impl PluginRegistry {
    /// Create a new plugin registry
    pub fn new() -> Self {
        Self {
            plugins: Arc::new(RwLock::new(HashMap::new())),
        }
    }

    /// Register a plugin
    pub fn register<P: Plugin + 'static>(&self, plugin: P) -> Result<()> {
        let name = plugin.name().to_string();
        let mut plugins = self.plugins.write().map_err(|_| {
            Error::PluginError("Failed to acquire write lock on plugin registry".to_string())
        })?;

        if plugins.contains_key(&name) {
            return Err(Error::PluginError(format!(
                "Plugin '{}' is already registered",
                name
            )));
        }

        plugins.insert(name, Arc::new(plugin));
        Ok(())
    }

    /// Get a plugin by name
    pub fn get(&self, name: &str) -> Result<Arc<dyn Plugin>> {
        let plugins = self.plugins.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on plugin registry".to_string())
        })?;

        plugins
            .get(name)
            .cloned()
            .ok_or_else(|| Error::PluginError(format!("Plugin '{}' not found", name)))
    }

    /// List all registered plugins
    pub fn list(&self) -> Result<Vec<String>> {
        let plugins = self.plugins.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on plugin registry".to_string())
        })?;

        Ok(plugins.keys().cloned().collect())
    }

    /// Get plugin metadata
    pub fn metadata(&self, name: &str) -> Result<PluginMetadata> {
        let plugin = self.get(name)?;
        Ok(PluginMetadata {
            name: plugin.name().to_string(),
            version: plugin.version().to_string(),
            description: plugin.description().to_string(),
            ..Default::default()
        })
    }
}

impl Default for PluginRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// Solver factory trait for creating solver instances
pub trait SolverFactory: Send + Sync + 'static {
    /// Solver type produced by this factory
    type Solver: SimulationPlugin;

    /// Create a new solver instance
    fn create_solver(&self, params: SolverParams) -> Result<Self::Solver>;

    /// Get factory metadata
    fn metadata(&self) -> FactoryMetadata {
        FactoryMetadata::default()
    }
}

/// Parameters for solver creation
#[derive(Debug, Clone)]
pub struct SolverParams {
    /// Solver-specific parameters
    pub params: HashMap<String, serde_json::Value>,
}

impl SolverParams {
    /// Create new solver parameters
    pub fn new() -> Self {
        Self {
            params: HashMap::new(),
        }
    }

    /// Add a parameter
    pub fn with<T: serde::Serialize>(mut self, key: impl Into<String>, value: T) -> Result<Self> {
        let value = serde_json::to_value(value)
            .map_err(|e| Error::SerializationError(e.to_string()))?;
        self.params.insert(key.into(), value);
        Ok(self)
    }

    /// Get a parameter
    pub fn get<T: serde::de::DeserializeOwned>(&self, key: &str) -> Result<T> {
        let value = self
            .params
            .get(key)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Parameter '{}' not found", key)))?;

        serde_json::from_value(value.clone())
            .map_err(|e| Error::SerializationError(e.to_string()))
    }
}

impl Default for SolverParams {
    fn default() -> Self {
        Self::new()
    }
}

/// Factory metadata
#[derive(Debug, Clone, Default)]
pub struct FactoryMetadata {
    /// Factory name
    pub name: String,
    /// Supported solver types
    pub solver_types: Vec<String>,
    /// Factory description
    pub description: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    struct TestPlugin {
        name: String,
    }

    impl Plugin for TestPlugin {
        fn name(&self) -> &str {
            &self.name
        }

        fn version(&self) -> &str {
            "1.0.0"
        }

        fn description(&self) -> &str {
            "Test plugin"
        }

        fn as_any(&self) -> &dyn Any {
            self
        }
    }

    #[test]
    fn test_plugin_registry() {
        let registry = PluginRegistry::new();
        let plugin = TestPlugin {
            name: "test".to_string(),
        };

        // Register plugin
        registry.register(plugin).unwrap();

        // Get plugin
        let retrieved = registry.get("test").unwrap();
        assert_eq!(retrieved.name(), "test");

        // List plugins
        let list = registry.list().unwrap();
        assert_eq!(list.len(), 1);
        assert!(list.contains(&"test".to_string()));
    }

    #[test]
    fn test_solver_params() {
        let params = SolverParams::new()
            .with("tolerance", 1e-6)
            .unwrap()
            .with("max_iterations", 1000)
            .unwrap();

        let tolerance: f64 = params.get("tolerance").unwrap();
        let max_iter: i32 = params.get("max_iterations").unwrap();

        assert_eq!(tolerance, 1e-6);
        assert_eq!(max_iter, 1000);
    }
}