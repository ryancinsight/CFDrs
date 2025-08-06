//! Plugin system for extensible CFD solvers.
//!
//! This module provides a flexible plugin architecture following SOLID principles,
//! enabling runtime extensibility while maintaining type safety and zero-cost abstractions.

use crate::{Error, Result};
use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

/// Trait for simulation plugins following the Interface Segregation Principle
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

/// Plugin lifecycle hooks for advanced control
pub trait PluginLifecycle: SimulationPlugin {
    /// Called before simulation starts
    fn on_start(&self, _state: &mut Self::State) -> Result<()> {
        Ok(())
    }

    /// Called after simulation ends
    fn on_stop(&self, _state: &Self::State) -> Result<()> {
        Ok(())
    }

    /// Called on simulation reset
    fn on_reset(&self, state: &mut Self::State, config: Self::Config) -> Result<()> {
        *state = self.initialize(config)?;
        Ok(())
    }
}

/// Plugin composition trait for combining multiple plugins
pub trait ComposablePlugin: SimulationPlugin {
    /// Compose with another plugin
    fn compose<P: SimulationPlugin>(self, other: P) -> ComposedPlugin<Self, P>
    where
        Self: Sized,
    {
        ComposedPlugin::new(self, other)
    }
}

/// Composed plugin combining two plugins
pub struct ComposedPlugin<P1, P2> {
    #[allow(dead_code)]
    plugin1: P1,
    #[allow(dead_code)]
    plugin2: P2,
}

impl<P1, P2> ComposedPlugin<P1, P2> {
    /// Create a new composed plugin
    pub fn new(plugin1: P1, plugin2: P2) -> Self {
        Self { plugin1, plugin2 }
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

    /// Get type ID for downcasting
    fn type_id(&self) -> TypeId;

    /// As any for downcasting
    fn as_any(&self) -> &dyn Any;

    /// As any mut for downcasting
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

/// Plugin metadata
#[derive(Debug, Clone, Default)]
pub struct PluginMetadata {
    /// Plugin name
    pub name: String,
    /// Plugin version following semver
    pub version: String,
    /// Plugin description
    pub description: String,
    /// Plugin author
    pub author: Option<String>,
    /// Plugin license
    pub license: Option<String>,
    /// Plugin capabilities
    pub capabilities: Vec<String>,
}

/// Plugin registry for managing available plugins
#[derive(Clone)]
pub struct PluginRegistry {
    plugins: Arc<RwLock<HashMap<String, Arc<dyn Plugin>>>>,
    factories: Arc<RwLock<HashMap<TypeId, Arc<dyn Any + Send + Sync>>>>,
}

impl PluginRegistry {
    /// Create a new plugin registry
    pub fn new() -> Self {
        Self {
            plugins: Arc::new(RwLock::new(HashMap::new())),
            factories: Arc::new(RwLock::new(HashMap::new())),
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

    /// Register a solver factory
    pub fn register_factory<F: SolverFactory + 'static>(&self, factory: F) -> Result<()> {
        let type_id = TypeId::of::<F::Solver>();
        let mut factories = self.factories.write().map_err(|_| {
            Error::PluginError("Failed to acquire write lock on factory registry".to_string())
        })?;

        factories.insert(type_id, Arc::new(factory));
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

    /// Get a solver factory by type
    pub fn get_factory<S: SimulationPlugin + 'static>(&self) -> Result<Arc<dyn Any + Send + Sync>> {
        let type_id = TypeId::of::<S>();
        let factories = self.factories.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on factory registry".to_string())
        })?;

        factories
            .get(&type_id)
            .cloned()
            .ok_or_else(|| Error::PluginError("Factory not found for solver type".to_string()))
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

    /// Iterate over all plugins
    pub fn iter(&self) -> Result<impl Iterator<Item = (String, Arc<dyn Plugin>)>> {
        let plugins = self.plugins.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on plugin registry".to_string())
        })?;

        Ok(plugins.clone().into_iter())
    }

    /// Filter plugins by capability
    pub fn filter_by_capability(&self, capability: &str) -> Result<Vec<String>> {
        let plugins = self.plugins.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on plugin registry".to_string())
        })?;

        Ok(plugins
            .iter()
            .filter_map(|(name, _)| {
                self.metadata(name).ok().and_then(|meta| {
                    if meta.capabilities.contains(&capability.to_string()) {
                        Some(name.clone())
                    } else {
                        None
                    }
                })
            })
            .collect())
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

    /// Validate parameters before solver creation
    fn validate_params(&self, _params: &SolverParams) -> Result<()> {
        // Default implementation does no validation
        Ok(())
    }
}

/// Builder pattern for solver parameters
#[derive(Debug, Clone, Default)]
pub struct SolverParams {
    /// Solver-specific parameters
    params: HashMap<String, serde_json::Value>,
}

impl SolverParams {
    /// Create new solver parameters
    pub fn new() -> Self {
        Self {
            params: HashMap::new(),
        }
    }

    /// Builder method to add a parameter
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

    /// Get optional parameter
    pub fn get_optional<T: serde::de::DeserializeOwned>(&self, key: &str) -> Result<Option<T>> {
        match self.params.get(key) {
            Some(value) => serde_json::from_value(value.clone())
                .map(Some)
                .map_err(|e| Error::SerializationError(e.to_string())),
            None => Ok(None),
        }
    }

    /// Check if parameter exists
    pub fn contains(&self, key: &str) -> bool {
        self.params.contains_key(key)
    }

    /// Iterate over all parameters
    pub fn iter(&self) -> impl Iterator<Item = (&String, &serde_json::Value)> {
        self.params.iter()
    }

    /// Merge with another set of parameters
    pub fn merge(mut self, other: Self) -> Self {
        self.params.extend(other.params);
        self
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
    /// Required parameters
    pub required_params: Vec<String>,
    /// Optional parameters with defaults
    pub optional_params: HashMap<String, serde_json::Value>,
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

        fn type_id(&self) -> TypeId {
            TypeId::of::<Self>()
        }

        fn as_any(&self) -> &dyn Any {
            self
        }

        fn as_any_mut(&mut self) -> &mut dyn Any {
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

        // Test optional parameter
        let optional: Option<String> = params.get_optional("missing").unwrap();
        assert!(optional.is_none());
    }

    #[test]
    fn test_params_builder_pattern() {
        let params = SolverParams::new()
            .with("a", 1)
            .unwrap()
            .with("b", 2.0)
            .unwrap()
            .with("c", "test")
            .unwrap();

        assert_eq!(params.get::<i32>("a").unwrap(), 1);
        assert_eq!(params.get::<f64>("b").unwrap(), 2.0);
        assert_eq!(params.get::<String>("c").unwrap(), "test");
    }
}