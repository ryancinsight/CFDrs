//! Tests for the plugin system.

use super::*;
use std::sync::Arc;

/// Test plugin implementation
struct TestPlugin {
    name: String,
    dependencies: Vec<String>,
}

impl Plugin for TestPlugin {
    fn name(&self) -> &str {
        &self.name
    }

    fn version(&self) -> &'static str {
        "1.0.0"
    }

    fn description(&self) -> &'static str {
        "Test plugin"
    }

    fn dependencies(&self) -> Vec<&str> {
        self.dependencies.iter().map(std::string::String::as_str).collect()
    }
}

#[test]
fn test_plugin_registration() {
    let mut registry = PluginRegistry::new();

    let plugin = Arc::new(TestPlugin {
        name: "test".to_string(),
        dependencies: vec![],
    });

    assert!(registry.register(plugin).is_ok());
    assert!(registry.get("test").is_some());
    assert_eq!(registry.list(), vec!["test"]);
}

#[test]
fn test_dependency_resolution() {
    let mut registry = PluginRegistry::new();

    // Register plugins with dependencies
    let plugin_a = Arc::new(TestPlugin {
        name: "A".to_string(),
        dependencies: vec![],
    });

    let plugin_b = Arc::new(TestPlugin {
        name: "B".to_string(),
        dependencies: vec!["A".to_string()],
    });

    let plugin_c = Arc::new(TestPlugin {
        name: "C".to_string(),
        dependencies: vec!["A".to_string(), "B".to_string()],
    });

    // Register in any order - should work
    assert!(registry.register(plugin_a).is_ok());
    assert!(registry.register(plugin_b).is_ok());
    assert!(registry.register(plugin_c).is_ok());

    // Check load order
    let load_order = registry.get_load_order();
    let a_idx = load_order
        .iter()
        .position(|x| x == "A")
        .expect("Plugin A not found in load order");
    let b_idx = load_order
        .iter()
        .position(|x| x == "B")
        .expect("Plugin B not found in load order");
    let c_idx = load_order
        .iter()
        .position(|x| x == "C")
        .expect("Plugin C not found in load order");

    assert!(a_idx < b_idx);
    assert!(b_idx < c_idx);
}

#[test]
fn test_circular_dependency_detection() {
    use dependency::DependencyResolver;

    let mut resolver = DependencyResolver::new();

    // Create circular dependency: A -> B -> C -> A
    assert!(resolver
        .add_plugin("A".to_string(), vec!["B".to_string()])
        .is_ok());
    assert!(resolver
        .add_plugin("B".to_string(), vec!["C".to_string()])
        .is_ok());

    // This should fail due to circular dependency
    assert!(resolver
        .add_plugin("C".to_string(), vec!["A".to_string()])
        .is_err());
}

#[test]
fn test_missing_dependency() {
    let mut registry = PluginRegistry::new();

    let plugin = Arc::new(TestPlugin {
        name: "test".to_string(),
        dependencies: vec!["nonexistent".to_string()],
    });

    // Should fail due to missing dependency
    assert!(registry.register(plugin).is_err());
}

#[test]
fn test_health_monitoring() {
    let mut registry = PluginRegistry::new();

    let plugin = Arc::new(TestPlugin {
        name: "test".to_string(),
        dependencies: vec![],
    });

    registry
        .register(plugin)
        .expect("Failed to register plugin");

    // Check initial health (should be Unknown)
    assert!(matches!(
        registry.get_health("test"),
        Some(PluginHealthStatus::Unknown)
    ));

    // Update health
    registry.update_health(
        "test",
        PluginHealthStatus::Degraded("Test issue".to_string()),
    );

    assert!(matches!(
        registry.get_health("test"),
        Some(PluginHealthStatus::Degraded(_))
    ));

    // Check system summary
    let summary = registry.get_system_health();
    assert_eq!(summary.degraded_plugins, 1);
    assert!(matches!(summary.status, SystemStatus::Degraded));
}

#[test]
fn test_plugin_removal() {
    let mut registry = PluginRegistry::new();

    let plugin_a = Arc::new(TestPlugin {
        name: "A".to_string(),
        dependencies: vec![],
    });

    let plugin_b = Arc::new(TestPlugin {
        name: "B".to_string(),
        dependencies: vec!["A".to_string()],
    });

    registry
        .register(plugin_a)
        .expect("Failed to register plugin A");
    registry
        .register(plugin_b)
        .expect("Failed to register plugin B");

    // Cannot remove A because B depends on it
    assert!(registry.unregister("A").is_err());

    // Can remove B
    assert!(registry.unregister("B").is_ok());

    // Now can remove A
    assert!(registry.unregister("A").is_ok());

    assert_eq!(registry.list().len(), 0);
}
