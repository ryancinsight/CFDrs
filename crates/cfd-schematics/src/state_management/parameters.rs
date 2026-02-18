//! Configurable parameter types and metadata
//!
//! This module provides the core configurable parameter type that supports
//! validation, adaptation, dependencies, and metadata tracking.

use crate::state_management::{
    adaptive::{AdaptationError, AdaptiveParameter, ChannelGenerationContext},
    constraints::ParameterConstraints,
    errors::{ParameterError, ParameterResult},
};
use std::collections::HashMap;
use std::fmt::Debug;

/// Metadata associated with a parameter
#[derive(Debug, Clone)]
pub struct ParameterMetadata {
    /// Human-readable name
    pub name: String,

    /// Description of the parameter's purpose
    pub description: String,

    /// Units of measurement (if applicable)
    pub units: Option<String>,

    /// Category or group this parameter belongs to
    pub category: String,

    /// Whether this parameter can be modified at runtime
    pub is_mutable: bool,

    /// Whether this parameter affects other parameters
    pub affects_others: bool,

    /// Version when this parameter was introduced
    pub version: String,

    /// Additional custom metadata
    pub custom: HashMap<String, String>,
}

impl ParameterMetadata {
    /// Create new parameter metadata
    #[must_use]
    pub fn new(name: &str, description: &str, category: &str) -> Self {
        Self {
            name: name.to_string(),
            description: description.to_string(),
            units: None,
            category: category.to_string(),
            is_mutable: true,
            affects_others: false,
            version: "1.0.0".to_string(),
            custom: HashMap::new(),
        }
    }

    /// Set units for this parameter
    #[must_use]
    pub fn with_units(mut self, units: &str) -> Self {
        self.units = Some(units.to_string());
        self
    }

    /// Mark parameter as immutable
    #[must_use]
    pub const fn immutable(mut self) -> Self {
        self.is_mutable = false;
        self
    }

    /// Mark parameter as affecting others
    #[must_use]
    pub const fn affects_others(mut self) -> Self {
        self.affects_others = true;
        self
    }

    /// Set version information
    #[must_use]
    pub fn with_version(mut self, version: &str) -> Self {
        self.version = version.to_string();
        self
    }

    /// Add custom metadata
    #[must_use]
    pub fn with_custom(mut self, key: &str, value: &str) -> Self {
        self.custom.insert(key.to_string(), value.to_string());
        self
    }
}

/// Generic configurable parameter with validation and adaptation
#[derive(Debug)]
pub struct ConfigurableParameter<T>
where
    T: Clone + Debug + PartialEq + PartialOrd + 'static,
{
    /// Current parameter value
    value: T,

    /// Default value for reset operations
    default: T,

    /// Validation constraints
    constraints: ParameterConstraints<T>,

    /// Adaptive behavior (if any)
    adaptive_behavior: Option<Box<dyn AdaptiveParameter<T, ChannelGenerationContext>>>,

    /// Parameter dependencies
    dependencies: Vec<String>,

    /// Parameter metadata
    metadata: ParameterMetadata,

    /// Whether the parameter has been modified from default
    is_modified: bool,

    /// History of value changes (for debugging/auditing)
    change_history: Vec<(T, String)>, // (value, reason)
}

impl<T> ConfigurableParameter<T>
where
    T: Clone + Debug + PartialEq + PartialOrd + 'static,
{
    /// Create a new configurable parameter
    pub fn new(
        value: T,
        constraints: ParameterConstraints<T>,
        metadata: ParameterMetadata,
    ) -> Self {
        let default = value.clone();
        Self {
            value: value.clone(),
            default,
            constraints,
            adaptive_behavior: None,
            dependencies: Vec::new(),
            metadata,
            is_modified: false,
            change_history: vec![(value, "initial".to_string())],
        }
    }

    /// Add adaptive behavior to this parameter
    pub fn with_adaptive_behavior<A>(mut self, behavior: A) -> Self
    where
        A: AdaptiveParameter<T, ChannelGenerationContext> + 'static,
    {
        self.adaptive_behavior = Some(Box::new(behavior));
        self
    }

    /// Add dependencies to this parameter
    pub fn with_dependencies(mut self, deps: Vec<String>) -> Self {
        self.dependencies = deps;
        self
    }

    /// Get current value with optional context-based adaptation
    ///
    /// Uses fallback to base value if adaptation fails
    pub fn get_value(&self, context: Option<&ChannelGenerationContext>) -> T {
        match (&self.adaptive_behavior, context) {
            (Some(behavior), Some(ctx)) if behavior.is_adaptive() => {
                // Use the new error-handling adapt method with fallback
                if let Ok(adapted) = behavior.adapt(self.value.clone(), ctx) {
                    adapted
                } else {
                    #[cfg(debug_assertions)]
                    eprintln!(
                        "Warning: Adaptation failed for parameter '{}', using base value",
                        self.metadata.name
                    );
                    self.value.clone()
                }
            }
            _ => self.value.clone(),
        }
    }

    /// Get current value with context-based adaptation, returning Result
    ///
    /// This method provides full error information if adaptation fails
    pub fn try_get_value(
        &self,
        context: Option<&ChannelGenerationContext>,
    ) -> Result<T, AdaptationError> {
        match (&self.adaptive_behavior, context) {
            (Some(behavior), Some(ctx)) if behavior.is_adaptive() => {
                behavior.adapt(self.value.clone(), ctx)
            }
            _ => Ok(self.value.clone()),
        }
    }

    /// Get the raw value without adaptation
    pub const fn get_raw_value(&self) -> &T {
        &self.value
    }

    /// Set parameter value with validation
    pub fn set_value(&mut self, value: T, reason: &str) -> ParameterResult<()> {
        // Check if parameter is mutable
        if !self.metadata.is_mutable {
            return Err(ParameterError::read_only(&self.metadata.name));
        }

        // Validate the new value
        self.constraints.validate(&value)?;

        // Update value and tracking
        self.value = value.clone();
        self.is_modified = self.value != self.default;
        self.change_history.push((value, reason.to_string()));

        // Limit history size to prevent memory growth
        if self.change_history.len() > 100 {
            self.change_history.remove(0);
        }

        Ok(())
    }

    /// Validate current value
    pub fn validate(&self) -> ParameterResult<()> {
        self.constraints.validate(&self.value)
    }

    /// Reset to default value
    pub fn reset(&mut self, reason: &str) -> ParameterResult<()> {
        if !self.metadata.is_mutable {
            return Err(ParameterError::read_only(&self.metadata.name));
        }

        let default_value = self.default.clone();
        self.set_value(default_value, reason)
    }

    /// Check if parameter has been modified from default
    pub const fn is_modified(&self) -> bool {
        self.is_modified
    }

    /// Get parameter metadata
    pub const fn metadata(&self) -> &ParameterMetadata {
        &self.metadata
    }

    /// Get parameter dependencies
    pub fn dependencies(&self) -> &[String] {
        &self.dependencies
    }

    /// Get constraint description
    pub fn constraint_description(&self) -> String {
        self.constraints.description()
    }

    /// Get change history
    pub fn change_history(&self) -> &[(T, String)] {
        &self.change_history
    }

    /// Get change history (alias for compatibility)
    pub fn get_change_history(&self) -> &[(T, String)] {
        &self.change_history
    }

    /// Check if parameter has adaptive behavior
    pub fn is_adaptive(&self) -> bool {
        self.adaptive_behavior
            .as_ref()
            .is_some_and(|b| b.is_adaptive())
    }

    /// Get default value
    pub const fn default_value(&self) -> &T {
        &self.default
    }
}
