//! Error types for the state management system
//!
//! This module defines comprehensive error types for parameter management,
//! validation, and state synchronization operations.

use std::fmt;
use thiserror::Error;

/// Result type for state management operations
pub type StateManagementResult<T> = Result<T, StateManagementError>;

/// Result type for parameter operations
pub type ParameterResult<T> = Result<T, ParameterError>;

/// Main error type for state management operations
#[derive(Error, Debug)]
pub enum StateManagementError {
    /// Parameter-related errors
    #[error("Parameter error: {0}")]
    Parameter(#[from] ParameterError),

    /// Validation errors
    #[error("Validation error: {0}")]
    Validation(#[from] ValidationError),

    /// Registry errors
    #[error("Registry error: {0}")]
    Registry(#[from] RegistryError),

    /// Constraint errors
    #[error("Constraint error: {0}")]
    Constraint(#[from] ConstraintError),

    /// Dependency resolution errors
    #[error("Dependency error: {0}")]
    Dependency(#[from] DependencyError),

    /// Serialization/deserialization errors
    #[error("Serialization error: {0}")]
    Serialization(String),

    /// Configuration errors
    #[error("Configuration error: {0}")]
    Configuration(String),
}

/// Errors related to parameter operations
#[derive(Error, Debug, Clone)]
pub enum ParameterError {
    /// Parameter not found
    #[error("Parameter '{name}' not found in domain '{domain}'")]
    NotFound { name: String, domain: String },

    /// Invalid parameter value
    #[error("Invalid value for parameter '{name}': {value:?}. {constraint}")]
    InvalidValue {
        name: String,
        value: String,
        constraint: String,
    },

    /// Parameter type mismatch
    #[error("Type mismatch for parameter '{name}': expected {expected}, got {actual}")]
    TypeMismatch {
        name: String,
        expected: String,
        actual: String,
    },

    /// Parameter is read-only
    #[error("Parameter '{name}' is read-only and cannot be modified")]
    ReadOnly { name: String },

    /// Parameter dependency not satisfied
    #[error("Parameter '{name}' dependency '{dependency}' not satisfied")]
    DependencyNotSatisfied { name: String, dependency: String },

    /// Circular dependency detected
    #[error("Circular dependency detected involving parameter '{name}'")]
    CircularDependency { name: String },

    /// Adaptation failed
    #[error("Failed to adapt parameter '{name}': {reason}")]
    AdaptationFailed { name: String, reason: String },
}

/// Errors related to validation operations
#[derive(Error, Debug)]
pub enum ValidationError {
    /// Validation rule failed
    #[error("Validation failed for '{field}': {message}")]
    RuleFailed { field: String, message: String },

    /// Multiple validation failures
    #[error("Multiple validation failures: {failures:?}")]
    Multiple { failures: Vec<ValidationError> },

    /// Custom validation error
    #[error("Custom validation error: {message}")]
    Custom { message: String },

    /// Constraint violation
    #[error("Constraint violation: {constraint}")]
    ConstraintViolation { constraint: String },
}

/// Errors related to parameter registry operations
#[derive(Error, Debug)]
pub enum RegistryError {
    /// Manager not found
    #[error("Parameter manager for domain '{domain}' not found")]
    ManagerNotFound { domain: String },

    /// Manager already registered
    #[error("Parameter manager for domain '{domain}' already registered")]
    ManagerAlreadyRegistered { domain: String },

    /// Registry is locked
    #[error("Registry is locked and cannot be modified")]
    RegistryLocked,

    /// Initialization failed
    #[error("Registry initialization failed: {reason}")]
    InitializationFailed { reason: String },

    /// Update conflict
    #[error("Update conflict for parameter '{parameter}': {reason}")]
    UpdateConflict { parameter: String, reason: String },
}

/// Errors related to parameter constraints
#[derive(Error, Debug)]
pub enum ConstraintError {
    /// Range constraint violation
    #[error("Value {value:?} is outside allowed range [{min:?}, {max:?}]")]
    RangeViolation {
        value: String,
        min: String,
        max: String,
    },

    /// Set constraint violation
    #[error("Value {value:?} is not in allowed set: {allowed:?}")]
    SetViolation { value: String, allowed: Vec<String> },

    /// Custom constraint violation
    #[error("Custom constraint violation: {message}")]
    CustomViolation { message: String },

    /// Constraint composition error
    #[error("Constraint composition error: {reason}")]
    CompositionError { reason: String },
}

/// Errors related to parameter dependencies
#[derive(Error, Debug)]
pub enum DependencyError {
    /// Dependency cycle detected
    #[error("Dependency cycle detected: {cycle:?}")]
    CycleDetected { cycle: Vec<String> },

    /// Missing dependency
    #[error("Missing dependency '{dependency}' for parameter '{parameter}'")]
    MissingDependency {
        parameter: String,
        dependency: String,
    },

    /// Dependency resolution failed
    #[error("Failed to resolve dependencies for parameter '{parameter}': {reason}")]
    ResolutionFailed { parameter: String, reason: String },

    /// Invalid dependency graph
    #[error("Invalid dependency graph: {reason}")]
    InvalidGraph { reason: String },
}

impl From<ConstraintError> for ParameterError {
    fn from(error: ConstraintError) -> Self {
        Self::InvalidValue {
            name: "unknown".to_string(),
            value: "unknown".to_string(),
            constraint: error.to_string(),
        }
    }
}

impl From<ValidationError> for ParameterError {
    fn from(error: ValidationError) -> Self {
        Self::InvalidValue {
            name: "unknown".to_string(),
            value: "unknown".to_string(),
            constraint: error.to_string(),
        }
    }
}

impl ParameterError {
    /// Create a not found error
    #[must_use]
    pub fn not_found(name: &str, domain: &str) -> Self {
        Self::NotFound {
            name: name.to_string(),
            domain: domain.to_string(),
        }
    }

    /// Create an invalid value error
    pub fn invalid_value(name: &str, value: &dyn fmt::Debug, constraint: &str) -> Self {
        Self::InvalidValue {
            name: name.to_string(),
            value: format!("{value:?}"),
            constraint: constraint.to_string(),
        }
    }

    /// Create a type mismatch error
    #[must_use]
    pub fn type_mismatch(name: &str, expected: &str, actual: &str) -> Self {
        Self::TypeMismatch {
            name: name.to_string(),
            expected: expected.to_string(),
            actual: actual.to_string(),
        }
    }

    /// Create a read-only error
    #[must_use]
    pub fn read_only(name: &str) -> Self {
        Self::ReadOnly {
            name: name.to_string(),
        }
    }

    /// Create a dependency not satisfied error
    #[must_use]
    pub fn dependency_not_satisfied(name: &str, dependency: &str) -> Self {
        Self::DependencyNotSatisfied {
            name: name.to_string(),
            dependency: dependency.to_string(),
        }
    }

    /// Create a circular dependency error
    #[must_use]
    pub fn circular_dependency(name: &str) -> Self {
        Self::CircularDependency {
            name: name.to_string(),
        }
    }

    /// Create an adaptation failed error
    #[must_use]
    pub fn adaptation_failed(name: &str, reason: &str) -> Self {
        Self::AdaptationFailed {
            name: name.to_string(),
            reason: reason.to_string(),
        }
    }
}

impl ValidationError {
    /// Create a rule failed error
    #[must_use]
    pub fn rule_failed(field: &str, message: &str) -> Self {
        Self::RuleFailed {
            field: field.to_string(),
            message: message.to_string(),
        }
    }

    /// Create a multiple failures error
    #[must_use]
    pub const fn multiple(failures: Vec<Self>) -> Self {
        Self::Multiple { failures }
    }

    /// Create a custom error
    #[must_use]
    pub fn custom(message: &str) -> Self {
        Self::Custom {
            message: message.to_string(),
        }
    }

    /// Create a constraint violation error
    #[must_use]
    pub fn constraint_violation(constraint: &str) -> Self {
        Self::ConstraintViolation {
            constraint: constraint.to_string(),
        }
    }
}

impl ConstraintError {
    /// Create a range violation error
    pub fn range_violation(
        value: &dyn fmt::Debug,
        min: &dyn fmt::Debug,
        max: &dyn fmt::Debug,
    ) -> Self {
        Self::RangeViolation {
            value: format!("{value:?}"),
            min: format!("{min:?}"),
            max: format!("{max:?}"),
        }
    }

    /// Create a set violation error
    pub fn set_violation(value: &dyn fmt::Debug, allowed: &[&dyn fmt::Debug]) -> Self {
        Self::SetViolation {
            value: format!("{value:?}"),
            allowed: allowed.iter().map(|v| format!("{v:?}")).collect(),
        }
    }

    /// Create a custom violation error
    #[must_use]
    pub fn custom_violation(message: &str) -> Self {
        Self::CustomViolation {
            message: message.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parameter_error_creation() {
        let error = ParameterError::not_found("amplitude", "serpentine");
        assert!(error.to_string().contains("amplitude"));
        assert!(error.to_string().contains("serpentine"));
    }

    #[test]
    fn test_validation_error_creation() {
        let error = ValidationError::rule_failed("wavelength", "must be positive");
        assert!(error.to_string().contains("wavelength"));
        assert!(error.to_string().contains("must be positive"));
    }

    #[test]
    fn test_constraint_error_creation() {
        let error = ConstraintError::range_violation(&-1.0, &0.0, &10.0);
        assert!(error.to_string().contains("-1"));
        assert!(error.to_string().contains("0"));
        assert!(error.to_string().contains("10"));
    }
}
