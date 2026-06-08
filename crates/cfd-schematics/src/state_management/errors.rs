//! Error types for the state management system.
//!
//! All error types are centralized in `cfd_core::error`. This module provides
//! re-exports and backward-compatible result type aliases.

use std::fmt;
use cfd_core::error::{
    ConstraintErrorKind, DependencyErrorKind, Error, ParameterErrorKind, RegistryErrorKind,
    Result, ValidationErrorKind,
};

/// Result type for state management operations
pub type StateManagementResult<T> = Result<T>;

/// Result type for parameter operations
pub type ParameterResult<T> = std::result::Result<T, ParameterErrorKind>;

// Re-export Kind types for direct construction
pub use cfd_core::error::{ConstraintErrorKind as ConstraintError, DependencyErrorKind as DependencyError, ParameterErrorKind as ParameterError, RegistryErrorKind as RegistryError, ValidationErrorKind as ValidationError};

// Re-export the main error type
pub use cfd_core::error::Error as StateManagementError;

/// Convenience constructors for parameter errors
pub trait ParameterErrorExt {
    /// Create a not-found error
    fn not_found(name: &str, domain: &str) -> Self;
    /// Create an invalid-value error
    fn invalid_value(name: &str, value: &dyn fmt::Debug, constraint: &str) -> Self;
    /// Create a type-mismatch error
    fn type_mismatch(name: &str, expected: &str, actual: &str) -> Self;
    /// Create a read-only error
    fn read_only(name: &str) -> Self;
    /// Create a dependency-not-satisfied error
    fn dependency_not_satisfied(name: &str, dependency: &str) -> Self;
    /// Create a circular-dependency error
    fn circular_dependency(name: &str) -> Self;
    /// Create an adaptation-failed error
    fn adaptation_failed(name: &str, reason: &str) -> Self;
}

impl ParameterErrorExt for Error {
    fn not_found(name: &str, domain: &str) -> Self {
        Error::Parameter(ParameterErrorKind::NotFound {
            name: name.to_string(),
            domain: domain.to_string(),
        })
    }

    fn invalid_value(name: &str, value: &dyn fmt::Debug, constraint: &str) -> Self {
        Error::Parameter(ParameterErrorKind::InvalidValue {
            name: name.to_string(),
            value: format!("{value:?}"),
            constraint: constraint.to_string(),
        })
    }

    fn type_mismatch(name: &str, expected: &str, actual: &str) -> Self {
        Error::Parameter(ParameterErrorKind::TypeMismatch {
            name: name.to_string(),
            expected: expected.to_string(),
            actual: actual.to_string(),
        })
    }

    fn read_only(name: &str) -> Self {
        Error::Parameter(ParameterErrorKind::ReadOnly {
            name: name.to_string(),
        })
    }

    fn dependency_not_satisfied(name: &str, dependency: &str) -> Self {
        Error::Parameter(ParameterErrorKind::DependencyNotSatisfied {
            name: name.to_string(),
            dependency: dependency.to_string(),
        })
    }

    fn circular_dependency(name: &str) -> Self {
        Error::Parameter(ParameterErrorKind::CircularDependency {
            name: name.to_string(),
        })
    }

    fn adaptation_failed(name: &str, reason: &str) -> Self {
        Error::Parameter(ParameterErrorKind::AdaptationFailed {
            name: name.to_string(),
            reason: reason.to_string(),
        })
    }
}

/// Convenience constructors for validation errors
pub trait ValidationErrorExt {
    /// Create a rule-failed error
    fn rule_failed(field: &str, message: &str) -> Self;
    /// Create a multiple-failures error
    fn multiple(failures: Vec<ValidationErrorKind>) -> Self;
    /// Create a custom error
    fn custom(message: &str) -> Self;
    /// Create a constraint-violation error
    fn constraint_violation(constraint: &str) -> Self;
}

impl ValidationErrorExt for Error {
    fn rule_failed(field: &str, message: &str) -> Self {
        Error::Validation(ValidationErrorKind::RuleFailed {
            field: field.to_string(),
            message: message.to_string(),
        })
    }

    fn multiple(failures: Vec<ValidationErrorKind>) -> Self {
        Error::Validation(ValidationErrorKind::Multiple { failures })
    }

    fn custom(message: &str) -> Self {
        Error::Validation(ValidationErrorKind::Custom {
            message: message.to_string(),
        })
    }

    fn constraint_violation(constraint: &str) -> Self {
        Error::Validation(ValidationErrorKind::ConstraintViolation {
            constraint: constraint.to_string(),
        })
    }
}

/// Convenience constructors for constraint errors
pub trait ConstraintErrorExt {
    /// Create a range-violation error
    fn range_violation(value: &dyn fmt::Debug, min: &dyn fmt::Debug, max: &dyn fmt::Debug) -> Self;
    /// Create a set-violation error
    fn set_violation(value: &dyn fmt::Debug, allowed: &[&dyn fmt::Debug]) -> Self;
    /// Create a custom-violation error
    fn custom_violation(message: &str) -> Self;
}

impl ConstraintErrorExt for Error {
    fn range_violation(value: &dyn fmt::Debug, min: &dyn fmt::Debug, max: &dyn fmt::Debug) -> Self {
        Error::Constraint(ConstraintErrorKind::RangeViolation {
            value: format!("{value:?}"),
            min: format!("{min:?}"),
            max: format!("{max:?}"),
        })
    }

    fn set_violation(value: &dyn fmt::Debug, allowed: &[&dyn fmt::Debug]) -> Self {
        Error::Constraint(ConstraintErrorKind::SetViolation {
            value: format!("{value:?}"),
            allowed: allowed.iter().map(|v| format!("{v:?}")).collect(),
        })
    }

    fn custom_violation(message: &str) -> Self {
        Error::Constraint(ConstraintErrorKind::CustomViolation {
            message: message.to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_core::error::ParameterErrorKind;

    #[test]
    fn test_parameter_error_creation() {
        let error = Error::not_found("amplitude", "serpentine");
        assert!(error.to_string().contains("amplitude"));
        assert!(error.to_string().contains("serpentine"));
    }

    #[test]
    fn test_validation_error_creation() {
        let error = Error::rule_failed("wavelength", "must be positive");
        assert!(error.to_string().contains("wavelength"));
        assert!(error.to_string().contains("must be positive"));
    }

    #[test]
    fn test_constraint_error_creation() {
        let error = Error::range_violation(&-1.0, &0.0, &10.0);
        assert!(error.to_string().contains("-1"));
        assert!(error.to_string().contains("0"));
        assert!(error.to_string().contains("10"));
    }
}
