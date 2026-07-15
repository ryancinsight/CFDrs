//! Error types for the state management system.
//!
//! All error types are centralized in `cfd_core::error`. This module provides
//! re-exports and backward-compatible result type aliases.

use cfd_core::error::{ParameterErrorKind, Result};

/// Result type for state management operations
pub type StateManagementResult<T> = Result<T>;

/// Result type for parameter operations
pub type ParameterResult<T> = std::result::Result<T, ParameterErrorKind>;

// Re-export Kind types for direct construction via inherent methods
pub use cfd_core::error::{
    ConstraintErrorKind as ConstraintError, DependencyErrorKind as DependencyError,
    ParameterErrorKind as ParameterError, RegistryErrorKind as RegistryError,
    ValidationErrorKind as ValidationError,
};

// Re-export the main error type
pub use cfd_core::error::Error as StateManagementError;

#[cfg(test)]
mod tests {
    use cfd_core::error::{ConstraintErrorKind, ParameterErrorKind, ValidationErrorKind};

    #[test]
    fn test_parameter_error_creation() {
        let error = ParameterErrorKind::not_found("amplitude", "serpentine");
        assert!(error.to_string().contains("amplitude"));
        assert!(error.to_string().contains("serpentine"));
    }

    #[test]
    fn test_validation_error_creation() {
        let error = ValidationErrorKind::rule_failed("wavelength", "must be positive");
        assert!(error.to_string().contains("wavelength"));
        assert!(error.to_string().contains("must be positive"));
    }

    #[test]
    fn test_constraint_error_creation() {
        let error = ConstraintErrorKind::range_violation(&-1.0, &0.0, &10.0);
        assert!(error.to_string().contains("-1"));
        assert!(error.to_string().contains("0"));
        assert!(error.to_string().contains("10"));
    }
}
