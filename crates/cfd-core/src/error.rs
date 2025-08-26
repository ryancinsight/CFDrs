//! Error types and result aliases for the CFD simulation suite.
//!
//! This module provides a comprehensive error handling system with:
//! - Structured error types for specific failure modes
//! - Automatic conversion from external error types
//! - Context extension trait for adding error context

use std::fmt;
use thiserror::Error;

/// Core error type for CFD operations
#[derive(Debug, Error)]
pub enum Error {
    /// Invalid input parameters
    #[error("Invalid input: {0}")]
    InvalidInput(String),

    /// Invalid configuration
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),

    /// Numerical computation error
    #[error("Numerical error: {0}")]
    Numerical(NumericalErrorKind),

    /// Convergence failure
    #[error("Convergence failed: {0}")]
    Convergence(ConvergenceErrorKind),

    /// Plugin-related errors
    #[error("Plugin error: {0}")]
    Plugin(PluginErrorKind),

    /// Solver-specific errors
    #[error("Solver error: {0}")]
    Solver(String),

    /// I/O errors
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Dimension mismatch
    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch {
        /// Expected dimension
        expected: usize,
        /// Actual dimension
        actual: usize,
    },

    /// Index out of bounds
    #[error("Index out of bounds: {index} >= {size}")]
    IndexOutOfBounds {
        /// Index that was accessed
        index: usize,
        /// Size of the collection
        size: usize,
    },

    /// Not implemented
    #[error("Not implemented: {0}")]
    NotImplemented(String),

    /// Generic error with context
    #[error("{context}: {source}")]
    WithContext {
        /// Context description
        context: String,
        /// Underlying error
        #[source]
        source: Box<Error>,
    },
}

/// Plugin error variants
#[derive(Debug, Clone)]
pub enum PluginErrorKind {
    /// Plugin not found
    NotFound {
        /// Name of the plugin
        name: String,
    },
    /// Plugin already registered
    AlreadyRegistered {
        /// Name of the plugin
        name: String,
    },
    /// Plugin initialization failed
    InitializationFailed {
        /// Name of the plugin
        name: String,
        /// Reason for failure
        reason: String,
    },
    /// Plugin execution failed
    ExecutionFailed {
        /// Name of the plugin
        name: String,
        /// Reason for failure
        reason: String,
    },
    /// Dependency not satisfied
    DependencyNotSatisfied {
        /// Plugin name
        plugin: String,
        /// Missing dependency
        dependency: String,
    },
    /// Circular dependency detected
    CircularDependency {
        /// Chain of dependencies forming the cycle
        chain: Vec<String>,
    },
    /// Invalid plugin configuration
    InvalidConfiguration {
        /// Name of the plugin
        name: String,
        /// Reason for invalid configuration
        reason: String,
    },
}

impl fmt::Display for PluginErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFound { name } => write!(f, "Plugin '{}' not found", name),
            Self::AlreadyRegistered { name } => write!(f, "Plugin '{}' already registered", name),
            Self::InitializationFailed { name, reason } => {
                write!(f, "Plugin '{}' initialization failed: {}", name, reason)
            }
            Self::ExecutionFailed { name, reason } => {
                write!(f, "Plugin '{}' execution failed: {}", name, reason)
            }
            Self::DependencyNotSatisfied { plugin, dependency } => write!(
                f,
                "Plugin '{}' dependency '{}' not satisfied",
                plugin, dependency
            ),
            Self::CircularDependency { chain } => {
                write!(f, "Circular dependency detected: {}", chain.join(" -> "))
            }
            Self::InvalidConfiguration { name, reason } => {
                write!(f, "Invalid configuration for plugin '{}': {}", name, reason)
            }
        }
    }
}

/// Numerical error variants
#[derive(Debug, Clone)]
pub enum NumericalErrorKind {
    /// Division by zero
    DivisionByZero,
    /// Invalid value (NaN or Inf)
    InvalidValue { value: String },
    /// Underflow occurred
    Underflow { value: f64 },
    /// Overflow occurred
    Overflow { value: f64 },
    /// Matrix is singular
    SingularMatrix,
    /// Matrix is not positive definite
    NotPositiveDefinite,
    /// Invalid tolerance
    InvalidTolerance { tolerance: f64 },
    /// Insufficient precision
    InsufficientPrecision { achieved: f64, required: f64 },
}

impl fmt::Display for NumericalErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DivisionByZero => write!(f, "Division by zero"),
            Self::InvalidValue { value } => write!(f, "Invalid numerical value: {}", value),
            Self::Underflow { value } => write!(f, "Numerical underflow: {:.2e}", value),
            Self::Overflow { value } => write!(f, "Numerical overflow: {:.2e}", value),
            Self::SingularMatrix => write!(f, "Matrix is singular"),
            Self::NotPositiveDefinite => write!(f, "Matrix is not positive definite"),
            Self::InvalidTolerance { tolerance } => {
                write!(f, "Invalid tolerance: {:.2e}", tolerance)
            }
            Self::InsufficientPrecision { achieved, required } => write!(
                f,
                "Insufficient precision: achieved {:.2e}, required {:.2e}",
                achieved, required
            ),
        }
    }
}

/// Convergence error variants
#[derive(Debug, Clone)]
pub enum ConvergenceErrorKind {
    /// Maximum iterations exceeded
    MaxIterationsExceeded { max: usize },
    /// Residual did not decrease
    StagnatedResidual { residual: f64 },
    /// Solution diverged
    Diverged { norm: f64 },
    /// NaN or Inf detected
    InvalidValue,
}

impl fmt::Display for ConvergenceErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaxIterationsExceeded { max } => {
                write!(f, "Maximum iterations ({}) exceeded", max)
            }
            Self::StagnatedResidual { residual } => {
                write!(f, "Residual stagnated at {:.2e}", residual)
            }
            Self::Diverged { norm } => write!(f, "Solution diverged with norm {:.2e}", norm),
            Self::InvalidValue => write!(f, "Invalid value (NaN or Inf) detected"),
        }
    }
}

/// Result type alias for CFD operations
pub type Result<T> = std::result::Result<T, Error>;

/// Extension trait for adding context to errors
pub trait ErrorContext<T> {
    /// Add context to an error
    fn context(self, msg: impl Into<String>) -> Result<T>;

    /// Add context with a closure (lazy evaluation)
    fn with_context<F>(self, f: F) -> Result<T>
    where
        F: FnOnce() -> String;
}

impl<T> ErrorContext<T> for Result<T> {
    fn context(self, msg: impl Into<String>) -> Result<T> {
        self.map_err(|e| Error::WithContext {
            context: msg.into(),
            source: Box::new(e),
        })
    }

    fn with_context<F>(self, f: F) -> Result<T>
    where
        F: FnOnce() -> String,
    {
        self.map_err(|e| Error::WithContext {
            context: f(),
            source: Box::new(e),
        })
    }
}

/// Helper function to convert Option to Result
pub fn require<T>(opt: Option<T>, msg: impl Into<String>) -> Result<T> {
    opt.ok_or_else(|| Error::InvalidInput(msg.into()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_context() {
        let result: Result<()> = Err(Error::InvalidInput("test".into()));
        let with_context = result.context("Additional context");
        assert!(with_context.is_err());
        let error_msg = format!("{}", with_context.unwrap_err());
        assert!(error_msg.contains("Additional context"));
    }

    #[test]
    fn test_require() {
        let some_value = Some(42);
        assert_eq!(require(some_value, "missing").unwrap(), 42);

        let none_value: Option<i32> = None;
        assert!(require(none_value, "missing").is_err());
    }
}
