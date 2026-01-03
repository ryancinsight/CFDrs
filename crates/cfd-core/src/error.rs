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

    /// GPU compute error
    #[error("GPU compute error: {0}")]
    GpuCompute(String),

    /// Convergence failure
    #[error("Convergence failed: {0}")]
    Convergence(ConvergenceErrorKind),

    /// Plugin-related errors
    #[error("Plugin error: {0}")]
    Plugin(PluginErrorKind),

    /// Solver-specific errors
    #[error("Solver error: {0}")]
    Solver(String),

    /// Numeric conversion error
    #[error("Conversion error: {0}")]
    ConversionError(String),

    /// I/O errors
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Serialization/deserialization errors
    #[error("Serialization error: {0}")]
    Serialization(String),

    /// Unsupported operation
    #[error("Unsupported operation: {0}")]
    UnsupportedOperation(String),

    /// Physical invariant violation
    #[error("Physical invariant violation: {0}")]
    PhysicsViolation(String),

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
    /// Dependency missing
    DependencyMissing {
        /// Plugin name
        plugin: String,
        /// Missing dependency
        dependency: String,
    },
    /// Plugin has dependents
    HasDependents {
        /// Plugin name
        plugin: String,
        /// List of dependent plugins
        dependents: Vec<String>,
    },
}

impl fmt::Display for PluginErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFound { name } => write!(f, "Plugin '{name}' not found"),
            Self::AlreadyRegistered { name } => write!(f, "Plugin '{name}' already registered"),
            Self::InitializationFailed { name, reason } => {
                write!(f, "Plugin '{name}' initialization failed: {reason}")
            }
            Self::ExecutionFailed { name, reason } => {
                write!(f, "Plugin '{name}' execution failed: {reason}")
            }
            Self::DependencyNotSatisfied { plugin, dependency } => write!(
                f,
                "Plugin '{plugin}' dependency '{dependency}' not satisfied"
            ),
            Self::CircularDependency { chain } => {
                write!(f, "Circular dependency detected: {}", chain.join(" -> "))
            }
            Self::InvalidConfiguration { name, reason } => {
                write!(f, "Invalid configuration for plugin '{name}': {reason}")
            }
            Self::DependencyMissing { plugin, dependency } => {
                write!(f, "Plugin '{plugin}' missing dependency '{dependency}'")
            }
            Self::HasDependents { plugin, dependents } => {
                write!(
                    f,
                    "Plugin '{}' has dependents: {}",
                    plugin,
                    dependents.join(", ")
                )
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
    InvalidValue {
        /// String representation of the invalid value
        value: String,
    },
    /// Underflow occurred
    Underflow {
        /// Value that caused underflow
        value: f64,
    },
    /// Overflow occurred
    Overflow {
        /// Value that caused overflow
        value: f64,
    },
    /// Matrix is singular
    SingularMatrix,
    /// Matrix is not positive definite
    NotPositiveDefinite,
    /// Invalid tolerance
    InvalidTolerance {
        /// The invalid tolerance value
        tolerance: f64,
    },
    /// Insufficient precision
    InsufficientPrecision {
        /// Precision achieved
        achieved: f64,
        /// Precision required
        required: f64,
    },
}

impl fmt::Display for NumericalErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DivisionByZero => write!(f, "Division by zero"),
            Self::InvalidValue { value } => write!(f, "Invalid numerical value: {value}"),
            Self::Underflow { value } => write!(f, "Numerical underflow: {value:.2e}"),
            Self::Overflow { value } => write!(f, "Numerical overflow: {value:.2e}"),
            Self::SingularMatrix => write!(f, "Matrix is singular"),
            Self::NotPositiveDefinite => write!(f, "Matrix is not positive definite"),
            Self::InvalidTolerance { tolerance } => {
                write!(f, "Invalid tolerance: {tolerance:.2e}")
            }
            Self::InsufficientPrecision { achieved, required } => write!(
                f,
                "Insufficient precision: achieved {achieved:.2e}, required {required:.2e}"
            ),
        }
    }
}

/// Convergence error variants
#[derive(Debug, Clone)]
pub enum ConvergenceErrorKind {
    /// Maximum iterations exceeded
    MaxIterationsExceeded {
        /// Maximum iteration limit
        max: usize,
    },
    /// Residual did not decrease
    StagnatedResidual {
        /// Current residual value
        residual: f64,
    },
    /// Solution diverged
    Diverged {
        /// Norm of the diverging solution
        norm: f64,
    },
    /// NaN or Inf detected
    InvalidValue,
    /// Algorithm breakdown (e.g., zero inner product in Krylov methods)
    Breakdown,
}

impl fmt::Display for ConvergenceErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaxIterationsExceeded { max } => {
                write!(f, "Maximum iterations ({max}) exceeded")
            }
            Self::StagnatedResidual { residual } => {
                write!(f, "Residual stagnated at {residual:.2e}")
            }
            Self::Diverged { norm } => write!(f, "Solution diverged with norm {norm:.2e}"),
            Self::InvalidValue => write!(f, "Invalid value (NaN or Inf) detected"),
            Self::Breakdown => write!(f, "Algorithm breakdown (zero inner product)"),
        }
    }
}

/// Result type alias for CFD operations
pub type Result<T> = std::result::Result<T, Error>;

/// Extension trait for adding context to errors
pub trait ErrorContext<T> {
    /// Add context to an error
    ///
    /// # Errors
    /// Returns the original error with added context
    fn context(self, msg: impl Into<String>) -> Result<T>;

    /// Add context with a closure (lazy evaluation)
    ///
    /// # Errors
    /// Returns the original error with context added via the closure
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
///
/// # Errors
/// Returns `Error::InvalidInput` if the option is None
pub fn require<T>(opt: Option<T>, msg: impl Into<String>) -> Result<T> {
    opt.ok_or_else(|| Error::InvalidInput(msg.into()))
}

impl From<&str> for Error {
    fn from(msg: &str) -> Self {
        Error::InvalidInput(msg.to_string())
    }
}

impl From<String> for Error {
    fn from(msg: String) -> Self {
        Error::InvalidInput(msg)
    }
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
