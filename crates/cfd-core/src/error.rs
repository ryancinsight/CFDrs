//! Error types and result aliases for the CFD simulation suite.
//!
//! This module provides a comprehensive error handling system with:
//! - Structured error types for specific failure modes
//! - Automatic conversion from external error types
//! - Context extension trait for adding error context

use thiserror::Error;

/// Specific types of numerical errors that can occur in CFD computations
#[derive(Error, Debug, Clone, PartialEq)]
pub enum NumericalErrorKind {
    /// Division by zero encountered
    #[error("Division by zero")]
    DivisionByZero,
    
    /// Matrix is singular or ill-conditioned
    #[error("Matrix is singular or ill-conditioned")]
    SingularMatrix,
    
    /// Invalid floating point operation (NaN or Inf)
    #[error("Invalid floating point operation (NaN or Inf)")]
    InvalidFpOperation,
    
    /// Value out of valid range
    #[error("Value out of valid range: {0}")]
    OutOfRange(String),
    
    /// Underflow in computation
    #[error("Numerical underflow")]
    Underflow,
    
    /// Overflow in computation
    #[error("Numerical overflow")]
    Overflow,
    
    /// Loss of precision
    #[error("Loss of precision: {0}")]
    PrecisionLoss(String),
}

/// Specific types of convergence failures
#[derive(Error, Debug, Clone, PartialEq)]
pub enum ConvergenceErrorKind {
    /// Maximum iterations exceeded
    #[error("Maximum iterations ({max}) exceeded")]
    MaxIterationsExceeded { max: usize },
    
    /// Divergence detected
    #[error("Solution diverging: residual = {residual:.3e}")]
    Divergence { residual: f64 },
    
    /// Stagnation detected
    #[error("Solution stagnated: change = {change:.3e}")]
    Stagnation { change: f64 },
    
    /// Oscillation detected
    #[error("Solution oscillating")]
    Oscillation,
}

/// Specific types of mesh/geometry errors
#[derive(Error, Debug, Clone, PartialEq)]
pub enum MeshErrorKind {
    /// Invalid topology
    #[error("Invalid topology: {0}")]
    InvalidTopology(String),
    
    /// Degenerate element
    #[error("Degenerate element at index {index}")]
    DegenerateElement { index: usize },
    
    /// Non-manifold geometry
    #[error("Non-manifold geometry detected")]
    NonManifold,
    
    /// Insufficient quality
    #[error("Mesh quality below threshold: {quality} < {threshold}")]
    InsufficientQuality { quality: f64, threshold: f64 },
}

/// Specific types of plugin errors
#[derive(Error, Debug, Clone, PartialEq)]
pub enum PluginErrorKind {
    /// Plugin not found
    #[error("Plugin not found: {name}")]
    NotFound { name: String },
    
    /// Plugin already registered
    #[error("Plugin already registered: {0}")]
    AlreadyRegistered(String),
    
    /// Plugin initialization failed
    #[error("Plugin initialization failed: {reason}")]
    InitializationFailed { reason: String },
    
    /// Plugin incompatible version
    #[error("Plugin version {version} incompatible with required {required}")]
    IncompatibleVersion { version: String, required: String },
    
    /// Plugin execution error
    #[error("Plugin execution error: {0}")]
    ExecutionError(String),
    
    /// Missing dependency
    #[error("Plugin {plugin} requires missing dependency: {dependency}")]
    MissingDependency { plugin: String, dependency: String },
    
    /// Circular dependency detected
    #[error("Circular dependency detected involving plugin: {0}")]
    CircularDependency(String),
    
    /// Plugin is in use by another plugin
    #[error("Cannot remove plugin {plugin} as it is used by {used_by}")]
    InUse { plugin: String, used_by: String },
}

/// Main error type for CFD operations
#[derive(Error, Debug)]
pub enum Error {
    /// Invalid configuration provided
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),

    /// Invalid input provided
    #[error("Invalid input: {0}")]
    InvalidInput(String),

    /// Invalid state encountered
    #[error("Invalid state: {0}")]
    InvalidState(String),

    /// Timeout error
    #[error("Timeout after {seconds} seconds")]
    Timeout { seconds: u64 },

    /// Mesh or geometry error
    #[error("Mesh error: {0}")]
    Mesh(#[from] MeshErrorKind),

    /// Numerical error during computation
    #[error("Numerical error: {0}")]
    Numerical(#[from] NumericalErrorKind),

    /// Convergence failure
    #[error("Convergence failure: {0}")]
    Convergence(#[from] ConvergenceErrorKind),

    /// Plugin-related error
    #[error("Plugin error: {0}")]
    Plugin(#[from] PluginErrorKind),

    /// I/O error from the standard library
    #[error(transparent)]
    Io(#[from] std::io::Error),
    
    /// JSON serialization error from serde_json
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),
    
    /// CSV processing error
    #[error("CSV error: {0}")]
    CsvError(String),
    
    /// Serialization error
    #[error("Serialization error: {0}")]
    SerializationError(String),
    
    /// IO error variant for compatibility
    #[error("IO error: {0}")]
    IoError(std::io::Error),

    /// Not implemented
    #[error("Not implemented: {0}")]
    NotImplemented(String),

    /// Error with additional context
    #[error("{message}")]
    WithContext {
        /// Context message
        message: String,
        /// Source error
        #[source]
        source: Box<Error>,
    },
    
    /// Generic error for external errors not covered above
    /// This should be used sparingly and only for truly unexpected errors
    #[error("External error: {0}")]
    External(String),
}

/// Convenience type alias for Results in this crate
pub type Result<T> = std::result::Result<T, Error>;

impl Error {
    /// Add context to an error
    #[must_use]
    pub fn context<S: Into<String>>(self, context: S) -> Self {
        Self::WithContext {
            message: context.into(),
            source: Box::new(self),
        }
    }
    
    /// Create a numerical division by zero error
    pub fn division_by_zero() -> Self {
        Self::Numerical(NumericalErrorKind::DivisionByZero)
    }
    
    /// Create a singular matrix error
    pub fn singular_matrix() -> Self {
        Self::Numerical(NumericalErrorKind::SingularMatrix)
    }
    
    /// Create a max iterations exceeded error
    pub fn max_iterations_exceeded(max: usize) -> Self {
        Self::Convergence(ConvergenceErrorKind::MaxIterationsExceeded { max })
    }
    
    /// Create a divergence error
    pub fn divergence(residual: f64) -> Self {
        Self::Convergence(ConvergenceErrorKind::Divergence { residual })
    }
}

/// Extension trait for adding context to results
pub trait Context<T> {
    /// Adds a string context to an error
    fn context<S: Into<String>>(self, context: S) -> Result<T>;
    
    /// Adds a lazy context to an error (only evaluated if there's an error)
    fn with_context<S, F>(self, f: F) -> Result<T>
    where
        S: Into<String>,
        F: FnOnce() -> S;
}

impl<T, E> Context<T> for std::result::Result<T, E>
where
    E: Into<Error>,
{
    fn context<S: Into<String>>(self, context: S) -> Result<T> {
        self.map_err(|error| Error::WithContext {
            message: context.into(),
            source: Box::new(error.into()),
        })
    }
    
    fn with_context<S, F>(self, f: F) -> Result<T>
    where
        S: Into<String>,
        F: FnOnce() -> S,
    {
        self.map_err(|error| Error::WithContext {
            message: f().into(),
            source: Box::new(error.into()),
        })
    }
}

// Convenience constructors for common error patterns
impl Error {
    /// Create an invalid configuration error
    pub fn invalid_config<S: Into<String>>(msg: S) -> Self {
        Self::InvalidConfiguration(msg.into())
    }
    
    /// Create an invalid input error
    pub fn invalid_input<S: Into<String>>(msg: S) -> Self {
        Self::InvalidInput(msg.into())
    }
    
    /// Check if this is a numerical error
    pub fn is_numerical(&self) -> bool {
        matches!(self, Self::Numerical(_))
    }
    
    /// Check if this is a convergence error
    pub fn is_convergence(&self) -> bool {
        matches!(self, Self::Convergence(_))
    }
    
    /// Check if this error is recoverable (e.g., convergence issues that might be fixed with different parameters)
    pub fn is_recoverable(&self) -> bool {
        matches!(
            self,
            Self::Convergence(
                ConvergenceErrorKind::MaxIterationsExceeded { .. } | ConvergenceErrorKind::Stagnation { .. }
            ) | Self::Timeout { .. }
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_numerical_error_kinds() {
        let err = Error::division_by_zero();
        assert!(err.is_numerical());
        assert!(!err.is_convergence());
        
        let err = Error::singular_matrix();
        assert!(err.is_numerical());
    }
    
    #[test]
    fn test_convergence_error_kinds() {
        let err = Error::max_iterations_exceeded(1000);
        assert!(err.is_convergence());
        assert!(err.is_recoverable());
        
        let err = Error::divergence(1e10);
        assert!(err.is_convergence());
        assert!(!err.is_recoverable());
    }
    
    #[test]
    fn test_context_trait() {
        use std::fs;
        
        // This would fail but demonstrates the API
        let result: Result<String> = fs::read_to_string("nonexistent.txt")
            .context("Failed to read configuration file");
        
        assert!(result.is_err());
        if let Err(e) = result {
            let error_string = e.to_string();
            assert!(error_string.contains("Failed to read configuration file"));
        }
    }
    
    #[test]
    fn test_with_context_lazy() {
        let expensive_context = || {
            // This would only be evaluated if there's an error
            format!("Context computed at {}", std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("CRITICAL: Add proper error handling")
                .as_secs())
        };
        
        let result: Result<i32> = Err(Error::invalid_input("test"))
            .with_context(expensive_context);
        
        assert!(result.is_err());
    }
    
    #[test]
    fn test_error_conversion() {
        // Test that std::io::Error converts automatically
        fn read_file() -> Result<String> {
            std::fs::read_to_string("nonexistent.txt")?;
            Ok(String::new())
        }
        
        assert!(read_file().is_err());
        
        // Test that serde_json::Error converts automatically
        fn parse_json() -> Result<serde_json::Value> {
            let invalid_json = "{ invalid json }";
            Ok(serde_json::from_str(invalid_json)?)
        }
        
        assert!(parse_json().is_err());
    }
}