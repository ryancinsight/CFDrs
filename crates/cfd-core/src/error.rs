//! Error types and result aliases for the CFD simulation suite.
//!
//! This module provides a comprehensive error handling system with:
//! - Structured error types for specific failure modes
//! - Automatic conversion from external error types
//! - Context extension trait for adding error context

use thiserror::Error;
use std::fmt;

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
    Numerical(String),
    
    /// Convergence failure
    #[error("Convergence failed: {0}")]
    Convergence(ConvergenceErrorKind),
    
    /// Solver-specific errors
    #[error("Solver error: {0}")]
    Solver(String),
    
    /// I/O errors
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    
    /// Dimension mismatch
    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },
    
    /// Index out of bounds
    #[error("Index out of bounds: {index} >= {size}")]
    IndexOutOfBounds { index: usize, size: usize },
    
    /// Not implemented
    #[error("Not implemented: {0}")]
    NotImplemented(String),
    
    /// Generic error with context
    #[error("{context}: {source}")]
    WithContext {
        context: String,
        #[source]
        source: Box<Error>,
    },
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
            Self::MaxIterationsExceeded { max } => 
                write!(f, "Maximum iterations ({}) exceeded", max),
            Self::StagnatedResidual { residual } => 
                write!(f, "Residual stagnated at {:.2e}", residual),
            Self::Diverged { norm } => 
                write!(f, "Solution diverged with norm {:.2e}", norm),
            Self::InvalidValue => 
                write!(f, "Invalid value (NaN or Inf) detected"),
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