//! Professional error handling for CFD computations
//!
//! This module provides structured, machine-readable errors following
//! Rust best practices and leveraging thiserror properly.

use thiserror::Error;
use std::time::Duration;

/// Configuration-related errors with structured information
#[derive(Error, Debug, Clone, PartialEq)]
pub enum ConfigurationError {
    #[error("Missing required field '{field}' in configuration")]
    MissingField { field: String },
    
    #[error("Invalid value for field '{field}': {reason}")]
    InvalidValue { field: String, reason: String },
    
    #[error("Incompatible configuration: {reason}")]
    Incompatible { reason: String },
    
    #[error("Configuration file not found: {path}")]
    FileNotFound { path: String },
    
    #[error("Failed to parse configuration: {reason}")]
    ParseError { reason: String },
}

/// Input validation errors with specific details
#[derive(Error, Debug, Clone, PartialEq)]
pub enum InputError {
    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },
    
    #[error("Invalid range: [{min}, {max}]")]
    InvalidRange { min: f64, max: f64 },
    
    #[error("Negative value not allowed for {parameter}")]
    NegativeValue { parameter: String },
    
    #[error("Value {value} out of bounds: must be in [{min}, {max}]")]
    OutOfBounds { value: f64, min: f64, max: f64 },
    
    #[error("Empty input provided for {parameter}")]
    EmptyInput { parameter: String },
}

/// Numerical computation errors
#[derive(Error, Debug, Clone, PartialEq)]
pub enum NumericalError {
    #[error("Division by zero in {operation}")]
    DivisionByZero { operation: String },
    
    #[error("Overflow in {operation}")]
    Overflow { operation: String },
    
    #[error("Underflow in {operation}")]
    Underflow { operation: String },
    
    #[error("Invalid floating point operation: {details}")]
    InvalidOperation { details: String },
    
    #[error("Matrix is singular with condition number {cond:.3e}")]
    SingularMatrix { cond: f64 },
    
    #[error("Matrix dimensions incompatible: ({rows1}, {cols1}) vs ({rows2}, {cols2})")]
    MatrixDimensionMismatch {
        rows1: usize,
        cols1: usize,
        rows2: usize,
        cols2: usize,
    },
}

/// Convergence-related errors with recoverability information
#[derive(Error, Debug, Clone, PartialEq)]
pub enum ConvergenceError {
    #[error("Maximum iterations ({max}) exceeded with residual {residual:.3e}")]
    MaxIterationsExceeded { max: usize, residual: f64 },
    
    #[error("Solution diverging: residual increased from {old:.3e} to {new:.3e}")]
    Divergence { old: f64, new: f64 },
    
    #[error("Solution stagnated: change {change:.3e} below threshold {threshold:.3e}")]
    Stagnation { change: f64, threshold: f64 },
    
    #[error("NaN detected in {location} at iteration {iteration}")]
    NaNDetected { location: String, iteration: usize },
}

impl ConvergenceError {
    /// Determine if this convergence error is potentially recoverable
    pub fn is_recoverable(&self) -> bool {
        match self {
            Self::MaxIterationsExceeded { .. } => true,  // Can increase iterations
            Self::Stagnation { .. } => true,              // Can adjust tolerance
            Self::Divergence { .. } => false,             // Fundamental problem
            Self::NaNDetected { .. } => false,            // Numerical instability
        }
    }
}

/// Mesh-related errors
#[derive(Error, Debug, Clone, PartialEq)]
pub enum MeshError {
    #[error("Invalid mesh topology: {reason}")]
    InvalidTopology { reason: String },
    
    #[error("Degenerate element at index {index}: {reason}")]
    DegenerateElement { index: usize, reason: String },
    
    #[error("Mesh quality below threshold: {metric} = {value:.3e} < {threshold:.3e}")]
    PoorQuality {
        metric: String,
        value: f64,
        threshold: f64,
    },
    
    #[error("Incompatible mesh dimensions: {reason}")]
    IncompatibleDimensions { reason: String },
}

/// Solver-specific errors
#[derive(Error, Debug, Clone, PartialEq)]
pub enum SolverError {
    #[error("Linear solver failed")]
    LinearSolver(#[from] LinearSolverError),
    
    #[error("Nonlinear solver failed")]
    NonlinearSolver(#[from] NonlinearSolverError),
    
    #[error("Time integration failed at t = {time:.3e}")]
    TimeIntegration { time: f64, source: Box<Error> },
    
    #[error("Preconditioner setup failed: {reason}")]
    PreconditionerSetup { reason: String },
}

/// Linear solver specific errors
#[derive(Error, Debug, Clone, PartialEq)]
pub enum LinearSolverError {
    #[error("Matrix is singular or nearly singular")]
    SingularMatrix,
    
    #[error("Insufficient memory for factorization")]
    InsufficientMemory,
    
    #[error("Convergence failed after {iterations} iterations")]
    ConvergenceFailed { iterations: usize },
    
    #[error("Preconditioner application failed")]
    PreconditionerFailed,
}

/// Nonlinear solver specific errors
#[derive(Error, Debug, Clone, PartialEq)]
pub enum NonlinearSolverError {
    #[error("Newton method failed to converge")]
    NewtonConvergenceFailed,
    
    #[error("Line search failed to find acceptable step")]
    LineSearchFailed,
    
    #[error("Jacobian computation failed")]
    JacobianComputationFailed,
}

/// I/O related errors
#[derive(Error, Debug)]
pub enum IOError {
    #[error("Failed to read file {path}")]
    ReadError {
        path: String,
        #[source]
        source: std::io::Error,
    },
    
    #[error("Failed to write file {path}")]
    WriteError {
        path: String,
        #[source]
        source: std::io::Error,
    },
    
    #[error("Invalid file format for {path}: expected {expected}, got {actual}")]
    InvalidFormat {
        path: String,
        expected: String,
        actual: String,
    },
    
    #[error("Parsing error in file {path} at line {line}")]
    ParseError {
        path: String,
        line: usize,
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}

/// Main error type for the CFD library
#[derive(Error, Debug)]
pub enum Error {
    #[error("Configuration error")]
    Configuration(#[from] ConfigurationError),
    
    #[error("Invalid input")]
    Input(#[from] InputError),
    
    #[error("Numerical error")]
    Numerical(#[from] NumericalError),
    
    #[error("Convergence error")]
    Convergence(#[from] ConvergenceError),
    
    #[error("Mesh error")]
    Mesh(#[from] MeshError),
    
    #[error("Solver error")]
    Solver(#[from] SolverError),
    
    #[error("I/O error")]
    IO(#[from] IOError),
    
    #[error("Operation timed out after {duration:?}")]
    Timeout { duration: Duration },
    
    #[error("Feature not implemented: {feature}")]
    NotImplemented { feature: String },
    
    #[error("External library error: {library}")]
    External {
        library: String,
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}

impl Error {
    /// Check if this error is potentially recoverable
    pub fn is_recoverable(&self) -> bool {
        match self {
            Self::Convergence(e) => e.is_recoverable(),
            Self::Timeout { .. } => true,  // Can retry with longer timeout
            Self::Input(_) => false,       // Bad input won't fix itself
            Self::Configuration(_) => false, // Configuration must be fixed
            Self::Numerical(_) => false,   // Numerical issues are fundamental
            Self::Mesh(_) => false,        // Mesh must be regenerated
            Self::Solver(e) => match e {
                SolverError::LinearSolver(LinearSolverError::ConvergenceFailed { .. }) => true,
                _ => false,
            },
            Self::IO(_) => false,          // I/O errors need external fixes
            Self::NotImplemented { .. } => false,
            Self::External { .. } => false,
        }
    }
    
    /// Get a user-friendly description of how to potentially fix this error
    pub fn recovery_hint(&self) -> Option<&'static str> {
        match self {
            Self::Convergence(ConvergenceError::MaxIterationsExceeded { .. }) => {
                Some("Try increasing max_iterations in the solver configuration")
            }
            Self::Convergence(ConvergenceError::Stagnation { .. }) => {
                Some("Try relaxing the convergence tolerance or using a different solver")
            }
            Self::Timeout { .. } => {
                Some("Try increasing the timeout duration or optimizing the computation")
            }
            Self::Solver(SolverError::LinearSolver(LinearSolverError::ConvergenceFailed { .. })) => {
                Some("Try using a different preconditioner or increasing iterations")
            }
            _ => None,
        }
    }
}

/// Result type alias for the CFD library
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_structured_errors() {
        // Configuration errors are now structured
        let err = ConfigurationError::MissingField {
            field: "tolerance".into(),
        };
        assert_eq!(
            err.to_string(),
            "Missing required field 'tolerance' in configuration"
        );
        
        // Input errors provide specific information
        let err = InputError::DimensionMismatch {
            expected: 3,
            actual: 2,
        };
        assert_eq!(
            err.to_string(),
            "Dimension mismatch: expected 3, got 2"
        );
    }
    
    #[test]
    fn test_recoverability() {
        // Recoverable error
        let err = ConvergenceError::MaxIterationsExceeded {
            max: 100,
            residual: 1e-3,
        };
        assert!(err.is_recoverable());
        
        // Non-recoverable error
        let err = ConvergenceError::NaNDetected {
            location: "velocity field".into(),
            iteration: 42,
        };
        assert!(!err.is_recoverable());
    }
    
    #[test]
    fn test_error_conversion() {
        // Natural conversion using From trait
        let convergence_err = ConvergenceError::Divergence {
            old: 1e-3,
            new: 1e-2,
        };
        
        let main_err: Error = convergence_err.into();
        assert!(!main_err.is_recoverable());
    }
    
    #[test]
    fn test_recovery_hints() {
        let err = Error::Convergence(ConvergenceError::MaxIterationsExceeded {
            max: 100,
            residual: 1e-3,
        });
        
        assert!(err.recovery_hint().is_some());
        assert!(err.recovery_hint().unwrap().contains("max_iterations"));
    }
}