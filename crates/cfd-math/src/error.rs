//! Error types for mathematical operations

use std::error::Error;
use std::fmt;

/// Math-specific error types
#[derive(Debug, Clone)]
pub enum MathError {
    /// Dimension mismatch in vector/matrix operations
    DimensionMismatch,
    /// Invalid input parameters
    InvalidInput(String),
    /// Numerical overflow
    Overflow,
    /// Numerical underflow
    Underflow,
    /// Convergence failure
    ConvergenceFailed,
}

impl fmt::Display for MathError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MathError::DimensionMismatch => write!(f, "Dimension mismatch"),
            MathError::InvalidInput(msg) => write!(f, "Invalid input: {msg}"),
            MathError::Overflow => write!(f, "Numerical overflow"),
            MathError::Underflow => write!(f, "Numerical underflow"),
            MathError::ConvergenceFailed => write!(f, "Convergence failed"),
        }
    }
}

impl Error for MathError {}

impl From<MathError> for cfd_core::error::Error {
    fn from(err: MathError) -> Self {
        use cfd_core::error::{ConvergenceErrorKind, NumericalErrorKind};
        match err {
            MathError::DimensionMismatch => {
                cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string())
            }
            MathError::Overflow => {
                cfd_core::error::Error::Numerical(NumericalErrorKind::Overflow {
                    value: f64::INFINITY,
                })
            }
            MathError::Underflow => {
                cfd_core::error::Error::Numerical(NumericalErrorKind::Underflow { value: 0.0 })
            }
            MathError::ConvergenceFailed => {
                cfd_core::error::Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded {
                    max: 1000,
                })
            }
            MathError::InvalidInput(msg) => cfd_core::error::Error::InvalidInput(msg),
        }
    }
}

/// Result type for math operations
pub type Result<T> = std::result::Result<T, cfd_core::error::Error>;
