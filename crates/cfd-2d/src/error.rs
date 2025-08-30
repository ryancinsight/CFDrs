//! Error types for 2D CFD operations

use std::fmt;

/// Result type for 2D CFD operations
pub type Result<T> = std::result::Result<T, Error>;

/// Error types for 2D CFD operations
#[derive(Debug)]
pub enum Error {
    /// Invalid grid dimensions
    InvalidDimensions(String),
    /// Numerical convergence failure
    ConvergenceFailure(String),
    /// Invalid configuration
    InvalidConfiguration(String),
    /// Core error
    CoreError(cfd_core::error::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDimensions(msg) => write!(f, "Invalid dimensions: {}", msg),
            Self::ConvergenceFailure(msg) => write!(f, "Convergence failure: {}", msg),
            Self::InvalidConfiguration(msg) => write!(f, "Invalid configuration: {}", msg),
            Self::CoreError(e) => write!(f, "Core error: {}", e),
        }
    }
}

impl std::error::Error for Error {}

impl From<cfd_core::error::Error> for Error {
    fn from(err: cfd_core::error::Error) -> Self {
        Self::CoreError(err)
    }
}

impl From<String> for Error {
    fn from(msg: String) -> Self {
        Self::InvalidConfiguration(msg)
    }
}

impl From<&str> for Error {
    fn from(msg: &str) -> Self {
        Self::InvalidConfiguration(msg.to_string())
    }
}
