//! Extended error types for I/O operations.
//!
//! This module extends the core error types with I/O-specific errors.

use thiserror::Error;

/// I/O-specific error types
#[derive(Error, Debug)]
pub enum IoError {
    /// CSV error from the csv crate
    #[error("CSV error: {0}")]
    Csv(#[from] csv::Error),
    
    /// HDF5 error (when feature is enabled)
    #[cfg(feature = "hdf5")]
    #[error("HDF5 error: {0}")]
    Hdf5(#[from] hdf5::Error),
    
    /// File format error
    #[error("Unsupported file format: {0}")]
    UnsupportedFormat(String),
    
    /// Data validation error
    #[error("Data validation failed: {0}")]
    ValidationError(String),
    
    /// Missing required field
    #[error("Missing required field: {0}")]
    MissingField(String),
    
    /// Core error
    #[error(transparent)]
    Core(#[from] cfd_core::error::Error),
}

/// Result type for I/O operations
pub type Result<T> = std::result::Result<T, IoError>;

// Allow conversion from IoError to core Error for compatibility
impl From<IoError> for cfd_core::error::Error {
    fn from(err: IoError) -> Self {
        match err {
            IoError::Core(e) => e,
            other => cfd_core::error::Error::External(other.to_string()),
        }
    }
}

// Re-export the Context trait from core
pub use cfd_core::error::Context;

/// Extension methods for `IoError`
impl IoError {
    /// Create an unsupported format error
    pub fn unsupported_format<S: Into<String>>(format: S) -> Self {
        Self::UnsupportedFormat(format.into())
    }
    
    /// Create a validation error
    pub fn validation<S: Into<String>>(msg: S) -> Self {
        Self::ValidationError(msg.into())
    }
    
    /// Create a missing field error
    pub fn missing_field<S: Into<String>>(field: S) -> Self {
        Self::MissingField(field.into())
    }
}