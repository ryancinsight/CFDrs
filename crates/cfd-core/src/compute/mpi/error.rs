//! MPI-specific error types and handling

use crate::error::{Error, Result};
use std::fmt;

/// MPI-specific errors
#[derive(Debug, thiserror::Error)]
pub enum MpiError {
    /// MPI initialization failed
    #[error("MPI initialization failed: {0}")]
    InitializationError(String),

    /// MPI communication error
    #[error("MPI communication error: {0}")]
    CommunicationError(String),

    /// Invalid rank specification
    #[error("Invalid MPI rank: {rank}, expected 0 <= rank < {size}")]
    InvalidRank { rank: i32, size: i32 },

    /// Domain decomposition error
    #[error("Domain decomposition error: {0}")]
    DecompositionError(String),

    /// Load balancing failed
    #[error("Load balancing failed: {0}")]
    LoadBalancingError(String),

    /// Ghost cell exchange error
    #[error("Ghost cell exchange error: {0}")]
    GhostCellError(String),

    /// Data distribution error
    #[error("Data distribution error: {0}")]
    DistributionError(String),

    /// Parallel I/O error
    #[error("Parallel I/O error: {0}")]
    IoError(String),

    /// MPI not available
    #[error("MPI not available: {0}")]
    NotAvailable(String),
}

impl From<MpiError> for Error {
    fn from(err: MpiError) -> Self {
        Error::ExternalError(err.to_string())
    }
}

/// Result type for MPI operations
pub type MpiResult<T> = std::result::Result<T, MpiError>;
