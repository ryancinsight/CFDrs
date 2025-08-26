//! Error types for CSG operations

use thiserror::Error;
/// Error types for CSG operations
#[derive(Debug, Error)]
pub enum CsgError {
    /// Invalid mesh for CSG operation
    #[error("Invalid mesh: {0}")]
    InvalidMesh(String),
    /// CSG operation failed
    #[error("CSG operation failed: {0}")]
    OperationFailed(String),
    /// Invalid geometry parameters
    #[error("Invalid geometry parameters: {0}")]
    InvalidParameters(String),
    /// STL export error
    #[error("STL export failed: {0}")]
    ExportError(String),
}
