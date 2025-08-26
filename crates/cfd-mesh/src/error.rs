//! Error types for mesh operations

use std::fmt;
/// Mesh operation errors
#[derive(Debug, Clone)]
pub enum MeshError {
    /// Invalid mesh structure
    InvalidMesh(String),
    /// Refinement failed
    RefinementFailed(String),
    /// Grid generation error
    GridError(String),
}
impl fmt::Display for MeshError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MeshError::InvalidMesh(msg) => write!(f, "Invalid mesh: {msg}"),
            MeshError::RefinementFailed(msg) => write!(f, "Refinement failed: {msg}"),
            MeshError::GridError(msg) => write!(f, "Grid error: {msg}"),
        }
    }
impl std::error::Error for MeshError {}
/// Result type for mesh operations
pub type Result<T> = std::result::Result<T, MeshError>;
