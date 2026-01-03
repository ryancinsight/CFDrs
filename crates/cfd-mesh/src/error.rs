//! Error types for mesh operations
//!
//! This module provides error types for mesh topology and operations.
//! Errors are designed to provide detailed context for debugging mesh issues.

use thiserror::Error;

/// Mesh operation errors
#[derive(Debug, Clone, Error)]
pub enum MeshError {
    /// Invalid mesh structure
    #[error("Invalid mesh: {0}")]
    InvalidMesh(String),
    /// Refinement failed
    #[error("Refinement failed: {0}")]
    RefinementFailed(String),
    /// Feature not implemented
    #[error("Not implemented: {0}")]
    NotImplemented(String),
    /// Grid generation error
    #[error("Grid error: {0}")]
    GridError(String),
    /// Face topology error
    #[error(transparent)]
    Face(#[from] FaceError),
}

/// Face topology errors
///
/// Errors specific to face construction and validation.
/// A face is a 2-simplex requiring at least 3 distinct vertices.
///
/// # Mathematical Background
///
/// A face in a mesh represents a 2-dimensional simplex embedded in ℝ³.
/// The simplex condition requires:
/// - Minimum 3 vertices (triangular face)
/// - All vertices must be distinct (non-degenerate)
/// - Vertices should be non-collinear (positive area)
///
/// # References
///
/// - Edelsbrunner, H. "Geometry and Topology for Mesh Generation" (2001)
/// - Shewchuk, J. "Triangle: Engineering a 2D Quality Mesh Generator" (1996)
#[derive(Debug, Clone, Error, PartialEq, Eq)]
pub enum FaceError {
    /// Face has fewer than the required minimum vertices
    ///
    /// A valid face requires at least 3 vertices to form a 2-simplex.
    #[error("Insufficient vertices: got {count}, minimum required is {minimum}")]
    InsufficientVertices {
        /// Number of vertices provided
        count: usize,
        /// Minimum required (always 3 for faces)
        minimum: usize,
    },

    /// Face contains duplicate vertex indices
    ///
    /// Repeated vertices create degenerate geometry with zero area,
    /// causing undefined behavior in geometric calculations.
    #[error("Duplicate vertex {vertex} found at positions {positions:?}")]
    DuplicateVertex {
        /// The duplicated vertex index
        vertex: usize,
        /// Positions in the vertex list where duplicates occur
        positions: Vec<usize>,
    },

    /// Vertex index exceeds valid range
    ///
    /// Used when validating face against a mesh's vertex count.
    #[error("Invalid vertex index {vertex}: must be less than {max_valid}")]
    InvalidVertexIndex {
        /// The invalid vertex index
        vertex: usize,
        /// Maximum valid index (exclusive)
        max_valid: usize,
    },
}

/// Result type for mesh operations
pub type Result<T> = std::result::Result<T, MeshError>;

/// Result type specifically for face operations
pub type FaceResult<T> = std::result::Result<T, FaceError>;
