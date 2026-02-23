//! Error taxonomy for mesh operations.

use thiserror::Error;

use crate::core::index::{VertexId, FaceId, EdgeId};
use crate::core::scalar::Point3r;

/// Alias for `Result<T, MeshError>`.
pub type MeshResult<T> = Result<T, MeshError>;

/// Comprehensive error type covering all mesh operations.
#[derive(Error, Debug)]
pub enum MeshError {
    // ── Topology ──────────────────────────────────────────────
    /// A face references a vertex that does not exist.
    #[error("face {face} references invalid vertex {vertex}")]
    InvalidVertexRef {
        /// The face containing the bad reference.
        face: FaceId,
        /// The invalid vertex index.
        vertex: VertexId,
    },

    /// An edge is non-manifold (shared by != 2 faces).
    #[error("edge {edge} is non-manifold: shared by {count} faces (expected 2)")]
    NonManifoldEdge {
        /// The problematic edge.
        edge: EdgeId,
        /// Actual face count.
        count: usize,
    },

    /// A vertex is non-manifold (fan not a single disk).
    #[error("vertex {vertex} is non-manifold")]
    NonManifoldVertex {
        /// The problematic vertex.
        vertex: VertexId,
    },

    /// Mesh has inconsistent winding order.
    #[error("inconsistent winding at face {face}")]
    InconsistentWinding {
        /// The face with wrong orientation.
        face: FaceId,
    },

    /// Mesh is not connected (multiple components).
    #[error("mesh has {count} disconnected components (expected 1)")]
    Disconnected {
        /// Number of components.
        count: usize,
    },

    // ── Geometry ──────────────────────────────────────────────
    /// Degenerate face (zero area).
    #[error("face {face} is degenerate (area ≈ 0)")]
    DegenerateFace {
        /// The zero-area face.
        face: FaceId,
    },

    /// Self-intersection detected.
    #[error("self-intersection between faces {a} and {b}")]
    SelfIntersection {
        /// First intersecting face.
        a: FaceId,
        /// Second intersecting face.
        b: FaceId,
    },

    /// A coordinate is NaN or Inf.
    #[error("invalid coordinate at vertex {vertex}: {point:?}")]
    InvalidCoordinate {
        /// Vertex with bad coordinates.
        vertex: VertexId,
        /// The bad point.
        point: Point3r,
    },

    // ── Watertight ────────────────────────────────────────────
    /// Mesh is not watertight (boundary edges exist).
    #[error("mesh has {count} boundary edges — not watertight")]
    NotWatertight {
        /// Number of boundary edges.
        count: usize,
    },

    /// Gap detected between regions.
    #[error("gap of width {width:.6} mm between regions")]
    GapDetected {
        /// Gap width in model units.
        width: f64,
    },

    // ── CSG ───────────────────────────────────────────────────
    /// BSP tree construction failed.
    #[error("BSP construction failed: {reason}")]
    BspError {
        /// Reason for failure.
        reason: String,
    },

    /// Boolean operation produced empty result.
    #[error("boolean {op} produced empty mesh")]
    EmptyBooleanResult {
        /// Operation name.
        op: String,
    },

    // ── Quality ───────────────────────────────────────────────
    /// Mesh fails quality threshold.
    #[error("quality score {score:.3} below threshold {threshold:.3}")]
    QualityBelowThreshold {
        /// Achieved score.
        score: f64,
        /// Required minimum.
        threshold: f64,
    },

    // ── Channel ───────────────────────────────────────────────
    /// Invalid channel geometry.
    #[error("channel error: {message}")]
    ChannelError {
        /// Description.
        message: String,
    },

    // ── IO ────────────────────────────────────────────────────
    /// IO failure.
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// JSON parsing failure.
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),

    /// Generic message.
    #[error("{0}")]
    Other(String),
}
