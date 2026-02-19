//! Error types for the cfd-optim SDT design optimizer.

use thiserror::Error;

/// All errors that can arise during candidate generation, physics evaluation,
/// and ranking within the SDT optimizer.
#[derive(Debug, Error)]
pub enum OptimError {
    /// The candidate space is empty — every combination was filtered out.
    #[error("no design candidates were generated; check parameter bounds")]
    EmptyCandidates,

    /// Fewer ranked results exist than the number requested.
    #[error("only {available} candidates available but {requested} were requested")]
    InsufficientCandidates { requested: usize, available: usize },

    /// A 1-D physics model returned an error for the given candidate.
    #[error("physics evaluation failed for candidate '{id}': {reason}")]
    PhysicsError { id: String, reason: String },

    /// A parameter value was outside the valid physical range.
    #[error("invalid parameter: {0}")]
    InvalidParameter(String),

    // ── Legacy variants (backward compatibility with original API) ──────────

    /// Legacy: no scenario results were provided.
    #[error("no scenario results were provided")]
    EmptyInput,

    /// Legacy: invalid top_k value.
    #[error("invalid top_k value: {0}")]
    InvalidTopK(usize),

    // Legacy: mesh export failed.
    // #[cfg(feature = "mesh-export")]
    // #[error("mesh export failed for scenario {scenario_id}: {message}")]
    // MeshExportFailed { scenario_id: String, message: String },
}
