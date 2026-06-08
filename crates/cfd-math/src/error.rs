//! Error types for mathematical operations

/// Result type for math operations
///
/// Uses the unified `cfd_core::error::Error` hierarchy for all error variants.
pub type Result<T> = std::result::Result<T, cfd_core::error::Error>;
