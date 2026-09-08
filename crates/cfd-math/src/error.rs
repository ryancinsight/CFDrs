//! Error types for mathematical operations.
//!
//! The `From<LetoError>` bridge is defined in `cfd-core` so that any crate
//! depending on `cfd-core` can use `?` to propagate `leto::LetoError` through
//! functions returning `cfd_core::error::Result`.

/// Result type for math operations (uses `cfd_core::error::Error`).
pub type Result<T> = std::result::Result<T, cfd_core::error::Error>;
