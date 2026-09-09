//! Error types for CFDrs file-format and checkpoint I/O.

/// I/O result alias.
pub type Result<T> = std::result::Result<T, Error>;

/// Error type for CFDrs file-format and checkpoint I/O.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Input data violates an I/O format contract.
    #[error("Invalid input: {0}")]
    InvalidInput(String),

    /// Filesystem or stream I/O failed.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Serialization or backend encoding failed.
    #[error("Serialization error: {0}")]
    Serialization(String),
}

impl From<&str> for Error {
    fn from(message: &str) -> Self {
        Self::InvalidInput(message.to_string())
    }
}

impl From<String> for Error {
    fn from(message: String) -> Self {
        Self::InvalidInput(message)
    }
}
