//! Error types and result aliases for the CFD simulation suite.

use thiserror::Error;

/// Main error type for CFD operations
#[derive(Error, Debug)]
pub enum Error {
    /// Invalid configuration provided
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),

    /// Invalid mesh or geometry
    #[error("Invalid mesh: {0}")]
    InvalidMesh(String),

    /// Numerical error during computation
    #[error("Numerical error: {0}")]
    NumericalError(String),

    /// Convergence failure
    #[error("Convergence failure: {0}")]
    ConvergenceFailure(String),

    /// Plugin-related error
    #[error("Plugin error: {0}")]
    PluginError(String),

    /// I/O error
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    
    /// CSV error
    #[error("CSV error: {0}")]
    CsvError(String),

    /// Serialization error
    #[error("Serialization error: {0}")]
    SerializationError(String),

    /// Not implemented
    #[error("Not implemented: {0}")]
    NotImplemented(String),

    /// Error with context
    #[error("{message}: {source}")]
    WithContext {
        /// Context message
        message: String,
        /// Source error
        #[source]
        source: Box<Error>,
    },
}

/// Convenience type alias for Results in this crate
pub type Result<T> = std::result::Result<T, Error>;

impl Error {
    /// Add context to an error
    pub fn context<S: Into<String>>(self, context: S) -> Self {
        Self::WithContext {
            message: context.into(),
            source: Box::new(self),
        }
    }
}