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
    #[error("Failed to converge after {iterations} iterations (residual: {residual})")]
    ConvergenceFailure {
        /// Number of iterations performed
        iterations: usize,
        /// Final residual value
        residual: f64,
    },

    /// Plugin-related error
    #[error("Plugin error: {0}")]
    PluginError(String),

    /// I/O error
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),

    /// Serialization error
    #[error("Serialization error: {0}")]
    SerializationError(String),

    /// Generic error with context
    #[error("{context}: {source}")]
    WithContext {
        /// Context description
        context: String,
        /// Underlying error
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}

/// Convenience type alias for Results in this crate
pub type Result<T> = std::result::Result<T, Error>;

impl Error {
    /// Add context to an error
    pub fn context<S: Into<String>>(self, context: S) -> Self {
        Self::WithContext {
            context: context.into(),
            source: Box::new(self),
        }
    }
}