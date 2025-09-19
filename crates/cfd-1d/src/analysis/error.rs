//! Error types for analysis operations

use cfd_core::error::Error as CoreError;
use std::fmt;

/// Errors that can occur during resistance calculation
#[derive(Debug)]
pub enum ResistanceCalculationError {
    /// Hydraulic diameter is missing for a component
    MissingHydraulicDiameter,
    /// Resistance model calculation failed
    ModelError(String),
    /// Invalid flow conditions
    InvalidFlowConditions(String),
    /// Numerical computation error
    NumericalError(String),
}

impl fmt::Display for ResistanceCalculationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingHydraulicDiameter => {
                write!(f, "Hydraulic diameter is missing for component")
            }
            Self::ModelError(msg) => write!(f, "Resistance model calculation failed: {msg}"),
            Self::InvalidFlowConditions(msg) => write!(f, "Invalid flow conditions: {msg}"),
            Self::NumericalError(msg) => write!(f, "Numerical computation error: {msg}"),
        }
    }
}

impl std::error::Error for ResistanceCalculationError {}

impl From<ResistanceCalculationError> for CoreError {
    fn from(err: ResistanceCalculationError) -> Self {
        CoreError::InvalidInput(err.to_string())
    }
}
