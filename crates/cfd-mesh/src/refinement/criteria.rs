//! Refinement criteria definitions

use crate::topology::{Cell, Vertex};
use nalgebra::RealField;
use thiserror::Error;

/// Refinement errors
#[derive(Debug, Error)]
pub enum RefinementError {
    /// Invalid mesh structure provided
    #[error("Invalid mesh: {0}")]
    InvalidMesh(String),
    /// Maximum refinement level reached
    #[error("Refinement limit reached: level {0}")]
    MaxLevelReached(usize),
    /// Minimum cell size reached
    #[error("Cell too small: size {0}")]
    MinSizeReached(f64),
    /// Invalid refinement criteria specified
    #[error("Invalid refinement criteria: {0}")]
    InvalidCriteria(String),
}

/// Refinement criteria for adaptive mesh refinement
pub enum RefinementCriterion<T: RealField + Copy> {
    /// Refine based on solution gradient
    Gradient { 
        /// Solution field values
        field: Vec<T>, 
        /// Gradient threshold for refinement
        threshold: T 
    },
    /// Refine based on error estimate
    Error { 
        /// Error field values
        error_field: Vec<T>, 
        /// Error threshold for refinement
        threshold: T 
    },
    /// Refine based on geometric features
    Geometric {
        /// Curvature threshold for geometric refinement
        curvature_threshold: T,
        /// Feature angle threshold in radians
        feature_angle: T,
    },
    /// Custom refinement function
    Custom(Box<dyn Fn(&Cell, &[Vertex<T>]) -> bool + Send + Sync>),
}

impl<T: RealField + Copy> std::fmt::Debug for RefinementCriterion<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Gradient { field, threshold } => f
                .debug_struct("Gradient")
                .field("field_len", &field.len())
                .field("threshold", threshold)
                .finish(),
            Self::Error {
                error_field,
                threshold,
            } => f
                .debug_struct("Error")
                .field("error_field_len", &error_field.len())
                .field("threshold", threshold)
                .finish(),
            Self::Geometric {
                curvature_threshold,
                feature_angle,
            } => f
                .debug_struct("Geometric")
                .field("curvature_threshold", curvature_threshold)
                .field("feature_angle", feature_angle)
                .finish(),
            Self::Custom(_) => f.debug_struct("Custom").finish(),
        }
    }
}
