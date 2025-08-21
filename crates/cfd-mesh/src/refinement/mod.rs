//! Mesh refinement and adaptation module
//!
//! This module provides adaptive mesh refinement (AMR) capabilities for improving
//! solution accuracy in regions of high gradients or errors.

pub mod criteria;
pub mod strategies;
pub mod uniform;
pub mod adaptive;
pub mod quality;
pub mod operations;

// Re-export main types
pub use criteria::{RefinementCriterion, RefinementError};
pub use strategies::{RefinementStrategy, RefinementResult};
pub use uniform::UniformRefinement;
pub use adaptive::AdaptiveRefinement;
pub use quality::QualityBasedRefinement;
pub use operations::RefinementOperations;

// Named constants for refinement parameters
pub const DEFAULT_MAX_REFINEMENT_LEVEL: usize = 5;
pub const DEFAULT_MIN_CELL_SIZE: f64 = 1e-6;
pub const DEFAULT_ERROR_THRESHOLD: f64 = 1e-3;
pub const DEFAULT_GRADIENT_THRESHOLD: f64 = 0.1;