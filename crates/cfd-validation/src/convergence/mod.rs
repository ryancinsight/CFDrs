//! Convergence analysis module for CFD validation
//!
//! Provides tools for grid convergence studies, Richardson extrapolation,
//! and convergence rate analysis following ASME V&V 20-2009 standards.

mod study;
mod richardson;
mod analysis;
mod criteria;

pub use study::ConvergenceStudy;
pub use richardson::RichardsonExtrapolation;
pub use analysis::{ConvergenceAnalysis, ConvergenceOrder};
pub use criteria::{ConvergenceStatus, ConvergenceCriterion, GridConvergenceIndex};

// Re-export core functionality
pub use study::compute_convergence_rate;