//! Convergence analysis module for CFD validation
//!
//! Provides tools for grid convergence studies, Richardson extrapolation,
//! and convergence rate analysis following ASME V&V 20-2009 standards.

mod analysis;
mod criteria;
mod richardson;
mod study;

pub use analysis::{ConvergenceAnalysis, ConvergenceOrder};
pub use criteria::{
    ConvergenceCriterion, ConvergenceMonitor, ConvergenceStatus, GridConvergenceIndex,
};
pub use richardson::{richardson_extrapolate, RichardsonExtrapolation};
pub use study::ConvergenceStudy;

// Re-export core functionality
pub use study::compute_convergence_rate;
