//! Convection scheme correction computations
//!
//! This module contains helper functions for computing high-order convection corrections
//! used in deferred correction approaches (Patankar 1980 §5.4.3).

mod quick;
mod tvd;

pub use quick::{compute_quick_correction_x, compute_quick_correction_y};
pub use tvd::{compute_tvd_correction_x, compute_tvd_correction_y};
