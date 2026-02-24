//! Canonical CSG module surface (application-layer path).

pub mod arrangement;
pub mod boolean;
pub mod broad_phase;
pub mod clip;
pub mod coplanar;
pub mod corefine;
pub mod intersect;
pub mod reconstruct;

pub use boolean::{csg_boolean, csg_boolean_indexed, BooleanOp, CsgNode};
pub use intersect::IntersectionType;

use thiserror::Error;

/// Error type for CSG operations.
#[derive(Error, Debug)]
pub enum CsgError {
    /// BSP construction failed.
    #[error("BSP construction failed: {0}")]
    BspError(String),
    /// Boolean operation produced empty result.
    #[error("boolean operation produced empty mesh: {0}")]
    EmptyResult(String),
    /// Generic CSG error.
    #[error("CSG error: {0}")]
    Other(String),
}
