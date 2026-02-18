//! CSG boolean operations on indexed meshes.
//!
//! Adapted from the BSP-tree approach used in `csgrs`, but operating
//! on indexed vertices/faces rather than duplicated polygon soups.
//!
//! Feature-gated behind `csg`.

pub mod classify;
pub mod split;
pub mod bsp;
pub mod boolean;

pub use boolean::{BooleanOp, csg_boolean};

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
