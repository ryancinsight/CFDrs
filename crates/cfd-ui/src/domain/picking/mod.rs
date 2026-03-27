//! Picking domain — types for GPU and CPU-based entity selection.

pub mod color_id;
pub mod ray;

use crate::domain::scene::selection::{SelectionGranularity, SelectionTarget};

/// A request to pick an entity at a screen-space location.
#[derive(Clone, Copy, Debug)]
pub struct PickQuery {
    /// Screen-space X coordinate (pixels from left).
    pub x: u32,
    /// Screen-space Y coordinate (pixels from top).
    pub y: u32,
    /// The selection granularity to use for interpreting the pick.
    pub granularity: SelectionGranularity,
}

/// The result of a successful pick operation.
#[derive(Clone, Copy, Debug)]
pub struct PickResult {
    /// The entity that was picked.
    pub target: SelectionTarget,
}
