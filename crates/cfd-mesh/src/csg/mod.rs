//! CSG Boolean operations on indexed meshes.
//!
//! ## Pipeline
//!
//! ```text
//! Input: IndexedMesh A, IndexedMesh B
//!   │
//!   ├─ broad_phase  → candidate triangle pairs (AABB overlap)
//!   ├─ intersect    → exact triangle-triangle intersection (orient_3d)
//!   ├─ clip         → Sutherland–Hodgman polygon clipping
//!   ├─ reconstruct  → fragment merging into a fresh IndexedMesh
//!   └─ boolean      → Union / Intersection / Difference
//! ```
//!
//! ## Quick Start
//!
//! For an `IndexedMesh`-level API:
//! ```rust,ignore
//! use cfd_mesh::csg::{CsgNode, BooleanOp};
//!
//! let result = CsgNode::Difference {
//!     left:  Box::new(CsgNode::Leaf(cube_mesh)),
//!     right: Box::new(CsgNode::Leaf(sphere_mesh)),
//! }.evaluate().expect("CSG evaluation");
//! ```
//!
//! For a lower-level face-soup API (backward compatible):
//! ```rust,ignore
//! use cfd_mesh::csg::{BooleanOp, csg_boolean};
//! let faces = csg_boolean(BooleanOp::Union, &faces_a, &faces_b, &mut pool)?;
//! ```

// ── New exact-pipeline modules ────────────────────────────────────────────────

/// AABB broad phase — candidate triangle pair search.
pub mod broad_phase;
/// Exact triangle–triangle intersection (Shewchuk predicates).
pub mod intersect;
/// Sutherland–Hodgman polygon clipping with exact predicates.
pub mod clip;
/// Fragment mesh reconstruction after Boolean classification.
pub mod reconstruct;

// ── Legacy BSP modules (kept for backward compatibility) ──────────────────────

pub mod classify;
pub mod split;
pub mod bsp;
pub mod boolean;

// ── Re-exports ────────────────────────────────────────────────────────────────

pub use boolean::{BooleanOp, CsgNode, csg_boolean, csg_boolean_indexed};
pub use intersect::IntersectionType;

// ── Error type ────────────────────────────────────────────────────────────────

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
