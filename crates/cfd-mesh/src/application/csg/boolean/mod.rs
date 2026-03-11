//! Boolean operations: union, difference, intersection.
//!
//! This module exposes two API levels:
//!
//! ## Low-level face-soup API (backward compatible)
//!
//! [`csg_boolean`] operates on `Vec<FaceData>` + a shared `VertexPool`.
//! It preserves the historical API shape but now dispatches across:
//! - coplanar 2-D Boolean handling for flat meshes,
//! - containment fast-paths for nested/disjoint solids,
//! - arrangement/corefine reconstruction for intersecting meshes.
//!
//! ## High-level `IndexedMesh` API (new in Phase 8)
//!
//! [`csg_boolean`] takes two `IndexedMesh` objects, merges their
//! vertex pools, runs the Boolean pipeline, and reconstructs a fresh mesh.
//!
//! [`CsgNode`] provides a composable CSG tree that evaluates lazily.

pub mod containment;
pub mod indexed;
pub mod operations;
pub mod tree;
mod union_strategy;

pub use indexed::{
    csg_boolean, csg_boolean_indexed, csg_boolean_indexed_best_effort,
    csg_boolean_indexed_tolerant, csg_boolean_nary, csg_boolean_nary_union,
};
pub use operations::BooleanOp;
pub use tree::CsgNode;
