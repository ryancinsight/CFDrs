//! Unified polygon clipping engine — 3-D half-space and 2-D Boolean operations.
//!
//! # Why two sub-modules?
//!
//! The two clipping scenarios require fundamentally different algorithms:
//!
//! | Scenario | Sub-module | Algorithm | Predicate |
//! |----------|------------|-----------|-----------|
//! | Non-coplanar 3-D faces | [`halfspace`] | Sutherland-Hodgman half-space | `orient_3d` |
//! | Coplanar 2-D faces | [`polygon2d`] | CDT-backed Boolean algebra | 2-D cross product |
//!
//! When two faces are coplanar, `orient_3d` returns `Degenerate` for every
//! query point (signed tet volume = 0), making the 3-D half-space clipper
//! return false "inside" for all inputs.  The [`polygon2d`] path projects
//! geometry to the shared plane and uses exact 2-D predicates instead.
//!
//! # Public surface
//!
//! All callers use `csg::clip::*` — there is no need to know which sub-module
//! provides a given function:
//!
//! **3-D functions** (from [`halfspace`]):
//! - [`clip_polygon_to_halfplane`]
//! - [`clip_triangle_to_halfplane`]
//! - [`fan_triangulate`]
//!
//! **2-D functions** (from [`polygon2d`]):
//! - [`boolean_clip`] + [`ClipOp`] — canonical CDT-backed Boolean SPI
//! - [`clip_polygon_to_triangle`] — specialized convex triangle clip
//! - [`split_polygon_outside_triangle`] — complement decomposition
//! - [`polygon_area`]

/// 3-D Sutherland-Hodgman half-space clipping and fan triangulation.
pub mod halfspace;
/// 2-D polygon Boolean operations (CDT-backed production path + test oracles).
pub(crate) mod polygon2d;

// ── Flat re-exports — callers use `csg::clip::*` exclusively ─────────────────

pub use halfspace::{clip_polygon_to_halfplane, clip_triangle_to_halfplane, fan_triangulate};
pub use polygon2d::{
    boolean_clip, clip_polygon_to_triangle, polygon_area, split_polygon_outside_triangle, ClipOp,
};
