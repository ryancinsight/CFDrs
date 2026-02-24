//! Coplanar 2-D Boolean operations for flat surface meshes.
//!
//! When both operand meshes are flat (all triangles share a common plane),
//! the standard 3-D intersection/classification pipelines break down because
//! `orient_3d` degenerates for all coplanar query points.
//!
//! ## Algorithm — Canonical 2-D clipping (via `clip` module)
//!
//! For each source triangle T:
//!
//! 1. **Project** vertices to the shared 2-D plane via `PlaneBasis::project`.
//! 2. **Fast-path — disjoint**: if T's AABB has no overlap with any opposing
//!    AABB, T is fully disjoint from the opposing mesh → emit (outside) or
//!    skip (inside).
//! 3. **Fast-path — fully inside/outside**: test centroid + 3 vertices against
//!    the AABB-candidate triangles using 2-D point-in-triangle.  If unanimous
//!    inside → emit (inside) or skip (outside); if unanimous outside → skip
//!    (inside) or emit (outside).
//! 4. **Exact boundary clipping**: for each AABB-overlapping opposing triangle
//!    B_i, compute `clip_polygon_to_triangle(T, B_i)` via the `clip` module's
//!    canonical CDT-backed pipeline plus specialized convex helpers.
//!    - **Intersection**: emit each `T ∩ B_i` fragment (disjoint since B-tris
//!      tile the disk without interior overlap).
//!    - **Outside**: progressive subtraction — start with `remaining = {T}`;
//!      for each B_i remove `inside(B_i)` pieces from every remaining poly.
//!
//! ## Why 2-D clipping is necessary
//!
//! `clip_polygon_to_halfplane` from `csg::clip` uses `orient_3d` for the
//! inside/outside decision.  For points in the shared flat plane, every query
//! returns `Degenerate` (signed tet volume = 0), so all points appear "inside"
//! every half-space — making the 3-D clipper useless for this case.
//!
//! The canonical `clip::polygon2d` backend uses the 2-D cross product sign to
//! classify query points (exact for planar geometry) and a CDT-backed
//! constrained arrangement for robust polygon Boolean reconstruction.
//! Specialized convex helpers are retained for triangle-local
//! clipping/decomposition inside this coplanar pipeline.
//!
//! ## Output quality
//!
//! Boundary vertices are computed as exact 2-D edge-edge intersections and
//! lifted back to 3-D, then registered via `VertexPool::insert_or_weld` →
//! shared seam vertices across adjacent A-triangles produce a manifold,
//! smooth boundary curve with no staircase artefacts.

pub mod basis;
pub mod geometry2d;
pub mod operations;

pub(crate) use basis::{detect_flat_plane, PlaneBasis};
pub(crate) use operations::boolean_coplanar;
