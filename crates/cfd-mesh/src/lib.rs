//! # cfd-mesh
//!
//! State-of-the-art watertight CFD mesh generation for millifluidic devices.
//!
//! ## Quick Start — new GhostCell half-edge API
//!
//! ```rust,ignore
//! use cfd_mesh::{with_mesh, welding::SnappingGrid};
//!
//! let vol = with_mesh(|mut mesh, mut token| {
//!     let a = mesh.add_vertex([0.0, 0.0, 0.0], &token);
//!     let b = mesh.add_vertex([1.0, 0.0, 0.0], &token);
//!     let c = mesh.add_vertex([0.0, 1.0, 0.0], &token);
//!     let d = mesh.add_vertex([0.0, 0.0, 1.0], &token);
//!     mesh.add_triangle(a, b, c, &mut token).unwrap();
//!     mesh.add_triangle(a, c, d, &mut token).unwrap();
//!     mesh.add_triangle(a, d, b, &mut token).unwrap();
//!     mesh.add_triangle(b, d, c, &mut token).unwrap();
//!     mesh.signed_volume(&token)
//! });
//! ```
//!
//! ## Module Overview
//!
//! | Module | Contents |
//! |---|---|
//! | [`mesh`] | `HalfEdgeMesh`, `Mesh<T>`, `IndexedMesh`, `MeshBuilder` |
//! | [`topology`] | Half-edge structures, boundary patches, element types |
//! | [`geometry`] | Exact predicates, AABB, plane, builders |
//! | [`welding`] | 26-neighbor `SnappingGrid`, `SpatialHashGrid`, `MeshWelder` |
//! | [`storage`] | `VertexPool`, `FaceStore`, `EdgeStore`, `SlotPool` |
//! | [`watertight`] | Manifold check, Euler characteristic, repair |
//! | [`permission`] | `GhostToken`, `GhostCell`, `PermissionedArena` |
//! | [`core`] | Scalar types, indices (`VertexKey`, `VertexId`, …), errors |
//! | [`csg`] | BSP-tree + BVH boolean operations (feature-gated) |

#![allow(unused)]
#![allow(missing_docs)]

pub mod core;
pub mod permission;
pub mod storage;
pub mod topology;
pub mod geometry;
pub mod hierarchy;
pub mod grid;
pub mod mesh;
pub mod welding;
pub mod quality;
pub mod watertight;
pub mod channel;
pub mod io;

#[cfg(feature = "csg")]
pub mod csg;

// ── Primary re-exports ────────────────────────────────────────────────────────

/// Create an empty `HalfEdgeMesh` + `GhostToken` inside a branded closure.
pub use mesh::with_mesh;

/// The GhostCell-permissioned half-edge mesh (new API).
pub use mesh::HalfEdgeMesh;

/// Legacy generic FEM/FVM mesh.
pub use mesh::Mesh;

/// Legacy watertight-first indexed surface mesh.
pub use mesh::IndexedMesh;

/// Ergonomic builder for `IndexedMesh`.
pub use mesh::MeshBuilder;

// ── Convenience re-exports ────────────────────────────────────────────────────

/// Named CFD boundary patch (Inlet / Outlet / Wall / Symmetry / Periodic).
pub use topology::halfedge::BoundaryPatch;

/// CFD boundary patch type discriminant.
pub use topology::halfedge::PatchType;

/// Exact Shewchuk orientation result.
pub use geometry::Orientation;
