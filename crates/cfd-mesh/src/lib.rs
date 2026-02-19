//! # cfd-mesh
//!
//! State-of-the-art watertight CFD mesh generation for millifluidic devices.
//!
//! This crate provides three mesh representations (one modern, two legacy),
//! a complete Boolean CSG pipeline, spatial-hash vertex welding, exact
//! geometric predicates, manifold/watertight checking, and OpenFOAM-compatible
//! I/O — all targeting millimetre-scale microfluidic channel geometries.
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
//! assert!(vol > 0.0, "outward-oriented tetrahedron has positive volume");
//! ```
//!
//! ## Quick Start — IndexedMesh builder
//!
//! ```rust,ignore
//! use cfd_mesh::{MeshBuilder, core::scalar::Point3r};
//!
//! let mesh = MeshBuilder::new()
//!     .add_triangle_vertex_positions(/* ... */)
//!     .build();
//! assert!(mesh.is_watertight());
//! ```
//!
//! ## Architecture Diagram
//!
//! ```text
//! ┌─ cfd-mesh crate ─────────────────────────────────────────────────────────┐
//! │                                                                           │
//! │  Entry points                                                             │
//! │  ┌──────────────────┐    ┌───────────────────┐    ┌──────────────────┐   │
//! │  │  with_mesh(f)    │    │   MeshBuilder     │    │  csg_boolean_*   │   │
//! │  │  HalfEdgeMesh    │    │   IndexedMesh     │    │  CsgNode tree    │   │
//! │  └────────┬─────────┘    └────────┬──────────┘    └────────┬─────────┘   │
//! │           │                       │                         │             │
//! │  ┌────────▼─────────┐    ┌────────▼──────────┐   ┌────────▼──────────┐  │
//! │  │ permission/      │    │ storage/           │   │ csg/              │  │
//! │  │  GhostToken      │    │  VertexPool        │   │  BspTree          │  │
//! │  │  GhostCell       │    │  FaceStore         │   │  BvhTree          │  │
//! │  └────────┬─────────┘    │  EdgeStore         │   │  boolean pipeline │  │
//! │           │              └────────┬──────────┘   └───────────────────┘  │
//! │  ┌────────▼─────────┐            │                                       │
//! │  │ topology/        │    ┌────────▼──────────┐                           │
//! │  │  halfedge kernel │    │ geometry/          │                           │
//! │  │  BoundaryPatch   │    │  exact predicates  │                           │
//! │  │  ElementType     │    │  AABB, Plane, NURBS│                           │
//! │  └──────────────────┘    └───────────────────┘                           │
//! │                                                                           │
//! │  Cross-cutting: welding/  watertight/  quality/  io/  core/              │
//! └───────────────────────────────────────────────────────────────────────────┘
//! ```
//!
//! ## Module Overview
//!
//! | Module | Contents |
//! |--------|---------|
//! | [`mesh`] | `HalfEdgeMesh`, `Mesh<T>`, `IndexedMesh`, `MeshBuilder` |
//! | [`topology`] | Half-edge structures, boundary patches, element types |
//! | [`geometry`] | Exact predicates, AABB, plane, NURBS, builders |
//! | [`welding`] | 26-neighbor `SnappingGrid`, `SpatialHashGrid`, `MeshWelder` |
//! | [`storage`] | `VertexPool`, `FaceStore`, `EdgeStore`, `SlotPool` |
//! | [`watertight`] | Manifold check, Euler characteristic, repair |
//! | [`quality`] | Triangle quality metrics and validation reports |
//! | [`permission`] | `GhostToken`, `GhostCell`, `PermissionedArena` |
//! | [`core`] | Scalar types, indices (`VertexKey`, `VertexId`, …), errors |
//! | [`io`] | STL and VTK mesh I/O |
//! | [`csg`] | BSP-tree + BVH boolean operations |
//!
//! ## Invariants
//!
//! The following mesh invariants are enforced at all API boundaries:
//!
//! 1. **Manifold half-edge**: every interior edge is shared by exactly 2 faces;
//!    `twin(twin(he)) == he` and `next(prev(he)) == he`.
//! 2. **Spatial deduplication**: `VertexPool` and `SnappingGrid` guarantee that
//!    no two vertex positions are closer than `TOLERANCE` apart.
//! 3. **Watertight closure**: `IndexedMesh::is_watertight()` verifies zero
//!    boundary edges (every edge has exactly 2 adjacent faces).
//! 4. **Generational key safety**: `VertexKey` / `FaceKey` values are valid
//!    only within the mesh that created them; stale keys return `None`.

// Phase 9: Rustdoc enforcement. All public items in user-facing modules must
// have documentation. Internal modules (channel, grid, hierarchy, io internals,
// permission internals) suppress the lint with targeted #[allow(missing_docs)].
#![warn(missing_docs)]


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

pub mod csg;

// ── Primary re-exports ────────────────────────────────────────────────────────

/// Create an empty `HalfEdgeMesh` + `GhostToken` inside a branded closure.
pub use mesh::with_mesh;

/// The GhostCell-permissioned half-edge mesh (new API).
pub use mesh::HalfEdgeMesh;

/// Legacy generic FEM/FVM mesh — retained only for volume tools (`grid.rs`, `hex_to_tet.rs`).
#[allow(deprecated)]
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

/// Analytic mesh primitives (26 builders from tetrahedron to truncated icosahedron).
pub use geometry::primitives;

/// Primitive builder re-exports for ergonomic top-level access.
pub use geometry::{
    Tetrahedron, Cube, UvSphere, Cylinder, Cone, Torus, LinearSweep, RevolutionSweep,
    Octahedron, Icosahedron, Ellipsoid, Frustum, Capsule, Pipe, Elbow,
    BiconcaveDisk, SphericalShell, StadiumPrism, Dodecahedron, GeodesicSphere,
    HelixSweep, RoundedCube, Cuboctahedron, Pyramid, Antiprism, TruncatedIcosahedron,
    Disk,
};
