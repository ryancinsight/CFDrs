//! # cfd-mesh
//!
//! State-of-the-art watertight CFD mesh generation for millifluidic devices.
//!
//! This crate provides three mesh representations (one modern, two legacy),
//! a complete Boolean CSG pipeline, spatial-hash vertex welding, exact
//! geometric predicates, manifold/watertight checking, and OpenFOAM-compatible
//! I/O — all targeting millimetre-scale microfluidic channel geometries.
//!
//! ## Quick Start
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
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
// Mesh geometry code routinely casts index→float for vertex positions.
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_possible_wrap)]
// Mathematical variable names (i,j,k,x,y,z,n,t,u,v) are standard.
#![allow(clippy::similar_names)]
#![allow(clippy::many_single_char_names)]
// Mesh/geometry functions often need many geometric parameters.
#![allow(clippy::too_many_arguments)]
#![allow(clippy::too_many_lines)]
// Numeric closures improve readability in computational pipelines.
#![allow(clippy::redundant_closure_for_method_calls)]
// Internal helper definitions placed near their call site.
#![allow(clippy::items_after_statements)]
// Domain-specific naming matches geometric convention.
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::doc_markdown)]
// Error/panic doc sections deferred for internal mesh operations.
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]
// Mesh builder patterns return Self.
#![allow(clippy::return_self_not_must_use)]
#![allow(clippy::must_use_candidate)]
// Exact floating-point comparisons used in topology classification.
#![allow(clippy::float_cmp)]
#![allow(clippy::needless_pass_by_value)]
#![allow(clippy::unreadable_literal)]
#![allow(clippy::manual_let_else)]
#![allow(clippy::unnecessary_wraps)]
#![allow(clippy::match_same_arms)]
#![allow(clippy::useless_conversion)]
#![allow(clippy::inline_always)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::implicit_hasher)]
#![allow(clippy::new_without_default)]
#![allow(clippy::trivially_copy_pass_by_ref)]
#![allow(clippy::ptr_arg)]
#![allow(clippy::format_push_string)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::empty_line_after_doc_comments)]
#![allow(clippy::len_without_is_empty)]
#![allow(clippy::new_ret_no_self)]
#![allow(clippy::unused_self)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::manual_clamp)]
#![allow(clippy::duplicated_attributes)]
#![allow(clippy::same_item_push)]
#![allow(clippy::doc_overindented_list_items)]

pub mod application;
pub mod domain;
pub mod infrastructure;

/// Legacy watertight-first indexed surface mesh.
pub use domain::mesh::IndexedMesh;

/// Ergonomic builder for `IndexedMesh`.
pub use domain::mesh::MeshBuilder;

// ── Convenience re-exports ────────────────────────────────────────────────────

/// Normal-orientation analysis report for `IndexedMesh` surfaces.
pub use application::quality::{analyze_normals, NormalAnalysis};

/// Named CFD boundary patch (Inlet / Outlet / Wall / Symmetry / Periodic).
pub use domain::topology::halfedge::BoundaryPatch;

/// CFD boundary patch type discriminant.
pub use domain::topology::halfedge::PatchType;

/// Exact Shewchuk orientation result.
pub use domain::geometry::Orientation;

/// Analytic mesh primitives (26 builders from tetrahedron to truncated icosahedron).
pub use domain::geometry::primitives;

/// Primitive builder re-exports for ergonomic top-level access.
pub use domain::geometry::{
    Antiprism, BiconcaveDisk, Capsule, Cone, Cube, Cuboctahedron, Cylinder, Disk, Dodecahedron,
    Elbow, Ellipsoid, Frustum, GeodesicSphere, GyroidSphere, HelixSweep, Icosahedron, LinearSweep,
    Octahedron, Pipe, Pyramid, RevolutionSweep, RoundedCube, SchwarzDSphere, SchwarzPSphere,
    SerpentineTube, SphericalShell, StadiumPrism, Tetrahedron, Torus, TruncatedIcosahedron, UvSphere,
    // TPMS expansion: Neovius, Lidinoid, I-WP, Split P, FRD, Fischer-Koch C(Y)
    FischerKochCySphere, FrdSphere, IwpSphere, LidinoidSphere, NeoviusSphere, SplitPSphere,
};

/// Application-level channel builders.
pub use application::channel::{
    ChannelPath, JunctionType, ChannelProfile, SubstrateBuilder, SweepMesher,
    BranchingMeshBuilder, SerpentineMeshBuilder, VenturiMeshBuilder
};
