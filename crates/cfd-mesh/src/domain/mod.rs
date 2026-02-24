//! Domain layer — pure value objects, domain logic, no infrastructure dependencies.
//!
//! This layer contains the foundational types and algorithms that express the
//! mathematical domain of computational mesh geometry.  It has zero dependency
//! on `infrastructure/` (I/O, storage backends, permission systems).
//!
//! ## Sub-layers
//!
//! | Module | Responsibility |
//! |--------|---------------|
//! | [`core`] | Scalar types, indices, error taxonomy, constants |
//! | [`geometry`] | Exact predicates, AABB, plane, NURBS, 26 primitive builders |
//! | [`topology`] | Half-edge kernel, manifold predicates, adjacency, orientation |
//! | [`mesh`] | `IndexedMesh`, `HalfEdgeMesh`, `MeshBuilder`, legacy `Mesh<T>` |
//! | [`grid`] | Structured volume grid |

pub mod core;
pub mod geometry;
pub mod grid;
pub mod mesh;
pub mod topology;
