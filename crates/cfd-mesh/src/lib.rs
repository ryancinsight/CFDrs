//! # cfd-mesh
//!
//! CFD mesh generation for millifluidic devices.
//!
//! Provides:
//! - **`mesh`** — Generic `Mesh<T>` type for FEM/FVM solvers
//! - **`topology`** — `Cell`, `Face`, `Vertex<T>`, `ElementType`
//! - **`geometry`** — Mesh builders for Venturi, Serpentine, Branching geometries
//! - **`hierarchy`** — Mesh refinement (hex-to-tet conversion, P2 promotion)
//! - **`grid`** — Structured grid builder
//! - **`core`** — Scalar types, indices, error taxonomy
//! - **`csg`** — BSP-tree boolean operations (feature-gated)

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

/// Re-export the primary mesh types at crate root.
pub use mesh::Mesh;
pub use mesh::IndexedMesh;
pub use mesh::MeshBuilder;
