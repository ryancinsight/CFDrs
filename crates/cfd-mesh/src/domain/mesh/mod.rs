//! # Mesh Types
//!
//! This module provides the consolidated mesh type:
//!
//! - **[`IndexedMesh`]** — the single source of truth (SSOT) surface and volume mesh,
//!   combining `VertexPool` (spatial-hash dedup), `FaceStore`, `EdgeStore`,
//!   `AttributeStore`, `Vec<Cell>`, and `GhostCell` topology.

#[allow(missing_docs)]
pub mod indexed;
pub use indexed::{IndexedMesh, MeshBuilder};
