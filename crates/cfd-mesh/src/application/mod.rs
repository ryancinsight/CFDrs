//! Application layer — use-case services that compose domain types with
//! infrastructure to deliver CFD mesh operations.
//!
//! Each module here implements a single, cohesive use-case (SRP) and depends
//! only on `domain/` and `infrastructure/` — never the reverse.
//!
//! ## Sub-layers
//!
//! | Module | Responsibility |
//! |--------|---------------|
//! | [`csg`] | Complete 5-phase Boolean CSG pipeline (union, intersection, difference) |
//! | [`channel`] | Millifluidic channel geometry builders (sweep, junction, profile) |
//! | [`quality`] | Triangle quality metrics and mesh validation reports |
//! | [`watertight`] | Manifold checking, Euler characteristic, sealing, repair |
//! | [`welding`] | Vertex deduplication (spatial hash + union-find) |
//! | [`hierarchy`] | Mesh format conversions (hex ↔ tet) |

pub mod channel;
pub mod csg;
pub mod delaunay;
pub mod hierarchy;
pub mod quality;
pub mod watertight;
pub mod welding;
