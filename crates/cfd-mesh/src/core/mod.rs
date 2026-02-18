//! Core types: scalar abstraction, indices, errors, constants.
//!
//! Single source of truth for all fundamental types used across the crate.

pub mod scalar;
pub mod index;
pub mod error;
pub mod constants;

pub use scalar::Real;
pub use index::{VertexId, FaceId, EdgeId, HalfEdgeId, RegionId};
pub use error::{MeshError, MeshResult};
