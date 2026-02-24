//! Core types: scalar abstraction, indices, errors, constants.
//!
//! Single source of truth for all fundamental types used across the crate.

pub mod constants;
pub mod error;
pub mod index;
pub mod scalar;

pub use error::{MeshError, MeshResult};
pub use index::{EdgeId, FaceId, HalfEdgeId, RegionId, VertexId};
pub use scalar::{Real, Scalar};
