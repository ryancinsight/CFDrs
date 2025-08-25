//! CSG (Constructive Solid Geometry) operations
//!
//! This module provides comprehensive CSG mesh generation and boolean operations
//! using the csgrs crate, enabling creation of complex geometries through
//! constructive solid geometry techniques.

// mod primitives; // TODO: Implement when needed
// mod operations; // TODO: Implement when needed
// mod geometry; // TODO: Implement when needed
mod errors;
mod transform;

pub use primitives::{Primitives, PrimitiveBuilder};
pub use operations::{BooleanOperator, BooleanOperation};
pub use geometry::{CsgGeometry, GeometryBuilder};
pub use errors::CsgError;
pub use transform::{Transform, TransformBuilder};

// Re-export the main operator for backward compatibility
pub use operations::CsgOperator;