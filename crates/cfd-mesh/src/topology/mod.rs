//! Mesh topology components: vertices, edges, faces, and cells

mod cell;
mod edge;
mod face;
mod vertex;

pub use cell::Cell;
pub use edge::Edge;
pub use face::Face;
pub use vertex::Vertex;

// Re-export ElementType for convenience
pub use cfd_core::domains::mesh_operations::ElementType;
