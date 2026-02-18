//! Indexed data storage for mesh elements.
//!
//! Each element kind has a dedicated store that manages allocation, indexing,
//! and optional spatial lookup. Vertex deduplication uses a spatial hash grid
//! instead of the O(n²) linear scan in csgrs.

pub mod pool;
pub mod vertex_pool;
pub mod face_store;
pub mod edge_store;
pub mod attribute;

pub use pool::Pool;
pub use vertex_pool::VertexPool;
pub use face_store::FaceStore;
pub use edge_store::EdgeStore;
pub use attribute::AttributeStore;
