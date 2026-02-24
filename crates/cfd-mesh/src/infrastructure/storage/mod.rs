//! Indexed data storage for mesh elements.
//!
//! Each element kind has a dedicated store that manages allocation, indexing,
//! and optional spatial lookup. Vertex deduplication uses a spatial hash grid
//! instead of the O(n²) linear scan in csgrs.
//!
//! ## Storage tiers
//!
//! | Module | Backing | Used by |
//! |---|---|---|
//! | [`pool`] | `Vec<T>` | legacy generic pool |
//! | [`vertex_pool`] | `Vec<T>` + spatial hash | `IndexedMesh` |
//! | [`face_store`] | `Vec<FaceData>` | `IndexedMesh` |
//! | [`edge_store`] | `Vec<EdgeData>` + hash map | `IndexedMesh` |
//! | [`slotmap_pool`] | `SlotMap<K, V>` | `HalfEdgeMesh<'id>` |

pub mod attribute;
pub mod edge_store;
pub mod face_store;
pub mod pool;
pub mod slotmap_pool;
pub mod vertex_pool;

pub use attribute::AttributeStore;
pub use edge_store::EdgeStore;
pub use face_store::FaceStore;
pub use pool::Pool;
pub use slotmap_pool::{GhostSlotPool, SlotPool};
pub use vertex_pool::VertexPool;
