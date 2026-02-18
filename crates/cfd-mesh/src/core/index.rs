//! # Mesh Element Keys
//!
//! Two independent index systems coexist:
//!
//! ## New: SlotMap generational keys (`*Key` types)
//!
//! Used by [`crate::mesh::HalfEdgeMesh`] and [`crate::mesh::with_mesh`].
//! Each key is an 8-byte generational index; stale keys (after element
//! deletion) return `None` rather than silently aliasing a new element.
//!
//! ```text
//! SlotMap<VertexKey, GhostCell<'id, VertexData>>
//!   ├── VertexKey(gen=1, idx=0)  →  { pos: (0,0,0) }
//!   └── VertexKey(gen=2, idx=2)  →  { pos: (0,1,0) }   (slot recycled once)
//! ```
//!
//! ## Legacy: u32 newtype indices (`*Id` types)
//!
//! Used by [`crate::storage::VertexPool`], [`crate::storage::FaceStore`], and
//! [`crate::storage::EdgeStore`] — all of which are contiguous `Vec`-backed
//! stores indexed by a plain `u32`.  These are the types used by
//! [`crate::mesh::IndexedMesh`].
//!
//! # Key Invariant
//! `*Key` values are only valid within the `HalfEdgeMesh<'id>` that created
//! them.  Mixing keys from different meshes is a logic error (returns `None`,
//! never undefined behaviour).

use std::fmt;
use slotmap::new_key_type;

// ── SlotMap generational keys (HalfEdgeMesh) ─────────────────────────────────

new_key_type! {
    /// Key for a vertex in a [`crate::mesh::HalfEdgeMesh`].
    ///
    /// # Invariant
    /// Valid only within the mesh that produced it.
    /// Stale keys (after vertex deletion) return `None` on lookup.
    pub struct VertexKey;
}

new_key_type! {
    /// Key for a directed half-edge in a [`crate::mesh::HalfEdgeMesh`].
    ///
    /// # Invariant
    /// `twin(twin(he)) == he` — the twin relationship is an involution.
    /// See [`crate::topology::halfedge`] for the full invariant set.
    pub struct HalfEdgeKey;
}

new_key_type! {
    /// Key for a face (polygon) in a [`crate::mesh::HalfEdgeMesh`].
    pub struct FaceKey;
}

new_key_type! {
    /// Key for a named CFD boundary patch in a [`crate::mesh::HalfEdgeMesh`].
    pub struct PatchKey;
}

new_key_type! {
    /// Key for a mesh region (channel segment, junction) in a
    /// [`crate::mesh::HalfEdgeMesh`].
    pub struct RegionKey;
}

// ── Legacy u32 index types (IndexedMesh / VertexPool / FaceStore / EdgeStore) ─

/// Strongly-typed vertex index for [`crate::storage::VertexPool`].
///
/// A plain `u32` index into a contiguous vertex array.  Use [`VertexKey`] for
/// the new [`crate::mesh::HalfEdgeMesh`]-based API.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VertexId(pub u32);

impl VertexId {
    /// Create from a raw `u32` index.
    #[inline] pub fn new(raw: u32) -> Self { VertexId(raw) }
    /// Create from a `usize` index (panics on overflow in debug builds).
    #[inline] pub fn from_usize(n: usize) -> Self { VertexId(n as u32) }
    /// Return the raw `u32` index.
    #[inline] pub fn raw(self) -> u32 { self.0 }
    /// Return as `usize`.
    #[inline] pub fn as_usize(self) -> usize { self.0 as usize }
}

impl From<usize> for VertexId {
    #[inline] fn from(n: usize) -> Self { VertexId(n as u32) }
}

impl fmt::Display for VertexId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}", self.0) }
}

/// Strongly-typed face index for [`crate::storage::FaceStore`].
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct FaceId(pub u32);

impl FaceId {
    /// Create from a `usize` index.
    #[inline] pub fn from_usize(n: usize) -> Self { FaceId(n as u32) }
    /// Return as `usize`.
    #[inline] pub fn as_usize(self) -> usize { self.0 as usize }
}

impl From<usize> for FaceId {
    #[inline] fn from(n: usize) -> Self { FaceId(n as u32) }
}

impl fmt::Display for FaceId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}", self.0) }
}

/// Strongly-typed edge index for [`crate::storage::EdgeStore`].
///
/// Indexes the flattened edge list in `EdgeStore`; distinct from
/// [`HalfEdgeKey`] which indexes the directed half-edges of `HalfEdgeMesh`.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct EdgeId(pub u32);

impl EdgeId {
    /// Create from a `usize` index.
    #[inline] pub fn from_usize(n: usize) -> Self { EdgeId(n as u32) }
    /// Return as `usize`.
    #[inline] pub fn as_usize(self) -> usize { self.0 as usize }
}

impl From<usize> for EdgeId {
    #[inline] fn from(n: usize) -> Self { EdgeId(n as u32) }
}

impl fmt::Display for EdgeId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}", self.0) }
}

/// Strongly-typed region index for [`crate::storage::FaceStore`].
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct RegionId(pub u32);

impl RegionId {
    /// Sentinel value indicating "no region assigned".
    pub const INVALID: Self = Self(u32::MAX);
    /// Create from a `usize` index.
    #[inline] pub fn from_usize(n: usize) -> Self { RegionId(n as u32) }
    /// Return as `usize`.
    #[inline] pub fn as_usize(self) -> usize { self.0 as usize }
}

impl From<usize> for RegionId {
    #[inline] fn from(n: usize) -> Self { RegionId(n as u32) }
}

impl fmt::Display for RegionId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}", self.0) }
}

// ── Legacy alias for HalfEdgeId ───────────────────────────────────────────────

/// Legacy alias for directed half-edge within a half-edge mesh.
/// In the new `HalfEdgeMesh` system this maps to [`HalfEdgeKey`].
pub type HalfEdgeId = HalfEdgeKey;

#[cfg(test)]
mod tests {
    use super::*;
    use slotmap::SlotMap;

    #[test]
    fn vertex_key_generational_safety() {
        let mut map: SlotMap<VertexKey, u32> = SlotMap::with_key();
        let k = map.insert(42);
        assert_eq!(map[k], 42);
        map.remove(k);
        // Stale key returns None — no panic, no silent aliasing
        assert!(map.get(k).is_none());
    }

    #[test]
    fn vertex_id_roundtrip() {
        let id = VertexId::new(7);
        assert_eq!(id.raw(), 7);
        assert_eq!(id.as_usize(), 7);
        assert_eq!(VertexId::from_usize(42).as_usize(), 42);
    }

    #[test]
    fn region_id_invalid_sentinel() {
        assert_ne!(RegionId::INVALID, RegionId::from_usize(0));
        assert_eq!(RegionId::INVALID.0, u32::MAX);
    }
}
