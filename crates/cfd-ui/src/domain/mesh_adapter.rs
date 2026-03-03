//! Mesh adapter — converts `IndexedMesh<f64>` to GPU-ready vertex/index arrays.
//!
//! # Theorem — Bijection of Index Mapping
//!
//! The adapter maps each `VertexId` to a contiguous `u32` index via `iter()`
//! enumeration order, which is the insertion order of the `VertexPool`. Because
//! `VertexPool` assigns monotonically increasing `VertexId` values and never
//! reuses IDs, the mapping `VertexId(n) -> n` is a bijection for all vertices
//! in the pool.  ∎

use cfd_mesh::domain::core::index::VertexId;
use cfd_mesh::IndexedMesh;

/// A GPU-ready vertex with position, normal, and region tag.
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct MeshVertex {
    /// World-space position.
    pub position: [f32; 3],
    /// Surface normal (may be zero if not computed).
    pub normal: [f32; 3],
    /// Region identifier for color-coding boundary patches.
    pub region_id: u32,
    /// Per-vertex scalar field value (simulation results, 0.0 if unused).
    pub field_value: f32,
}

/// Result of converting an `IndexedMesh` to GPU-ready arrays.
pub struct GpuMeshData {
    /// Interleaved vertex data.
    pub vertices: Vec<MeshVertex>,
    /// Triangle index data (3 indices per face).
    pub indices: Vec<u32>,
}

/// Convert an `IndexedMesh<f64>` into GPU-ready vertex and index arrays.
///
/// Vertices are extracted in pool iteration order. Face indices are remapped
/// from `VertexId` to contiguous `u32` positions in the resulting vertex array.
pub fn convert_mesh(mesh: &IndexedMesh<f64>) -> GpuMeshData {
    let mut vertices = Vec::with_capacity(mesh.vertex_count());
    let mut id_to_index = Vec::new();

    // Build the vertex array and ID-to-index lookup.
    // VertexId(n) uses n as inner u32; the pool iterates in insertion order.
    for (vid, vdata) in mesh.vertices.iter() {
        let idx = vertices.len() as u32;
        let needed = vid.0 as usize + 1;
        if id_to_index.len() < needed {
            id_to_index.resize(needed, 0u32);
        }
        id_to_index[vid.0 as usize] = idx;

        vertices.push(MeshVertex {
            position: [
                vdata.position.x as f32,
                vdata.position.y as f32,
                vdata.position.z as f32,
            ],
            normal: [
                vdata.normal.x as f32,
                vdata.normal.y as f32,
                vdata.normal.z as f32,
            ],
            region_id: 0,
            field_value: 0.0,
        });
    }

    // Build the index array, assigning region IDs from faces.
    let mut indices = Vec::with_capacity(mesh.face_count() * 3);
    for face in mesh.faces.iter() {
        let [v0, v1, v2] = face.vertices;
        let i0 = lookup_index(&id_to_index, v0);
        let i1 = lookup_index(&id_to_index, v1);
        let i2 = lookup_index(&id_to_index, v2);

        // Tag each vertex with the region from its face.
        let region = face.region.0;
        vertices[i0 as usize].region_id = region;
        vertices[i1 as usize].region_id = region;
        vertices[i2 as usize].region_id = region;

        indices.push(i0);
        indices.push(i1);
        indices.push(i2);
    }

    GpuMeshData { vertices, indices }
}

/// Look up contiguous index for a `VertexId`.
fn lookup_index(id_to_index: &[u32], vid: VertexId) -> u32 {
    id_to_index
        .get(vid.0 as usize)
        .copied()
        .unwrap_or(0)
}
