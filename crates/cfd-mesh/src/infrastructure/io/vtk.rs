//! VTK unstructured grid export.
//!
//! Writes a VTK legacy ASCII file suitable for ParaView visualization.

use std::io::Write;
use std::collections::HashMap;

use crate::core::scalar::Real;
use crate::core::error::{MeshError, MeshResult};
use crate::mesh::{HalfEdgeMesh, IndexedMesh};
use crate::permission::GhostToken;
use crate::storage::face_store::FaceStore;
use crate::storage::vertex_pool::VertexPool;

/// Write an indexed mesh as a VTK legacy ASCII unstructured grid.
pub fn write_vtk<W: Write>(
    writer: &mut W,
    vertex_pool: &VertexPool,
    face_store: &FaceStore,
) -> MeshResult<()> {
    let n_verts = vertex_pool.len();
    let n_faces = face_store.len();

    writeln!(writer, "# vtk DataFile Version 3.0").map_err(MeshError::Io)?;
    writeln!(writer, "cfd-mesh output").map_err(MeshError::Io)?;
    writeln!(writer, "ASCII").map_err(MeshError::Io)?;
    writeln!(writer, "DATASET UNSTRUCTURED_GRID").map_err(MeshError::Io)?;

    // Points
    writeln!(writer, "POINTS {n_verts} double").map_err(MeshError::Io)?;
    for i in 0..n_verts {
        let vid = crate::core::index::VertexId::new(i as u32);
        let p = vertex_pool.position(vid);
        writeln!(writer, "{} {} {}", p.x, p.y, p.z).map_err(MeshError::Io)?;
    }

    // Cells (triangles, VTK type 5)
    let cell_data_size = n_faces * 4; // 3 vertices + count per face
    writeln!(writer, "CELLS {n_faces} {cell_data_size}").map_err(MeshError::Io)?;
    for (_, face) in face_store.iter_enumerated() {
        writeln!(
            writer,
            "3 {} {} {}",
            face.vertices[0].raw(),
            face.vertices[1].raw(),
            face.vertices[2].raw()
        )
        .map_err(MeshError::Io)?;
    }

    // Cell types
    writeln!(writer, "CELL_TYPES {n_faces}").map_err(MeshError::Io)?;
    for _ in 0..n_faces {
        writeln!(writer, "5").map_err(MeshError::Io)?; // VTK_TRIANGLE
    }

    Ok(())
}

/// Write an [`IndexedMesh`] as a VTK legacy ASCII unstructured grid (convenience).
pub fn write_vtk_indexed<W: Write>(writer: &mut W, mesh: &IndexedMesh) -> MeshResult<()> {
    write_vtk(writer, &mesh.vertices, &mesh.faces)
}

/// Write a [`HalfEdgeMesh`] as a VTK legacy ASCII unstructured grid.
///
/// Vertex keys are mapped to sequential indices 0..n in iteration order.
/// Only triangular faces (those with exactly 3 vertices) are emitted.
/// The VTK cell type is always `5` (VTK_TRIANGLE).
pub fn write_vtk_he<'id, W: Write>(
    writer: &mut W,
    mesh: &HalfEdgeMesh<'id>,
    token: &GhostToken<'id>,
) -> MeshResult<()> {
    // Build a sequential index for every vertex key.
    let vertex_keys: Vec<_> = mesh.vertex_keys().collect();
    let n_verts = vertex_keys.len();
    let vertex_index: HashMap<_, usize> =
        vertex_keys.iter().enumerate().map(|(i, &k)| (k, i)).collect();

    let face_keys: Vec<_> = mesh.face_keys().collect();
    let n_faces = face_keys.len();

    writeln!(writer, "# vtk DataFile Version 3.0").map_err(MeshError::Io)?;
    writeln!(writer, "cfd-mesh output").map_err(MeshError::Io)?;
    writeln!(writer, "ASCII").map_err(MeshError::Io)?;
    writeln!(writer, "DATASET UNSTRUCTURED_GRID").map_err(MeshError::Io)?;

    // Points
    writeln!(writer, "POINTS {n_verts} double").map_err(MeshError::Io)?;
    for &vk in &vertex_keys {
        let p = mesh.vertex_pos(vk, token)
            .unwrap_or_else(|| crate::core::scalar::Point3r::new(0.0, 0.0, 0.0));
        writeln!(writer, "{} {} {}", p.x, p.y, p.z).map_err(MeshError::Io)?;
    }

    // Cells (connectivity list: count + indices per row)
    let cell_data_size = n_faces * 4;
    writeln!(writer, "CELLS {n_faces} {cell_data_size}").map_err(MeshError::Io)?;
    for &fk in &face_keys {
        let verts = mesh.face_vertices(fk, token);
        if let [v0, v1, v2] = verts.as_slice() {
            let i0 = vertex_index.get(v0).copied().unwrap_or(0);
            let i1 = vertex_index.get(v1).copied().unwrap_or(0);
            let i2 = vertex_index.get(v2).copied().unwrap_or(0);
            writeln!(writer, "3 {i0} {i1} {i2}").map_err(MeshError::Io)?;
        }
    }

    // Cell types
    writeln!(writer, "CELL_TYPES {n_faces}").map_err(MeshError::Io)?;
    for _ in 0..n_faces {
        writeln!(writer, "5").map_err(MeshError::Io)?; // VTK_TRIANGLE
    }

    Ok(())
}
