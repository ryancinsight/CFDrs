//! VTK unstructured grid export.
//!
//! Writes a VTK legacy ASCII file suitable for ParaView visualization.

use std::io::Write;

use crate::core::scalar::Real;
use crate::core::error::{MeshError, MeshResult};
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
