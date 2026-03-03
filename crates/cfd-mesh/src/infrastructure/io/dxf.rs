//! DXF export for CNC fabrication workflows.
//!
//! Exports mesh faces as `3DFACE` entities in DXF R12 format. No external
//! dependencies are required — the DXF ASCII format is simple enough to
//! emit directly.

use std::io::Write;

use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::mesh::IndexedMesh;

/// Write an [`IndexedMesh`] as a DXF file with `3DFACE` entities.
///
/// Uses the minimal DXF R12 ASCII format, compatible with AutoCAD, LibreCAD,
/// and most CNC software.
pub fn write_dxf<W: Write>(writer: &mut W, mesh: &IndexedMesh) -> MeshResult<()> {
    // HEADER section (minimal)
    w(writer, "0")?;
    w(writer, "SECTION")?;
    w(writer, "2")?;
    w(writer, "HEADER")?;
    w(writer, "0")?;
    w(writer, "ENDSEC")?;

    // ENTITIES section
    w(writer, "0")?;
    w(writer, "SECTION")?;
    w(writer, "2")?;
    w(writer, "ENTITIES")?;

    for (_fid, face) in mesh.faces.iter_enumerated() {
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);

        w(writer, "0")?;
        w(writer, "3DFACE")?;
        w(writer, "8")?;
        w(writer, "0")?; // Layer 0

        // First vertex (group codes 10, 20, 30)
        write_point(writer, 10, 20, 30, a.x, a.y, a.z)?;
        // Second vertex (group codes 11, 21, 31)
        write_point(writer, 11, 21, 31, b.x, b.y, b.z)?;
        // Third vertex (group codes 12, 22, 32)
        write_point(writer, 12, 22, 32, c.x, c.y, c.z)?;
        // Fourth vertex (same as third for triangles)
        write_point(writer, 13, 23, 33, c.x, c.y, c.z)?;
    }

    w(writer, "0")?;
    w(writer, "ENDSEC")?;

    // EOF
    w(writer, "0")?;
    w(writer, "EOF")?;

    Ok(())
}

fn w<W: Write>(writer: &mut W, s: &str) -> MeshResult<()> {
    writeln!(writer, "{s}").map_err(MeshError::Io)
}

fn write_point<W: Write>(
    writer: &mut W,
    gx: u16,
    gy: u16,
    gz: u16,
    x: f64,
    y: f64,
    z: f64,
) -> MeshResult<()> {
    writeln!(writer, "{gx}").map_err(MeshError::Io)?;
    writeln!(writer, "{x}").map_err(MeshError::Io)?;
    writeln!(writer, "{gy}").map_err(MeshError::Io)?;
    writeln!(writer, "{y}").map_err(MeshError::Io)?;
    writeln!(writer, "{gz}").map_err(MeshError::Io)?;
    writeln!(writer, "{z}").map_err(MeshError::Io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;

    #[test]
    fn dxf_contains_3dface() {
        let mut mesh = IndexedMesh::new();
        let v0 = mesh.add_vertex_pos(Point3r::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3r::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3r::new(0.0, 1.0, 0.0));
        mesh.add_face(v0, v1, v2);

        let mut buf = Vec::new();
        write_dxf(&mut buf, &mesh).unwrap();

        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("3DFACE"));
        assert!(text.contains("EOF"));
    }
}
