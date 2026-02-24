//! STL import and export.
//!
//! Supports both ASCII and binary STL formats.

use std::io::{BufRead, BufReader, Read, Write};


use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::domain::mesh::IndexedMesh;
use crate::infrastructure::storage::face_store::{FaceData, FaceStore};
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Write an indexed mesh as ASCII STL.
pub fn write_ascii_stl<W: Write>(
    writer: &mut W,
    name: &str,
    vertex_pool: &VertexPool,
    face_store: &FaceStore,
) -> MeshResult<()> {
    writeln!(writer, "solid {name}").map_err(MeshError::Io)?;

    for (_, face) in face_store.iter_enumerated() {
        let a = vertex_pool.position(face.vertices[0]);
        let b = vertex_pool.position(face.vertices[1]);
        let c = vertex_pool.position(face.vertices[2]);

        let normal = crate::domain::geometry::normal::triangle_normal(&a, &b, &c)
            .unwrap_or_else(|| Vector3r::z());

        writeln!(
            writer,
            "  facet normal {} {} {}",
            normal.x, normal.y, normal.z
        )
        .map_err(MeshError::Io)?;
        writeln!(writer, "    outer loop").map_err(MeshError::Io)?;
        for p in [&a, &b, &c] {
            writeln!(writer, "      vertex {} {} {}", p.x, p.y, p.z).map_err(MeshError::Io)?;
        }
        writeln!(writer, "    endloop").map_err(MeshError::Io)?;
        writeln!(writer, "  endfacet").map_err(MeshError::Io)?;
    }

    writeln!(writer, "endsolid {name}").map_err(MeshError::Io)?;
    Ok(())
}

/// Write an indexed mesh as binary STL.
pub fn write_binary_stl<W: Write>(
    writer: &mut W,
    vertex_pool: &VertexPool,
    face_store: &FaceStore,
) -> MeshResult<()> {
    // 80-byte header
    let header = [0u8; 80];
    writer.write_all(&header).map_err(MeshError::Io)?;

    // Number of triangles
    let n_triangles = face_store.len() as u32;
    writer
        .write_all(&n_triangles.to_le_bytes())
        .map_err(MeshError::Io)?;

    for (_, face) in face_store.iter_enumerated() {
        let a = vertex_pool.position(face.vertices[0]);
        let b = vertex_pool.position(face.vertices[1]);
        let c = vertex_pool.position(face.vertices[2]);

        let normal = crate::domain::geometry::normal::triangle_normal(&a, &b, &c)
            .unwrap_or_else(|| Vector3r::z());

        // Normal (3 × f32)
        write_f32(writer, normal.x as f32)?;
        write_f32(writer, normal.y as f32)?;
        write_f32(writer, normal.z as f32)?;

        // Vertices (3 × 3 × f32)
        for p in [&a, &b, &c] {
            write_f32(writer, p.x as f32)?;
            write_f32(writer, p.y as f32)?;
            write_f32(writer, p.z as f32)?;
        }

        // Attribute byte count
        writer
            .write_all(&0u16.to_le_bytes())
            .map_err(MeshError::Io)?;
    }

    Ok(())
}

/// Read an ASCII STL into the vertex pool and face store.
pub fn read_ascii_stl<R: Read>(
    reader: R,
    vertex_pool: &mut VertexPool,
    face_store: &mut FaceStore,
    region: RegionId,
) -> MeshResult<usize> {
    let buf = BufReader::new(reader);
    let mut count = 0usize;
    let mut verts: Vec<Point3r> = Vec::with_capacity(3);

    for line in buf.lines() {
        let line = line.map_err(MeshError::Io)?;
        let trimmed = line.trim();

        if trimmed.starts_with("vertex") {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 4 {
                let x: Real = parts[1]
                    .parse()
                    .map_err(|_| MeshError::Other("bad vertex x".to_string()))?;
                let y: Real = parts[2]
                    .parse()
                    .map_err(|_| MeshError::Other("bad vertex y".to_string()))?;
                let z: Real = parts[3]
                    .parse()
                    .map_err(|_| MeshError::Other("bad vertex z".to_string()))?;
                verts.push(Point3r::new(x, y, z));
            }
        }

        if trimmed.starts_with("endfacet") && verts.len() == 3 {
            let normal =
                crate::domain::geometry::normal::triangle_normal(&verts[0], &verts[1], &verts[2])
                    .unwrap_or_else(|| Vector3r::z());

            let v0 = vertex_pool.insert_or_weld(verts[0], normal);
            let v1 = vertex_pool.insert_or_weld(verts[1], normal);
            let v2 = vertex_pool.insert_or_weld(verts[2], normal);

            face_store.push(FaceData {
                vertices: [v0, v1, v2],
                region,
            });

            count += 1;
            verts.clear();
        }
    }

    Ok(count)
}

/// Write a single f32 in little-endian.
fn write_f32<W: Write>(w: &mut W, v: f32) -> MeshResult<()> {
    w.write_all(&v.to_le_bytes()).map_err(MeshError::Io)
}

// =============================================================================
//  Low-level binary STL reader (to VertexPool + FaceStore)
// =============================================================================

/// Read a binary STL into the vertex pool and face store.
///
/// Binary STL format: 80-byte header, u32 triangle count, then for each
/// triangle: 12-byte normal, 3 × 12-byte vertices, 2-byte attribute count.
/// Vertex normals are recomputed from face geometry rather than read from the
/// file (the spec does not require them to be correct).
pub fn read_binary_stl<R: Read>(
    reader: R,
    vertex_pool: &mut VertexPool,
    face_store: &mut FaceStore,
    region: RegionId,
) -> MeshResult<usize> {
    let mut r = BufReader::new(reader);
    let mut header = [0u8; 80];
    r.read_exact(&mut header).map_err(MeshError::Io)?;
    let mut count_bytes = [0u8; 4];
    r.read_exact(&mut count_bytes).map_err(MeshError::Io)?;
    let n = u32::from_le_bytes(count_bytes) as usize;

    for _ in 0..n {
        // Skip the stored normal (12 bytes) — we recompute it.
        let mut skip = [0u8; 12];
        r.read_exact(&mut skip).map_err(MeshError::Io)?;

        let mut verts = [Point3r::new(0.0, 0.0, 0.0); 3];
        for vert in &mut verts {
            let mut vbuf = [0u8; 12];
            r.read_exact(&mut vbuf).map_err(MeshError::Io)?;
            let x = f32::from_le_bytes([vbuf[0], vbuf[1], vbuf[2], vbuf[3]]) as Real;
            let y = f32::from_le_bytes([vbuf[4], vbuf[5], vbuf[6], vbuf[7]]) as Real;
            let z = f32::from_le_bytes([vbuf[8], vbuf[9], vbuf[10], vbuf[11]]) as Real;
            *vert = Point3r::new(x, y, z);
        }
        // Skip attribute byte count (2 bytes).
        let mut attr = [0u8; 2];
        r.read_exact(&mut attr).map_err(MeshError::Io)?;

        let normal =
            crate::domain::geometry::normal::triangle_normal(&verts[0], &verts[1], &verts[2])
                .unwrap_or_else(|| Vector3r::z());
        let v0 = vertex_pool.insert_or_weld(verts[0], normal);
        let v1 = vertex_pool.insert_or_weld(verts[1], normal);
        let v2 = vertex_pool.insert_or_weld(verts[2], normal);
        face_store.push(FaceData {
            vertices: [v0, v1, v2],
            region,
        });
    }
    Ok(n)
}

// =============================================================================
//  High-level IndexedMesh helpers
// =============================================================================

/// Read an STL file (auto-detecting ASCII vs binary) into a new [`IndexedMesh`].
///
/// Detection is based on the binary record-size invariant:
/// `file_bytes == 84 + triangle_count * 50`.  Any file that satisfies this
/// is parsed as binary; everything else is attempted as ASCII.
/// All faces are tagged with `RegionId(0)` (wall).
pub fn read_stl<R: Read>(reader: R) -> MeshResult<IndexedMesh> {
    let mut data = Vec::new();
    // read_to_end needs the Read trait in scope — it is via `use std::io::Read`.
    BufReader::new(reader)
        .read_to_end(&mut data)
        .map_err(MeshError::Io)?;

    let region = RegionId::from_usize(0);
    let mut mesh = IndexedMesh::new();

    let is_binary = data.len() >= 84
        && data.len()
            == 84 + u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize * 50;

    if is_binary {
        read_binary_stl(
            std::io::Cursor::new(data),
            &mut mesh.vertices,
            &mut mesh.faces,
            region,
        )?;
    } else {
        read_ascii_stl(
            std::io::Cursor::new(data),
            &mut mesh.vertices,
            &mut mesh.faces,
            region,
        )?;
    }
    Ok(mesh)
}

/// Write an [`IndexedMesh`] as ASCII STL (convenience wrapper).
pub fn write_stl_ascii<W: Write>(writer: &mut W, name: &str, mesh: &IndexedMesh) -> MeshResult<()> {
    write_ascii_stl(writer, name, &mesh.vertices, &mesh.faces)
}

/// Write an [`IndexedMesh`] as binary STL (convenience wrapper).
pub fn write_stl_binary<W: Write>(writer: &mut W, mesh: &IndexedMesh) -> MeshResult<()> {
    write_binary_stl(writer, &mesh.vertices, &mesh.faces)
}

// =============================================================================
//  Fuzz entry point
// =============================================================================

/// Fuzz entry point for STL parsing.
///
/// Accepts arbitrary bytes and attempts to parse them as STL.  This function
/// must **never panic** — all errors are returned as `Err`.  Suitable as the
/// inner body of a `cargo-fuzz` target.
///
/// # Example (in a fuzz target)
/// ```rust,ignore
/// #![no_main]
/// libfuzzer_sys::fuzz_target!(|data: &[u8]| {
///     let _ = cfd_mesh::infrastructure::io::stl::fuzz_read_stl(data);
/// });
/// ```
pub fn fuzz_read_stl(data: &[u8]) -> MeshResult<IndexedMesh> {
    read_stl(std::io::Cursor::new(data))
}

// =============================================================================
//  Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;


    // ── ASCII round-trip ──────────────────────────────────────────────────

    #[test]
    fn ascii_stl_round_trip_indexed() {
        let mut mesh = IndexedMesh::new();
        let v0 = mesh.add_vertex_pos(Point3r::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3r::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3r::new(0.0, 1.0, 0.0));
        mesh.add_face(v0, v1, v2);

        let mut buf = Vec::new();
        write_stl_ascii(&mut buf, "test", &mesh).unwrap();

        let mesh2 = read_stl(std::io::Cursor::new(&buf)).unwrap();
        assert_eq!(mesh2.face_count(), 1);
        assert_eq!(mesh2.vertex_count(), 3);
    }

    // ── Binary round-trip ─────────────────────────────────────────────────

    #[test]
    fn binary_stl_round_trip_indexed() {
        let mut mesh = IndexedMesh::new();
        let v0 = mesh.add_vertex_pos(Point3r::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3r::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3r::new(0.0, 1.0, 0.0));
        mesh.add_face(v0, v1, v2);

        let mut buf = Vec::new();
        write_stl_binary(&mut buf, &mesh).unwrap();

        let mesh2 = read_stl(std::io::Cursor::new(&buf)).unwrap();
        assert_eq!(mesh2.face_count(), 1);
        assert_eq!(mesh2.vertex_count(), 3);
    }

    // ── Fuzz entry point never panics ─────────────────────────────────────

    #[test]
    fn fuzz_target_handles_empty_input() {
        let result = fuzz_read_stl(b"");
        // May succeed (empty mesh) or return an error — must not panic.
        let _ = result;
    }

    #[test]
    fn fuzz_target_handles_truncated_binary() {
        let result = fuzz_read_stl(&[0u8; 84]);
        let _ = result;
    }
}
