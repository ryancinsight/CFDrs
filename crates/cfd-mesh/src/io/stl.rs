//! STL import and export.
//!
//! Supports both ASCII and binary STL formats.

use std::io::{Read, Write, BufRead, BufReader};

use crate::core::index::RegionId;
use crate::core::scalar::{Real, Point3r, Vector3r};
use crate::core::error::{MeshError, MeshResult};
use crate::storage::face_store::{FaceData, FaceStore};
use crate::storage::vertex_pool::VertexPool;

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

        let normal = crate::geometry::normal::triangle_normal(&a, &b, &c)
            .unwrap_or_else(|| Vector3r::z());

        writeln!(writer, "  facet normal {} {} {}", normal.x, normal.y, normal.z)
            .map_err(MeshError::Io)?;
        writeln!(writer, "    outer loop").map_err(MeshError::Io)?;
        for p in [&a, &b, &c] {
            writeln!(writer, "      vertex {} {} {}", p.x, p.y, p.z)
                .map_err(MeshError::Io)?;
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
    writer
        .write_all(&header)
        .map_err(MeshError::Io)?;

    // Number of triangles
    let n_triangles = face_store.len() as u32;
    writer
        .write_all(&n_triangles.to_le_bytes())
        .map_err(MeshError::Io)?;

    for (_, face) in face_store.iter_enumerated() {
        let a = vertex_pool.position(face.vertices[0]);
        let b = vertex_pool.position(face.vertices[1]);
        let c = vertex_pool.position(face.vertices[2]);

        let normal = crate::geometry::normal::triangle_normal(&a, &b, &c)
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
            let normal = crate::geometry::normal::triangle_normal(&verts[0], &verts[1], &verts[2])
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
    w.write_all(&v.to_le_bytes())
        .map_err(MeshError::Io)
}
