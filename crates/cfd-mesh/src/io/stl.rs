//! STL import and export.
//!
//! Supports both ASCII and binary STL formats.

use std::io::{Read, Write, BufRead, BufReader};

use hashbrown::HashMap;

use crate::core::index::{RegionId, VertexKey};
use crate::core::scalar::{Real, Point3r, Vector3r};
use crate::core::error::{MeshError, MeshResult};
use crate::mesh::{HalfEdgeMesh, IndexedMesh};
use crate::permission::GhostToken;
use crate::storage::face_store::{FaceData, FaceStore};
use crate::storage::vertex_pool::VertexPool;
use crate::welding::snap::SnappingGrid;

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

        let normal = crate::geometry::normal::triangle_normal(&verts[0], &verts[1], &verts[2])
            .unwrap_or_else(|| Vector3r::z());
        let v0 = vertex_pool.insert_or_weld(verts[0], normal);
        let v1 = vertex_pool.insert_or_weld(verts[1], normal);
        let v2 = vertex_pool.insert_or_weld(verts[2], normal);
        face_store.push(FaceData { vertices: [v0, v1, v2], region });
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
    BufReader::new(reader).read_to_end(&mut data).map_err(MeshError::Io)?;

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
//  HalfEdgeMesh STL reader
// =============================================================================

/// Read an STL file directly into a [`HalfEdgeMesh`].
///
/// Vertices are deduplicated using a millifluidic [`SnappingGrid`] (ε = 1 μm)
/// before being added to the mesh, ensuring manifold topology.  Non-manifold
/// triangles that would violate the half-edge invariants (e.g. duplicate
/// directed edges) are silently skipped.
///
/// # Example
/// ```rust,ignore
/// use cfd_mesh::mesh::with_mesh;
/// use cfd_mesh::io::stl::read_stl_he;
///
/// let ascii = "solid s\n  facet normal 0 0 1\n    outer loop\n\
///              vertex 0 0 0\n      vertex 1 0 0\n      vertex 0 1 0\n\
///              endloop\n  endfacet\nendsolid s\n";
/// let faces = with_mesh(|mut mesh, mut token| {
///     read_stl_he(ascii.as_bytes(), &mut mesh, &mut token).unwrap();
///     mesh.face_count()
/// });
/// assert_eq!(faces, 1);
/// ```
pub fn read_stl_he<'id, R: Read>(
    reader: R,
    mesh: &mut HalfEdgeMesh<'id>,
    token: &mut GhostToken<'id>,
) -> MeshResult<usize> {
    let mut data = Vec::new();
    BufReader::new(reader).read_to_end(&mut data).map_err(MeshError::Io)?;

    let is_binary = data.len() >= 84
        && data.len()
            == 84 + u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize * 50;

    let mut snap = SnappingGrid::millifluidic();
    let mut key_map: HashMap<u32, VertexKey> = HashMap::new();

    if is_binary {
        read_binary_into_he(&data, mesh, token, &mut snap, &mut key_map)
    } else {
        read_ascii_into_he(&data, mesh, token, &mut snap, &mut key_map)
    }
}

fn read_binary_into_he<'id>(
    data: &[u8],
    mesh: &mut HalfEdgeMesh<'id>,
    token: &mut GhostToken<'id>,
    snap: &mut SnappingGrid,
    key_map: &mut HashMap<u32, VertexKey>,
) -> MeshResult<usize> {
    // Safety: caller guarantees data.len() == 84 + n*50
    let n = u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize;
    let mut count = 0usize;
    for i in 0..n {
        let base = 84 + i * 50;
        // Skip stored normal (bytes base..base+12), read three vertices.
        let mut verts = [Point3r::new(0.0, 0.0, 0.0); 3];
        for (j, vert) in verts.iter_mut().enumerate() {
            let off = base + 12 + j * 12;
            let x = f32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as Real;
            let y = f32::from_le_bytes(data[off + 4..off + 8].try_into().unwrap()) as Real;
            let z = f32::from_le_bytes(data[off + 8..off + 12].try_into().unwrap()) as Real;
            *vert = Point3r::new(x, y, z);
        }
        let vk = [
            snap.insert_or_weld_he(verts[0], key_map, mesh, token),
            snap.insert_or_weld_he(verts[1], key_map, mesh, token),
            snap.insert_or_weld_he(verts[2], key_map, mesh, token),
        ];
        // Skip degenerate triangles (snapped to fewer than 3 unique vertices).
        if vk[0] == vk[1] || vk[1] == vk[2] || vk[0] == vk[2] {
            continue;
        }
        if mesh.add_triangle(vk[0], vk[1], vk[2], token).is_ok() {
            count += 1;
        }
    }
    Ok(count)
}

fn read_ascii_into_he<'id>(
    data: &[u8],
    mesh: &mut HalfEdgeMesh<'id>,
    token: &mut GhostToken<'id>,
    snap: &mut SnappingGrid,
    key_map: &mut HashMap<u32, VertexKey>,
) -> MeshResult<usize> {
    let text = std::str::from_utf8(data)
        .map_err(|_| MeshError::Other("non-UTF-8 STL data".into()))?;
    let mut count = 0usize;
    let mut verts: Vec<Point3r> = Vec::with_capacity(3);

    for line in text.lines() {
        let t = line.trim();
        if t.starts_with("vertex ") {
            let p: Vec<&str> = t.split_whitespace().collect();
            if p.len() >= 4 {
                let x: Real = p[1].parse().map_err(|_| MeshError::Other("bad x".into()))?;
                let y: Real = p[2].parse().map_err(|_| MeshError::Other("bad y".into()))?;
                let z: Real = p[3].parse().map_err(|_| MeshError::Other("bad z".into()))?;
                verts.push(Point3r::new(x, y, z));
            }
        }
        if t.starts_with("endfacet") && verts.len() == 3 {
            let vk = [
                snap.insert_or_weld_he(verts[0], key_map, mesh, token),
                snap.insert_or_weld_he(verts[1], key_map, mesh, token),
                snap.insert_or_weld_he(verts[2], key_map, mesh, token),
            ];
            verts.clear();
            if vk[0] == vk[1] || vk[1] == vk[2] || vk[0] == vk[2] {
                continue;
            }
            if mesh.add_triangle(vk[0], vk[1], vk[2], token).is_ok() {
                count += 1;
            }
        }
    }
    Ok(count)
}

// =============================================================================
//  HalfEdgeMesh STL writers
// =============================================================================

/// Write a [`HalfEdgeMesh`] as ASCII STL.
///
/// Face normals are recomputed from vertex positions.
pub fn write_stl_ascii_he<'id, W: Write>(
    writer: &mut W,
    name: &str,
    mesh: &HalfEdgeMesh<'id>,
    token: &GhostToken<'id>,
) -> MeshResult<()> {
    writeln!(writer, "solid {name}").map_err(MeshError::Io)?;
    for fk in mesh.face_keys() {
        let verts = mesh.face_vertices(fk, token);
        if verts.len() != 3 {
            continue;
        }
        let pa = mesh.vertex_pos(verts[0], token).unwrap_or_else(|| Point3r::new(0.0, 0.0, 0.0));
        let pb = mesh.vertex_pos(verts[1], token).unwrap_or_else(|| Point3r::new(0.0, 0.0, 0.0));
        let pc = mesh.vertex_pos(verts[2], token).unwrap_or_else(|| Point3r::new(0.0, 0.0, 0.0));
        let n = crate::geometry::normal::triangle_normal(&pa, &pb, &pc)
            .unwrap_or_else(|| Vector3r::z());
        writeln!(writer, "  facet normal {} {} {}", n.x, n.y, n.z).map_err(MeshError::Io)?;
        writeln!(writer, "    outer loop").map_err(MeshError::Io)?;
        for p in [&pa, &pb, &pc] {
            writeln!(writer, "      vertex {} {} {}", p.x, p.y, p.z)
                .map_err(MeshError::Io)?;
        }
        writeln!(writer, "    endloop").map_err(MeshError::Io)?;
        writeln!(writer, "  endfacet").map_err(MeshError::Io)?;
    }
    writeln!(writer, "endsolid {name}").map_err(MeshError::Io)?;
    Ok(())
}

/// Write a [`HalfEdgeMesh`] as binary STL.
///
/// Face normals are recomputed from vertex positions.
pub fn write_stl_binary_he<'id, W: Write>(
    writer: &mut W,
    mesh: &HalfEdgeMesh<'id>,
    token: &GhostToken<'id>,
) -> MeshResult<()> {
    let face_keys: Vec<_> = mesh.face_keys().collect();
    let n = face_keys.len() as u32;
    writer.write_all(&[0u8; 80]).map_err(MeshError::Io)?;
    writer.write_all(&n.to_le_bytes()).map_err(MeshError::Io)?;
    for fk in face_keys {
        let verts = mesh.face_vertices(fk, token);
        if let [v0, v1, v2] = verts.as_slice() {
            let pa = mesh.vertex_pos(*v0, token).unwrap_or_else(|| Point3r::new(0.0, 0.0, 0.0));
            let pb = mesh.vertex_pos(*v1, token).unwrap_or_else(|| Point3r::new(0.0, 0.0, 0.0));
            let pc = mesh.vertex_pos(*v2, token).unwrap_or_else(|| Point3r::new(0.0, 0.0, 0.0));
            let nrm = crate::geometry::normal::triangle_normal(&pa, &pb, &pc)
                .unwrap_or_else(|| Vector3r::z());
            write_f32(writer, nrm.x as f32)?;
            write_f32(writer, nrm.y as f32)?;
            write_f32(writer, nrm.z as f32)?;
            for p in [&pa, &pb, &pc] {
                write_f32(writer, p.x as f32)?;
                write_f32(writer, p.y as f32)?;
                write_f32(writer, p.z as f32)?;
            }
            writer.write_all(&0u16.to_le_bytes()).map_err(MeshError::Io)?;
        }
    }
    Ok(())
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
///     let _ = cfd_mesh::io::stl::fuzz_read_stl(data);
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
    use crate::mesh::with_mesh;

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

    // ── HalfEdgeMesh round-trip ───────────────────────────────────────────

    #[test]
    fn half_edge_ascii_stl_round_trip() {
        let ascii = "solid s\n  facet normal 0 0 1\n    outer loop\n      vertex 0 0 0\n      vertex 1 0 0\n      vertex 0 1 0\n    endloop\n  endfacet\nendsolid s\n";

        let face_count = with_mesh(|mut mesh, mut token| {
            let n = read_stl_he(ascii.as_bytes(), &mut mesh, &mut token).unwrap();
            assert_eq!(n, 1);

            // Write back as ASCII
            let mut out = Vec::new();
            write_stl_ascii_he(&mut out, "s", &mesh, &token).unwrap();
            let text = String::from_utf8(out).unwrap();
            assert!(text.contains("solid s"));
            assert!(text.contains("vertex"));

            mesh.face_count()
        });
        assert_eq!(face_count, 1);
    }

    #[test]
    fn half_edge_binary_stl_round_trip() {
        // Build a one-triangle IndexedMesh, write as binary, then read into HalfEdgeMesh.
        let mut indexed = IndexedMesh::new();
        let v0 = indexed.add_vertex_pos(Point3r::new(0.0, 0.0, 0.0));
        let v1 = indexed.add_vertex_pos(Point3r::new(1.0, 0.0, 0.0));
        let v2 = indexed.add_vertex_pos(Point3r::new(0.0, 1.0, 0.0));
        indexed.add_face(v0, v1, v2);
        let mut buf = Vec::new();
        write_stl_binary(&mut buf, &indexed).unwrap();

        let face_count = with_mesh(|mut mesh, mut token| {
            read_stl_he(buf.as_slice(), &mut mesh, &mut token).unwrap();
            mesh.face_count()
        });
        assert_eq!(face_count, 1);
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
