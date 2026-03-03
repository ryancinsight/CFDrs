//! Wavefront OBJ import and export.
//!
//! Supports triangle meshes with vertex positions and normals.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};

use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

// =============================================================================
//  Export
// =============================================================================

/// Write an [`IndexedMesh`] as Wavefront OBJ.
///
/// Emits `v` (position), `vn` (normal), and `f` (face) records.
/// OBJ uses 1-based indexing.
pub fn write_obj<W: Write>(writer: &mut W, mesh: &IndexedMesh) -> MeshResult<()> {
    writeln!(writer, "# OBJ exported by cfd-mesh").map_err(MeshError::Io)?;

    // Build a contiguous index map: VertexId -> 0-based index.
    let mut id_to_idx: HashMap<crate::domain::core::index::VertexId, usize> = HashMap::new();
    let mut idx = 0usize;

    // Emit vertices and normals in insertion order.
    for (vid, vdata) in mesh.vertices.iter() {
        id_to_idx.insert(vid, idx);
        idx += 1;
        let p = &vdata.position;
        writeln!(writer, "v {} {} {}", p.x, p.y, p.z).map_err(MeshError::Io)?;
    }

    for (_vid, vdata) in mesh.vertices.iter() {
        let n = &vdata.normal;
        writeln!(writer, "vn {} {} {}", n.x, n.y, n.z).map_err(MeshError::Io)?;
    }

    // Emit faces (1-indexed).
    for (_fid, face) in mesh.faces.iter_enumerated() {
        let i0 = id_to_idx[&face.vertices[0]] + 1;
        let i1 = id_to_idx[&face.vertices[1]] + 1;
        let i2 = id_to_idx[&face.vertices[2]] + 1;
        writeln!(writer, "f {i0}//{i0} {i1}//{i1} {i2}//{i2}").map_err(MeshError::Io)?;
    }

    Ok(())
}

// =============================================================================
//  Import
// =============================================================================

/// Read a Wavefront OBJ file into a new [`IndexedMesh`].
///
/// Supports `v`, `vn`, and `f` records. Face specifications may be:
/// - `f v1 v2 v3` (position only)
/// - `f v1//vn1 v2//vn2 v3//vn3` (position + normal)
/// - `f v1/vt1/vn1 ...` (position + texcoord + normal, texcoord ignored)
///
/// Polygonal faces with more than 3 vertices are fan-triangulated.
pub fn read_obj<R: Read>(reader: R) -> MeshResult<IndexedMesh> {
    let buf = BufReader::new(reader);
    let mut positions: Vec<Point3r> = Vec::new();
    let mut normals: Vec<Vector3r> = Vec::new();
    let mut mesh = IndexedMesh::new();
    let region = RegionId::from_usize(0);

    for line in buf.lines() {
        let line = line.map_err(MeshError::Io)?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        match parts[0] {
            "v" if parts.len() >= 4 => {
                let x = parse_f64(parts[1])?;
                let y = parse_f64(parts[2])?;
                let z = parse_f64(parts[3])?;
                positions.push(Point3r::new(x, y, z));
            }
            "vn" if parts.len() >= 4 => {
                let x = parse_f64(parts[1])?;
                let y = parse_f64(parts[2])?;
                let z = parse_f64(parts[3])?;
                normals.push(Vector3r::new(x, y, z));
            }
            "f" if parts.len() >= 4 => {
                // Parse face vertex indices.
                let verts: Vec<(usize, Option<usize>)> = parts[1..]
                    .iter()
                    .map(|s| parse_face_vertex(s))
                    .collect::<MeshResult<_>>()?;

                // Convert to mesh vertex IDs.
                let ids: Vec<_> = verts
                    .iter()
                    .map(|&(pi, ni)| {
                        let pos = positions[pi];
                        let normal = ni
                            .and_then(|i| normals.get(i).copied())
                            .unwrap_or_else(Vector3r::zeros);
                        mesh.add_vertex(pos, normal)
                    })
                    .collect();

                // Fan-triangulate polygonal faces.
                for i in 1..ids.len() - 1 {
                    mesh.add_face_with_region(ids[0], ids[i], ids[i + 1], region);
                }
            }
            _ => {} // Skip unsupported records (vt, mtllib, usemtl, o, g, s, etc.)
        }
    }

    Ok(mesh)
}

/// Parse a single f64 value, returning a `MeshError` on failure.
fn parse_f64(s: &str) -> MeshResult<f64> {
    s.parse::<f64>()
        .map_err(|_| MeshError::Other(format!("invalid number: {s}")))
}

/// Parse a face vertex specification like `v`, `v//vn`, or `v/vt/vn`.
/// Returns `(position_index, optional_normal_index)`, converted to 0-based.
fn parse_face_vertex(s: &str) -> MeshResult<(usize, Option<usize>)> {
    let parts: Vec<&str> = s.split('/').collect();

    let pos_idx: usize = parts[0]
        .parse::<usize>()
        .map_err(|_| MeshError::Other(format!("invalid face vertex: {s}")))?
        .checked_sub(1)
        .ok_or_else(|| MeshError::Other(format!("face vertex index 0 is invalid: {s}")))?;

    let normal_idx = if parts.len() == 3 && !parts[2].is_empty() {
        Some(
            parts[2]
                .parse::<usize>()
                .map_err(|_| MeshError::Other(format!("invalid normal index: {s}")))?
                .checked_sub(1)
                .ok_or_else(|| MeshError::Other(format!("normal index 0 is invalid: {s}")))?,
        )
    } else {
        None
    };

    Ok((pos_idx, normal_idx))
}

// =============================================================================
//  Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn obj_round_trip() {
        let mut mesh = IndexedMesh::new();
        let v0 = mesh.add_vertex_pos(Point3r::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3r::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3r::new(0.0, 1.0, 0.0));
        mesh.add_face(v0, v1, v2);

        let mut buf = Vec::new();
        write_obj(&mut buf, &mesh).unwrap();

        let mesh2 = read_obj(std::io::Cursor::new(&buf)).unwrap();
        assert_eq!(mesh2.vertex_count(), 3);
        assert_eq!(mesh2.face_count(), 1);
    }

    #[test]
    fn obj_quad_fan_triangulated() {
        let obj = b"v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1 2 3 4\n";
        let mesh = read_obj(std::io::Cursor::new(&obj[..])).unwrap();
        // A quad should produce 2 triangles via fan triangulation.
        assert_eq!(mesh.face_count(), 2);
    }
}
