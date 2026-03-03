//! Stanford PLY import and export.
//!
//! Supports ASCII PLY with vertex positions, normals, and triangular faces.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};

use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

// =============================================================================
//  Export
// =============================================================================

/// Write an [`IndexedMesh`] as ASCII PLY.
pub fn write_ply<W: Write>(writer: &mut W, mesh: &IndexedMesh) -> MeshResult<()> {
    let vertex_count = mesh.vertex_count();
    let face_count = mesh.face_count();

    // Header
    writeln!(writer, "ply").map_err(MeshError::Io)?;
    writeln!(writer, "format ascii 1.0").map_err(MeshError::Io)?;
    writeln!(writer, "comment exported by cfd-mesh").map_err(MeshError::Io)?;
    writeln!(writer, "element vertex {vertex_count}").map_err(MeshError::Io)?;
    writeln!(writer, "property float x").map_err(MeshError::Io)?;
    writeln!(writer, "property float y").map_err(MeshError::Io)?;
    writeln!(writer, "property float z").map_err(MeshError::Io)?;
    writeln!(writer, "property float nx").map_err(MeshError::Io)?;
    writeln!(writer, "property float ny").map_err(MeshError::Io)?;
    writeln!(writer, "property float nz").map_err(MeshError::Io)?;
    writeln!(writer, "element face {face_count}").map_err(MeshError::Io)?;
    writeln!(writer, "property list uchar int vertex_indices").map_err(MeshError::Io)?;
    writeln!(writer, "end_header").map_err(MeshError::Io)?;

    // Build contiguous index map.
    let mut id_to_idx: HashMap<crate::domain::core::index::VertexId, usize> = HashMap::new();
    let mut idx = 0usize;

    // Vertex data
    for (vid, vdata) in mesh.vertices.iter() {
        id_to_idx.insert(vid, idx);
        idx += 1;
        let p = &vdata.position;
        let n = &vdata.normal;
        writeln!(writer, "{} {} {} {} {} {}", p.x, p.y, p.z, n.x, n.y, n.z)
            .map_err(MeshError::Io)?;
    }

    // Face data
    for (_fid, face) in mesh.faces.iter_enumerated() {
        let i0 = id_to_idx[&face.vertices[0]];
        let i1 = id_to_idx[&face.vertices[1]];
        let i2 = id_to_idx[&face.vertices[2]];
        writeln!(writer, "3 {i0} {i1} {i2}").map_err(MeshError::Io)?;
    }

    Ok(())
}

// =============================================================================
//  Import
// =============================================================================

/// Read an ASCII PLY file into a new [`IndexedMesh`].
pub fn read_ply<R: Read>(reader: R) -> MeshResult<IndexedMesh> {
    let buf = BufReader::new(reader);
    let mut lines = buf.lines();

    // Parse header.
    let mut vertex_count = 0usize;
    let mut face_count = 0usize;
    let mut in_vertex_props = false;
    let mut normal_prop_count = 0u8;

    // Read the "ply" magic.
    let magic = next_line(&mut lines)?;
    if magic.trim() != "ply" {
        return Err(MeshError::Other("not a PLY file".to_owned()));
    }

    loop {
        let line = next_line(&mut lines)?;
        let trimmed = line.trim();

        if trimmed == "end_header" {
            break;
        }
        if trimmed.starts_with("element vertex") {
            vertex_count = trimmed
                .split_whitespace()
                .nth(2)
                .ok_or_else(|| MeshError::Other("bad vertex count".to_owned()))?
                .parse()
                .map_err(|_| MeshError::Other("bad vertex count".to_owned()))?;
            in_vertex_props = true;
        } else if trimmed.starts_with("element face") {
            face_count = trimmed
                .split_whitespace()
                .nth(2)
                .ok_or_else(|| MeshError::Other("bad face count".to_owned()))?
                .parse()
                .map_err(|_| MeshError::Other("bad face count".to_owned()))?;
            in_vertex_props = false;
        } else if trimmed.starts_with("element") {
            in_vertex_props = false;
        }

        if in_vertex_props {
            if trimmed.contains(" nx") || trimmed.contains(" ny") || trimmed.contains(" nz") {
                normal_prop_count += 1;
            }
        }
    }

    let has_normals = normal_prop_count >= 3;

    // Read vertex data.
    let mut positions = Vec::with_capacity(vertex_count);
    let mut normals_vec = Vec::with_capacity(if has_normals { vertex_count } else { 0 });

    for _ in 0..vertex_count {
        let line = next_line(&mut lines)?;
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        if parts.len() < 3 {
            return Err(MeshError::Other("vertex line too short".to_owned()));
        }
        let x = parse_f64(parts[0])?;
        let y = parse_f64(parts[1])?;
        let z = parse_f64(parts[2])?;
        positions.push(Point3r::new(x, y, z));

        if has_normals && parts.len() >= 6 {
            let nx = parse_f64(parts[3])?;
            let ny = parse_f64(parts[4])?;
            let nz = parse_f64(parts[5])?;
            normals_vec.push(Vector3r::new(nx, ny, nz));
        }
    }

    // Build indexed mesh from vertices.
    let mut mesh = IndexedMesh::new();
    let region = RegionId::from_usize(0);

    let vertex_ids: Vec<_> = (0..positions.len())
        .map(|i| {
            let n = normals_vec
                .get(i)
                .copied()
                .unwrap_or_else(Vector3r::zeros);
            mesh.add_vertex(positions[i], n)
        })
        .collect();

    // Read face data.
    for _ in 0..face_count {
        let line = next_line(&mut lines)?;
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }
        let n_verts: usize = parts[0]
            .parse()
            .map_err(|_| MeshError::Other("bad face vertex count".to_owned()))?;
        if parts.len() < n_verts + 1 {
            return Err(MeshError::Other("face line too short".to_owned()));
        }

        let face_verts: Vec<usize> = parts[1..=n_verts]
            .iter()
            .map(|s| {
                s.parse::<usize>()
                    .map_err(|_| MeshError::Other(format!("bad face index: {s}")))
            })
            .collect::<MeshResult<_>>()?;

        // Fan-triangulate.
        for i in 1..face_verts.len() - 1 {
            mesh.add_face_with_region(
                vertex_ids[face_verts[0]],
                vertex_ids[face_verts[i]],
                vertex_ids[face_verts[i + 1]],
                region,
            );
        }
    }

    Ok(mesh)
}

fn next_line(
    lines: &mut std::io::Lines<BufReader<impl Read>>,
) -> MeshResult<String> {
    lines
        .next()
        .ok_or_else(|| MeshError::Other("unexpected end of PLY file".to_owned()))?
        .map_err(MeshError::Io)
}

fn parse_f64(s: &str) -> MeshResult<f64> {
    s.parse::<f64>()
        .map_err(|_| MeshError::Other(format!("invalid number: {s}")))
}

// =============================================================================
//  Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ply_round_trip() {
        let mut mesh = IndexedMesh::new();
        let v0 = mesh.add_vertex_pos(Point3r::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3r::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3r::new(0.0, 1.0, 0.0));
        mesh.add_face(v0, v1, v2);

        let mut buf = Vec::new();
        write_ply(&mut buf, &mesh).unwrap();

        let mesh2 = read_ply(std::io::Cursor::new(&buf)).unwrap();
        assert_eq!(mesh2.vertex_count(), 3);
        assert_eq!(mesh2.face_count(), 1);
    }
}
