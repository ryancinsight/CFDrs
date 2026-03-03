//! glTF 2.0 Binary (GLB) export.
//!
//! Produces a single `.glb` file containing the mesh geometry in a binary
//! buffer. No external dependencies are required — the JSON and binary
//! payloads are assembled manually per the glTF 2.0 specification.

use std::collections::HashMap;
use std::io::Write;

use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::mesh::IndexedMesh;

/// Write an [`IndexedMesh`] as a glTF 2.0 Binary (`.glb`) file.
pub fn write_glb<W: Write>(writer: &mut W, mesh: &IndexedMesh) -> MeshResult<()> {
    // Build contiguous vertex data (position f32x3 + normal f32x3 = 24 bytes/vertex).
    let mut id_to_idx: HashMap<crate::domain::core::index::VertexId, u32> = HashMap::new();
    let mut positions: Vec<[f32; 3]> = Vec::new();
    let mut normals: Vec<[f32; 3]> = Vec::new();
    let mut min_pos = [f32::MAX; 3];
    let mut max_pos = [f32::MIN; 3];

    for (idx, (vid, vdata)) in mesh.vertices.iter().enumerate() {
        id_to_idx.insert(vid, idx as u32);
        let p = [
            vdata.position.x as f32,
            vdata.position.y as f32,
            vdata.position.z as f32,
        ];
        let n = [
            vdata.normal.x as f32,
            vdata.normal.y as f32,
            vdata.normal.z as f32,
        ];
        for i in 0..3 {
            min_pos[i] = min_pos[i].min(p[i]);
            max_pos[i] = max_pos[i].max(p[i]);
        }
        positions.push(p);
        normals.push(n);
    }

    // Build index buffer (u32).
    let mut indices: Vec<u32> = Vec::new();
    for (_fid, face) in mesh.faces.iter_enumerated() {
        indices.push(id_to_idx[&face.vertices[0]]);
        indices.push(id_to_idx[&face.vertices[1]]);
        indices.push(id_to_idx[&face.vertices[2]]);
    }

    let vertex_count = positions.len();
    let index_count = indices.len();

    // Build binary buffer: indices, then positions, then normals.
    let indices_byte_len = index_count * 4;
    let positions_byte_len = vertex_count * 12;
    let normals_byte_len = vertex_count * 12;
    let total_buffer_len = indices_byte_len + positions_byte_len + normals_byte_len;

    let mut bin_buf: Vec<u8> = Vec::with_capacity(total_buffer_len);

    // Indices
    for &i in &indices {
        bin_buf.extend_from_slice(&i.to_le_bytes());
    }
    // Positions
    for p in &positions {
        for &v in p {
            bin_buf.extend_from_slice(&v.to_le_bytes());
        }
    }
    // Normals
    for n in &normals {
        for &v in n {
            bin_buf.extend_from_slice(&v.to_le_bytes());
        }
    }

    // Pad binary buffer to 4-byte alignment.
    while !bin_buf.len().is_multiple_of(4) {
        bin_buf.push(0);
    }
    let padded_bin_len = bin_buf.len();

    // Build the JSON string.
    let json_str = format!(
        r#"{{"asset":{{"version":"2.0","generator":"cfd-mesh"}},"scene":0,"scenes":[{{"nodes":[0]}}],"nodes":[{{"mesh":0}}],"meshes":[{{"primitives":[{{"attributes":{{"POSITION":1,"NORMAL":2}},"indices":0}}]}}],"accessors":[{{"bufferView":0,"componentType":5125,"count":{index_count},"type":"SCALAR"}},{{"bufferView":1,"componentType":5126,"count":{vertex_count},"type":"VEC3","min":[{},{},{}],"max":[{},{},{}]}},{{"bufferView":2,"componentType":5126,"count":{vertex_count},"type":"VEC3"}}],"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":{indices_byte_len},"target":34963}},{{"buffer":0,"byteOffset":{indices_byte_len},"byteLength":{positions_byte_len},"target":34962}},{{"buffer":0,"byteOffset":{},"byteLength":{normals_byte_len},"target":34962}}],"buffers":[{{"byteLength":{total_buffer_len}}}]}}"#,
        min_pos[0],
        min_pos[1],
        min_pos[2],
        max_pos[0],
        max_pos[1],
        max_pos[2],
        indices_byte_len + positions_byte_len,
    );

    // Pad JSON to 4-byte alignment with spaces (GLB spec requirement).
    let mut json_bytes = json_str.into_bytes();
    while json_bytes.len() % 4 != 0 {
        json_bytes.push(b' ');
    }
    let padded_json_len = json_bytes.len();

    // GLB header: magic + version + total length.
    let total_len = 12 + 8 + padded_json_len + 8 + padded_bin_len;

    // GLB magic
    writer.write_all(b"glTF").map_err(MeshError::Io)?;
    // Version 2
    writer
        .write_all(&2u32.to_le_bytes())
        .map_err(MeshError::Io)?;
    // Total file length
    writer
        .write_all(&(total_len as u32).to_le_bytes())
        .map_err(MeshError::Io)?;

    // JSON chunk header
    writer
        .write_all(&(padded_json_len as u32).to_le_bytes())
        .map_err(MeshError::Io)?;
    writer
        .write_all(&0x4E4F534Au32.to_le_bytes())
        .map_err(MeshError::Io)?; // "JSON"
    writer.write_all(&json_bytes).map_err(MeshError::Io)?;

    // BIN chunk header
    writer
        .write_all(&(padded_bin_len as u32).to_le_bytes())
        .map_err(MeshError::Io)?;
    writer
        .write_all(&0x004E4942u32.to_le_bytes())
        .map_err(MeshError::Io)?; // "BIN\0"
    writer.write_all(&bin_buf).map_err(MeshError::Io)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;

    #[test]
    fn glb_write_produces_valid_header() {
        let mut mesh = IndexedMesh::new();
        let v0 = mesh.add_vertex_pos(Point3r::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3r::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3r::new(0.0, 1.0, 0.0));
        mesh.add_face(v0, v1, v2);

        let mut buf = Vec::new();
        write_glb(&mut buf, &mesh).unwrap();

        // Check GLB magic.
        assert_eq!(&buf[0..4], b"glTF");
        // Check version 2.
        assert_eq!(u32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]), 2);
        // Check total length matches buffer.
        let total = u32::from_le_bytes([buf[8], buf[9], buf[10], buf[11]]) as usize;
        assert_eq!(total, buf.len());
    }
}
