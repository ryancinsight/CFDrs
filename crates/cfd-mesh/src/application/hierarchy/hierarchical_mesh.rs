//! Hierarchical mesh refinement — P1 to P2 (quadratic) promotion.
//!
//! The Taylor-Hood element (P2/P1) requires mid-edge nodes in addition to
//! corner nodes.  `P2MeshConverter::convert_to_p2` inserts these edge midpoints.
//!
//! This module is a **volume/FEM tool** — intentional `Mesh<T>` usage for P1/P2
//! nodal topology.

// Volume tool: Mesh<T> is the correct type here.
#![allow(deprecated)]

use crate::domain::mesh::IndexedMesh;
use nalgebra::Point3;
use crate::domain::core::scalar::Scalar;
use crate::domain::core::index::VertexId;
use std::collections::HashMap;

/// Converts a P1 (linear) mesh to a P2 (quadratic) mesh by inserting
/// mid-edge nodes on every tetrahedron edge.
pub struct P2MeshConverter;

impl P2MeshConverter {
    /// Promote a linear mesh to P2 (quadratic) elements.
    ///
    /// Each edge of every tetrahedron gets a new mid-point node.
    /// The returned mesh has the same cells but double (roughly) the vertices.
    pub fn convert_to_p2<T: Scalar>(mesh: &IndexedMesh<T>) -> IndexedMesh<T> {
        let mut out = mesh.clone();

        // Insert midpoint nodes for every unique face edge.
        let mut edge_mid: HashMap<(VertexId, VertexId), VertexId> = HashMap::new();

        for (_, face) in mesh.faces.iter_enumerated() {
            let verts = &face.vertices;
            let n = verts.len();
            for i in 0..n {
                let a = verts[i];
                let b = verts[(i + 1) % n];
                let key = if a.as_usize() < b.as_usize() { (a, b) } else { (b, a) };
                if !edge_mid.contains_key(&key) {
                    let va = mesh.vertices.position(a);
                    let vb = mesh.vertices.position(b);
                    let two = T::one() + T::one();
                    let mid = Point3::new(
                        (va.x + vb.x) / two,
                        (va.y + vb.y) / two,
                        (va.z + vb.z) / two,
                    );
                    let mid_idx = out.add_vertex_pos(mid);
                    edge_mid.insert(key, mid_idx);
                }
            }
        }

        out
    }
}
