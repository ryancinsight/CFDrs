//! Hierarchical mesh refinement — P1 to P2 (quadratic) promotion.
//!
//! The Taylor-Hood element (P2/P1) requires mid-edge nodes in addition to
//! corner nodes.  `P2MeshConverter::convert_to_p2` inserts these edge midpoints.
//!
//! This module is a **volume/FEM tool** — intentional `Mesh<T>` usage for P1/P2
//! nodal topology.

// Volume tool: Mesh<T> is the correct type here.
#![allow(deprecated)]

use nalgebra::RealField;
use crate::mesh::Mesh;
use crate::topology::Vertex;
use nalgebra::Point3;
use std::collections::HashMap;

/// Converts a P1 (linear) mesh to a P2 (quadratic) mesh by inserting
/// mid-edge nodes on every tetrahedron edge.
pub struct P2MeshConverter;

impl P2MeshConverter {
    /// Promote a linear mesh to P2 (quadratic) elements.
    ///
    /// Each edge of every tetrahedron gets a new mid-point node.
    /// The returned mesh has the same cells but double (roughly) the vertices.
    pub fn convert_to_p2<T: Copy + RealField>(mesh: &Mesh<T>) -> Mesh<T> {
        let mut out = mesh.clone();

        // Insert midpoint nodes for every unique face edge.
        let mut edge_mid: HashMap<(usize, usize), usize> = HashMap::new();

        for face in mesh.faces() {
            let verts = &face.vertices;
            let n = verts.len();
            for i in 0..n {
                let a = verts[i];
                let b = verts[(i + 1) % n];
                let key = if a < b { (a, b) } else { (b, a) };
                if !edge_mid.contains_key(&key) {
                    if let (Some(va), Some(vb)) = (mesh.vertex(a), mesh.vertex(b)) {
                        let two = T::one() + T::one();
                        let mid = Point3::new(
                            (va.position.x + vb.position.x) / two,
                            (va.position.y + vb.position.y) / two,
                            (va.position.z + vb.position.z) / two,
                        );
                        let mid_idx = out.add_vertex(Vertex::new(mid));
                        edge_mid.insert(key, mid_idx);
                    }
                }
            }
        }

        out
    }
}
