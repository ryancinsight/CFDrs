//! Hierarchical mesh refinement — P1 to P2 (quadratic) promotion.
//!
//! Subdivides each triangle into 4 sub-triangles by inserting mid-edge
//! nodes. This is standard 1:4 refinement used for P2 element support
//! in Taylor-Hood (P2/P1) finite elements.

use crate::domain::mesh::IndexedMesh;
use crate::infrastructure::storage::face_store::FaceData;
use nalgebra::Point3;
use crate::domain::core::scalar::Scalar;
use crate::domain::core::index::VertexId;
use std::collections::HashMap;

/// Converts a P1 (linear) mesh to a refined mesh by 1:4 triangle subdivision.
///
/// Each triangle edge gets a midpoint vertex, and the original triangle is
/// replaced by four sub-triangles: three corner triangles plus one central
/// triangle formed by the three midpoints.
pub struct P2MeshConverter;

impl P2MeshConverter {
    /// Refine a mesh by 1:4 subdivision of every triangle.
    ///
    /// Each edge of every triangle gets a new mid-point node.
    /// Each original triangle is replaced by 4 sub-triangles.
    pub fn convert_to_p2<T: Scalar>(mesh: &IndexedMesh<T>) -> IndexedMesh<T> {
        let mut out = mesh.clone();

        // Insert midpoint nodes for every unique edge.
        let mut edge_mid: HashMap<(VertexId, VertexId), VertexId> = HashMap::new();

        // Collect new faces to replace the originals.
        let mut new_faces: Vec<FaceData> = Vec::new();

        for (_f_id, face) in mesh.faces.iter_enumerated() {
            let [v0, v1, v2] = face.vertices;
            let region = face.region;

            // Get or create midpoints for each edge.
            let m01 = Self::get_or_create_midpoint(&mut edge_mid, &mut out, v0, v1);
            let m12 = Self::get_or_create_midpoint(&mut edge_mid, &mut out, v1, v2);
            let m20 = Self::get_or_create_midpoint(&mut edge_mid, &mut out, v2, v0);

            // 1:4 subdivision: 3 corner triangles + 1 central triangle.
            new_faces.push(FaceData::new(v0, m01, m20, region));
            new_faces.push(FaceData::new(m01, v1, m12, region));
            new_faces.push(FaceData::new(m20, m12, v2, region));
            new_faces.push(FaceData::new(m01, m12, m20, region));
        }

        // Replace all faces with the subdivided set.
        out.faces.clear();
        for face in new_faces {
            out.faces.push(face);
        }

        out
    }

    /// Get or create the midpoint vertex for an edge.
    fn get_or_create_midpoint<T: Scalar>(
        edge_mid: &mut HashMap<(VertexId, VertexId), VertexId>,
        mesh: &mut IndexedMesh<T>,
        a: VertexId,
        b: VertexId,
    ) -> VertexId {
        let key = if a.as_usize() < b.as_usize() { (a, b) } else { (b, a) };
        *edge_mid.entry(key).or_insert_with(|| {
            let va = mesh.vertices.position(a);
            let vb = mesh.vertices.position(b);
            let two = T::one() + T::one();
            let mid = Point3::new(
                (va.x + vb.x) / two,
                (va.y + vb.y) / two,
                (va.z + vb.z) / two,
            );
            mesh.add_vertex_pos(mid)
        })
    }
}
