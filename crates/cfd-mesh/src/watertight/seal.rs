//! Boundary sealing: close holes in an otherwise-manifold mesh.

use crate::core::index::{VertexId, FaceId, RegionId};
use crate::core::scalar::{Point3r, Vector3r};
use crate::storage::edge_store::EdgeStore;
use crate::storage::face_store::{FaceData, FaceStore};
use crate::storage::vertex_pool::VertexPool;

/// Seal boundary loops by fan triangulation from the centroid.
///
/// For each connected boundary loop, insert a centroid vertex and create
/// triangles from each boundary edge to the centroid.
///
/// Returns the number of faces added.
pub fn seal_boundary_loops(
    vertex_pool: &mut VertexPool,
    face_store: &mut FaceStore,
    edge_store: &EdgeStore,
    region: RegionId,
) -> usize {
    let boundary = edge_store.boundary_edges();
    if boundary.is_empty() {
        return 0;
    }

    // Collect boundary edges as directed pairs
    // Collect boundary edges as directed pairs
    let mut boundary_pairs: Vec<(VertexId, VertexId)> = Vec::new();
    for &eid in &boundary {
        let edge = edge_store.get(eid);
        // For boundary edges (valence 1), the single adjacent face determines the winding.
        let face_id = edge.faces[0];
        let face = face_store.get(face_id);
        
        let (v0, v1) = edge.vertices;
        
        // Check if the face uses the edge as (v0, v1) or (v1, v0).
        // The boundary loop should run in the opposite direction to seal it.
        // Face edges are: (f.v[0], f.v[1]), (f.v[1], f.v[2]), (f.v[2], f.v[0])
        let mut is_forward = false;
        let [a, b, c] = face.vertices;
        if (a == v0 && b == v1) || (b == v0 && c == v1) || (c == v0 && a == v1) {
            is_forward = true;
        }

        if is_forward {
            // Face has (v0 -> v1). Boundary loop must be (v1 -> v0).
            boundary_pairs.push((v1, v0));
        } else {
            // Face has (v1 -> v0). Boundary loop must be (v0 -> v1).
            boundary_pairs.push((v0, v1));
        }
    }

    // Find connected loops
    let loops = extract_boundary_loops(&boundary_pairs);

    let mut faces_added = 0;
    for boundary_loop in &loops {
        if boundary_loop.len() < 3 {
            continue;
        }

        // Compute centroid
        let mut centroid = Point3r::origin();
        for &vid in boundary_loop {
            centroid.coords += vertex_pool.position(vid).coords;
        }
        centroid.coords /= boundary_loop.len() as crate::core::scalar::Real;

        let centroid_id = vertex_pool.insert_or_weld(centroid, Vector3r::zeros());

        // Fan triangulate
        for i in 0..boundary_loop.len() {
            let j = (i + 1) % boundary_loop.len();
            face_store.push(FaceData {
                vertices: [boundary_loop[i], boundary_loop[j], centroid_id],
                region,
            });
            faces_added += 1;
        }
    }

    faces_added
}

/// Extract connected loops from a set of directed edges.
fn extract_boundary_loops(edges: &[(VertexId, VertexId)]) -> Vec<Vec<VertexId>> {
    use hashbrown::HashMap;

    let mut next_map: HashMap<VertexId, VertexId> = HashMap::new();
    for &(a, b) in edges {
        next_map.insert(a, b);
    }

    let mut loops = Vec::new();
    let mut visited = hashbrown::HashSet::new();

    for &(start, _) in edges {
        if visited.contains(&start) {
            continue;
        }

        let mut loop_verts = Vec::new();
        let mut current = start;

        loop {
            if visited.contains(&current) {
                break;
            }
            visited.insert(current);
            loop_verts.push(current);

            match next_map.get(&current) {
                Some(&next) => current = next,
                None => break,
            }
        }

        if loop_verts.len() >= 3 {
            loops.push(loop_verts);
        }
    }

    loops
}
