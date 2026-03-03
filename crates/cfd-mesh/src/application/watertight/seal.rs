//! Boundary sealing: close holes in an otherwise-manifold mesh.

use crate::domain::core::index::{RegionId, VertexId};
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::infrastructure::storage::edge_store::EdgeStore;
use crate::infrastructure::storage::face_store::{FaceData, FaceStore};
use crate::infrastructure::storage::vertex_pool::VertexPool;

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
        centroid.coords /= boundary_loop.len() as crate::domain::core::scalar::Real;

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
    use hashbrown::{HashMap, HashSet};

    // Build adjacency that preserves all outgoing boundary links.
    let mut adj: HashMap<VertexId, Vec<VertexId>> = HashMap::new();
    for &(a, b) in edges {
        adj.entry(a).or_default().push(b);
    }
    for nexts in adj.values_mut() {
        nexts.sort();
    }

    let mut used: HashSet<(VertexId, VertexId)> = HashSet::new();
    let mut loops = Vec::new();
    let mut starts: Vec<VertexId> = adj.keys().copied().collect();
    starts.sort();

    for start in starts {
        let successors = match adj.get(&start) {
            Some(s) => s.clone(),
            None => continue,
        };
        for first_next in successors {
            if used.contains(&(start, first_next)) {
                continue;
            }
            let mut path: Vec<VertexId> = vec![start, first_next];
            used.insert((start, first_next));
            let mut cur = first_next;
            let mut closed = false;

            loop {
                if path.len() > 4096 {
                    break;
                }
                let nexts = match adj.get(&cur) {
                    Some(s) => s,
                    None => break,
                };
                let mut found = false;
                for &n in nexts {
                    if used.contains(&(cur, n)) {
                        continue;
                    }
                    used.insert((cur, n));
                    if n == start {
                        closed = true;
                        found = true;
                        break;
                    }
                    // Split figure-8 style inner cycles into separate simple loops.
                    if let Some(pos) = path.iter().position(|&v| v == n) {
                        let inner = path[pos..].to_vec();
                        if inner.len() >= 3 {
                            loops.push(inner);
                        }
                        path.truncate(pos + 1);
                        cur = n;
                        found = true;
                        break;
                    }
                    path.push(n);
                    cur = n;
                    found = true;
                    break;
                }
                if !found || closed {
                    break;
                }
            }

            if closed && path.len() >= 3 {
                loops.push(path);
            }
        }
    }

    loops
}
