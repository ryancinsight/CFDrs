//! Orthographic projector — projects 3D mesh geometry into 2D line drawings.
//!
//! # Theorem — Orthographic Projection Preserves Parallel Lines
//!
//! An orthographic projection P: R^3 -> R^2 defined by dropping one coordinate
//! (after rotation into view space) maps parallel line segments to parallel
//! line segments.  **Proof sketch**: If lines L1 and L2 have direction d in R^3,
//! then P(L1) and P(L2) have direction P(d) in R^2.  Since P is linear,
//! P(d) is the same for both, so the projected lines are parallel.  ∎

use cfd_mesh::IndexedMesh;
use nalgebra::{Matrix3, Point2, Point3, Vector3};

/// Edges classified by visibility after projection.
pub struct ProjectedEdges {
    /// Edges visible from the view direction.
    pub visible: Vec<(Point2<f64>, Point2<f64>)>,
    /// Edges hidden behind other geometry.
    pub hidden: Vec<(Point2<f64>, Point2<f64>)>,
    /// Silhouette edges (boundary between front-facing and back-facing triangles).
    pub silhouette: Vec<(Point2<f64>, Point2<f64>)>,
}

/// Projects meshes into 2D engineering drawing views.
pub struct OrthographicProjector {
    /// View direction (world space, from camera toward model).
    view_dir: Vector3<f64>,
    /// Up direction in view space.
    _up_dir: Vector3<f64>,
    /// Rotation matrix from world space to view space.
    rotation: Matrix3<f64>,
}

impl OrthographicProjector {
    /// Create a projector with the given view and up directions.
    #[must_use]
    pub fn new(view_dir: Vector3<f64>, up_dir: Vector3<f64>) -> Self {
        let forward = view_dir.normalize();
        let right = forward.cross(&up_dir).normalize();
        let up = right.cross(&forward).normalize();
        // Rotation maps world coords to (right, up, -forward) view coords.
        let rotation = Matrix3::from_rows(&[right.transpose(), up.transpose(), (-forward).transpose()]);
        Self { view_dir: forward, _up_dir: up, rotation }
    }

    /// Project a 3D point to 2D view coordinates.
    #[must_use]
    pub fn project_point(&self, p: &Point3<f64>) -> Point2<f64> {
        let v = self.rotation * p.coords;
        Point2::new(v.x, v.y)
    }

    /// Project all edges of a mesh and classify them as visible, hidden, or silhouette.
    ///
    /// Uses a simplified visibility heuristic: edges whose face normals both
    /// face away from the view direction are classified as hidden.
    #[must_use]
    pub fn project_mesh(&self, mesh: &IndexedMesh<f64>) -> ProjectedEdges {
        let mut visible = Vec::new();
        let mut hidden = Vec::new();
        let mut silhouette = Vec::new();

        // Precompute face normals and front-facing flags.
        let face_data: Vec<(bool, [Point3<f64>; 3])> = mesh
            .faces
            .iter()
            .map(|face| {
                let [v0, v1, v2] = face.vertices;
                let p0 = *mesh.vertices.position(v0);
                let p1 = *mesh.vertices.position(v1);
                let p2 = *mesh.vertices.position(v2);
                let normal = (p1 - p0).cross(&(p2 - p0));
                let front_facing = normal.dot(&(-self.view_dir)) > 0.0;
                (front_facing, [p0, p1, p2])
            })
            .collect();

        // Build edge -> face adjacency.
        let mut edge_faces: std::collections::HashMap<(usize, usize), Vec<usize>> =
            std::collections::HashMap::new();
        for (fi, face) in mesh.faces.iter().enumerate() {
            let verts = face.vertices.map(|v| v.0 as usize);
            for &(a, b) in &[(verts[0], verts[1]), (verts[1], verts[2]), (verts[2], verts[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                edge_faces.entry(key).or_default().push(fi);
            }
        }

        // Classify each unique edge.
        for (&(a, b), faces) in &edge_faces {
            let pa = *mesh.vertices.position(cfd_mesh::domain::core::index::VertexId(a as u32));
            let pb = *mesh.vertices.position(cfd_mesh::domain::core::index::VertexId(b as u32));
            let p2d_a = self.project_point(&pa);
            let p2d_b = self.project_point(&pb);

            if faces.len() == 2 {
                let f0_front = face_data[faces[0]].0;
                let f1_front = face_data[faces[1]].0;
                if f0_front != f1_front {
                    silhouette.push((p2d_a, p2d_b));
                } else if f0_front {
                    visible.push((p2d_a, p2d_b));
                } else {
                    hidden.push((p2d_a, p2d_b));
                }
            } else if faces.len() == 1 {
                // Boundary edge — always visible.
                if face_data[faces[0]].0 {
                    visible.push((p2d_a, p2d_b));
                } else {
                    silhouette.push((p2d_a, p2d_b));
                }
            }
        }

        ProjectedEdges { visible, hidden, silhouette }
    }
}
