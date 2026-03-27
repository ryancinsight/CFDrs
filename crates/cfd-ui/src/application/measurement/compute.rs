//! Measurement computation — geometric queries on meshes and points.
//!
//! Reuses cfd-mesh geometry primitives for face area and mesh volume.

use cfd_mesh::domain::core::index::FaceId;
use cfd_mesh::IndexedMesh;
use nalgebra::{Point3, Vector3};

/// Euclidean distance between two points.
#[must_use]
pub fn compute_point_to_point(a: &Point3<f64>, b: &Point3<f64>) -> f64 {
    (b - a).norm()
}

/// Length of an edge defined by two vertices.
#[must_use]
pub fn compute_edge_length(v0: &Point3<f64>, v1: &Point3<f64>) -> f64 {
    (v1 - v0).norm()
}

/// Dihedral angle between two faces of a mesh, in degrees.
///
/// Computes the angle between the face normals. Returns `None` if either
/// face index is out of range or either face is degenerate.
#[must_use]
pub fn compute_dihedral_angle(
    mesh: &IndexedMesh<f64>,
    face_a: u32,
    face_b: u32,
) -> Option<f64> {
    let na = face_normal(mesh, face_a)?;
    let nb = face_normal(mesh, face_b)?;
    let cos_angle = na.dot(&nb).clamp(-1.0, 1.0);
    Some(cos_angle.acos().to_degrees())
}

/// Area of a single triangular face.
///
/// Uses the cross-product magnitude formula: `area = 0.5 * |e1 x e2|`.
#[must_use]
pub fn compute_face_area(mesh: &IndexedMesh<f64>, face_idx: u32) -> Option<f64> {
    let (v0, v1, v2) = face_vertices(mesh, face_idx)?;
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    Some(0.5 * e1.cross(&e2).norm())
}

/// Total signed volume of a closed mesh using the divergence theorem.
///
/// # Theorem — Divergence Theorem Volume
///
/// For a closed oriented surface, the enclosed volume is:
/// `V = (1/6) * sum_over_faces(v0 . (v1 x v2))`
/// for each triangular face `(v0, v1, v2)` in CCW winding order.  ∎
#[must_use]
pub fn compute_mesh_volume(mesh: &IndexedMesh<f64>) -> f64 {
    let mut volume: f64 = 0.0;
    for face in mesh.faces.iter() {
        let [vi0, vi1, vi2] = face.vertices;
        let d0 = mesh.vertices.get(vi0);
        let d1 = mesh.vertices.get(vi1);
        let d2 = mesh.vertices.get(vi2);
        let v0 = d0.position.coords;
        let v1 = d1.position.coords;
        let v2 = d2.position.coords;
        volume += v0.dot(&v1.cross(&v2));
    }
    (volume / 6.0).abs()
}

/// Center of a triangular face (centroid).
#[must_use]
pub fn face_centroid(mesh: &IndexedMesh<f64>, face_idx: u32) -> Option<Point3<f64>> {
    let (v0, v1, v2) = face_vertices(mesh, face_idx)?;
    Some(Point3::from((v0.coords + v1.coords + v2.coords) / 3.0))
}

/// Unit normal of a triangular face.
fn face_normal(mesh: &IndexedMesh<f64>, face_idx: u32) -> Option<Vector3<f64>> {
    let (v0, v1, v2) = face_vertices(mesh, face_idx)?;
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    let n = e1.cross(&e2);
    let len = n.norm();
    if len < 1e-12 {
        return None; // Degenerate face.
    }
    Some(n / len)
}

/// Extract the three vertex positions of a face by index.
fn face_vertices(
    mesh: &IndexedMesh<f64>,
    face_idx: u32,
) -> Option<(Point3<f64>, Point3<f64>, Point3<f64>)> {
    let fid = FaceId::from_usize(face_idx as usize);
    if face_idx as usize >= mesh.faces.len() {
        return None;
    }
    let face = mesh.faces.get(fid);
    let [vi0, vi1, vi2] = face.vertices;
    let d0 = mesh.vertices.get(vi0);
    let d1 = mesh.vertices.get(vi1);
    let d2 = mesh.vertices.get(vi2);
    Some((d0.position, d1.position, d2.position))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn point_to_point_distance() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(3.0, 4.0, 0.0);
        assert_relative_eq!(compute_point_to_point(&a, &b), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn edge_length_is_distance() {
        let v0 = Point3::new(1.0, 0.0, 0.0);
        let v1 = Point3::new(1.0, 3.0, 4.0);
        assert_relative_eq!(compute_edge_length(&v0, &v1), 5.0, epsilon = 1e-12);
    }
}
