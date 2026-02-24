//! Tests for conversion to IndexedMesh.

use crate::application::delaunay::convert::indexed_mesh::to_indexed_mesh;
use crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation;

// ── IndexedMesh conversion ────────────────────────────────────────────────

#[test]
fn to_indexed_mesh_basic() {
    let dt = DelaunayTriangulation::from_points(&[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]);

    let mesh = to_indexed_mesh(&dt);
    assert_eq!(mesh.vertex_count(), 3);
    assert_eq!(mesh.face_count(), 1);
}

#[test]
fn to_indexed_mesh_square() {
    let dt = DelaunayTriangulation::from_points(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);

    let mesh = to_indexed_mesh(&dt);
    assert_eq!(mesh.vertex_count(), 4);
    assert_eq!(mesh.face_count(), 2);
}

#[test]
fn to_indexed_mesh_preserves_vertex_count() {
    let pts: Vec<(f64, f64)> = (0..10)
        .map(|i| {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / 10.0;
            (angle.cos(), angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    let mesh = to_indexed_mesh(&dt);

    assert_eq!(mesh.vertex_count(), 10);
    assert_eq!(mesh.face_count(), dt.triangle_count());
}
