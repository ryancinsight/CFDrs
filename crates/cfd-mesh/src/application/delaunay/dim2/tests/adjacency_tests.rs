//! Tests for the triangle-triangle adjacency bookkeeping module.
//!
//! `adjacency.rs` provides link, find_edge, edge_opposite_vertex,
//! neighbor, is_hull_triangle, and verify_symmetry.  These tests validate
//! correctness of each operation and the symmetric-adjacency invariant.

use crate::application::delaunay::dim2::pslg::vertex::PslgVertexId;
use crate::application::delaunay::dim2::triangulation::adjacency::Adjacency;
use crate::application::delaunay::dim2::triangulation::bowyer_watson::DelaunayTriangulation;
use crate::application::delaunay::dim2::triangulation::triangle::{Triangle, TriangleId, GHOST_TRIANGLE};

// ── link / link_one ───────────────────────────────────────────────────────

/// Verify that `link(t1, e1, t2, e2)` sets both directions symmetrically.
#[test]
fn link_sets_symmetric_adjacency() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let mut tris = vec![
        Triangle::new(v(0), v(1), v(2)),
        Triangle::new(v(1), v(3), v(2)),
    ];
    let t0 = TriangleId::new(0);
    let t1 = TriangleId::new(1);

    Adjacency::link(&mut tris, t0, 0, t1, 2);

    assert_eq!(tris[0].adj[0], t1);
    assert_eq!(tris[1].adj[2], t0);
}

/// Verify that `link_one` sets only the specified side, leaving the other
/// unchanged (e.g. at GHOST_TRIANGLE).
#[test]
fn link_one_sets_single_direction() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let mut tris = vec![Triangle::new(v(0), v(1), v(2))];
    let t0 = TriangleId::new(0);

    // All adjacencies start as GHOST_TRIANGLE.
    assert_eq!(tris[0].adj[0], GHOST_TRIANGLE);
    assert_eq!(tris[0].adj[1], GHOST_TRIANGLE);
    assert_eq!(tris[0].adj[2], GHOST_TRIANGLE);

    Adjacency::link_one(&mut tris, t0, 1, TriangleId::new(42));

    assert_eq!(tris[0].adj[1], TriangleId::new(42));
    // Other edges unchanged.
    assert_eq!(tris[0].adj[0], GHOST_TRIANGLE);
    assert_eq!(tris[0].adj[2], GHOST_TRIANGLE);
}

// ── find_edge ─────────────────────────────────────────────────────────────

/// `find_edge` returns the correct local edge index.
#[test]
fn find_edge_returns_correct_index() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let mut tris = vec![
        Triangle::new(v(0), v(1), v(2)),
        Triangle::new(v(1), v(3), v(2)),
    ];
    let t0 = TriangleId::new(0);
    let t1 = TriangleId::new(1);

    Adjacency::link(&mut tris, t0, 0, t1, 2);

    assert_eq!(Adjacency::find_edge(&tris, t0, t1), Some(0));
    assert_eq!(Adjacency::find_edge(&tris, t1, t0), Some(2));
}

/// `find_edge` returns `None` for non-adjacent triangles.
#[test]
fn find_edge_returns_none_for_non_adjacent() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let tris = vec![
        Triangle::new(v(0), v(1), v(2)),
        Triangle::new(v(3), v(4), v(5)),
    ];
    let t0 = TriangleId::new(0);
    let t1 = TriangleId::new(1);

    assert_eq!(Adjacency::find_edge(&tris, t0, t1), None);
}

// ── edge_opposite_vertex ──────────────────────────────────────────────────

/// `edge_opposite_vertex` maps each vertex to the expected edge index.
#[test]
fn edge_opposite_vertex_per_vertex() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let tris = vec![Triangle::new(v(0), v(1), v(2))];
    let t0 = TriangleId::new(0);

    // v[0] is opposite edge 0.
    assert_eq!(Adjacency::edge_opposite_vertex(&tris, t0, v(0)), Some(0));
    // v[1] is opposite edge 1.
    assert_eq!(Adjacency::edge_opposite_vertex(&tris, t0, v(1)), Some(1));
    // v[2] is opposite edge 2.
    assert_eq!(Adjacency::edge_opposite_vertex(&tris, t0, v(2)), Some(2));
    // Non-existent vertex.
    assert_eq!(Adjacency::edge_opposite_vertex(&tris, t0, v(99)), None);
}

// ── neighbor ──────────────────────────────────────────────────────────────

/// `neighbor` returns the adjacent triangle across the specified edge.
#[test]
fn neighbor_returns_adjacent_triangle() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let mut tris = vec![
        Triangle::new(v(0), v(1), v(2)),
        Triangle::new(v(1), v(3), v(2)),
    ];
    let t0 = TriangleId::new(0);
    let t1 = TriangleId::new(1);

    Adjacency::link(&mut tris, t0, 0, t1, 2);

    assert_eq!(Adjacency::neighbor(&tris, t0, 0), t1);
    assert_eq!(Adjacency::neighbor(&tris, t1, 2), t0);
    // Unlinked edge → GHOST.
    assert_eq!(Adjacency::neighbor(&tris, t0, 1), GHOST_TRIANGLE);
}

// ── is_hull_triangle ──────────────────────────────────────────────────────

/// A triangle with all adjacencies set is NOT a hull triangle.
#[test]
fn fully_linked_triangle_is_not_hull() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let mut tris = vec![
        Triangle::new(v(0), v(1), v(2)),
        Triangle::new(v(1), v(3), v(2)),
        Triangle::new(v(2), v(3), v(4)),
        Triangle::new(v(0), v(2), v(4)),
    ];
    let t0 = TriangleId::new(0);
    let t1 = TriangleId::new(1);
    let t2 = TriangleId::new(2);
    let t3 = TriangleId::new(3);

    Adjacency::link(&mut tris, t0, 0, t1, 2);
    Adjacency::link(&mut tris, t0, 1, t3, 0);
    Adjacency::link(&mut tris, t0, 2, t2, 1);

    assert!(!Adjacency::is_hull_triangle(&tris, t0));
}

/// A triangle with at least one GHOST_TRIANGLE neighbor IS a hull triangle.
#[test]
fn partial_ghost_is_hull() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let tris = vec![Triangle::new(v(0), v(1), v(2))];
    let t0 = TriangleId::new(0);

    // All adjacencies default to GHOST_TRIANGLE.
    assert!(Adjacency::is_hull_triangle(&tris, t0));
}

// ── verify_symmetry ───────────────────────────────────────────────────────

/// verify_symmetry succeeds on a properly linked pair.
#[test]
fn verify_symmetry_on_valid_pair() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let mut tris = vec![
        Triangle::new(v(0), v(1), v(2)),
        Triangle::new(v(1), v(3), v(2)),
    ];
    let t0 = TriangleId::new(0);
    let t1 = TriangleId::new(1);

    Adjacency::link(&mut tris, t0, 0, t1, 2);

    assert!(Adjacency::verify_symmetry(&tris));
}

/// verify_symmetry detects a one-sided (asymmetric) adjacency.
#[test]
fn verify_symmetry_detects_asymmetry() {
    let v = |i: usize| PslgVertexId::from_usize(i);
    let mut tris = vec![
        Triangle::new(v(0), v(1), v(2)),
        Triangle::new(v(1), v(3), v(2)),
    ];
    let t0 = TriangleId::new(0);
    let t1 = TriangleId::new(1);

    // Deliberately set only one direction.
    Adjacency::link_one(&mut tris, t0, 0, t1);

    assert!(!Adjacency::verify_symmetry(&tris));
}

/// verify_symmetry holds on a real Delaunay triangulation.
#[test]
fn verify_symmetry_on_triangulation() {
    let dt = DelaunayTriangulation::from_points(&[
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (0.5, 0.5),
    ]);
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── Integration: adjacency after BW insertion ─────────────────────────────

/// After incremental insertion, every alive triangle has symmetric adjacency
/// and no GHOST_TRIANGLE neighbor points to a dead triangle.
#[test]
fn adjacency_integrity_after_incremental_insertion() {
    let pts: Vec<(f64, f64)> = (0..50)
        .map(|i| {
            let angle = 2.0 * std::f64::consts::PI * f64::from(i) / 50.0;
            (angle.cos(), angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);

    // Symmetric adjacency.
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));

    // Every neighbor of a live triangle is either GHOST or another live triangle.
    for (tid, tri) in dt.all_alive_triangles() {
        for e in 0..3 {
            let nbr = tri.adj[e];
            if nbr != GHOST_TRIANGLE {
                assert!(
                    dt.triangle(nbr).alive,
                    "triangle {tid} edge {e} points to dead neighbor {nbr}"
                );
            }
        }
    }
}

/// After inserting on a regular grid (co-circular stress), adjacency is still
/// symmetric.
#[test]
fn adjacency_integrity_grid_stress() {
    let mut pts = Vec::new();
    for i in 0..15 {
        for j in 0..15 {
            pts.push((f64::from(i), f64::from(j)));
        }
    }
    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    assert!(dt.is_delaunay());
}
