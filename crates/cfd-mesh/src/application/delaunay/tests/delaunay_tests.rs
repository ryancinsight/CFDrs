//! Tests for the core Bowyer-Watson Delaunay triangulation.

use crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation;
use crate::application::delaunay::triangulation::triangle::GHOST_TRIANGLE;
use std::f64::consts::PI;

// ── Basic construction ────────────────────────────────────────────────────

#[test]
fn single_triangle() {
    let dt = DelaunayTriangulation::from_points(&[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]);
    assert_eq!(dt.vertex_count(), 3);
    assert_eq!(dt.triangle_count(), 1);
    assert!(dt.is_delaunay());
}

#[test]
fn four_points_square() {
    let dt = DelaunayTriangulation::from_points(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
    assert_eq!(dt.vertex_count(), 4);
    assert_eq!(dt.triangle_count(), 2);
    assert!(dt.is_delaunay());
}

#[test]
fn five_points_pentagon() {
    let pts: Vec<(f64, f64)> = (0..5)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / 5.0;
            (angle.cos(), angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 5);
    assert!(dt.triangle_count() >= 3);
    assert!(dt.is_delaunay());
}

// ── Delaunay property verification ────────────────────────────────────────

/// Theorem: Empty-circumcircle property.
/// For every interior triangle, no other vertex lies strictly inside its
/// circumcircle.
#[test]
fn delaunay_empty_circumcircle_random_20() {
    let points: Vec<(f64, f64)> = [
        (0.1, 0.2),
        (0.9, 0.3),
        (0.5, 0.8),
        (0.2, 0.7),
        (0.8, 0.9),
        (0.15, 0.5),
        (0.85, 0.1),
        (0.4, 0.6),
        (0.6, 0.4),
        (0.3, 0.15),
        (0.7, 0.65),
        (0.45, 0.35),
        (0.55, 0.75),
        (0.25, 0.45),
        (0.75, 0.55),
        (0.35, 0.85),
        (0.65, 0.25),
        (0.05, 0.95),
        (0.95, 0.05),
        (0.5, 0.5),
    ]
    .into();
    let dt = DelaunayTriangulation::from_points(&points);
    assert!(dt.is_delaunay(), "Empty-circumcircle property violated");
}

#[test]
fn delaunay_grid_4x4() {
    let mut pts = Vec::new();
    for i in 0..4 {
        for j in 0..4 {
            pts.push((i as f64, j as f64));
        }
    }
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 16);
    assert!(dt.is_delaunay());
}

// ── Collinear points ──────────────────────────────────────────────────────

#[test]
fn collinear_points() {
    let dt = DelaunayTriangulation::from_points(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]);
    // Collinear — degenerate, but should not panic.
    assert_eq!(dt.vertex_count(), 4);
}

// ── Duplicate points ──────────────────────────────────────────────────────

#[test]
fn duplicate_points_skipped() {
    let dt = DelaunayTriangulation::from_points(&[
        (0.0, 0.0),
        (1.0, 0.0),
        (0.5, 1.0),
        (0.0, 0.0), // duplicate
    ]);
    // Should handle gracefully (duplicate is skipped).
    assert!(dt.triangle_count() >= 1);
}

// ── Adjacency invariant ──────────────────────────────────────────────────

/// Theorem: Adjacency symmetry.
/// For every triangle t with neighbor n via edge e, triangle n has t as
/// its neighbor via the shared edge (reversal).
#[test]
fn adjacency_symmetry() {
    let pts: Vec<(f64, f64)> = (0..10)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / 10.0;
            (angle.cos(), angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);

    for (tid, tri) in dt.all_alive_triangles() {
        for edge in 0..3 {
            let nbr = tri.adj[edge];
            if nbr == GHOST_TRIANGLE {
                continue;
            }
            let nbr_tri = dt.triangle(nbr);
            assert!(nbr_tri.alive, "Neighbor {:?} of {:?} is dead", nbr, tid);
            let back = nbr_tri.shared_edge(tid);
            assert!(
                back.is_some(),
                "Adjacency asymmetry: {:?} → {:?} but {:?} does not point back",
                tid,
                nbr,
                nbr
            );
        }
    }
}

// ── Triangle count bounds ─────────────────────────────────────────────────

/// For n non-collinear points in general position: F ≈ 2n - 2 - h
/// where h is the number of convex-hull vertices.
#[test]
fn triangle_count_bound() {
    let pts: Vec<(f64, f64)> = (0..20)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / 20.0;
            let r = 1.0 + 0.3 * (i as f64 * 0.7).sin();
            (r * angle.cos(), r * angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    let n = dt.vertex_count();
    let f = dt.triangle_count();

    // Euler bound: F ≤ 2n - 5 (for n ≥ 3 with no collinear)
    assert!(
        f <= 2 * n,
        "Triangle count {} exceeds 2n = {} for n = {}",
        f,
        2 * n,
        n
    );
}

// ── Larger stress test ────────────────────────────────────────────────────

#[test]
fn stress_100_points() {
    // Deterministic pseudo-random via LCG.
    let mut rng = 42_u64;
    let mut pts = Vec::with_capacity(100);
    for _ in 0..100 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let x = (rng >> 33) as f64 / (1u64 << 31) as f64;
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let y = (rng >> 33) as f64 / (1u64 << 31) as f64;
        pts.push((x, y));
    }
    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(
        dt.is_delaunay(),
        "Delaunay violated for 100-point stress test"
    );
    assert_eq!(dt.vertex_count(), 100);
}

// ── Interior triangle iteration ───────────────────────────────────────────

#[test]
fn interior_triangles_exclude_super() {
    let dt = DelaunayTriangulation::from_points(&[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]);
    let interior_count = dt.interior_triangles().count();
    let raw_count = dt.triangle_count_raw();
    // Interior count should be less than raw (some touch super-triangle)
    assert!(
        interior_count <= raw_count,
        "Interior {} > raw {}",
        interior_count,
        raw_count
    );
    assert_eq!(interior_count, dt.triangle_count());
}
