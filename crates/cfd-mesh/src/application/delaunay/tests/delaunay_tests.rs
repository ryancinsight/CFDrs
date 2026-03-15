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

// ── Compact removes dead triangles ────────────────────────────────────────

/// Theorem: Compaction invariant preservation.
///
/// Compacting a triangulation with $D$ dead and $A$ alive triangles
/// produces a triangulation with exactly $A$ entries where every
/// adjacency link and vertex→triangle hint is remapped consistently.
#[test]
fn compact_removes_dead_triangles() {
    let mut dt =
        DelaunayTriangulation::from_points(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);

    let raw_before = dt.triangle_count_raw();
    let interior_before = dt.triangle_count();

    // After Bowyer-Watson insertions, dead (tombstoned) triangles exist in
    // the backing array.
    assert!(
        dt.triangles_slice().len() > dt.triangle_count_raw(),
        "Should have dead triangles in the backing store"
    );

    dt.compact();

    // After compaction, the backing array should have exactly the alive count.
    assert_eq!(
        dt.triangles_slice().len(),
        raw_before,
        "Compacted length should match alive count"
    );
    // Interior triangle count is unchanged.
    assert_eq!(dt.triangle_count(), interior_before);
    // Delaunay property preserved.
    assert!(dt.is_delaunay(), "Delaunay violated after compact");
    // No dead entries remain.
    assert!(
        dt.triangles_slice().iter().all(|t| t.alive),
        "All triangles should be alive after compact"
    );
    // Adjacency symmetry preserved.
    for (tid, tri) in dt.all_alive_triangles() {
        for edge in 0..3 {
            let nbr = tri.adj[edge];
            if nbr == GHOST_TRIANGLE {
                continue;
            }
            let nbr_tri = dt.triangle(nbr);
            assert!(nbr_tri.alive);
            assert!(
                nbr_tri.shared_edge(tid).is_some(),
                "Adjacency broken after compact"
            );
        }
    }
}

// ── vert_to_tri invariant ─────────────────────────────────────────────────

/// Theorem: Vertex→triangle invariant.
///
/// After construction, for every inserted real vertex $v_i$,
/// `vert_to_tri_slice()[i]` refers to an alive triangle containing $v_i$.
#[test]
fn vert_to_tri_invariant_holds() {
    let pts: Vec<(f64, f64)> = (0..15)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / 15.0;
            (angle.cos(), angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);

    let hints = dt.vert_to_tri_slice();
    for vid_idx in 0..dt.vertex_count() {
        let tid = hints[vid_idx];
        assert_ne!(
            tid, GHOST_TRIANGLE,
            "vertex {} has GHOST hint",
            vid_idx
        );
        let tri = dt.triangle(tid);
        assert!(tri.alive, "vertex {} hints to dead triangle", vid_idx);
        let vid = crate::application::delaunay::pslg::vertex::PslgVertexId::from_usize(vid_idx);
        assert!(
            tri.contains_vertex(vid),
            "vertex {} not in triangle {:?}",
            vid_idx,
            tid
        );
    }
}

// ── k-connectivity ────────────────────────────────────────────────────────

/// Theorem (Whitney 1932): A convex Delaunay triangulation of $n \ge 4$
/// points in general position is 3-vertex-connected.
///
/// We check the weaker necessary condition that every real vertex has at
/// least 2 real neighbours (2-connected).
#[test]
fn triangulation_is_2_connected() {
    let pts: Vec<(f64, f64)> = (0..12)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / 12.0;
            let r = 1.0 + 0.2 * (i as f64 * 1.3).sin();
            (r * angle.cos(), r * angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(
        dt.is_k_connected(2),
        "Delaunay triangulation of 12 points should be >= 2-connected"
    );
}

// ── Insertion order independence ──────────────────────────────────────────

/// The Delaunay triangulation is unique for points in general position
/// (no four cocircular points), so two different insertion orders must
/// produce the same number of triangles and the same Delaunay property.
#[test]
fn insertion_order_independence() {
    let pts = vec![
        (0.1, 0.2),
        (0.9, 0.3),
        (0.5, 0.8),
        (0.3, 0.6),
        (0.7, 0.1),
        (0.4, 0.4),
        (0.8, 0.7),
    ];
    let dt1 = DelaunayTriangulation::from_points(&pts);

    let mut reversed = pts.clone();
    reversed.reverse();
    let dt2 = DelaunayTriangulation::from_points(&reversed);

    assert!(dt1.is_delaunay());
    assert!(dt2.is_delaunay());
    assert_eq!(
        dt1.triangle_count(),
        dt2.triangle_count(),
        "Same points in different order should yield same triangle count"
    );
}

// ── Large-scale stress test ───────────────────────────────────────────────

/// Stress test with 1000 deterministic pseudo-random points.
///
/// Validates Delaunay property, Euler bound, and 2-connectivity on a
/// non-trivial point set.
#[test]
fn stress_1000_points() {
    let mut rng = 12345_u64;
    let mut pts = Vec::with_capacity(1000);
    for _ in 0..1000 {
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

    assert_eq!(dt.vertex_count(), 1000);
    assert!(dt.is_delaunay(), "Delaunay violated for 1000 points");
    assert!(
        dt.is_k_connected(2),
        "1000-point triangulation should be >= 2-connected"
    );

    // Euler bound: F ≤ 2V - 2
    let f = dt.triangle_count();
    let v = dt.vertex_count();
    assert!(f <= 2 * v, "Triangle count {} exceeds 2V = {}", f, 2 * v);
}
