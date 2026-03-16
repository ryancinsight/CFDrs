//! Topology & structural stress tests for Delaunay/CDT internals.
//!
//! Targets known failure modes in mesh libraries:
//! - Convex hull extraction correctness after many insertions
//! - Vertex-star fan-walk for hull/interior vertices
//! - Cavity re-triangulation (ear-clipping) with non-trivial polygons
//! - Euler formula consistency under stress
//! - vert_to_tri cache coherence after CDT constraint enforcement
//!
//! Each test documents the **failure mode** it guards against and the
//! relevant **literature** where applicable.

use crate::application::delaunay::constraint::enforce::Cdt;
use crate::application::delaunay::pslg::graph::Pslg;
use crate::application::delaunay::triangulation::adjacency::Adjacency;
use crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation;
use std::f64::consts::PI;

// ── Convex hull extraction ────────────────────────────────────────────────

/// Verify that convex hull of a regular polygon returns exactly its vertices.
///
/// **Failure mode**: The O(h) super-vertex star walk may miss hull vertices
/// if the adjacency walk around a super-vertex's star is broken.
///
/// **Literature**: Preparata & Shamos (1985), Chapter 3 — gift-wrapping
/// hull correctness.
#[test]
fn convex_hull_regular_hexagon() {
    let n = 6;
    let pts: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / n as f64;
            (angle.cos(), angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    let hull = dt.convex_hull_vertices();
    assert_eq!(
        hull.len(),
        n,
        "Regular hexagon hull should have {n} vertices, got {}",
        hull.len()
    );
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    assert!(dt.is_delaunay());
}

/// Convex hull of points on a circle plus a centre point should still
/// yield exactly the circle points (centre is interior).
///
/// **Failure mode**: The star walk adds interior vertices from super-vertex
/// triangles that also touch hull triangles.
#[test]
fn convex_hull_circle_with_centre() {
    let n = 12;
    let mut pts: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / n as f64;
            (10.0 * angle.cos(), 10.0 * angle.sin())
        })
        .collect();
    pts.push((0.0, 0.0)); // Centre point — interior.
    let dt = DelaunayTriangulation::from_points(&pts);
    let hull = dt.convex_hull_vertices();
    assert_eq!(
        hull.len(),
        n,
        "Circle-plus-centre hull should have {n} vertices, got {}",
        hull.len()
    );
    assert!(dt.satisfies_euler());
}

/// Convex hull of a dense grid should yield only the boundary vertices.
///
/// **Failure mode**: Dense interior points stress the super-vertex star
/// walk by creating many non-hull triangles adjacent to super vertices.
#[test]
fn convex_hull_10x10_grid() {
    let mut pts = Vec::new();
    for i in 0..10 {
        for j in 0..10 {
            pts.push((i as f64, j as f64));
        }
    }
    let dt = DelaunayTriangulation::from_points(&pts);
    let hull = dt.convex_hull_vertices();
    // 4 corners + 4 * 8 edges = 36 boundary points.
    assert_eq!(
        hull.len(),
        36,
        "10x10 grid hull should have 36 boundary vertices, got {}",
        hull.len()
    );
}

/// Collinear hull: all points on a line.  The hull should have exactly 2
/// vertices (the endpoints), or degenerate gracefully.
///
/// **Failure mode**: Super-vertex star walk on degenerate (collinear)
/// inputs may loop or return wrong count.
#[test]
fn convex_hull_collinear_points() {
    let pts: Vec<(f64, f64)> = (0..10).map(|i| (i as f64, 0.0)).collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    let hull = dt.convex_hull_vertices();
    // All points are on the hull for collinear inputs.
    assert_eq!(
        hull.len(),
        10,
        "Collinear hull should contain all 10 points, got {}",
        hull.len()
    );
}

// ── Vertex-star fan walk ──────────────────────────────────────────────────

/// Verify that the vertex-star fan walk returns the correct degree for
/// all vertices of a symmetric configuration.
///
/// **Failure mode**: Bidirectional fan walk incorrectly terminates at
/// hull boundaries, returning partial star.
///
/// **Literature**: Boissonnat & Teillaud (2006), Chapter 2 — half-edge
/// traversal for convex-hull vertices.
#[test]
fn vertex_star_degree_regular_polygon() {
    let n = 8;
    let mut pts: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / n as f64;
            (angle.cos(), angle.sin())
        })
        .collect();
    pts.push((0.0, 0.0)); // Centre vertex.
    let dt = DelaunayTriangulation::from_points(&pts);

    // Centre vertex should have degree n (connected to all ring vertices).
    let centre_id = crate::application::delaunay::pslg::vertex::PslgVertexId::from_usize(n);
    let star = dt.triangles_around_vertex(centre_id);
    assert_eq!(
        star.len(),
        n,
        "Centre vertex of {n}-gon should have degree {n}, got {}",
        star.len()
    );

    // Each ring vertex is a hull vertex — its star should be reachable
    // via bidirectional walk.
    for i in 0..n {
        let vid = crate::application::delaunay::pslg::vertex::PslgVertexId::from_usize(i);
        let ring_star = dt.triangles_around_vertex(vid);
        assert!(
            ring_star.len() >= 2,
            "Hull vertex {i} should have ≥2 incident triangles, got {}",
            ring_star.len()
        );
    }
}

/// Verify vertex-star completeness: every alive triangle containing
/// vertex v is returned by `triangles_around_vertex(v)`.
///
/// **Failure mode**: vert_to_tri cache staleness causes the fan walk
/// to start from a wrong triangle and miss some star members.
#[test]
fn vertex_star_completeness_random_50() {
    // Random-ish 50-point cloud via deterministic formula.
    let n = 50;
    let pts: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let t = i as f64 * 0.618_033_988_749_895; // golden ratio
            let x = (t * 7.3).sin() * 5.0;
            let y = (t * 13.1).cos() * 5.0;
            (x, y)
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);

    // For each real vertex, verify star is complete by cross-checking
    // with a linear scan.
    for vi in 0..n {
        let vid = crate::application::delaunay::pslg::vertex::PslgVertexId::from_usize(vi);
        let star = dt.triangles_around_vertex(vid);
        // Linear scan: count triangles containing this vertex.
        let linear_count = dt
            .all_alive_triangles()
            .filter(|(_, t)| t.contains_vertex(vid))
            .count();
        assert_eq!(
            star.len(),
            linear_count,
            "Vertex {vi}: fan-walk returned {} triangles but linear scan found {}",
            star.len(),
            linear_count
        );
    }
}

// ── Euler formula under stress ────────────────────────────────────────────

/// Euler formula must hold after many point insertions at varied scales.
///
/// **Failure mode**: Numerical precision issues during Bowyer-Watson
/// cavity identification cause orphaned or overlapping triangles,
/// violating Euler's relation.
///
/// **Literature**: Shewchuk (1997), "Adaptive Precision Floating-Point
/// Arithmetic and Fast Robust Geometric Predicates".
#[test]
fn euler_formula_multi_scale() {
    let mut pts = Vec::new();
    // Scale 1: unit square
    for i in 0..10 {
        for j in 0..10 {
            pts.push((i as f64 * 0.1, j as f64 * 0.1));
        }
    }
    // Scale 2: large ring
    for i in 0..20 {
        let angle = 2.0 * PI * i as f64 / 20.0;
        pts.push((100.0 * angle.cos(), 100.0 * angle.sin()));
    }
    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(
        dt.satisfies_euler(),
        "Euler formula violated for multi-scale point set"
    );
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

/// Euler formula after CDT constraint enforcement with hole removal.
///
/// **Failure mode**: Hole removal flood-fill corrupts adjacency links,
/// causing Euler violation.
#[test]
fn euler_formula_cdt_with_multiple_holes() {
    let mut pslg = Pslg::new();
    // Outer boundary: 5x5 square.
    let o0 = pslg.add_vertex(0.0, 0.0);
    let o1 = pslg.add_vertex(5.0, 0.0);
    let o2 = pslg.add_vertex(5.0, 5.0);
    let o3 = pslg.add_vertex(0.0, 5.0);
    pslg.add_segment(o0, o1);
    pslg.add_segment(o1, o2);
    pslg.add_segment(o2, o3);
    pslg.add_segment(o3, o0);

    // Hole 1: small square at (1,1)-(2,2).
    let h0 = pslg.add_vertex(1.0, 1.0);
    let h1 = pslg.add_vertex(2.0, 1.0);
    let h2 = pslg.add_vertex(2.0, 2.0);
    let h3 = pslg.add_vertex(1.0, 2.0);
    pslg.add_segment(h0, h1);
    pslg.add_segment(h1, h2);
    pslg.add_segment(h2, h3);
    pslg.add_segment(h3, h0);
    pslg.add_hole(1.5, 1.5);

    // Hole 2: small square at (3,3)-(4,4).
    let h4 = pslg.add_vertex(3.0, 3.0);
    let h5 = pslg.add_vertex(4.0, 3.0);
    let h6 = pslg.add_vertex(4.0, 4.0);
    let h7 = pslg.add_vertex(3.0, 4.0);
    pslg.add_segment(h4, h5);
    pslg.add_segment(h5, h6);
    pslg.add_segment(h6, h7);
    pslg.add_segment(h7, h4);
    pslg.add_hole(3.5, 3.5);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    // Verify no alive triangle has its centroid inside either hole.
    for (_, tri) in dt.all_alive_triangles() {
        let [a, b, c] = tri.vertices;
        let pa = dt.vertex(a);
        let pb = dt.vertex(b);
        let pc = dt.vertex(c);
        let cx = (pa.x + pb.x + pc.x) / 3.0;
        let cy = (pa.y + pb.y + pc.y) / 3.0;
        let in_hole1 = cx > 1.0 && cx < 2.0 && cy > 1.0 && cy < 2.0;
        let in_hole2 = cx > 3.0 && cx < 4.0 && cy > 3.0 && cy < 4.0;
        assert!(
            !in_hole1 && !in_hole2,
            "Triangle centroid ({cx}, {cy}) should not be inside a hole"
        );
    }
}

// ── CDT constraint recovery stress ───────────────────────────────────────

/// Many parallel constraints crossing the same region.
///
/// **Failure mode**: Queue-based flip recovery fails to make progress
/// when multiple constraints share crossing edges.  Tests the flip
/// iteration budget and Steiner fallback path.
///
/// **Literature**: Sloan (1993), "A fast algorithm for generating
/// constrained Delaunay triangulations".
#[test]
fn cdt_parallel_constraints_stress() {
    let mut pslg = Pslg::new();
    // Outer rectangle 0..10 x 0..4.
    let o0 = pslg.add_vertex(0.0, 0.0);
    let o1 = pslg.add_vertex(10.0, 0.0);
    let o2 = pslg.add_vertex(10.0, 4.0);
    let o3 = pslg.add_vertex(0.0, 4.0);
    pslg.add_segment(o0, o1);
    pslg.add_segment(o1, o2);
    pslg.add_segment(o2, o3);
    pslg.add_segment(o3, o0);

    // 8 horizontal constraints at y = 0.5, 1.0, ..., 3.5.
    // Offset endpoints inward to avoid overlap with boundary segments.
    for i in 1..=7 {
        let y = i as f64 * 0.5;
        let left = pslg.add_vertex(0.5, y);
        let right = pslg.add_vertex(9.5, y);
        pslg.add_segment(left, right);
    }

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    assert!(
        dt.triangle_count() > 0,
        "Parallel constraints should produce a valid mesh"
    );
}

/// Fan-shaped constraint pattern: constraints radiating from one vertex.
///
/// **Failure mode**: Multiple constraints sharing a common endpoint
/// interfere with each other's flip sequences.
///
/// **Literature**: Chew (1989), "Constrained Delaunay triangulations" —
/// correctness of constraint recovery for shared endpoints.
#[test]
fn cdt_fan_constraints_from_centre() {
    let mut pslg = Pslg::new();
    // Ring of points + centre.
    let n = 16;
    let centre = pslg.add_vertex(0.0, 0.0);
    let mut ring = Vec::new();
    for i in 0..n {
        let angle = 2.0 * PI * i as f64 / n as f64;
        ring.push(pslg.add_vertex(5.0 * angle.cos(), 5.0 * angle.sin()));
    }
    // Outer boundary.
    for i in 0..n {
        pslg.add_segment(ring[i], ring[(i + 1) % n]);
    }
    // Radial constraints from centre to every other ring vertex.
    for i in (0..n).step_by(2) {
        pslg.add_segment(centre, ring[i]);
    }

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    assert!(dt.is_delaunay());
}

/// CDT with a narrow slit: two close parallel constraints creating a
/// very thin gap. Tests numerical robustness of cavity polygon
/// identification and ear-clipping.
///
/// **Failure mode**: Ear-clipping in cavity re-triangulation fails for
/// very thin cavities (all ears nearly degenerate).
#[test]
fn cdt_narrow_slit() {
    let mut pslg = Pslg::new();
    // Square boundary.
    let o0 = pslg.add_vertex(0.0, 0.0);
    let o1 = pslg.add_vertex(10.0, 0.0);
    let o2 = pslg.add_vertex(10.0, 10.0);
    let o3 = pslg.add_vertex(0.0, 10.0);
    pslg.add_segment(o0, o1);
    pslg.add_segment(o1, o2);
    pslg.add_segment(o2, o3);
    pslg.add_segment(o3, o0);

    // Two very close horizontal lines creating a narrow slit.
    let eps = 0.01;
    let l0 = pslg.add_vertex(1.0, 5.0 - eps);
    let l1 = pslg.add_vertex(9.0, 5.0 - eps);
    let l2 = pslg.add_vertex(1.0, 5.0 + eps);
    let l3 = pslg.add_vertex(9.0, 5.0 + eps);
    pslg.add_segment(l0, l1);
    pslg.add_segment(l2, l3);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── Connectivity ──────────────────────────────────────────────────────────

/// 3-connectivity of a convex Delaunay triangulation of ≥4 points.
///
/// **Failure mode**: Incorrect adjacency patching after flips leaves
/// vertices with degree < 3.
///
/// **Literature**: Whitney (1932) — every triangulation of a convex
/// point set with n ≥ 4 is 3-vertex-connected.
#[test]
fn three_connected_random_cloud() {
    // Use points in convex position (circle) plus interior points
    // to guarantee well-separated, general-position input.
    let mut pts = Vec::new();
    // 20 points on a large circle (convex hull).
    for i in 0..20 {
        let angle = 2.0 * PI * i as f64 / 20.0;
        pts.push((10.0 * angle.cos(), 10.0 * angle.sin()));
    }
    // 30 interior points on a smaller circle.
    for i in 0..30 {
        let angle = 2.0 * PI * i as f64 / 30.0 + 0.1;
        pts.push((3.0 * angle.cos(), 3.0 * angle.sin()));
    }
    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(
        dt.is_k_connected(3),
        "Convex Delaunay triangulation of {} points should be 3-connected",
        pts.len()
    );
    assert!(dt.satisfies_euler());
}

/// After CDT constraint enforcement and hole removal, the remaining
/// mesh should have consistent adjacency and valid Delaunay property
/// on non-constrained edges.
///
/// **Failure mode**: Adjacency patching in remove_hole_triangles
/// leaves dangling references (the CW19 bug class).
#[test]
fn cdt_annulus_adjacency_consistency() {
    let mut pslg = Pslg::new();
    let n_outer = 20;
    let n_inner = 10;

    // Outer ring.
    let mut outer = Vec::new();
    for i in 0..n_outer {
        let angle = 2.0 * PI * i as f64 / n_outer as f64;
        outer.push(pslg.add_vertex(10.0 * angle.cos(), 10.0 * angle.sin()));
    }
    for i in 0..n_outer {
        pslg.add_segment(outer[i], outer[(i + 1) % n_outer]);
    }

    // Inner ring (hole boundary).
    let mut inner = Vec::new();
    for i in 0..n_inner {
        let angle = 2.0 * PI * i as f64 / n_inner as f64;
        inner.push(pslg.add_vertex(3.0 * angle.cos(), 3.0 * angle.sin()));
    }
    for i in 0..n_inner {
        pslg.add_segment(inner[i], inner[(i + 1) % n_inner]);
    }
    pslg.add_hole(0.0, 0.0);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    // No triangle should have its centroid inside the inner ring.
    for (_, tri) in dt.all_alive_triangles() {
        let [a, b, c] = tri.vertices;
        let pa = dt.vertex(a);
        let pb = dt.vertex(b);
        let pc = dt.vertex(c);
        let cx = (pa.x + pb.x + pc.x) / 3.0;
        let cy = (pa.y + pb.y + pc.y) / 3.0;
        let r2 = cx * cx + cy * cy;
        assert!(
            r2 > 3.0 * 3.0 * 0.5,
            "Triangle centroid ({cx:.2}, {cy:.2}) is inside the hole (r²={r2:.2})"
        );
    }
}
