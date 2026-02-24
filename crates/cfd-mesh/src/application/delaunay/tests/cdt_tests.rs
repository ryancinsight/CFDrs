//! Tests for the Constrained Delaunay Triangulation (CDT).

use crate::application::delaunay::constraint::enforce::Cdt;
use crate::application::delaunay::pslg::graph::Pslg;
use crate::application::delaunay::pslg::graph::PslgValidationError;

// ── Simple triangle constraint ────────────────────────────────────────────

#[test]
fn triangle_cdt() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 0.0);
    let c = pslg.add_vertex(0.5, 1.0);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    assert!(dt.is_delaunay(), "CDT should satisfy Delaunay property");
    assert!(
        dt.triangle_count() >= 1,
        "CDT should have at least 1 triangle"
    );
}

// ── Rectangle with diagonal constraint ────────────────────────────────────

#[test]
fn rectangle_with_diagonal() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(2.0, 0.0);
    let c = pslg.add_vertex(2.0, 1.0);
    let d = pslg.add_vertex(0.0, 1.0);

    // Boundary segments.
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, d);
    pslg.add_segment(d, a);

    // Diagonal constraint.
    pslg.add_segment(a, c);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    // The diagonal (a, c) must exist as an edge.
    assert!(
        cdt.is_constrained(a, c),
        "Diagonal constraint should be enforced"
    );
    assert_eq!(dt.triangle_count(), 2, "Rectangle should have 2 triangles");
}

// ── Constrained edges are preserved ───────────────────────────────────────

#[test]
fn constraint_edges_present_in_triangulation() {
    let mut pslg = Pslg::new();
    let v0 = pslg.add_vertex(0.0, 0.0);
    let v1 = pslg.add_vertex(3.0, 0.0);
    let v2 = pslg.add_vertex(3.0, 3.0);
    let v3 = pslg.add_vertex(0.0, 3.0);
    let _v4 = pslg.add_vertex(1.5, 1.5); // interior point

    pslg.add_segment(v0, v1);
    pslg.add_segment(v1, v2);
    pslg.add_segment(v2, v3);
    pslg.add_segment(v3, v0);

    let cdt = Cdt::from_pslg(&pslg);

    // All boundary edges should be constrained.
    assert!(cdt.is_constrained(v0, v1));
    assert!(cdt.is_constrained(v1, v2));
    assert!(cdt.is_constrained(v2, v3));
    assert!(cdt.is_constrained(v3, v0));
}

// ── Hole removal ──────────────────────────────────────────────────────────

#[test]
fn square_with_hole() {
    let mut pslg = Pslg::new();

    // Outer square.
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(4.0, 0.0);
    let c = pslg.add_vertex(4.0, 4.0);
    let d = pslg.add_vertex(0.0, 4.0);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, d);
    pslg.add_segment(d, a);

    // Inner square (hole boundary).
    let e = pslg.add_vertex(1.0, 1.0);
    let f = pslg.add_vertex(3.0, 1.0);
    let g = pslg.add_vertex(3.0, 3.0);
    let h = pslg.add_vertex(1.0, 3.0);
    pslg.add_segment(e, f);
    pslg.add_segment(f, g);
    pslg.add_segment(g, h);
    pslg.add_segment(h, e);

    // Hole seed inside the inner square.
    pslg.add_hole(2.0, 2.0);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    // The hole region should have no triangles.
    // Every surviving triangle's centroid should be outside the inner square.
    for (_, tri) in dt.interior_triangles() {
        let v0 = dt.vertex(tri.vertices[0]);
        let v1 = dt.vertex(tri.vertices[1]);
        let v2 = dt.vertex(tri.vertices[2]);

        let cx = (v0.x + v1.x + v2.x) / 3.0;
        let cy = (v0.y + v1.y + v2.y) / 3.0;

        let in_hole = cx > 1.0 && cx < 3.0 && cy > 1.0 && cy < 3.0;
        assert!(
            !in_hole,
            "Triangle centroid ({}, {}) is inside the hole region",
            cx, cy
        );
    }
}

// ── CDT with many interior points ─────────────────────────────────────────

#[test]
fn cdt_with_interior_points() {
    let mut pslg = Pslg::new();

    // Outer boundary.
    let v0 = pslg.add_vertex(0.0, 0.0);
    let v1 = pslg.add_vertex(5.0, 0.0);
    let v2 = pslg.add_vertex(5.0, 5.0);
    let v3 = pslg.add_vertex(0.0, 5.0);
    pslg.add_segment(v0, v1);
    pslg.add_segment(v1, v2);
    pslg.add_segment(v2, v3);
    pslg.add_segment(v3, v0);

    // Interior points on a grid.
    for i in 1..5 {
        for j in 1..5 {
            pslg.add_vertex(i as f64, j as f64);
        }
    }

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    assert!(
        dt.is_delaunay(),
        "CDT with interior points should be Delaunay"
    );
    assert!(dt.triangle_count() > 0, "CDT should produce triangles");
}

#[test]
fn try_from_pslg_rejects_crossing_segments() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 1.0);
    let c = pslg.add_vertex(0.0, 1.0);
    let d = pslg.add_vertex(1.0, 0.0);
    pslg.add_segment(a, b);
    pslg.add_segment(c, d);

    let err = Cdt::try_from_pslg(&pslg)
        .err()
        .expect("crossing PSLG should be rejected");
    assert!(matches!(
        err,
        PslgValidationError::IntersectingSegments { .. }
    ));
}
