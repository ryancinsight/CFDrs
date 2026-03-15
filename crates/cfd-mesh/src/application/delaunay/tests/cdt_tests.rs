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

// ── Compact preserves constraints ─────────────────────────────────────────

/// Theorem: Compaction invariant preservation (CDT extension).
///
/// After compacting the underlying triangulation, every constrained edge
/// recorded in the CDT must still be present and marked constrained.
#[test]
fn cdt_compact_preserves_constraints() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(2.0, 0.0);
    let c = pslg.add_vertex(2.0, 2.0);
    let d = pslg.add_vertex(0.0, 2.0);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, d);
    pslg.add_segment(d, a);
    // Diagonal constraint.
    pslg.add_segment(a, c);

    let mut cdt = Cdt::from_pslg(&pslg);

    // Compact the underlying triangulation.
    cdt.triangulation_mut().compact();

    // All constraints must still be present and marked constrained.
    assert!(cdt.is_constrained(a, b), "bottom edge lost after compact");
    assert!(cdt.is_constrained(b, c), "right edge lost after compact");
    assert!(cdt.is_constrained(c, d), "top edge lost after compact");
    assert!(cdt.is_constrained(d, a), "left edge lost after compact");
    assert!(cdt.is_constrained(a, c), "diagonal lost after compact");

    // Delaunay property still holds.
    assert!(
        cdt.triangulation().is_delaunay(),
        "Delaunay violated after CDT compact"
    );
}

// ── Constraint recovery on internal segments ──────────────────────────────

/// Inserting a constraint that cuts across existing Delaunay edges forces
/// the CDT to flip or re-triangulate the affected region.  The constraint
/// must appear as an edge in the result.
#[test]
fn cdt_constraint_recovery_internal_segments() {
    let mut pslg = Pslg::new();
    // Unit square with interior point at (0.5, 0.5).
    let v0 = pslg.add_vertex(0.0, 0.0);
    let v1 = pslg.add_vertex(1.0, 0.0);
    let v2 = pslg.add_vertex(1.0, 1.0);
    let v3 = pslg.add_vertex(0.0, 1.0);
    let _v4 = pslg.add_vertex(0.5, 0.5);

    // Boundary.
    pslg.add_segment(v0, v1);
    pslg.add_segment(v1, v2);
    pslg.add_segment(v2, v3);
    pslg.add_segment(v3, v0);

    // Internal constraint cutting the domain.
    pslg.add_segment(v0, v2);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    // The diagonal (v0, v2) must exist and be constrained.
    assert!(
        cdt.is_constrained(v0, v2),
        "Internal diagonal constraint not recovered"
    );
    // Interior point should be present.
    assert!(dt.vertex_count() >= 5);
    // CDT should still be Delaunay (modulo constrained edges).
    assert!(dt.is_delaunay());
}

// ── Complex non-convex domain ─────────────────────────────────────────────

/// L-shaped non-convex domain — CDT must correctly triangulate both arms.
#[test]
fn cdt_l_shaped_domain() {
    let mut pslg = Pslg::new();
    //   (0,0)──(3,0)
    //     |      |
    //   (0,2)──(1,2)
    //           |
    //   (1,3)──(3,3)
    // Wait — L-shape:
    //   (0,0)→(2,0)→(2,1)→(1,1)→(1,2)→(0,2)→(0,0)
    let v0 = pslg.add_vertex(0.0, 0.0);
    let v1 = pslg.add_vertex(2.0, 0.0);
    let v2 = pslg.add_vertex(2.0, 1.0);
    let v3 = pslg.add_vertex(1.0, 1.0);
    let v4 = pslg.add_vertex(1.0, 2.0);
    let v5 = pslg.add_vertex(0.0, 2.0);
    pslg.add_segment(v0, v1);
    pslg.add_segment(v1, v2);
    pslg.add_segment(v2, v3);
    pslg.add_segment(v3, v4);
    pslg.add_segment(v4, v5);
    pslg.add_segment(v5, v0);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    // All boundary segments exist.
    assert!(cdt.is_constrained(v0, v1));
    assert!(cdt.is_constrained(v1, v2));
    assert!(cdt.is_constrained(v2, v3));
    assert!(cdt.is_constrained(v3, v4));
    assert!(cdt.is_constrained(v4, v5));
    assert!(cdt.is_constrained(v5, v0));

    // Must have at least 4 triangles for a 6-vertex L-shape.
    assert!(
        dt.triangle_count() >= 4,
        "L-shape should have >= 4 triangles, got {}",
        dt.triangle_count()
    );
    assert!(dt.is_delaunay());
}

// ── Constraint bidirectionality ───────────────────────────────────────────

/// Constrained status must be symmetric: `is_constrained(a, b) == is_constrained(b, a)`.
#[test]
fn constraint_edges_bidirectional() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 0.0);
    let c = pslg.add_vertex(0.5, 1.0);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    let cdt = Cdt::from_pslg(&pslg);

    assert_eq!(cdt.is_constrained(a, b), cdt.is_constrained(b, a));
    assert_eq!(cdt.is_constrained(b, c), cdt.is_constrained(c, b));
    assert_eq!(cdt.is_constrained(c, a), cdt.is_constrained(a, c));
}

// ── Scale-extreme CDT: micro scale ────────────────────────────────────────

/// Theorem: CDT correctness is scale-invariant.
///
/// Statement: A CDT built from coordinates uniformly scaled by a factor
/// $\alpha > 0$ produces the same combinatorial triangulation as the
/// un-scaled input.  In particular, the Delaunay property is preserved.
///
/// Proof sketch: The in-circle and orientation predicates are homogeneous
/// polynomials in the input coordinates.  Scaling all coordinates by $\alpha$
/// multiplies the predicate value by $\alpha^k$ for some $k > 0$, preserving
/// the sign.  Therefore scaling does not change any combinatorial decision.
/// QED.
#[test]
fn micro_scale_cdt() {
    let scale = 1e-15;
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0 * scale, 0.0 * scale);
    let b = pslg.add_vertex(1.0 * scale, 0.0 * scale);
    let c = pslg.add_vertex(0.5 * scale, 0.866 * scale);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    assert!(dt.is_delaunay(), "Micro-scale CDT should be Delaunay");
    assert!(
        dt.triangle_count() >= 1,
        "Micro-scale CDT should produce triangles"
    );
    assert!(cdt.is_constrained(a, b));
    assert!(cdt.is_constrained(b, c));
    assert!(cdt.is_constrained(c, a));
}

// ── Scale-extreme CDT: macro scale ────────────────────────────────────────

#[test]
fn macro_scale_cdt() {
    let scale = 1e12;
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0 * scale, 0.0 * scale);
    let b = pslg.add_vertex(1.0 * scale, 0.0 * scale);
    let c = pslg.add_vertex(0.5 * scale, 0.866 * scale);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    assert!(dt.is_delaunay(), "Macro-scale CDT should be Delaunay");
    assert!(
        dt.triangle_count() >= 1,
        "Macro-scale CDT should produce triangles"
    );
    assert!(cdt.is_constrained(a, b));
    assert!(cdt.is_constrained(b, c));
    assert!(cdt.is_constrained(c, a));
}

// ── Constraint forces Delaunay restoration ────────────────────────────────

/// Theorem: CDT Locally-Delaunay Property (Chew 1989, Lee & Lin 1986).
///
/// Statement: In a constrained Delaunay triangulation, every non-constrained
/// edge satisfies the empty-circumcircle (in-circle) criterion.
///
/// This test builds a square with 4 interior points and a diagonal
/// constraint.  The constraint cuts through the unconstrained Delaunay
/// triangulation, and `restore_delaunay_near_constraint` must flip the
/// surrounding edges so all non-constrained edges pass the in-circle test.
#[test]
fn cdt_constraint_restores_delaunay() {
    let mut pslg = Pslg::new();
    // 3 × 3 grid points.
    let v00 = pslg.add_vertex(0.0, 0.0);
    let _v10 = pslg.add_vertex(1.0, 0.0);
    let v20 = pslg.add_vertex(2.0, 0.0);
    let _v01 = pslg.add_vertex(0.0, 1.0);
    let _v11 = pslg.add_vertex(1.0, 1.0);
    let _v21 = pslg.add_vertex(2.0, 1.0);
    let v02 = pslg.add_vertex(0.0, 2.0);
    let _v12 = pslg.add_vertex(1.0, 2.0);
    let v22 = pslg.add_vertex(2.0, 2.0);

    // Boundary.
    pslg.add_segment(v00, v20);
    pslg.add_segment(v20, v22);
    pslg.add_segment(v22, v02);
    pslg.add_segment(v02, v00);

    // Diagonal constraint cutting across multiple Delaunay edges.
    pslg.add_segment(v00, v22);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    assert!(
        cdt.is_constrained(v00, v22),
        "Diagonal constraint must be enforced"
    );
    assert!(
        dt.is_delaunay(),
        "CDT with grid + diagonal should be Delaunay"
    );
    // All 9 points should be in the triangulation.
    assert!(dt.vertex_count() >= 9);
    // Verify other boundary constraints.
    assert!(cdt.is_constrained(v00, v20));
    assert!(cdt.is_constrained(v20, v22));
    assert!(cdt.is_constrained(v22, v02));
    assert!(cdt.is_constrained(v02, v00));
}

// ── PSLG rejects coincident vertices ──────────────────────────────────────

/// Two vertices at the same position should be rejected by validation.
#[test]
fn pslg_rejects_coincident_vertices() {
    let mut pslg = Pslg::new();
    pslg.add_vertex(1.0, 2.0);
    pslg.add_vertex(1.0, 2.0); // exact duplicate

    let err = pslg.validate().err().expect("should reject coincident verts");
    assert!(
        matches!(err, PslgValidationError::CoincidentVertices { .. }),
        "expected CoincidentVertices, got: {err}"
    );
}

/// Near-coincident vertices (within tolerance) should also be rejected.
#[test]
fn pslg_rejects_near_coincident_vertices() {
    let mut pslg = Pslg::new();
    pslg.add_vertex(1.0, 2.0);
    pslg.add_vertex(1.0 + 1e-16, 2.0 + 1e-16); // within tolerance

    let err = pslg.validate().err().expect("should reject near-coincident");
    assert!(
        matches!(err, PslgValidationError::CoincidentVertices { .. }),
        "expected CoincidentVertices, got: {err}"
    );
}

/// Vertices that are clearly distinct should pass validation.
#[test]
fn pslg_accepts_distinct_vertices() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 0.0);
    let c = pslg.add_vertex(0.5, 1.0);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    assert!(pslg.validate().is_ok(), "distinct vertices should pass");
}

// ── Adversarial near-degenerate constraint ────────────────────────────────

/// A very thin, long triangle with a constraint across its shortest edge.
/// Tests that scale-relative tolerances handle near-degenerate geometry.
#[test]
fn cdt_near_degenerate_thin_triangle() {
    let mut pslg = Pslg::new();
    // Very elongated triangle: base 1000, height ~0.001
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1000.0, 0.0);
    let c = pslg.add_vertex(500.0, 0.001);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    assert!(
        dt.triangle_count() >= 1,
        "Near-degenerate CDT should produce triangles"
    );
    assert!(cdt.is_constrained(a, b));
    assert!(cdt.is_constrained(b, c));
    assert!(cdt.is_constrained(c, a));
}
