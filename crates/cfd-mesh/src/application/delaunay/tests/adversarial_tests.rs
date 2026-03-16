//! Adversarial tests targeting known failure modes in mesh libraries.
//!
//! These tests exercise pathological geometries and edge cases that commonly
//! expose bugs in Delaunay triangulation, CDT constraint recovery, and
//! Ruppert refinement implementations.
//!
//! # Motivation
//!
//! Production mesh libraries (Triangle, CGAL, Gmsh, etc.) have documented
//! issues with:
//! - Co-circular (co-spherical in 3-D) point sets → degenerate incircle
//! - Collinear / near-collinear point sets → degenerate orient2d
//! - Near-coincident points → floating-point resolution collapse
//! - Very thin slivers → near-zero area, circumradius blowup
//! - Regular grids → systematic co-circularity
//! - Walk cycles → infinite loops in Lawson walk
//! - Constraint through existing DT edges → no-op vs. spurious flip
//! - Small input angles → Ruppert non-termination
//!
//! Each test is annotated with the failure mode it targets and a reference
//! to the relevant literature where applicable.

use crate::application::delaunay::constraint::enforce::Cdt;
use crate::application::delaunay::pslg::graph::Pslg;
use crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation;
use std::f64::consts::PI;

// ── Co-circular point sets ────────────────────────────────────────────────

/// **Failure mode**: 4 co-circular points produce a degenerate incircle test
/// (the test returns exactly zero).  Naive implementations crash or produce
/// non-Delaunay output.
///
/// **Literature**: Shewchuk (1997), "Adaptive Precision Floating-Point
/// Arithmetic and Fast Robust Geometric Predicates."
#[test]
fn cocircular_4_points_on_unit_circle() {
    // Four points on the unit circle at 0°, 90°, 180°, 270°.
    let pts: Vec<(f64, f64)> = (0..4)
        .map(|i| {
            let angle = PI / 2.0 * i as f64;
            (angle.cos(), angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 4);
    assert_eq!(dt.triangle_count(), 2, "4 co-circular points → exactly 2 triangles");
    assert!(dt.is_delaunay());
}

/// **Failure mode**: 8 co-circular points — every subset of 4 is co-circular.
/// Maximises the number of degenerate incircle tests during construction.
#[test]
fn cocircular_8_points_regular_octagon() {
    let pts: Vec<(f64, f64)> = (0..8)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / 8.0;
            (angle.cos(), angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 8);
    // Euler: T = 2V - 2 - h  (h = 0 for convex hull), so T = 2·8 - 2 = 14? No.
    // For a convex polygon with n vertices: T = n - 2.
    assert_eq!(dt.triangle_count(), 6, "convex 8-gon → 6 interior triangles");
    assert!(dt.is_delaunay());
}

/// **Failure mode**: 16 co-circular points stress both incircle and the walk.
/// With Hilbert ordering, consecutive insertions may land on near-degenerate
/// incircle quadrilaterals.
#[test]
fn cocircular_16_points_on_circle() {
    let pts: Vec<(f64, f64)> = (0..16)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / 16.0;
            (angle.cos(), angle.sin())
        })
        .collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 16);
    assert!(dt.triangle_count() >= 14, "at least n-2 triangles");
    assert!(dt.is_delaunay());
}

/// **Failure mode**: co-circular with a centre point.  The centre is
/// equidistant from all circumcircles — maximum degeneracy.
#[test]
fn cocircular_with_centre_point() {
    let mut pts: Vec<(f64, f64)> = (0..8)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / 8.0;
            (angle.cos(), angle.sin())
        })
        .collect();
    pts.push((0.0, 0.0)); // centre point
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 9);
    // With centre: each sector is a triangle → 8 triangles.
    assert_eq!(dt.triangle_count(), 8, "8 sectors from centre");
    assert!(dt.is_delaunay());
}

// ── Collinear point sets ──────────────────────────────────────────────────

/// **Failure mode**: all points collinear → zero-area triangulation.
/// Degenerate orient2d for every triple.  Some implementations crash
/// or produce invalid topology.
#[test]
fn collinear_horizontal_5_points() {
    let pts: Vec<(f64, f64)> = (0..5).map(|i| (i as f64, 0.0)).collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 5);
    // All collinear → no interior triangles (all degenerate or on convex hull).
    assert_eq!(dt.triangle_count(), 0, "collinear points → 0 interior triangles");
}

/// **Failure mode**: collinear diagonal — tests orient2d with non-axis-aligned
/// collinearity.
#[test]
fn collinear_diagonal_5_points() {
    let pts: Vec<(f64, f64)> = (0..5).map(|i| (i as f64, i as f64)).collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 5);
    assert_eq!(dt.triangle_count(), 0, "diagonal collinear → 0 interior triangles");
}

/// **Failure mode**: near-collinear points with one barely-off-axis vertex.
/// The orient2d result is near-zero → floating-point sensitivity.
#[test]
fn near_collinear_with_perturbation() {
    let pts = vec![
        (0.0, 0.0),
        (1.0, 0.0),
        (2.0, 0.0),
        (3.0, 0.0),
        (1.5, 0.05), // slightly above the line (circumradius ~22 < margin 60)
    ];
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 5);
    assert!(dt.is_delaunay());
    // Should produce at least 2 real triangles with the perturbed point.
    assert!(dt.triangle_count() >= 2);
}

// ── Near-coincident points ────────────────────────────────────────────────

/// **Failure mode**: two nearly-coincident points separated by machine
/// epsilon.  Tests that duplicate detection and orient2d remain correct.
#[test]
fn near_coincident_pair() {
    let eps = 1e-14;
    let pts = vec![
        (0.0, 0.0),
        (1.0, 0.0),
        (0.5, 1.0),
        (0.5 + eps, 1.0 + eps), // near-duplicate of (0.5, 1.0)
    ];
    let dt = DelaunayTriangulation::from_points(&pts);
    // Should not crash.  May merge or keep as separate — either is OK.
    assert!(dt.vertex_count() >= 3);
    assert!(dt.is_delaunay());
}

/// **Failure mode**: cluster of near-coincident points.
#[test]
fn near_coincident_cluster() {
    let centre = (0.5, 0.5);
    let eps = 1e-12;
    let mut pts = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    // Add 8 points clustered within eps of centre.
    for i in 0..8 {
        let angle = 2.0 * PI * i as f64 / 8.0;
        pts.push((centre.0 + eps * angle.cos(), centre.1 + eps * angle.sin()));
    }
    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(dt.is_delaunay());
}

// ── Thin slivers ──────────────────────────────────────────────────────────

/// **Failure mode**: thin triangle with aspect ratio ~100:1.
/// Tests circumradius computation stability and incircle robustness.
///
/// Note: the circumradius must fit within the super-triangle (margin = 20·dmax)
/// for the interior triangle to survive the Delaunay in-circle test against
/// super-triangle vertices.
#[test]
fn thin_sliver_triangle() {
    let pts = vec![
        (0.0, 0.0),
        (1.0, 0.0),
        (0.5, 0.01), // aspect ratio ~100:1, circumradius ~12.5 < margin 20
    ];
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 3);
    assert_eq!(dt.triangle_count(), 1);
    assert!(dt.is_delaunay());
}

/// **Failure mode**: needle triangle with high aspect ratio.
#[test]
fn needle_triangle() {
    let pts = vec![
        (0.0, 0.0),
        (100.0, 0.0),
        (50.0, 1.0), // aspect ratio ~100:1
    ];
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 3);
    assert_eq!(dt.triangle_count(), 1);
    assert!(dt.is_delaunay());
}

/// **Failure mode**: extremely thin triangle whose circumradius exceeds the
/// super-triangle bounding box.  Documents the known limitation: the BW
/// algorithm produces no interior triangles because all connect to
/// super-triangle vertices.
///
/// **Literature**: super-triangle–based BW implementations break when the
/// circumradius of the input greatly exceeds the super-triangle margin
/// (Shewchuk 1996).
#[test]
fn extreme_needle_documents_super_tri_limitation() {
    let pts = vec![
        (0.0, 0.0),
        (1e4, 0.0),
        (5000.0, 1e-4), // circumradius ~1.25e11 >> super-tri margin ~2e5
    ];
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 3);
    // Known limitation: circumradius exceeds super-triangle → 0 interior triangles.
    assert_eq!(dt.triangle_count(), 0);
    assert!(dt.is_delaunay());
}

// ── Regular grids ─────────────────────────────────────────────────────────

/// **Failure mode**: regular grid creates systematic co-circularity.
/// Every set of 4 adjacent grid points forms a rectangle (co-circular).
/// This maximises degenerate incircle tests.
///
/// **Literature**: Shewchuk (1996) notes that regular grids are a worst-case
/// for incremental Delaunay since every quad is degenerate.
#[test]
fn regular_grid_10x10() {
    let mut pts = Vec::with_capacity(100);
    for i in 0..10 {
        for j in 0..10 {
            pts.push((i as f64, j as f64));
        }
    }
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 100);
    // Grid: (n-1)×(n-1) squares × 2 triangles = 9×9×2 = 162.
    assert!(
        dt.triangle_count() >= 160 && dt.triangle_count() <= 165,
        "10×10 grid should have ~162 triangles, got {}",
        dt.triangle_count()
    );
    assert!(dt.is_delaunay());
}

/// **Failure mode**: larger grid stress test.
#[test]
fn regular_grid_20x20() {
    let mut pts = Vec::with_capacity(400);
    for i in 0..20 {
        for j in 0..20 {
            pts.push((i as f64, j as f64));
        }
    }
    let dt = DelaunayTriangulation::from_points(&pts);
    assert_eq!(dt.vertex_count(), 400);
    assert!(dt.is_delaunay());
}

// ── Stress tests ──────────────────────────────────────────────────────────

/// **Failure mode**: performance regression.  5000 uniformly distributed
/// points should complete in reasonable time with Hilbert ordering.
#[test]
fn stress_5000_uniform_random() {
    // Deterministic pseudo-random via LCG (no rand dependency).
    let mut rng: u64 = 0xDEAD_BEEF_CAFE_1234;
    let mut next_f64 = || -> f64 {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
        (rng >> 11) as f64 / (1u64 << 53) as f64
    };
    let pts: Vec<(f64, f64)> = (0..5000).map(|_| (next_f64(), next_f64())).collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(dt.vertex_count() <= 5000); // some may be deduplicated
    assert!(dt.is_delaunay());
    // Euler formula: T ≈ 2V - 2 - h for convex hull with h boundary vertices.
    assert!(dt.triangle_count() > 0);
}

// ── Edge cases ────────────────────────────────────────────────────────────

/// **Failure mode**: single point → no triangles.
#[test]
fn single_point() {
    let dt = DelaunayTriangulation::from_points(&[(42.0, 17.0)]);
    assert_eq!(dt.vertex_count(), 1);
    assert_eq!(dt.triangle_count(), 0);
}

/// **Failure mode**: two points → no triangles (degenerate 1-simplex).
#[test]
fn two_points() {
    let dt = DelaunayTriangulation::from_points(&[(0.0, 0.0), (1.0, 1.0)]);
    assert_eq!(dt.vertex_count(), 2);
    assert_eq!(dt.triangle_count(), 0);
}

/// **Failure mode**: empty input.
#[test]
fn empty_input() {
    let dt = DelaunayTriangulation::from_points(&[]);
    assert_eq!(dt.vertex_count(), 0);
    assert_eq!(dt.triangle_count(), 0);
}

/// **Failure mode**: all duplicate points → should reduce to 1 effective vertex.
#[test]
fn all_duplicates() {
    let pts: Vec<(f64, f64)> = (0..50).map(|_| (1.0, 2.0)).collect();
    let dt = DelaunayTriangulation::from_points(&pts);
    // Duplicates are skipped during insertion — only 1 effective vertex.
    assert_eq!(dt.triangle_count(), 0);
}

// ── Walk cycle adversarial ────────────────────────────────────────────────

/// **Failure mode**: Lawson walk can cycle between 2-3 triangles in
/// degenerate configurations.  This arrangement — several nearly-collinear
/// points with slightly offset perpendiculars — is designed to create
/// near-zero orient2d values that can cause the walk to oscillate.
///
/// **Literature**: Devillers, Pion, Teillaud (2002), "Walking in a
/// triangulation."
#[test]
fn walk_cycle_adversarial() {
    let mut pts = Vec::new();
    // Long thin strip with perpendicular perturbations.
    for i in 0..20 {
        let x = i as f64;
        pts.push((x, 0.0));
        pts.push((x, 1e-6));
    }
    // Query-like points that force the walk through the thin strip.
    pts.push((10.0, 5e-7));
    pts.push((15.0, 5e-7));
    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(dt.is_delaunay());
    assert!(dt.triangle_count() > 0);
}

// ── CDT adversarial ──────────────────────────────────────────────────────

/// **Failure mode**: constraint along an existing Delaunay edge.
/// The constraint enforcement should be a no-op (no flips needed), but
/// naive implementations may spuriously flip the edge.
#[test]
fn cdt_constraint_along_existing_dt_edge() {
    let mut pslg = Pslg::new();
    let v0 = pslg.add_vertex(0.0, 0.0);
    let v1 = pslg.add_vertex(1.0, 0.0);
    let v2 = pslg.add_vertex(0.5, 1.0);
    // Constraint matches an edge that Delaunay would naturally produce.
    pslg.add_segment(v0, v1);
    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert_eq!(dt.vertex_count(), 3);
    assert_eq!(dt.triangle_count(), 1);
    assert!(dt.is_delaunay());
    let _ = v2;
}

/// **Failure mode**: co-circular CDT with constraint through the centre.
/// Combines co-circularity degeneracy with constraint enforcement stress.
#[test]
fn cdt_cocircular_with_diameter_constraint() {
    let mut pslg = Pslg::new();
    let mut vids = Vec::new();
    for i in 0..8 {
        let angle = 2.0 * PI * i as f64 / 8.0;
        vids.push(pslg.add_vertex(angle.cos(), angle.sin()));
    }
    // Diameter constraint: vertex 0 to vertex 4 (opposites).
    pslg.add_segment(vids[0], vids[4]);
    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    // The constraint must be present in the output.
    assert!(
        cdt.constrained_edges().contains(&(vids[0], vids[4]))
            || cdt.constrained_edges().contains(&(vids[4], vids[0])),
        "diameter constraint should be preserved"
    );
}

/// **Failure mode**: CDT with constraints forming a complex polygon with
/// narrow angles.  Tests constraint recovery near-degenerate angles.
#[test]
fn cdt_narrow_angle_polygon() {
    let mut pslg = Pslg::new();
    // Star-shaped polygon with narrow tips.
    let n = 5;
    let mut vids = Vec::new();
    for i in 0..n {
        let angle = 2.0 * PI * i as f64 / n as f64;
        // Outer vertex (tip of star).
        let r_outer = 2.0;
        vids.push(pslg.add_vertex(r_outer * angle.cos(), r_outer * angle.sin()));
        // Inner vertex (notch between tips).
        let inner_angle = angle + PI / n as f64;
        let r_inner = 0.5;
        vids.push(pslg.add_vertex(r_inner * inner_angle.cos(), r_inner * inner_angle.sin()));
    }
    // Close the star polygon with segments.
    for i in 0..vids.len() {
        pslg.add_segment(vids[i], vids[(i + 1) % vids.len()]);
    }
    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(dt.triangle_count() > 0);
}

/// **Failure mode**: CDT with many crossing constraints that must be resolved
/// via Steiner points.  Tests `resolve_crossings` and constraint recovery
/// interaction.
#[test]
fn cdt_dense_crosshatch_constraints() {
    let mut pslg = Pslg::new();
    // 5×5 grid of vertices in the unit square.
    let mut grid = Vec::new();
    for i in 0..=4 {
        let mut row = Vec::new();
        for j in 0..=4 {
            row.push(pslg.add_vertex(i as f64 / 4.0, j as f64 / 4.0));
        }
        grid.push(row);
    }
    // Horizontal segments.
    for row in &grid {
        for w in row.windows(2) {
            pslg.add_segment(w[0], w[1]);
        }
    }
    // Vertical segments.
    for j in 0..5 {
        for i in 0..4 {
            pslg.add_segment(grid[i][j], grid[i + 1][j]);
        }
    }
    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    // 4×4 grid = 16 cells × 2 triangles = 32 expected.
    assert!(
        dt.triangle_count() >= 30,
        "4×4 grid CDT should have ~32 triangles, got {}",
        dt.triangle_count()
    );
}

/// **Failure mode**: CDT with a hole entirely inside the domain.
/// Tests that hole-seeding and flood-fill removal correctly identify
/// interior vs. exterior regions.
#[test]
fn cdt_hole_inside_domain() {
    let mut pslg = Pslg::new();
    // Outer square.
    let o0 = pslg.add_vertex(0.0, 0.0);
    let o1 = pslg.add_vertex(4.0, 0.0);
    let o2 = pslg.add_vertex(4.0, 4.0);
    let o3 = pslg.add_vertex(0.0, 4.0);
    pslg.add_segment(o0, o1);
    pslg.add_segment(o1, o2);
    pslg.add_segment(o2, o3);
    pslg.add_segment(o3, o0);

    // Inner square (hole).
    let h0 = pslg.add_vertex(1.0, 1.0);
    let h1 = pslg.add_vertex(3.0, 1.0);
    let h2 = pslg.add_vertex(3.0, 3.0);
    let h3 = pslg.add_vertex(1.0, 3.0);
    pslg.add_segment(h0, h1);
    pslg.add_segment(h1, h2);
    pslg.add_segment(h2, h3);
    pslg.add_segment(h3, h0);

    // Seed the hole.
    pslg.add_hole(2.0, 2.0);
    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    // Verify no triangles exist inside the hole.
    for (_, tri) in dt.interior_triangles() {
        let centroid_x: f64 = tri
            .vertices
            .iter()
            .map(|v| dt.vertex(*v).x)
            .sum::<f64>()
            / 3.0;
        let centroid_y: f64 = tri
            .vertices
            .iter()
            .map(|v| dt.vertex(*v).y)
            .sum::<f64>()
            / 3.0;
        let inside_hole = centroid_x > 1.1 && centroid_x < 2.9 && centroid_y > 1.1 && centroid_y < 2.9;
        assert!(
            !inside_hole,
            "triangle centroid ({centroid_x}, {centroid_y}) should not be inside hole"
        );
    }
}

// ── Refinement adversarial ────────────────────────────────────────────────

/// **Failure mode**: Ruppert refinement with acute input angles (< 20°).
/// Known to cause infinite Steiner point insertion loops in naive
/// implementations.  Robust implementations must handle encroachment
/// splitting of obtuse Steiner vertices.
///
/// **Literature**: Ruppert (1995), "A Delaunay Refinement Algorithm for
/// Quality 2-Dimensional Mesh Generation." Theorem 4: termination requires
/// minimum input angle ≥ ~20.7°.
#[test]
fn ruppert_acute_input_angle() {
    use crate::application::delaunay::refinement::ruppert::RuppertRefiner;

    let mut pslg = Pslg::new();
    // Isoceles triangle with very narrow top — ~10° at the apex.
    let v0 = pslg.add_vertex(0.0, 0.0);
    let v1 = pslg.add_vertex(2.0, 0.0);
    // Apex angle ~10°.
    let apex_half_angle = 5.0_f64.to_radians();
    let h = 1.0 / apex_half_angle.tan();
    let v2 = pslg.add_vertex(1.0, h);
    pslg.add_segment(v0, v1);
    pslg.add_segment(v1, v2);
    pslg.add_segment(v2, v0);

    let cdt = Cdt::from_pslg(&pslg);
    let mut refiner = RuppertRefiner::new(cdt);
    refiner.set_max_steiner(500);
    let _n_steiner = refiner.refine();
    // Should terminate (the Steiner limit prevents infinite loops).
    let dt = refiner.cdt().triangulation();
    assert!(dt.triangle_count() > 0);
    assert!(dt.is_delaunay());
}

/// **Failure mode**: Ruppert on a highly elongated rectangle.  Tests that
/// refinement correctly handles high aspect-ratio input domains.
#[test]
fn ruppert_elongated_rectangle() {
    use crate::application::delaunay::refinement::ruppert::RuppertRefiner;

    let mut pslg = Pslg::new();
    let v0 = pslg.add_vertex(0.0, 0.0);
    let v1 = pslg.add_vertex(10.0, 0.0);
    let v2 = pslg.add_vertex(10.0, 1.0);
    let v3 = pslg.add_vertex(0.0, 1.0);
    pslg.add_segment(v0, v1);
    pslg.add_segment(v1, v2);
    pslg.add_segment(v2, v3);
    pslg.add_segment(v3, v0);

    let cdt = Cdt::from_pslg(&pslg);
    let mut refiner = RuppertRefiner::new(cdt);
    refiner.set_max_area(0.5);
    refiner.set_max_steiner(10000);
    let _n_steiner = refiner.refine();
    let dt = refiner.cdt().triangulation();
    assert!(
        dt.triangle_count() > 10,
        "10×1 rectangle with max_area=0.5 should have many triangles, got {}",
        dt.triangle_count()
    );
    assert!(dt.is_delaunay());
}
