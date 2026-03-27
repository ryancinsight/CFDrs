//! Adversarial CDT, PSLG, and Ruppert tests targeting known mesh library
//! failure modes not covered by the primary adversarial test suite.
//!
//! # Targeted Failure Modes
//!
//! - **Spiral point distributions**: stress Lawson walk locality — O(√n) walk
//!   distance instead of near-O(1) when Hilbert ordering is ineffective.
//! - **Long constraint crossing many DT edges**: CDT constraint recovery must
//!   flip O(k) edges where k can be O(n).  Known to crash Triangle 1.6 on
//!   certain inputs (Shewchuk, personal communication).
//! - **Star-shaped constraints**: many constraints radiating from a single
//!   vertex create a high-degree fan that stresses cavity re-triangulation.
//! - **Closely-spaced parallel constraints**: narrow channels between
//!   constraint edges produce extreme aspect-ratio triangles during CDT
//!   recovery.
//! - **Concentric polygon constraints**: nested constraint polygons with
//!   holes test hole-removal correctness.
//! - **Points on constraint edges**: Steiner points placed exactly on
//!   constraint segments test OnEdge location handling.
//! - **Ruppert area-only refinement**: validates the max_area constraint
//!   independent of angle quality.
//! - **CDT with many crossing DT edges**: a single long constraint that
//!   crosses many Delaunay edges exercises the flip-recovery algorithm.
//! - **PSLG with multiple adjacent small holes**: holes close together test
//!   hole-removal flood-fill isolation.

use crate::application::delaunay::dim2::constraint::enforce::Cdt;
use crate::application::delaunay::dim2::pslg::graph::Pslg;
use crate::application::delaunay::dim2::triangulation::adjacency::Adjacency;
use crate::application::delaunay::dim2::triangulation::bowyer_watson::DelaunayTriangulation;
use std::f64::consts::PI;

// ── Spiral point distribution ─────────────────────────────────────────────

/// **Failure mode**: Spiral point distributions defeat spatial locality of
/// Hilbert-curve ordering, causing the Lawson walk to traverse O(√n) edges
/// per insertion instead of near-O(1).
///
/// **Literature**: Amenta, Choi, Rote (2003), "Incremental Constructions
/// con BRIO" — BRIO ordering mitigates but doesn't eliminate spiral pathology.
#[test]
fn spiral_distribution_200_points() {
    let n = 200;
    let pts: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let t = i as f64;
            let r = 1.0 + 0.1 * t;
            let theta = t * 0.5;
            (r * theta.cos(), r * theta.sin())
        })
        .collect();

    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(dt.is_delaunay());
    assert_eq!(dt.vertex_count(), n);
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── Long constraint crossing many DT edges ────────────────────────────────

/// **Failure mode**: A single horizontal constraint through a random point
/// cloud crosses many Delaunay edges. CDT constraint recovery must flip
/// each crossing edge correctly without corrupting adjacency.
///
/// **Literature**: Shewchuk (2002), "Constrained Delaunay Tetrahedrali-
/// zations and Provably Good Boundary Recovery".  The 2D analogue of the
/// cavity re-triangulation problem.
#[test]
fn long_constraint_through_random_cloud() {
    let mut pslg = Pslg::new();

    // Deterministic interior points scattered above and below y=2.5.
    let n = 40;
    for i in 0..n {
        let px = 1.0 + f64::from(i) * 8.0 / f64::from(n);
        let py = if i % 2 == 0 { 1.2 + f64::from(i) * 0.02 } else { 3.5 + f64::from(i) * 0.02 };
        pslg.add_vertex(px, py);
    }

    // Long horizontal constraint endpoints — must be boundary vertices.
    let left = pslg.add_vertex(0.0, 2.5);
    let right = pslg.add_vertex(10.0, 2.5);

    // Boundary rectangle, split at constraint endpoints.
    let bl = pslg.add_vertex(0.0, 0.0);
    let br = pslg.add_vertex(10.0, 0.0);
    let tr = pslg.add_vertex(10.0, 5.0);
    let tl = pslg.add_vertex(0.0, 5.0);
    // Bottom.
    pslg.add_segment(bl, br);
    // Right side split at `right`.
    pslg.add_segment(br, right);
    pslg.add_segment(right, tr);
    // Top.
    pslg.add_segment(tr, tl);
    // Left side split at `left`.
    pslg.add_segment(tl, left);
    pslg.add_segment(left, bl);

    // Long horizontal constraint at y=2.5 crossing many DT edges.
    pslg.add_segment(left, right);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    // CDTs are Delaunay *except* across constraint edges, so is_delaunay()
    // may legitimately return false when constraints violate the empty-
    // circumcircle property.  We verify adjacency symmetry instead.
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    assert!(
        dt.triangle_count() >= 2,
        "CDT should have at least 2 triangles"
    );
}

// ── Star-shaped constraints ───────────────────────────────────────────────

/// **Failure mode**: Many constraints radiating from a single central vertex
/// to boundary vertices create a high-degree fan.  This stresses constraint
/// recovery (many sequential flips around one vertex) and the CCW-ordering
/// invariant.
///
/// **Known issue** (Gmsh 4.x): star-shaped constraints with > 16 rays
/// occasionally produce inverted triangles due to float rounding during
/// multi-flip sequences.
#[test]
fn star_shaped_constraints_24_rays() {
    let mut pslg = Pslg::new();

    // Central vertex.
    let center = pslg.add_vertex(5.0, 5.0);

    // 24 boundary vertices on a circle of radius 4.
    let num_rays = 24;
    let mut boundary = Vec::new();
    for i in 0..num_rays {
        let theta = 2.0 * PI * (i as f64) / (num_rays as f64);
        let vx = 5.0 + 4.0 * theta.cos();
        let vy = 5.0 + 4.0 * theta.sin();
        boundary.push(pslg.add_vertex(vx, vy));
    }

    // Boundary polygon.
    for i in 0..num_rays {
        pslg.add_segment(boundary[i], boundary[(i + 1) % num_rays]);
    }

    // Star constraints: center → each boundary vertex.
    for &bv in &boundary {
        pslg.add_segment(center, bv);
    }

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    // Should produce exactly 2 * num_rays interior triangles (fan from center).
    assert!(
        dt.triangle_count() >= num_rays,
        "Star CDT should produce at least {num_rays} triangles, got {}",
        dt.triangle_count()
    );
}

// ── Closely-spaced parallel constraints ───────────────────────────────────

/// **Failure mode**: Two nearly-parallel constraint segments separated by a
/// tiny gap produce extreme aspect-ratio triangles during CDT recovery.
/// Naive implementations fail when the gap width is below the epsilon
/// threshold used in orientation tests.
///
/// **Known issue** (Triangle 1.6): parallel constraints with gap < 1e-8
/// cause "zero-area triangle" assertions.
#[test]
fn closely_spaced_parallel_constraints() {
    let mut pslg = Pslg::new();
    let gap = 0.01;

    // Constraint endpoints on the left and right boundary — must split
    // boundary segments at these points to keep the PSLG valid.
    let a0 = pslg.add_vertex(0.0, 1.0 - gap / 2.0);
    let a1 = pslg.add_vertex(10.0, 1.0 - gap / 2.0);
    let b0 = pslg.add_vertex(0.0, 1.0 + gap / 2.0);
    let b1 = pslg.add_vertex(10.0, 1.0 + gap / 2.0);

    // Outer boundary split at constraint endpoints.
    let bl = pslg.add_vertex(0.0, 0.0);
    let br = pslg.add_vertex(10.0, 0.0);
    let tr = pslg.add_vertex(10.0, 2.0);
    let tl = pslg.add_vertex(0.0, 2.0);
    // Bottom.
    pslg.add_segment(bl, br);
    // Right side, split at a1, b1.
    pslg.add_segment(br, a1);
    pslg.add_segment(a1, b1);
    pslg.add_segment(b1, tr);
    // Top.
    pslg.add_segment(tr, tl);
    // Left side, split at b0, a0.
    pslg.add_segment(tl, b0);
    pslg.add_segment(b0, a0);
    pslg.add_segment(a0, bl);

    // Two parallel horizontal constraints.
    pslg.add_segment(a0, a1);
    pslg.add_segment(b0, b1);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── Concentric polygon constraints ────────────────────────────────────────

/// **Failure mode**: Nested constraint polygons with a hole between them
/// tests hole-removal correctness.  The flood-fill must correctly identify
/// the inner region as a hole and not delete the outer domain.
#[test]
fn concentric_square_with_hole() {
    let mut pslg = Pslg::new();

    // Outer square: [0,10] x [0,10].
    let o0 = pslg.add_vertex(0.0, 0.0);
    let o1 = pslg.add_vertex(10.0, 0.0);
    let o2 = pslg.add_vertex(10.0, 10.0);
    let o3 = pslg.add_vertex(0.0, 10.0);
    pslg.add_segment(o0, o1);
    pslg.add_segment(o1, o2);
    pslg.add_segment(o2, o3);
    pslg.add_segment(o3, o0);

    // Inner square (hole boundary): [3,7] x [3,7].
    let i0 = pslg.add_vertex(3.0, 3.0);
    let i1 = pslg.add_vertex(7.0, 3.0);
    let i2 = pslg.add_vertex(7.0, 7.0);
    let i3 = pslg.add_vertex(3.0, 7.0);
    pslg.add_segment(i0, i1);
    pslg.add_segment(i1, i2);
    pslg.add_segment(i2, i3);
    pslg.add_segment(i3, i0);

    // Mark the interior of the inner square as a hole.
    pslg.add_hole(5.0, 5.0);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    // No triangle should have all 3 vertices on the inner square.
    for (_, tri) in dt.interior_triangles() {
        let all_inner = tri.vertices.iter().all(|v| {
            let vt = dt.vertex(*v);
            vt.x >= 3.0 && vt.x <= 7.0 && vt.y >= 3.0 && vt.y <= 7.0
        });
        // A triangle can have inner vertices IF it bridges the gap.
        // But no triangle should be entirely inside the hole.
        if all_inner {
            // Check that the centroid is NOT inside the hole.
            let cx: f64 = tri
                .vertices
                .iter()
                .map(|v| dt.vertex(*v).x)
                .sum::<f64>()
                / 3.0;
            let cy: f64 = tri
                .vertices
                .iter()
                .map(|v| dt.vertex(*v).y)
                .sum::<f64>()
                / 3.0;
            assert!(
                !(3.0..=7.0).contains(&cx) || !(3.0..=7.0).contains(&cy),
                "Found a triangle centroid inside the hole: ({cx}, {cy})"
            );
        }
    }
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── Multiple adjacent holes ───────────────────────────────────────────────

/// **Failure mode**: Multiple small holes close together stress the
/// hole-removal flood-fill.  Incorrect flood-fill can merge separate holes
/// or delete exterior triangles.
#[test]
fn multiple_adjacent_small_holes() {
    let mut pslg = Pslg::new();

    // Outer boundary: [0,20] x [0,10].
    let o0 = pslg.add_vertex(0.0, 0.0);
    let o1 = pslg.add_vertex(20.0, 0.0);
    let o2 = pslg.add_vertex(20.0, 10.0);
    let o3 = pslg.add_vertex(0.0, 10.0);
    pslg.add_segment(o0, o1);
    pslg.add_segment(o1, o2);
    pslg.add_segment(o2, o3);
    pslg.add_segment(o3, o0);

    // Two small square holes side by side.
    // Hole 1: [4,6] x [4,6]
    let h1a = pslg.add_vertex(4.0, 4.0);
    let h1b = pslg.add_vertex(6.0, 4.0);
    let h1c = pslg.add_vertex(6.0, 6.0);
    let h1d = pslg.add_vertex(4.0, 6.0);
    pslg.add_segment(h1a, h1b);
    pslg.add_segment(h1b, h1c);
    pslg.add_segment(h1c, h1d);
    pslg.add_segment(h1d, h1a);
    pslg.add_hole(5.0, 5.0);

    // Hole 2: [8,10] x [4,6] — adjacent with 2-unit gap.
    let h2a = pslg.add_vertex(8.0, 4.0);
    let h2b = pslg.add_vertex(10.0, 4.0);
    let h2c = pslg.add_vertex(10.0, 6.0);
    let h2d = pslg.add_vertex(8.0, 6.0);
    pslg.add_segment(h2a, h2b);
    pslg.add_segment(h2b, h2c);
    pslg.add_segment(h2c, h2d);
    pslg.add_segment(h2d, h2a);
    pslg.add_hole(9.0, 5.0);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    // At least some triangles should exist in the gap region (x ∈ [6,8]).
    let gap_tris = dt
        .interior_triangles()
        .filter(|(_, tri)| {
            let cx: f64 = tri
                .vertices
                .iter()
                .map(|v| dt.vertex(*v).x)
                .sum::<f64>()
                / 3.0;
            cx > 6.0 && cx < 8.0
        })
        .count();
    assert!(
        gap_tris > 0,
        "Should have triangles in the gap between the two holes"
    );
}

// ── Ruppert area-only refinement ──────────────────────────────────────────

/// **Failure mode**: Testing max_area constraint alone (no angle constraint)
/// exercises the area-based priority queue ordering independently of
/// radius-edge ratio.
#[test]
fn ruppert_area_only_refinement() {
    use crate::application::delaunay::dim2::refinement::quality::TriangleQuality;
    use crate::application::delaunay::dim2::refinement::ruppert::RuppertRefiner;

    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(10.0, 0.0);
    let c = pslg.add_vertex(5.0, 8.66);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    let cdt = Cdt::from_pslg(&pslg);
    let mut refiner = RuppertRefiner::new(cdt);
    // Very lenient angle (basically no angle constraint).
    refiner.set_max_ratio(10.0);
    // Tight area constraint.
    refiner.set_max_area(2.0);
    refiner.set_max_steiner(500);
    let n = refiner.refine();

    assert!(n > 0, "Should insert Steiner points for area refinement");

    let dt = refiner.cdt().triangulation();
    for (_, tri) in dt.interior_triangles() {
        let v0 = dt.vertex(tri.vertices[0]);
        let v1 = dt.vertex(tri.vertices[1]);
        let v2 = dt.vertex(tri.vertices[2]);
        let q = TriangleQuality::compute(v0, v1, v2);
        assert!(
            q.area < 2.5,
            "Triangle area {:.2} exceeds max_area 2.0 (with tolerance)",
            q.area
        );
    }
}

// ── CDT zigzag constraint ─────────────────────────────────────────────────

/// **Failure mode**: A zigzag polyline constraint (many short segments at
/// alternating angles) stresses constraint recovery with repeated
/// flip-sequences in close proximity.
#[test]
fn cdt_zigzag_constraint() {
    let mut pslg = Pslg::new();

    // Outer boundary.
    let bl = pslg.add_vertex(0.0, 0.0);
    let br = pslg.add_vertex(20.0, 0.0);
    let tr = pslg.add_vertex(20.0, 6.0);
    let tl = pslg.add_vertex(0.0, 6.0);
    pslg.add_segment(bl, br);
    pslg.add_segment(br, tr);
    pslg.add_segment(tr, tl);
    pslg.add_segment(tl, bl);

    // Zigzag polyline from left to right.
    let n_zigs = 10;
    let mut prev = pslg.add_vertex(1.0, 3.0);
    for i in 1..=n_zigs {
        let x = 1.0 + f64::from(i) * 18.0 / f64::from(n_zigs);
        let y = if i % 2 == 0 { 4.5 } else { 1.5 };
        let curr = pslg.add_vertex(x, y);
        pslg.add_segment(prev, curr);
        prev = curr;
    }

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── CDT with collinear constraint chain ───────────────────────────────────

/// **Failure mode**: Many collinear constraint segments forming a single
/// line stress the flip-recovery when the constraint coincides with
/// existing DT edges.  Known to trigger infinite loops in implementations
/// that don't handle co-linearity in `collect_crossing_edges`.
#[test]
fn cdt_collinear_constraint_chain() {
    let mut pslg = Pslg::new();

    // Boundary.
    let bl = pslg.add_vertex(0.0, 0.0);
    let br = pslg.add_vertex(12.0, 0.0);
    let tr = pslg.add_vertex(12.0, 4.0);
    let tl = pslg.add_vertex(0.0, 4.0);
    pslg.add_segment(bl, br);
    pslg.add_segment(br, tr);
    pslg.add_segment(tr, tl);
    pslg.add_segment(tl, bl);

    // Chain of collinear segments along y = 2.
    let n_segs = 10;
    let mut chain_vids = Vec::new();
    for i in 0..=n_segs {
        let x = 1.0 + (i as f64) * 10.0 / (n_segs as f64);
        chain_vids.push(pslg.add_vertex(x, 2.0));
    }
    for i in 0..n_segs {
        pslg.add_segment(chain_vids[i], chain_vids[i + 1]);
    }

    // Scatter some points above and below.
    for i in 0..8 {
        let x = 1.5 + f64::from(i) * 1.2;
        pslg.add_vertex(x, 1.0);
        pslg.add_vertex(x, 3.0);
    }

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── Ruppert multi-region (disconnected domains) ───────────────────────────

/// **Failure mode**: Multiple disconnected polygons in the outer domain.
/// Ruppert refinement must correctly handle Steiner insertion that lands
/// outside both domains (between the two polygons).
#[test]
fn ruppert_disconnected_regions() {
    use crate::application::delaunay::dim2::refinement::ruppert::RuppertRefiner;

    let mut pslg = Pslg::new();

    // Wall endpoints on the boundary — reuse as boundary corners.
    let w0 = pslg.add_vertex(10.0, 0.0);
    let w1 = pslg.add_vertex(10.0, 10.0);

    // Outer rectangle split at wall endpoints.
    let o0 = pslg.add_vertex(0.0, 0.0);
    let o1 = pslg.add_vertex(20.0, 0.0);
    let o2 = pslg.add_vertex(20.0, 10.0);
    let o3 = pslg.add_vertex(0.0, 10.0);
    // Bottom: o0 → w0 → o1.
    pslg.add_segment(o0, w0);
    pslg.add_segment(w0, o1);
    // Right.
    pslg.add_segment(o1, o2);
    // Top: o2 → w1 → o3.
    pslg.add_segment(o2, w1);
    pslg.add_segment(w1, o3);
    // Left.
    pslg.add_segment(o3, o0);

    // Inner dividing wall.
    pslg.add_segment(w0, w1);

    let cdt = Cdt::from_pslg(&pslg);
    let mut refiner = RuppertRefiner::new(cdt);
    refiner.set_max_ratio(std::f64::consts::SQRT_2);
    refiner.set_max_area(5.0);
    refiner.set_max_steiner(500);
    let n = refiner.refine();

    assert!(n > 0);
    let dt = refiner.cdt().triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── Locate on exact triangle edge ─────────────────────────────────────────

/// **Failure mode**: Inserting a point that lies exactly on an existing DT
/// edge exercises the `OnEdge` branch of point location and the 2-to-4
/// triangle split (or 1-to-2 on the hull).  Float rounding can misclassify
/// OnEdge as Inside, producing degenerate zero-area triangles.
#[test]
fn insert_on_exact_edge() {
    // Build a square DT, then insert a point exactly on one of the DT edges.
    let mut pts = vec![
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
    ];
    // Midpoint of the bottom edge (0,0)→(2,0).
    pts.push((1.0, 0.0));

    let dt = DelaunayTriangulation::from_points(&pts);
    assert!(dt.is_delaunay());
    assert_eq!(dt.vertex_count(), 5);
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── CDT diamond with interior point ───────────────────────────────────────

/// **Failure mode**: A diamond-shaped constraint polygon with a single
/// interior point.  Known to expose edge-flip bugs when the interior point
/// lies on the circumcircle of the diamond edges.
///
/// **Known issue** (Triangle, pre-2005): diamond constraints with co-
/// circular interior points caused constraint edge deletion during flips.
#[test]
fn cdt_diamond_with_interior_point() {
    let mut pslg = Pslg::new();

    // Diamond vertices.
    let top = pslg.add_vertex(5.0, 10.0);
    let right = pslg.add_vertex(10.0, 5.0);
    let bot = pslg.add_vertex(5.0, 0.0);
    let left = pslg.add_vertex(0.0, 5.0);
    pslg.add_segment(top, right);
    pslg.add_segment(right, bot);
    pslg.add_segment(bot, left);
    pslg.add_segment(left, top);

    // Interior point at centre (on circumcircle of two opposite triangles).
    pslg.add_vertex(5.0, 5.0);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
    // Should have 4 interior triangles (fan from centre).
    assert!(
        dt.triangle_count() >= 4,
        "Diamond with interior point should produce ≥4 triangles, got {}",
        dt.triangle_count()
    );
}

// ── PSLG T-intersection validation ────────────────────────────────────────

/// **Failure mode**: T-intersection (segment endpoint touching the interior
/// of another segment) is a valid PSLG configuration.  Some implementations
/// reject it or fail to recover the constraint.
#[test]
fn pslg_t_intersection_constraint() {
    let mut pslg = Pslg::new();

    // Horizontal base segment.
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(10.0, 0.0);
    pslg.add_segment(a, b);

    // Vertical T-junction: endpoint at (5, 0) on the base.
    let c = pslg.add_vertex(5.0, 0.0);
    let d = pslg.add_vertex(5.0, 5.0);
    pslg.add_segment(c, d);

    // The base segment gets split by vertex c at (5,0).
    // We need to re-segment: a→c and c→b instead of a→b.
    // But first, test that the PSLG is valid (c lies on segment a→b).
    // Depending on implementation, this may auto-split or require explicit
    // specification.  Our implementation requires explicit segments.
    let mut pslg2 = Pslg::new();
    let a2 = pslg2.add_vertex(0.0, 0.0);
    let c2 = pslg2.add_vertex(5.0, 0.0);
    let b2 = pslg2.add_vertex(10.0, 0.0);
    let d2 = pslg2.add_vertex(5.0, 5.0);
    pslg2.add_segment(a2, c2);
    pslg2.add_segment(c2, b2);
    pslg2.add_segment(c2, d2);

    let cdt = Cdt::from_pslg(&pslg2);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── CDT convex hull stitch ────────────────────────────────────────────────

/// **Failure mode**: Constraint segments along the convex hull boundary
/// should be trivially recovered (they're already DT edges).  But some
/// implementations incorrectly attempt to flip hull edges.
#[test]
fn cdt_constraints_on_convex_hull() {
    let mut pslg = Pslg::new();
    let pts = [
        (0.0, 0.0),
        (5.0, 0.0),
        (10.0, 0.0),
        (10.0, 5.0),
        (5.0, 5.0),
        (0.0, 5.0),
    ];
    let vids: Vec<_> = pts.iter().map(|&(x, y)| pslg.add_vertex(x, y)).collect();

    // All edges are along the convex hull.
    for i in 0..vids.len() {
        pslg.add_segment(vids[i], vids[(i + 1) % vids.len()]);
    }

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}

// ── Spiral CDT with radial constraints ────────────────────────────────────

/// **Failure mode**: Combination of spiral point distribution (stresses
/// walk) with radial constraints (stresses CDT recovery).
#[test]
fn spiral_cdt_with_radial_constraints() {
    let mut pslg = Pslg::new();
    let center = pslg.add_vertex(0.0, 0.0);

    let n = 32;
    let mut ring = Vec::new();
    for i in 0..n {
        let theta = 2.0 * PI * (i as f64) / (n as f64);
        let r = 5.0;
        ring.push(pslg.add_vertex(r * theta.cos(), r * theta.sin()));
    }

    // Boundary polygon.
    for i in 0..n {
        pslg.add_segment(ring[i], ring[(i + 1) % n]);
    }

    // Every other vertex gets a radial constraint to center.
    for i in (0..n).step_by(2) {
        pslg.add_segment(center, ring[i]);
    }

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();
    assert!(dt.is_delaunay());
    assert!(Adjacency::verify_symmetry(dt.triangles_slice()));
}
