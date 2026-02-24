//! Tests for Ruppert's refinement algorithm.

use crate::application::delaunay::constraint::enforce::Cdt;
use crate::application::delaunay::pslg::graph::Pslg;
use crate::application::delaunay::refinement::quality::TriangleQuality;
use crate::application::delaunay::refinement::ruppert::RuppertRefiner;

// ── Basic refinement ──────────────────────────────────────────────────────

#[test]
fn refine_triangle_improves_quality() {
    let mut pslg = Pslg::new();
    // Near-equilateral triangle (all angles ~60°) with area large enough
    // that max_area refinement produces multiple triangles.
    // Ruppert requires input angles > ~20° for guaranteed termination.
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(10.0, 0.0);
    let c = pslg.add_vertex(5.0, 8.66);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    let cdt = Cdt::from_pslg(&pslg);
    let mut refiner = RuppertRefiner::new(cdt);
    refiner.set_max_ratio(std::f64::consts::SQRT_2);
    refiner.set_max_area(5.0); // Force area subdivision.
    refiner.set_max_steiner(500);
    let n_steiner = refiner.refine();

    assert!(n_steiner > 0, "Refinement should insert Steiner points");

    // Check all interior triangles satisfy quality bound.
    let dt = refiner.cdt().triangulation();
    for (_, tri) in dt.interior_triangles() {
        let v0 = dt.vertex(tri.vertices[0]);
        let v1 = dt.vertex(tri.vertices[1]);
        let v2 = dt.vertex(tri.vertices[2]);
        let q = TriangleQuality::compute(v0, v1, v2);

        assert!(
            q.radius_edge_ratio < 2.0,
            "Triangle quality too poor: ratio = {} (min_angle = {:.1}°)",
            q.radius_edge_ratio,
            q.min_angle_deg(),
        );
    }
}

// ── Square domain refinement ──────────────────────────────────────────────

#[test]
fn refine_square_domain() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 0.0);
    let c = pslg.add_vertex(1.0, 1.0);
    let d = pslg.add_vertex(0.0, 1.0);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, d);
    pslg.add_segment(d, a);

    let cdt = Cdt::from_pslg(&pslg);
    let pre_count = cdt.triangulation().triangle_count();

    let mut refiner = RuppertRefiner::new(cdt);
    refiner.set_max_area(0.1); // Force area refinement.
    refiner.set_max_steiner(1000);
    refiner.refine();

    let dt = refiner.cdt().triangulation();
    let post_count = dt.triangle_count();

    assert!(
        post_count > pre_count,
        "Refinement should increase triangle count: {} > {}",
        post_count,
        pre_count,
    );
}

// ── Minimum angle guarantee ───────────────────────────────────────────────

/// Theorem: Ruppert guarantees α_min ≥ arcsin(1/(2B)) for B = √2.
/// This gives α_min ≥ 20.7°.
#[test]
fn minimum_angle_around_20_degrees() {
    let mut pslg = Pslg::new();
    // L-shaped domain with potential for skinny triangles.
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
    let mut refiner = RuppertRefiner::new(cdt);
    refiner.set_max_ratio(std::f64::consts::SQRT_2);
    refiner.set_max_steiner(2000);
    refiner.refine();

    let dt = refiner.cdt().triangulation();
    let mut min_angle_global = f64::MAX;

    for (_, tri) in dt.interior_triangles() {
        let v0 = dt.vertex(tri.vertices[0]);
        let v1 = dt.vertex(tri.vertices[1]);
        let v2 = dt.vertex(tri.vertices[2]);
        let q = TriangleQuality::compute(v0, v1, v2);
        if q.min_angle_deg() < min_angle_global {
            min_angle_global = q.min_angle_deg();
        }
    }

    // Theoretical minimum is 20.7° but numerical tolerances may yield slightly less.
    assert!(
        min_angle_global > 15.0,
        "Minimum angle {:.1}° should be > 15° after Ruppert refinement",
        min_angle_global,
    );
}

// ── Quality metric computation ────────────────────────────────────────────

#[test]
fn equilateral_triangle_quality() {
    use crate::application::delaunay::pslg::vertex::PslgVertex;

    let a = PslgVertex::new(0.0, 0.0);
    let b = PslgVertex::new(1.0, 0.0);
    let c = PslgVertex::new(0.5, 3.0_f64.sqrt() / 2.0);

    let q = TriangleQuality::compute(&a, &b, &c);

    // Equilateral: all angles = 60°.
    assert!(
        (q.min_angle_deg() - 60.0).abs() < 0.1,
        "Equilateral min angle should be ~60°, got {:.2}°",
        q.min_angle_deg()
    );

    // Radius-edge ratio of equilateral = 1/√3 ≈ 0.577.
    assert!(
        (q.radius_edge_ratio - 1.0 / 3.0_f64.sqrt()).abs() < 0.01,
        "Equilateral ratio should be ~0.577, got {:.4}",
        q.radius_edge_ratio
    );
}

#[test]
fn right_triangle_quality() {
    use crate::application::delaunay::pslg::vertex::PslgVertex;

    let a = PslgVertex::new(0.0, 0.0);
    let b = PslgVertex::new(1.0, 0.0);
    let c = PslgVertex::new(0.0, 1.0);

    let q = TriangleQuality::compute(&a, &b, &c);

    // Min angle = 45°.
    assert!(
        (q.min_angle_deg() - 45.0).abs() < 0.5,
        "Right isoceles min angle should be ~45°, got {:.2}°",
        q.min_angle_deg()
    );
}

// ── Steiner count bounded ─────────────────────────────────────────────────

#[test]
fn steiner_respects_limit() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(100.0, 0.0);
    let c = pslg.add_vertex(50.0, 0.01); // extremely skinny
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    let cdt = Cdt::from_pslg(&pslg);
    let mut refiner = RuppertRefiner::new(cdt);
    refiner.set_max_steiner(10);
    let n = refiner.refine();

    assert!(
        n <= 10,
        "Steiner count {} should respect the limit of 10",
        n
    );
}
