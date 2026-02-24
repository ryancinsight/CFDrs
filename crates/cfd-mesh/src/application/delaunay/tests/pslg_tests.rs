//! Tests for PSLG data structures.

use crate::application::delaunay::pslg::graph::Pslg;
use crate::application::delaunay::pslg::graph::PslgValidationError;
use crate::application::delaunay::pslg::segment::PslgSegment;
use crate::application::delaunay::pslg::vertex::{PslgVertex, PslgVertexId, GHOST_VERTEX};

// ── Vertex ────────────────────────────────────────────────────────────────

#[test]
fn vertex_creation() {
    let v = PslgVertex::new(3.14, 2.72);
    assert!((v.x - 3.14).abs() < 1e-15);
    assert!((v.y - 2.72).abs() < 1e-15);
}

#[test]
fn vertex_distance() {
    let a = PslgVertex::new(0.0, 0.0);
    let b = PslgVertex::new(3.0, 4.0);
    assert!((a.dist(&b) - 5.0).abs() < 1e-12);
    assert!((a.dist_sq(&b) - 25.0).abs() < 1e-12);
}

#[test]
fn vertex_midpoint() {
    let a = PslgVertex::new(0.0, 0.0);
    let b = PslgVertex::new(2.0, 4.0);
    let m = a.midpoint(&b);
    assert!((m.x - 1.0).abs() < 1e-15);
    assert!((m.y - 2.0).abs() < 1e-15);
}

#[test]
fn vertex_id_sentinel() {
    assert_eq!(GHOST_VERTEX.idx(), u32::MAX as usize);
}

// ── Segment ───────────────────────────────────────────────────────────────

#[test]
fn segment_canonical() {
    let s = PslgSegment::new(PslgVertexId::new(5), PslgVertexId::new(2));
    let (a, b) = s.canonical();
    assert!(a <= b, "canonical should order (min, max)");
}

// ── PSLG graph ────────────────────────────────────────────────────────────

#[test]
fn pslg_add_vertices_and_segments() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 0.0);
    let c = pslg.add_vertex(0.5, 1.0);

    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    pslg.add_segment(c, a);

    assert_eq!(pslg.vertices().len(), 3);
    assert_eq!(pslg.segments().len(), 3);
}

#[test]
fn pslg_bounding_box() {
    let mut pslg = Pslg::new();
    pslg.add_vertex(-5.0, -3.0);
    pslg.add_vertex(7.0, 11.0);
    pslg.add_vertex(2.0, 4.0);

    let (lo, hi) = pslg.bounding_box().expect("non-empty PSLG");
    assert!((lo.x - (-5.0)).abs() < 1e-15);
    assert!((lo.y - (-3.0)).abs() < 1e-15);
    assert!((hi.x - 7.0).abs() < 1e-15);
    assert!((hi.y - 11.0).abs() < 1e-15);
}

#[test]
fn pslg_add_hole() {
    let mut pslg = Pslg::new();
    pslg.add_hole(1.0, 1.0);
    assert_eq!(pslg.holes().len(), 1);
    assert!((pslg.holes()[0].x - 1.0).abs() < 1e-15);
}

#[test]
fn pslg_empty() {
    let pslg = Pslg::new();
    assert_eq!(pslg.vertices().len(), 0);
    assert_eq!(pslg.segments().len(), 0);
    assert_eq!(pslg.holes().len(), 0);
}

#[test]
fn pslg_validate_ok_shared_endpoint() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 0.0);
    let c = pslg.add_vertex(1.0, 1.0);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);
    assert!(pslg.validate().is_ok());
}

#[test]
fn pslg_validate_detects_duplicate_segments() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 0.0);
    pslg.add_segment(a, b);
    pslg.add_segment(b, a);

    let err = pslg
        .validate()
        .expect_err("expected duplicate segment error");
    assert!(matches!(err, PslgValidationError::DuplicateSegment { .. }));
}

#[test]
fn pslg_validate_detects_crossing_segments() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 1.0);
    let c = pslg.add_vertex(0.0, 1.0);
    let d = pslg.add_vertex(1.0, 0.0);

    pslg.add_segment(a, b);
    pslg.add_segment(c, d);

    let err = pslg.validate().expect_err("expected intersection error");
    assert!(matches!(
        err,
        PslgValidationError::IntersectingSegments { .. }
    ));
}

#[test]
fn pslg_validate_detects_degenerate_segment() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    // Bypass add_segment debug assertions to build an invalid PSLG explicitly.
    pslg.segments_mut_for_test_only()
        .push(PslgSegment::new(a, a));

    let err = pslg
        .validate()
        .expect_err("expected degenerate segment error");
    assert!(matches!(err, PslgValidationError::DegenerateSegment { .. }));
}
