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

// ── resolve_crossings ──────────────────────────────────────────────────────

#[test]
fn resolve_crossings_splits_crossing_pair() {
    let mut pslg = Pslg::new();
    // Two crossing diagonals: (0,0)→(1,1) and (0,1)→(1,0), cross at (0.5, 0.5).
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 1.0);
    let c = pslg.add_vertex(0.0, 1.0);
    let d = pslg.add_vertex(1.0, 0.0);
    pslg.add_segment(a, b);
    pslg.add_segment(c, d);

    assert!(pslg.validate().is_err(), "should be invalid before resolve");

    pslg.resolve_crossings();

    assert!(pslg.validate().is_ok(), "should be valid after resolve");
    assert_eq!(pslg.vertex_count(), 5, "4 original + 1 intersection vertex");
    assert_eq!(pslg.segment_count(), 4, "2 original → 4 sub-segments");

    // The intersection vertex was appended last: index 4.
    let xv = pslg.vertices()[4];
    assert!((xv.x - 0.5).abs() < 1e-12, "intersection x ≈ 0.5, got {}", xv.x);
    assert!((xv.y - 0.5).abs() < 1e-12, "intersection y ≈ 0.5, got {}", xv.y);
}

#[test]
fn resolve_crossings_noop_on_valid_pslg() {
    let mut pslg = Pslg::new();
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(1.0, 0.0);
    let c = pslg.add_vertex(1.0, 1.0);
    pslg.add_segment(a, b);
    pslg.add_segment(b, c);

    assert!(pslg.validate().is_ok());
    pslg.resolve_crossings();
    assert!(pslg.validate().is_ok(), "valid PSLG must stay valid");
    assert_eq!(pslg.vertex_count(), 3, "no new vertices on valid PSLG");
    assert_eq!(pslg.segment_count(), 2, "no new segments on valid PSLG");
}

#[test]
fn resolve_crossings_handles_multiple_crossings() {
    let mut pslg = Pslg::new();
    // A long horizontal segment crossed by two separate vertical segments.
    //   seg 0: (0,0)→(4,0)   long horizontal
    //   seg 1: (1,-1)→(1,1)  vertical at x=1, crosses seg 0 at (1, 0)
    //   seg 2: (3,-1)→(3,1)  vertical at x=3, crosses seg 0 at (3, 0)
    // The two crossings are at distinct points so no degenerate three-way
    // intersection arises.
    let a = pslg.add_vertex(0.0, 0.0);
    let b = pslg.add_vertex(4.0, 0.0);
    let c = pslg.add_vertex(1.0, -1.0);
    let d = pslg.add_vertex(1.0, 1.0);
    let e = pslg.add_vertex(3.0, -1.0);
    let f = pslg.add_vertex(3.0, 1.0);
    pslg.add_segment(a, b); // horizontal
    pslg.add_segment(c, d); // vertical at x=1
    pslg.add_segment(e, f); // vertical at x=3

    // At least one crossing before.
    assert!(pslg.validate().is_err());

    pslg.resolve_crossings();

    // After resolution the PSLG must be valid.
    assert!(pslg.validate().is_ok(), "all crossings should be resolved");
    // 6 original + 2 intersection vertices.
    assert_eq!(pslg.vertex_count(), 8);
    // 1 horizontal → 3 sub-segments; 2 verticals → 2×2 = 4 sub-segments = 7 total.
    assert_eq!(pslg.segment_count(), 7);
}
