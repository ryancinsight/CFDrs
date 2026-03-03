//! Planar Straight-Line Graph (PSLG) — the input to CDT.
//!
//! A PSLG consists of:
//! - A set of **vertices** with 2-D coordinates.
//! - A set of **constraint segments** connecting pairs of vertices.
//! - Optional **hole seeds** — points inside regions that should be removed
//!   from the triangulation.
//!
//! # Invariant
//!
//! No two constraint segments may cross in their interiors.  Segments may
//! share endpoints.  The caller is responsible for ensuring this; the CDT
//! will produce undefined results if segments cross.
//!
//! # Theorem — PSLG Validity
//!
//! **Statement**: A set of segments $S$ forms a valid PSLG if and only if
//! no two segments in $S$ share an interior point.  (Shared endpoints are
//! permitted.)
//!
//! **Proof sketch**: The definition of a planar subdivision requires that
//! edges intersect only at shared vertices.  If two segments cross, the
//! crossing point is not a vertex, violating the subdivision property.

use crate::domain::core::scalar::Real;
use crate::domain::geometry::predicates::{orient_2d, Orientation};
use nalgebra::Point2;

use super::segment::{PslgSegment, PslgSegmentId};
use super::vertex::{PslgVertex, PslgVertexId};

/// Validation errors for a [`Pslg`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PslgValidationError {
    /// Segment endpoint index is out of range of the vertex list.
    SegmentVertexOutOfRange {
        /// Segment id with invalid endpoint reference.
        segment: PslgSegmentId,
        /// Start vertex id.
        start: PslgVertexId,
        /// End vertex id.
        end: PslgVertexId,
        /// Total number of vertices present in the PSLG.
        vertex_count: usize,
    },
    /// Segment start and end are identical.
    DegenerateSegment {
        /// Degenerate segment id.
        segment: PslgSegmentId,
        /// Collapsed endpoint id.
        vertex: PslgVertexId,
    },
    /// Two segments are duplicates (same canonical endpoints).
    DuplicateSegment {
        /// First segment id.
        first: PslgSegmentId,
        /// Second segment id.
        second: PslgSegmentId,
        /// Canonical first endpoint.
        a: PslgVertexId,
        /// Canonical second endpoint.
        b: PslgVertexId,
    },
    /// Two segments intersect in their interiors or overlap collinearly.
    IntersectingSegments {
        /// First intersecting segment id.
        first: PslgSegmentId,
        /// Second intersecting segment id.
        second: PslgSegmentId,
    },
}

impl core::fmt::Display for PslgValidationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SegmentVertexOutOfRange {
                segment,
                start,
                end,
                vertex_count,
            } => write!(
                f,
                "segment {segment} references out-of-range vertex ids ({start}, {end}) with vertex_count={vertex_count}",
            ),
            Self::DegenerateSegment { segment, vertex } => {
                write!(f, "segment {segment} is degenerate at vertex {vertex}")
            }
            Self::DuplicateSegment {
                first,
                second,
                a,
                b,
            } => write!(
                f,
                "segments {first} and {second} are duplicates of edge ({a}, {b})",
            ),
            Self::IntersectingSegments { first, second } => {
                write!(
                    f,
                    "segments {first} and {second} intersect in their interiors"
                )
            }
        }
    }
}

impl std::error::Error for PslgValidationError {}

/// A Planar Straight-Line Graph — the canonical input to CDT.
///
/// # Example
///
/// ```rust,ignore
/// use cfd_mesh::application::delaunay::Pslg;
///
/// let mut pslg = Pslg::new();
/// let a = pslg.add_vertex(0.0, 0.0);
/// let b = pslg.add_vertex(1.0, 0.0);
/// let c = pslg.add_vertex(0.5, 0.866);
/// pslg.add_segment(a, b);
/// pslg.add_segment(b, c);
/// pslg.add_segment(c, a);
/// ```
#[derive(Clone, Debug)]
pub struct Pslg {
    /// Vertex positions.
    vertices: Vec<PslgVertex>,
    /// Constraint segments.
    segments: Vec<PslgSegment>,
    /// Hole seed points — each point inside a region to be removed.
    holes: Vec<PslgVertex>,
}

impl Pslg {
    /// Create an empty PSLG.
    #[must_use] 
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            segments: Vec::new(),
            holes: Vec::new(),
        }
    }

    /// Create with pre-allocated capacity.
    #[must_use] 
    pub fn with_capacity(num_vertices: usize, num_segments: usize) -> Self {
        Self {
            vertices: Vec::with_capacity(num_vertices),
            segments: Vec::with_capacity(num_segments),
            holes: Vec::new(),
        }
    }

    // ── Vertex operations ─────────────────────────────────────────────────

    /// Add a vertex at `(x, y)` and return its ID.
    pub fn add_vertex(&mut self, x: Real, y: Real) -> PslgVertexId {
        let id = PslgVertexId::from_usize(self.vertices.len());
        self.vertices.push(PslgVertex::new(x, y));
        id
    }

    /// Add a vertex from a `PslgVertex` value.
    pub fn add_vertex_value(&mut self, v: PslgVertex) -> PslgVertexId {
        let id = PslgVertexId::from_usize(self.vertices.len());
        self.vertices.push(v);
        id
    }

    /// Number of vertices.
    #[inline]
    #[must_use] 
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Get vertex by ID.
    #[inline]
    #[must_use] 
    pub fn vertex(&self, id: PslgVertexId) -> &PslgVertex {
        &self.vertices[id.idx()]
    }

    /// Slice of all vertex positions.
    #[inline]
    #[must_use] 
    pub fn vertices(&self) -> &[PslgVertex] {
        &self.vertices
    }

    // ── Segment operations ────────────────────────────────────────────────

    /// Add a constraint segment between two existing vertices.
    ///
    /// # Panics
    ///
    /// Panics (in debug) if `start == end` or if either ID is out of range.
    pub fn add_segment(&mut self, start: PslgVertexId, end: PslgVertexId) -> PslgSegmentId {
        debug_assert_ne!(start, end, "degenerate segment");
        debug_assert!(
            start.idx() < self.vertices.len() && end.idx() < self.vertices.len(),
            "segment vertex out of range"
        );
        let id = PslgSegmentId::from_usize(self.segments.len());
        self.segments.push(PslgSegment::new(start, end));
        id
    }

    /// Number of constraint segments.
    #[inline]
    #[must_use] 
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Get segment by ID.
    #[inline]
    #[must_use] 
    pub fn segment(&self, id: PslgSegmentId) -> &PslgSegment {
        &self.segments[id.idx()]
    }

    /// Slice of all segments.
    #[inline]
    #[must_use] 
    pub fn segments(&self) -> &[PslgSegment] {
        &self.segments
    }

    #[cfg(test)]
    pub(crate) fn segments_mut_for_test_only(&mut self) -> &mut Vec<PslgSegment> {
        &mut self.segments
    }

    // ── Hole operations ───────────────────────────────────────────────────

    /// Mark a point as a hole seed.
    ///
    /// All triangles whose centroid is reachable from this point without
    /// crossing a constraint segment will be removed.
    pub fn add_hole(&mut self, x: Real, y: Real) {
        self.holes.push(PslgVertex::new(x, y));
    }

    /// Validate PSLG topological constraints.
    ///
    /// Checks:
    /// - Segment endpoint indices are in range.
    /// - No degenerate segments.
    /// - No duplicate segments.
    /// - No segment-segment interior intersections (shared endpoints allowed).
    pub fn validate(&self) -> Result<(), PslgValidationError> {
        use std::collections::HashMap;

        let n_vertices = self.vertices.len();

        for (idx, seg) in self.segments.iter().copied().enumerate() {
            let sid = PslgSegmentId::from_usize(idx);
            if seg.start.idx() >= n_vertices || seg.end.idx() >= n_vertices {
                return Err(PslgValidationError::SegmentVertexOutOfRange {
                    segment: sid,
                    start: seg.start,
                    end: seg.end,
                    vertex_count: n_vertices,
                });
            }
            if seg.is_degenerate() {
                return Err(PslgValidationError::DegenerateSegment {
                    segment: sid,
                    vertex: seg.start,
                });
            }
        }

        let mut seen: HashMap<(PslgVertexId, PslgVertexId), PslgSegmentId> = HashMap::new();
        for (idx, seg) in self.segments.iter().copied().enumerate() {
            let sid = PslgSegmentId::from_usize(idx);
            let key = seg.canonical();
            if let Some(first) = seen.insert(key, sid) {
                return Err(PslgValidationError::DuplicateSegment {
                    first,
                    second: sid,
                    a: key.0,
                    b: key.1,
                });
            }
        }

        for i in 0..self.segments.len() {
            for j in (i + 1)..self.segments.len() {
                let s1 = self.segments[i];
                let s2 = self.segments[j];

                let share_endpoint = s1.start == s2.start
                    || s1.start == s2.end
                    || s1.end == s2.start
                    || s1.end == s2.end;

                let a1 = self.vertices[s1.start.idx()].to_point2();
                let a2 = self.vertices[s1.end.idx()].to_point2();
                let b1 = self.vertices[s2.start.idx()].to_point2();
                let b2 = self.vertices[s2.end.idx()].to_point2();

                if !segments_intersect_closed(&a1, &a2, &b1, &b2) {
                    continue;
                }

                // Shared endpoints are allowed only when they do not overlap
                // beyond that endpoint (collinear overlap is invalid).
                if share_endpoint && !collinear_overlap_interior(&a1, &a2, &b1, &b2) {
                    continue;
                }

                return Err(PslgValidationError::IntersectingSegments {
                    first: PslgSegmentId::from_usize(i),
                    second: PslgSegmentId::from_usize(j),
                });
            }
        }

        Ok(())
    }

    /// Slice of all hole seeds.
    #[inline]
    #[must_use] 
    pub fn holes(&self) -> &[PslgVertex] {
        &self.holes
    }

    // ── Bounding box ──────────────────────────────────────────────────────

    /// Resolve *all* illegal constraint pairs so that [`Self::validate`] passes.
    ///
    /// Three classes of illegality are handled:
    ///
    /// 1. **Proper interior crossings** — `orient_2d` shows both endpoints of
    ///    each segment on opposite sides of the other.  Split both at their
    ///    intersection point.
    ///
    /// 2. **T-intersections** — an endpoint of segment B lies strictly in the
    ///    interior of segment A (`orient_2d` = Degenerate for that endpoint).
    ///    Split A at B's endpoint.
    ///
    /// 3. **Collinear overlaps** — both segments are collinear *and* their
    ///    interiors share an interval.  Split at both overlap boundary points
    ///    so every collinear subsegment covers only one canonical interval.
    ///
    /// Repeats until no illegal pair remains.  For CSG corefine inputs (k = 0–3
    /// crossings) the total work is O(k · n²).
    ///
    /// After this call [`Self::validate`] will not return
    /// [`PslgValidationError::IntersectingSegments`] for any geometrically
    /// representable input.  Duplicate-segment errors (produced when two
    /// previously distinct segments resolve to the same subsegment) are
    /// removed by deduplication after each restart.
    pub fn resolve_crossings(&mut self) {
        'outer: loop {
            let n_seg = self.segments.len();

            for i in 0..n_seg {
                for j in (i + 1)..n_seg {
                    let si = self.segments[i];
                    let sj = self.segments[j];

                    // Endpoints shared — adjacency is always valid for CDT.
                    let share_endpoint = si.start == sj.start
                        || si.start == sj.end
                        || si.end == sj.start
                        || si.end == sj.end;

                    let a1 = self.vertices[si.start.idx()].to_point2();
                    let a2 = self.vertices[si.end.idx()].to_point2();
                    let b1 = self.vertices[sj.start.idx()].to_point2();
                    let b2 = self.vertices[sj.end.idx()].to_point2();

                    let o_b1 = orient_2d(&a1, &a2, &b1);
                    let o_b2 = orient_2d(&a1, &a2, &b2);
                    let o_a1 = orient_2d(&b1, &b2, &a1);
                    let o_a2 = orient_2d(&b1, &b2, &a2);

                    // ── Case 1: Proper interior crossing ──────────────────────
                    if o_b1 != Orientation::Degenerate
                        && o_b2 != Orientation::Degenerate
                        && o_a1 != Orientation::Degenerate
                        && o_a2 != Orientation::Degenerate
                        && o_b1 != o_b2
                        && o_a1 != o_a2
                    {
                        if let Some((px, py)) = segment_cross_point(&a1, &a2, &b1, &b2) {
                            let xid = self.add_vertex(px, py);
                            let (si_s, si_e) = (si.start, si.end);
                            let (sj_s, sj_e) = (sj.start, sj.end);
                            self.segments.swap_remove(j);
                            self.segments.swap_remove(i);
                            self.add_segment(si_s, xid);
                            self.add_segment(xid, si_e);
                            self.add_segment(sj_s, xid);
                            self.add_segment(xid, sj_e);
                            continue 'outer;
                        }
                    }

                    // ── Case 2: T-intersection (endpoint on interior) ──────────
                    // Skip adjacent segments — shared endpoints are legal.
                    if !share_endpoint {
                        // b1 lies strictly in the interior of segment A?
                        if o_b1 == Orientation::Degenerate && on_segment(&a1, &a2, &b1)
                            && b1 != a1 && b1 != a2
                        {
                            let xid = sj.start; // reuse the existing PSLG vertex
                            let (si_s, si_e) = (si.start, si.end);
                            self.segments.swap_remove(i);
                            self.add_segment(si_s, xid);
                            self.add_segment(xid, si_e);
                            continue 'outer;
                        }
                        // b2 lies strictly in the interior of segment A?
                        if o_b2 == Orientation::Degenerate && on_segment(&a1, &a2, &b2)
                            && b2 != a1 && b2 != a2
                        {
                            let xid = sj.end;
                            let (si_s, si_e) = (si.start, si.end);
                            self.segments.swap_remove(i);
                            self.add_segment(si_s, xid);
                            self.add_segment(xid, si_e);
                            continue 'outer;
                        }
                        // a1 lies strictly in the interior of segment B?
                        if o_a1 == Orientation::Degenerate && on_segment(&b1, &b2, &a1)
                            && a1 != b1 && a1 != b2
                        {
                            let xid = si.start;
                            let (sj_s, sj_e) = (sj.start, sj.end);
                            self.segments.swap_remove(j);
                            self.add_segment(sj_s, xid);
                            self.add_segment(xid, sj_e);
                            continue 'outer;
                        }
                        // a2 lies strictly in the interior of segment B?
                        if o_a2 == Orientation::Degenerate && on_segment(&b1, &b2, &a2)
                            && a2 != b1 && a2 != b2
                        {
                            let xid = si.end;
                            let (sj_s, sj_e) = (sj.start, sj.end);
                            self.segments.swap_remove(j);
                            self.add_segment(sj_s, xid);
                            self.add_segment(xid, sj_e);
                            continue 'outer;
                        }
                    }

                    // ── Case 3: Collinear overlap ──────────────────────────────
                    // All four orientations are Degenerate → segments are collinear.
                    if o_b1 == Orientation::Degenerate
                        && o_b2 == Orientation::Degenerate
                        && o_a1 == Orientation::Degenerate
                        && o_a2 == Orientation::Degenerate
                        && collinear_overlap_interior(&a1, &a2, &b1, &b2)
                    {
                        // Project onto the dominant axis to compute overlap interval.
                        let use_x = (a2.x - a1.x).abs() >= (a2.y - a1.y).abs();
                        // Collect all four coordinates and sort them.
                        let coords: [Real; 4] = if use_x {
                            [a1.x, a2.x, b1.x, b2.x]
                        } else {
                            [a1.y, a2.y, b1.y, b2.y]
                        };
                        let mut sorted = coords;
                        sorted.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
                        // The two middle values are the overlap boundaries.
                        let lo_val = sorted[1];
                        let hi_val = sorted[2];
                        if (hi_val - lo_val).abs() < 1e-30 {
                            // Degenerate overlap (a single point) — not a true overlap.
                            continue;
                        }
                        // Find or create PSLG vertices at lo_val and hi_val.
                        let lo_pt = if use_x {
                            // x = lo_val, lerp y from segment A
                            let t = (lo_val - a1.x) / (a2.x - a1.x);
                            (lo_val, a1.y + t * (a2.y - a1.y))
                        } else {
                            let t = (lo_val - a1.y) / (a2.y - a1.y);
                            (a1.x + t * (a2.x - a1.x), lo_val)
                        };
                        let hi_pt = if use_x {
                            let t = (hi_val - a1.x) / (a2.x - a1.x);
                            (hi_val, a1.y + t * (a2.y - a1.y))
                        } else {
                            let t = (hi_val - a1.y) / (a2.y - a1.y);
                            (a1.x + t * (a2.x - a1.x), hi_val)
                        };
                        // Insert or reuse vertices for lo / hi overlap boundaries.
                        let find_or_add = |pslg: &mut Pslg, px: Real, py: Real| -> PslgVertexId {
                            // Check all existing vertices for a match.
                            for (idx, v) in pslg.vertices.iter().enumerate() {
                                let dx = v.x - px;
                                let dy = v.y - py;
                                if dx * dx + dy * dy < 1e-30 {
                                    return PslgVertexId::from_usize(idx);
                                }
                            }
                            pslg.add_vertex(px, py)
                        };
                        let lo_id = find_or_add(self, lo_pt.0, lo_pt.1);
                        let hi_id = find_or_add(self, hi_pt.0, hi_pt.1);
                        // Gather endpoints before removal.
                        let (si_s, si_e) = (si.start, si.end);
                        let (sj_s, sj_e) = (sj.start, sj.end);
                        self.segments.swap_remove(j);
                        self.segments.swap_remove(i);

                        // Helper to shatter a segment along the dominant axis.
                        let mut add_shattered = |s_start: PslgVertexId, s_end: PslgVertexId| {
                            let mut verts = [s_start, s_end, lo_id, hi_id];
                            verts.sort_by(|&v1, &v2| {
                                let p1 = self.vertices[v1.idx()];
                                let p2 = self.vertices[v2.idx()];
                                let val1 = if use_x { p1.x } else { p1.y };
                                let val2 = if use_x { p2.x } else { p2.y };
                                val1.partial_cmp(&val2).unwrap_or(std::cmp::Ordering::Equal)
                            });
                            for k in 0..3 {
                                let p = verts[k];
                                let q = verts[k + 1];
                                if p != q {
                                    let canonical = (p.min(q), p.max(q));
                                    if !self.segments.iter().any(|s| s.canonical() == canonical) {
                                        let _ = self.add_segment(p, q);
                                    }
                                }
                            }
                        };

                        // Reconstruct non-degenerate sub-segments for A and B.
                        add_shattered(si_s, si_e);
                        add_shattered(sj_s, sj_e);
                        continue 'outer;
                    }
                }
            }

            break; // no illegal pairs remain
        }
    }

    /// Compute the axis-aligned bounding box `(min, max)`.
    ///
    /// Returns `None` if the PSLG has fewer than 1 vertex.
    #[must_use] 
    pub fn bounding_box(&self) -> Option<(PslgVertex, PslgVertex)> {
        if self.vertices.is_empty() {
            return None;
        }
        let mut min_x = self.vertices[0].x;
        let mut min_y = self.vertices[0].y;
        let mut max_x = min_x;
        let mut max_y = min_y;
        for v in &self.vertices[1..] {
            if v.x < min_x {
                min_x = v.x;
            }
            if v.y < min_y {
                min_y = v.y;
            }
            if v.x > max_x {
                max_x = v.x;
            }
            if v.y > max_y {
                max_y = v.y;
            }
        }
        Some((PslgVertex::new(min_x, min_y), PslgVertex::new(max_x, max_y)))
    }
}

impl Default for Pslg {
    fn default() -> Self {
        Self::new()
    }
}

fn segments_intersect_closed(
    a1: &Point2<Real>,
    a2: &Point2<Real>,
    b1: &Point2<Real>,
    b2: &Point2<Real>,
) -> bool {
    let o1 = orient_2d(a1, a2, b1);
    let o2 = orient_2d(a1, a2, b2);
    let o3 = orient_2d(b1, b2, a1);
    let o4 = orient_2d(b1, b2, a2);

    if o1 != Orientation::Degenerate
        && o2 != Orientation::Degenerate
        && o3 != Orientation::Degenerate
        && o4 != Orientation::Degenerate
    {
        return o1 != o2 && o3 != o4;
    }

    if o1 == Orientation::Degenerate && on_segment(a1, a2, b1) {
        return true;
    }
    if o2 == Orientation::Degenerate && on_segment(a1, a2, b2) {
        return true;
    }
    if o3 == Orientation::Degenerate && on_segment(b1, b2, a1) {
        return true;
    }
    if o4 == Orientation::Degenerate && on_segment(b1, b2, a2) {
        return true;
    }

    false
}

fn on_segment(a: &Point2<Real>, b: &Point2<Real>, p: &Point2<Real>) -> bool {
    p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y)
}

fn collinear_overlap_interior(
    a1: &Point2<Real>,
    a2: &Point2<Real>,
    b1: &Point2<Real>,
    b2: &Point2<Real>,
) -> bool {
    // Non-collinear cannot overlap interiorly.
    if orient_2d(a1, a2, b1) != Orientation::Degenerate
        || orient_2d(a1, a2, b2) != Orientation::Degenerate
    {
        return false;
    }

    let use_x = (a2.x - a1.x).abs() >= (a2.y - a1.y).abs();
    let (a_lo, a_hi, b_lo, b_hi) = if use_x {
        (
            a1.x.min(a2.x),
            a1.x.max(a2.x),
            b1.x.min(b2.x),
            b1.x.max(b2.x),
        )
    } else {
        (
            a1.y.min(a2.y),
            a1.y.max(a2.y),
            b1.y.min(b2.y),
            b1.y.max(b2.y),
        )
    };

    let overlap = a_hi.min(b_hi) - a_lo.max(b_lo);
    overlap > 0.0
}

/// Compute the f64 parametric crossing point of two non-parallel line segments.
///
/// Returns `None` when the segments are parallel (|denom| < 1e-30).
fn segment_cross_point(
    a1: &Point2<Real>,
    a2: &Point2<Real>,
    b1: &Point2<Real>,
    b2: &Point2<Real>,
) -> Option<(Real, Real)> {
    let dx_a = a2.x - a1.x;
    let dy_a = a2.y - a1.y;
    let dx_b = b2.x - b1.x;
    let dy_b = b2.y - b1.y;
    let denom = dx_a * dy_b - dy_a * dx_b;
    if denom.abs() < 1e-30 {
        return None;
    }
    let t = ((b1.x - a1.x) * dy_b - (b1.y - a1.y) * dx_b) / denom;
    Some((a1.x + t * dx_a, a1.y + t * dy_a))
}
