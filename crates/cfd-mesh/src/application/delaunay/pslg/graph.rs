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
                "segment {} references out-of-range vertex ids ({}, {}) with vertex_count={}",
                segment, start, end, vertex_count,
            ),
            Self::DegenerateSegment { segment, vertex } => {
                write!(f, "segment {} is degenerate at vertex {}", segment, vertex)
            }
            Self::DuplicateSegment {
                first,
                second,
                a,
                b,
            } => write!(
                f,
                "segments {} and {} are duplicates of edge ({}, {})",
                first, second, a, b,
            ),
            Self::IntersectingSegments { first, second } => {
                write!(
                    f,
                    "segments {} and {} intersect in their interiors",
                    first, second
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
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            segments: Vec::new(),
            holes: Vec::new(),
        }
    }

    /// Create with pre-allocated capacity.
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
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Get vertex by ID.
    #[inline]
    pub fn vertex(&self, id: PslgVertexId) -> &PslgVertex {
        &self.vertices[id.idx()]
    }

    /// Slice of all vertex positions.
    #[inline]
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
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Get segment by ID.
    #[inline]
    pub fn segment(&self, id: PslgSegmentId) -> &PslgSegment {
        &self.segments[id.idx()]
    }

    /// Slice of all segments.
    #[inline]
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
    pub fn holes(&self) -> &[PslgVertex] {
        &self.holes
    }

    // ── Bounding box ──────────────────────────────────────────────────────

    /// Resolve any interior-crossing segment pairs by iteratively finding a
    /// crossing pair, inserting their 2-D intersection as a new vertex, and
    /// replacing each crossing segment with two sub-segments that share the
    /// new vertex.  Repeats until no interior crossings remain.
    ///
    /// After this call [`Self::validate`] will not return
    /// [`PslgValidationError::IntersectingSegments`] for proper interior
    /// crossings.  Collinear-overlap errors, if any, are left unmodified.
    ///
    /// # Complexity
    ///
    /// O(k · n²) where k is the number of crossing pairs in the original
    /// PSLG.  For the CSG corefine use-case k is typically 0–3.
    pub fn resolve_crossings(&mut self) {
        'outer: loop {
            let n_seg = self.segments.len();

            // Scan for proper interior-crossing pairs — compute intersection
            // vertex and split both crossing segments into halves.
            for i in 0..n_seg {
                for j in (i + 1)..n_seg {
                    let si = self.segments[i];
                    let sj = self.segments[j];

                    // Adjacent segments (shared endpoint) are always valid.
                    if si.start == sj.start
                        || si.start == sj.end
                        || si.end == sj.start
                        || si.end == sj.end
                    {
                        continue;
                    }

                    let a1 = self.vertices[si.start.idx()].to_point2();
                    let a2 = self.vertices[si.end.idx()].to_point2();
                    let b1 = self.vertices[sj.start.idx()].to_point2();
                    let b2 = self.vertices[sj.end.idx()].to_point2();

                    if !segments_cross_interior(&a1, &a2, &b1, &b2) {
                        continue;
                    }

                    let (px, py) = match segment_cross_point(&a1, &a2, &b1, &b2) {
                        Some(p) => p,
                        // orient_2d says crossing but formula says parallel —
                        // geometrically impossible; skip defensively.
                        None => continue,
                    };

                    let xid = self.add_vertex(px, py);

                    // Record endpoints before swap_removes invalidate the copies.
                    let (si_s, si_e) = (si.start, si.end);
                    let (sj_s, sj_e) = (sj.start, sj.end);

                    // Remove higher index first so lower index stays valid.
                    self.segments.swap_remove(j);
                    self.segments.swap_remove(i);

                    // Four sub-segments sharing the new intersection vertex.
                    self.add_segment(si_s, xid);
                    self.add_segment(xid, si_e);
                    self.add_segment(sj_s, xid);
                    self.add_segment(xid, sj_e);

                    continue 'outer; // restart with updated segment list
                }
            }

            break; // no interior crossings remain
        }
    }

    /// Compute the axis-aligned bounding box `(min, max)`.
    ///
    /// Returns `None` if the PSLG has fewer than 1 vertex.
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

/// Strict proper-interior crossing test — does **not** report endpoint
/// touches or collinear overlaps.  The caller must already have filtered out
/// segment pairs that share an endpoint.
fn segments_cross_interior(
    a1: &Point2<Real>,
    a2: &Point2<Real>,
    b1: &Point2<Real>,
    b2: &Point2<Real>,
) -> bool {
    let o1 = orient_2d(a1, a2, b1);
    let o2 = orient_2d(a1, a2, b2);
    let o3 = orient_2d(b1, b2, a1);
    let o4 = orient_2d(b1, b2, a2);

    // Any degenerate orientation means at least one endpoint is collinear
    // with the opposite segment — not a proper interior crossing.
    if o1 == Orientation::Degenerate
        || o2 == Orientation::Degenerate
        || o3 == Orientation::Degenerate
        || o4 == Orientation::Degenerate
    {
        return false;
    }

    o1 != o2 && o3 != o4
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
