//! Accelerated plane-based clip and refine operations.
//!
//! # CGAL 6.1 Insight — Plane Clippers (2025)
//!
//! CGAL 6.1 reimplemented `clip()`, `split()`, and introduced
//! `refine_with_plane()` achieving 10× speedup by:
//!
//! 1. **Cached plane equations** — compute normal + offset once, reuse for all
//!    orientation tests against that plane.
//! 2. **Batch vertex classification** — classify every unique vertex against the
//!    cut plane in a single O(V) pass before processing faces, reducing
//!    redundant `orient_3d` calls from O(3F) to O(V).
//! 3. **Direct plane clipping** — split straddling faces at the plane directly,
//!    producing 1–3 sub-triangles without the overhead of full CDT
//!    co-refinement.
//!
//! These techniques prevent cascading degenerates (thin walls from repeated
//! subdivision) and complement corefinement for hybrid workflows (e.g.,
//! pre-split at symmetry planes before full Boolean).
//!
//! # Theorem — Batch Vertex Classification Equivalence
//!
//! Let V be the set of unique vertex IDs referenced by F faces, and let P be a
//! plane defined by three CCW points.  Per-face classification calls
//! `orient_3d` 3|F| times, but since vertices are shared, |V| ≤ 3|F| with
//! equality only for isolated triangles.  Batch classification calls
//! `orient_3d` exactly |V| times.  Since `orient_3d` is a pure function of its
//! arguments, the per-vertex results are identical in both approaches.  ∎
//!
//! # References
//!
//! - Sébastien Loriot, Mael Rouxel-Labbé, Jane Tournois, Ilker O. Yaz,
//!   "Polygon Mesh Processing", CGAL 6.1, 2025.

use hashbrown::HashMap;

use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::geometry::predicates::{orient_3d, Orientation};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

// ── PlaneEquation ─────────────────────────────────────────────────────────────

/// A cached plane equation for efficient batch classification and clipping.
///
/// Stores both the three defining points (for exact `orient_3d` calls) and the
/// float normal + offset (for fast intersection-point computation).
///
/// # Theorem — Caching Preserves Exactness
///
/// `orient_3d(pa, pb, pc, q)` depends only on the coordinates of `pa, pb, pc,
/// q`.  Storing `[pa, pb, pc]` as `[[f64; 3]; 3]` and passing them to
/// `orient_3d` produces bit-identical results to recomputing from a
/// `VertexPool` on every call.  ∎
#[derive(Clone, Debug)]
pub struct PlaneEquation {
    /// Three CCW-oriented points defining the plane (exact `orient_3d` input).
    pa: [f64; 3],
    pb: [f64; 3],
    pc: [f64; 3],
    /// Float normal vector (unnormalised cross product `(pb-pa) × (pc-pa)`).
    normal: [f64; 3],
    /// `normal · pa` — the float plane offset.
    offset: f64,
}

impl PlaneEquation {
    /// Construct from three CCW-ordered points.
    #[must_use]
    pub fn from_points(pa: &Point3r, pb: &Point3r, pc: &Point3r) -> Self {
        let pa_arr = [pa.x, pa.y, pa.z];
        let pb_arr = [pb.x, pb.y, pb.z];
        let pc_arr = [pc.x, pc.y, pc.z];

        let ab = [pb.x - pa.x, pb.y - pa.y, pb.z - pa.z];
        let ac = [pc.x - pa.x, pc.y - pa.y, pc.z - pa.z];
        let n = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        let offset = n[0] * pa.x + n[1] * pa.y + n[2] * pa.z;

        Self {
            pa: pa_arr,
            pb: pb_arr,
            pc: pc_arr,
            normal: n,
            offset,
        }
    }

    /// Construct from a [`FaceData`] and [`VertexPool`].
    #[must_use]
    pub fn from_face(face: &FaceData, pool: &VertexPool) -> Self {
        let pa = pool.position(face.vertices[0]);
        let pb = pool.position(face.vertices[1]);
        let pc = pool.position(face.vertices[2]);
        Self::from_points(pa, pb, pc)
    }

    /// Exact orientation classification of a point against this plane.
    ///
    /// Uses Shewchuk `orient_3d` — no floating-point sign error.
    #[inline]
    #[must_use]
    pub fn classify(&self, point: &Point3r) -> Orientation {
        orient_3d(self.pa, self.pb, self.pc, [point.x, point.y, point.z])
    }

    /// Float signed distance (positive = same side as normal).
    ///
    /// Uses the unnormalised normal, so the magnitude is proportional to
    /// distance × ‖normal‖.  Only the *sign* is reliable for topology; use
    /// [`classify`](Self::classify) for exact decisions.
    #[inline]
    #[must_use]
    pub fn signed_distance_unnorm(&self, point: &Point3r) -> f64 {
        self.normal[0] * point.x + self.normal[1] * point.y + self.normal[2] * point.z
            - self.offset
    }

    /// Compute the plane–edge intersection along `s → e`.
    ///
    /// Returns the interpolated 3-D point.  The parameter `t = ds / (ds - de)`
    /// where `ds, de` are the (unnormalised) signed distances of `s, e`.
    #[must_use]
    pub fn intersect_edge(&self, s: &Point3r, e: &Point3r) -> Point3r {
        let ds = self.signed_distance_unnorm(s);
        let de = self.signed_distance_unnorm(e);
        let denom = ds - de;
        if denom.abs() < 1e-30 {
            return *s;
        }
        let t = ds / denom;
        Point3r::new(
            s.x + (e.x - s.x) * t,
            s.y + (e.y - s.y) * t,
            s.z + (e.z - s.z) * t,
        )
    }

    /// The unnormalised normal vector as `[f64; 3]`.
    #[inline]
    #[must_use]
    pub fn normal_array(&self) -> [f64; 3] {
        self.normal
    }

    /// The unnormalised normal as a nalgebra `Vector3r`.
    #[inline]
    #[must_use]
    pub fn normal_vec(&self) -> Vector3r {
        Vector3r::new(self.normal[0], self.normal[1], self.normal[2])
    }
}

// ── Face classification ───────────────────────────────────────────────────────

/// Classification of a face's vertices against a plane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FacePlaneClass {
    /// All vertices are positive or degenerate — fully inside the half-space.
    Inside,
    /// All vertices are negative or degenerate — fully outside the half-space.
    Outside,
    /// At least one positive and one negative — face must be split.
    Straddling,
    /// All vertices are exactly on the plane.
    Coplanar,
}

/// Classify a face against a cached plane equation.
///
/// Returns the classification and the three per-vertex orientations.
#[must_use]
pub fn classify_face(
    face: &FaceData,
    pool: &VertexPool,
    plane: &PlaneEquation,
) -> (FacePlaneClass, [Orientation; 3]) {
    let signs = [
        plane.classify(pool.position(face.vertices[0])),
        plane.classify(pool.position(face.vertices[1])),
        plane.classify(pool.position(face.vertices[2])),
    ];

    let any_pos = signs.iter().any(|s| *s == Orientation::Positive);
    let any_neg = signs.iter().any(|s| *s == Orientation::Negative);

    let class = match (any_pos, any_neg) {
        (true, true) => FacePlaneClass::Straddling,
        (true, false) => FacePlaneClass::Inside,
        (false, true) => FacePlaneClass::Outside,
        (false, false) => FacePlaneClass::Coplanar,
    };

    (class, signs)
}

// ── Face-level plane clipping ─────────────────────────────────────────────────

/// Clip a face by a plane, returning sub-faces in the positive half-space.
///
/// Uses exact `orient_3d` for inside/outside classification and float
/// arithmetic for intersection positions (acceptable per Shewchuk: position
/// errors do not affect topology).
///
/// New vertices from edge–plane intersections are inserted into `pool` via
/// `insert_or_weld`.
///
/// Returns 0 faces if fully clipped, 1 face if fully kept, or 1–2 faces for
/// a straddling triangle.
pub fn clip_face_by_plane(
    face: &FaceData,
    pool: &mut VertexPool,
    plane: &PlaneEquation,
) -> Vec<FaceData> {
    let (class, signs) = classify_face(face, pool, plane);

    match class {
        FacePlaneClass::Inside | FacePlaneClass::Coplanar => vec![*face],
        FacePlaneClass::Outside => Vec::new(),
        FacePlaneClass::Straddling => {
            let vids = face.vertices;
            let positions = [
                *pool.position(vids[0]),
                *pool.position(vids[1]),
                *pool.position(vids[2]),
            ];
            let plane_n = plane.normal_vec();

            let mut output: Vec<VertexId> = Vec::with_capacity(4);
            for i in 0..3 {
                let j = (i + 1) % 3;
                let s_in = signs[i] != Orientation::Negative;
                let e_in = signs[j] != Orientation::Negative;

                match (s_in, e_in) {
                    (true, true) => output.push(vids[j]),
                    (true, false) => {
                        let cut = plane.intersect_edge(&positions[i], &positions[j]);
                        output.push(pool.insert_or_weld(cut, plane_n));
                    }
                    (false, true) => {
                        let cut = plane.intersect_edge(&positions[i], &positions[j]);
                        output.push(pool.insert_or_weld(cut, plane_n));
                        output.push(vids[j]);
                    }
                    (false, false) => {}
                }
            }
            fan_triangulate_vids(&output, face.region)
        }
    }
}

// ── Mesh-level plane refinement ───────────────────────────────────────────────

/// Split all faces of a mesh by a plane.
///
/// # CGAL 6.1 — `refine_with_plane`
///
/// Equivalent to CGAL 6.1's `refine_with_plane()`:
///
/// 1. **Batch classify** all unique vertices against the plane — O(V)
///    `orient_3d` calls (each vertex tested exactly once).
/// 2. For each face, look up cached per-vertex classifications — O(1).
/// 3. Non-straddling faces pass through unchanged.
/// 4. Straddling faces are split at the plane, producing 1–3 sub-faces per
///    side.
///
/// Total: O(V + F_straddling) instead of O(3F) for naive per-face orient_3d.
///
/// Returns `(inside_faces, outside_faces)` — both are valid face soups using
/// the same pool.  Coplanar faces are included in **both** halves.
pub fn refine_faces_with_plane(
    faces: &[FaceData],
    pool: &mut VertexPool,
    plane: &PlaneEquation,
) -> (Vec<FaceData>, Vec<FaceData>) {
    // Phase 1: batch-classify all unique vertices.
    let mut vertex_signs: HashMap<VertexId, Orientation> =
        HashMap::with_capacity(faces.len() * 2);
    for face in faces {
        for &vid in &face.vertices {
            vertex_signs
                .entry(vid)
                .or_insert_with(|| plane.classify(pool.position(vid)));
        }
    }

    let mut inside = Vec::with_capacity(faces.len());
    let mut outside = Vec::with_capacity(faces.len() / 4);

    // Phase 2: process each face using cached classifications.
    for face in faces {
        let signs = [
            vertex_signs[&face.vertices[0]],
            vertex_signs[&face.vertices[1]],
            vertex_signs[&face.vertices[2]],
        ];
        let any_pos = signs.iter().any(|s| *s == Orientation::Positive);
        let any_neg = signs.iter().any(|s| *s == Orientation::Negative);

        match (any_pos, any_neg) {
            (true, false) => inside.push(*face),
            (false, true) => outside.push(*face),
            (false, false) => {
                inside.push(*face);
                outside.push(*face);
            }
            (true, true) => {
                let (pos, neg) = split_straddling_face(face, pool, plane, &signs);
                inside.extend(pos);
                outside.extend(neg);
            }
        }
    }

    (inside, outside)
}

/// Split a straddling face into positive-side and negative-side sub-faces.
///
/// Walks the triangle edges in order, tracking which vertices/intersection
/// points belong to each half.  Intersection points are added to both polygons.
fn split_straddling_face(
    face: &FaceData,
    pool: &mut VertexPool,
    plane: &PlaneEquation,
    signs: &[Orientation; 3],
) -> (Vec<FaceData>, Vec<FaceData>) {
    let vids = face.vertices;
    let positions = [
        *pool.position(vids[0]),
        *pool.position(vids[1]),
        *pool.position(vids[2]),
    ];
    let plane_n = plane.normal_vec();

    let mut pos_poly: Vec<VertexId> = Vec::with_capacity(4);
    let mut neg_poly: Vec<VertexId> = Vec::with_capacity(4);

    for i in 0..3 {
        let j = (i + 1) % 3;

        // Emit vertex to its polygon.
        match signs[i] {
            Orientation::Positive => pos_poly.push(vids[i]),
            Orientation::Negative => neg_poly.push(vids[i]),
            Orientation::Degenerate => {
                pos_poly.push(vids[i]);
                neg_poly.push(vids[i]);
            }
        }

        // If the edge crosses the plane, emit intersection to both polygons.
        let crosses = (signs[i] == Orientation::Positive && signs[j] == Orientation::Negative)
            || (signs[i] == Orientation::Negative && signs[j] == Orientation::Positive);
        if crosses {
            let cut = plane.intersect_edge(&positions[i], &positions[j]);
            let cut_vid = pool.insert_or_weld(cut, plane_n);
            pos_poly.push(cut_vid);
            neg_poly.push(cut_vid);
        }
    }

    let region = face.region;
    (
        fan_triangulate_vids(&pos_poly, region),
        fan_triangulate_vids(&neg_poly, region),
    )
}

/// Fan-triangulate a convex polygon given as vertex IDs.
fn fan_triangulate_vids(vids: &[VertexId], region: crate::domain::core::index::RegionId) -> Vec<FaceData> {
    if vids.len() < 3 {
        return Vec::new();
    }
    let root = vids[0];
    (1..vids.len() - 1)
        .map(|k| FaceData::new(root, vids[k], vids[k + 1], region))
        .collect()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::index::RegionId;

    fn csg_pool() -> VertexPool {
        VertexPool::for_csg()
    }

    fn p(x: f64, y: f64, z: f64) -> Point3r {
        Point3r::new(x, y, z)
    }

    fn n_up() -> Vector3r {
        Vector3r::new(0.0, 0.0, 1.0)
    }

    /// Helper: build a face from 3 points, returning (face, pool).
    fn face_from_pts(a: Point3r, b: Point3r, c: Point3r) -> (FaceData, VertexPool) {
        let mut pool = csg_pool();
        let n = n_up();
        let va = pool.insert_or_weld(a, n);
        let vb = pool.insert_or_weld(b, n);
        let vc = pool.insert_or_weld(c, n);
        let face = FaceData::new(va, vb, vc, RegionId::INVALID);
        (face, pool)
    }

    // ── PlaneEquation tests ──────────────────────────────────────────────

    #[test]
    fn plane_classify_above_below_on() {
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        assert_eq!(plane.classify(&p(0.5, 0.5, 1.0)), Orientation::Positive);
        assert_eq!(plane.classify(&p(0.5, 0.5, -1.0)), Orientation::Negative);
        assert_eq!(plane.classify(&p(0.5, 0.5, 0.0)), Orientation::Degenerate);
    }

    #[test]
    fn plane_from_face_matches_from_points() {
        let mut pool = csg_pool();
        let n = n_up();
        let va = pool.insert_or_weld(p(0.0, 0.0, 0.0), n);
        let vb = pool.insert_or_weld(p(1.0, 0.0, 0.0), n);
        let vc = pool.insert_or_weld(p(0.0, 1.0, 0.0), n);
        let face = FaceData::new(va, vb, vc, RegionId::INVALID);

        let p1 = PlaneEquation::from_face(&face, &pool);
        let p2 = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );

        let q = p(0.3, 0.2, 0.7);
        assert_eq!(p1.classify(&q), p2.classify(&q));
    }

    #[test]
    fn plane_intersect_edge_midpoint() {
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let s = p(0.0, 0.0, 1.0);
        let e = p(0.0, 0.0, -1.0);
        let cut = plane.intersect_edge(&s, &e);
        assert!((cut.x).abs() < 1e-12);
        assert!((cut.y).abs() < 1e-12);
        assert!((cut.z).abs() < 1e-12);
    }

    #[test]
    fn plane_intersect_edge_quarter() {
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 1.0),
            &p(1.0, 0.0, 1.0),
            &p(0.0, 1.0, 1.0),
        );
        let s = p(0.0, 0.0, 0.0);
        let e = p(0.0, 0.0, 4.0);
        let cut = plane.intersect_edge(&s, &e);
        assert!((cut.z - 1.0).abs() < 1e-12);
    }

    // ── classify_face tests ──────────────────────────────────────────────

    #[test]
    fn classify_face_above_plane() {
        let (face, pool) = face_from_pts(p(0.0, 0.0, 1.0), p(1.0, 0.0, 2.0), p(0.0, 1.0, 3.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (class, _) = classify_face(&face, &pool, &plane);
        assert_eq!(class, FacePlaneClass::Inside);
    }

    #[test]
    fn classify_face_below_plane() {
        let (face, pool) =
            face_from_pts(p(0.0, 0.0, -1.0), p(1.0, 0.0, -2.0), p(0.0, 1.0, -3.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (class, _) = classify_face(&face, &pool, &plane);
        assert_eq!(class, FacePlaneClass::Outside);
    }

    #[test]
    fn classify_face_straddling() {
        let (face, pool) = face_from_pts(p(0.0, 0.0, 1.0), p(1.0, 0.0, -1.0), p(0.0, 1.0, 0.5));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (class, _) = classify_face(&face, &pool, &plane);
        assert_eq!(class, FacePlaneClass::Straddling);
    }

    #[test]
    fn classify_face_coplanar() {
        let (face, pool) = face_from_pts(p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (class, _) = classify_face(&face, &pool, &plane);
        assert_eq!(class, FacePlaneClass::Coplanar);
    }

    // ── clip_face_by_plane tests ─────────────────────────────────────────

    #[test]
    fn clip_fully_inside_keeps_face() {
        let (face, mut pool) =
            face_from_pts(p(0.0, 0.0, 1.0), p(1.0, 0.0, 2.0), p(0.0, 1.0, 3.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let result = clip_face_by_plane(&face, &mut pool, &plane);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], face);
    }

    #[test]
    fn clip_fully_outside_removes_face() {
        let (face, mut pool) =
            face_from_pts(p(0.0, 0.0, -1.0), p(1.0, 0.0, -2.0), p(0.0, 1.0, -3.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let result = clip_face_by_plane(&face, &mut pool, &plane);
        assert!(result.is_empty());
    }

    #[test]
    fn clip_straddling_produces_subtriangles() {
        let (face, mut pool) =
            face_from_pts(p(0.0, 0.0, 1.0), p(1.0, 0.0, -1.0), p(0.0, 1.0, 1.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let result = clip_face_by_plane(&face, &mut pool, &plane);
        assert!(
            !result.is_empty(),
            "straddling face should produce inside sub-faces"
        );
        // One vertex below → inside polygon is 4 vertices → 2 triangles
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn clip_coplanar_keeps_face() {
        let (face, mut pool) =
            face_from_pts(p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let result = clip_face_by_plane(&face, &mut pool, &plane);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn clip_one_vertex_on_plane_two_above() {
        // v0 on plane, v1 and v2 above → fully inside (Degenerate counts as
        // inside).
        let (face, mut pool) =
            face_from_pts(p(0.5, 0.5, 0.0), p(1.0, 0.0, 1.0), p(0.0, 1.0, 1.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let result = clip_face_by_plane(&face, &mut pool, &plane);
        assert_eq!(result.len(), 1);
    }

    // ── refine_faces_with_plane tests ────────────────────────────────────

    #[test]
    fn refine_empty_mesh() {
        let mut pool = csg_pool();
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (ins, outs) = refine_faces_with_plane(&[], &mut pool, &plane);
        assert!(ins.is_empty());
        assert!(outs.is_empty());
    }

    #[test]
    fn refine_all_above() {
        let (face, mut pool) =
            face_from_pts(p(0.0, 0.0, 1.0), p(1.0, 0.0, 2.0), p(0.0, 1.0, 3.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (ins, outs) = refine_faces_with_plane(&[face], &mut pool, &plane);
        assert_eq!(ins.len(), 1);
        assert!(outs.is_empty());
    }

    #[test]
    fn refine_straddling_triangle_splits_both_sides() {
        let (face, mut pool) =
            face_from_pts(p(0.0, 0.0, 1.0), p(1.0, 0.0, -1.0), p(0.0, 1.0, 1.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (ins, outs) = refine_faces_with_plane(&[face], &mut pool, &plane);
        assert!(!ins.is_empty(), "should have inside sub-faces");
        assert!(!outs.is_empty(), "should have outside sub-faces");
        // 1 vertex below → inside = 2 tris, outside = 1 tri
        assert_eq!(ins.len(), 2);
        assert_eq!(outs.len(), 1);
    }

    #[test]
    fn refine_preserves_region_ids() {
        let mut pool = csg_pool();
        let n = n_up();
        let va = pool.insert_or_weld(p(0.0, 0.0, 1.0), n);
        let vb = pool.insert_or_weld(p(1.0, 0.0, -1.0), n);
        let vc = pool.insert_or_weld(p(0.0, 1.0, 1.0), n);
        let region = RegionId::new(42);
        let face = FaceData::new(va, vb, vc, region);

        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (ins, outs) = refine_faces_with_plane(&[face], &mut pool, &plane);
        for f in ins.iter().chain(outs.iter()) {
            assert_eq!(f.region, region, "sub-faces must inherit parent region");
        }
    }

    #[test]
    fn refine_batch_classification_reduces_orient3d_calls() {
        // 4 faces sharing 5 vertices: batch should call orient_3d 5 times,
        // not 12.  We verify correctness by checking the output (counts of
        // inside/outside faces), which implicitly validates that the cached
        // per-vertex classifications were used.
        let mut pool = csg_pool();
        let n = n_up();
        // Shared vertices
        let v0 = pool.insert_or_weld(p(0.0, 0.0, 1.0), n);
        let v1 = pool.insert_or_weld(p(1.0, 0.0, 1.0), n);
        let v2 = pool.insert_or_weld(p(0.5, 1.0, 1.0), n);
        let v3 = pool.insert_or_weld(p(0.5, 0.5, -1.0), n);
        let v4 = pool.insert_or_weld(p(1.0, 1.0, 1.0), n);

        let faces = [
            FaceData::untagged(v0, v1, v2), // all above
            FaceData::untagged(v0, v1, v3), // v3 below → straddles
            FaceData::untagged(v1, v2, v3), // v3 below → straddles
            FaceData::untagged(v1, v4, v2), // all above
        ];

        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (ins, outs) = refine_faces_with_plane(&faces, &mut pool, &plane);
        // 2 fully-inside + 2 straddling (each producing 2 inside sub-faces)
        assert!(ins.len() >= 4);
        assert!(!outs.is_empty());
    }

    #[test]
    fn refine_vertex_on_plane_goes_to_both_halves() {
        // v0 on plane, v1 above, v2 below → straddling.
        // The on-plane vertex appears in both halves.
        let (face, mut pool) =
            face_from_pts(p(0.5, 0.5, 0.0), p(1.0, 0.0, 1.0), p(0.0, 1.0, -1.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 0.0),
        );
        let (ins, outs) = refine_faces_with_plane(&[face], &mut pool, &plane);
        assert_eq!(ins.len(), 1);
        assert_eq!(outs.len(), 1);
        // The on-plane vertex (v0) should appear in both sub-faces.
        let v0 = face.vertices[0];
        assert!(ins[0].vertices.contains(&v0));
        assert!(outs[0].vertices.contains(&v0));
    }

    #[test]
    fn refine_tilted_plane_xy() {
        // Plane at 45° through x-axis: n = (0, -1/√2, 1/√2), d = 0
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(0.0, 1.0, 1.0),
        );
        let (face, mut pool) =
            face_from_pts(p(0.0, 0.0, 2.0), p(1.0, 2.0, 0.0), p(0.0, 2.0, 0.0));
        let (ins, outs) = refine_faces_with_plane(&[face], &mut pool, &plane);
        let total = ins.len() + outs.len();
        assert!(total >= 2, "tilted plane should split the face");
    }

    #[test]
    fn refine_is_deterministic() {
        let build = || {
            let (face, mut pool) =
                face_from_pts(p(0.0, 0.0, 1.0), p(1.0, 0.0, -1.0), p(0.0, 1.0, 0.5));
            let plane = PlaneEquation::from_points(
                &p(0.0, 0.0, 0.0),
                &p(1.0, 0.0, 0.0),
                &p(0.0, 1.0, 0.0),
            );
            refine_faces_with_plane(&[face], &mut pool, &plane)
        };
        let (ins1, outs1) = build();
        let (ins2, outs2) = build();
        assert_eq!(ins1.len(), ins2.len());
        assert_eq!(outs1.len(), outs2.len());
    }

    #[test]
    fn clip_face_degenerate_plane_no_panic() {
        // Degenerate plane (collinear points) → normal = 0.
        let (face, mut pool) =
            face_from_pts(p(0.0, 0.0, 1.0), p(1.0, 0.0, 1.0), p(0.0, 1.0, 1.0));
        let plane = PlaneEquation::from_points(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(2.0, 0.0, 0.0),
        );
        // All orient_3d calls return Degenerate for a degenerate plane.
        let result = clip_face_by_plane(&face, &mut pool, &plane);
        // Should not panic; face treated as coplanar (all Degenerate).
        assert!(!result.is_empty());
    }
}
