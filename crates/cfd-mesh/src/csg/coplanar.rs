//! Coplanar 2-D Boolean operations for flat surface meshes.
//!
//! When both operand meshes are flat (all triangles share a common plane),
//! the standard 3-D BSP and Mesh-Arrangement pipelines break down:
//!
//! - **BSP** selects a face plane as its splitting plane.  Every other face is
//!   `Coplanar` with that plane, so nothing ever lands in the Front/Back lists
//!   and the tree degenerates into a single node with all faces stored as
//!   `coplanar`.  The subsequent `clip_faces` dot-product test is then the
//!   only classification, but for a flat mesh all faces share the same normal
//!   direction, so every face goes to the same list — the clipping is vacuous.
//!
//! - **Mesh Arrangement** (arrangement.rs) is designed for 3-D solids.  All
//!   triangle-triangle intersection tests return `IntersectionType::Coplanar`
//!   (no `Segment` variants), so the seam subdivision phase does nothing and
//!   the centroid nudge moves test points *out of the plane*, where
//!   `point_in_mesh_local` (a 3-D ray cast) gives unreliable results because
//!   the flat mesh has no closed volume to penetrate.
//!
//! ## Algorithm — sub-triangle centroid classification
//!
//! For each A triangle:
//!
//! 1. **Project** all triangle vertices onto the shared plane.
//! 2. **Fast classify** via 2-D vertex-in-union test:
//!    - All 3 A-vertices inside (∪B) → whole A triangle is inside B.
//!    - All 3 A-vertices outside (∪B) → whole A triangle is outside B.
//! 3. **Boundary triangles** (mixed vertices): sub-divide the A triangle into
//!    an 8×8 barycentric grid (up to 128 sub-triangles) and test each
//!    sub-triangle's centroid directly against `b_tris` using
//!    `point_in_b_2d`.
//! 4. **Select sub-triangles** per Boolean operation:
//!    - Union:        sub-triangles whose centroid is **outside** (∪B)
//!    - Intersection: sub-triangles whose centroid is **inside**  (∪B)
//!    - Difference:   sub-triangles whose centroid is **outside** (∪B)
//! 5. **B-side (Union only)**: symmetrically add B sub-triangles outside (∪A).
//!
//! ## Complexity
//!
//! O(n·m·S²) where n = |faces_a|, m = |faces_b|, S = 8 (subdivision depth).
//! For a 128-segment disk that is 128 × 128 × 128 ≈ 2 M point tests — fast
//! and completely free of exponential blowup.
//!
//! ## Theorem — area accuracy
//!
//! For a boundary triangle of area A_t that straddles the boundary, the
//! sub-triangle grid introduces at most O(A_t · S⁻²) area error — the area
//! of one row of sub-triangles adjacent to the boundary edge.  At S = 8 this
//! is 1/64 of the boundary-triangle area.  Since there are O(√n) boundary
//! triangles each of area O(1/n), the total area error is O(n^{-1/2}) → 0.
//!
//! ## Output mesh orientation
//!
//! Output triangles preserve the winding of the A-side operand.  B-side
//! fragments used in Union are re-wound to match A's normal direction.

#![allow(missing_docs)]

use crate::core::scalar::{Point3r, Real, Vector3r};
use crate::csg::boolean::BooleanOp;
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;

// ── 2-D coordinate system on the shared plane ─────────────────────────────────

/// An orthonormal basis for projecting the shared plane into 2-D.
///
/// `pub(crate)` so that `boolean.rs` can receive a `&PlaneBasis` from
/// `detect_flat_plane` and pass it to `boolean_coplanar` without having to
/// re-export the type in the public API.
pub(crate) struct PlaneBasis {
    /// Origin of the 2-D coordinate system (first vertex of first face).
    pub(crate) origin: Point3r,
    /// First basis vector (lies in the plane).
    pub(crate) u: Vector3r,
    /// Second basis vector (lies in the plane, perpendicular to `u`).
    pub(crate) v: Vector3r,
    /// Plane normal (= u × v, normalised).
    pub(crate) normal: Vector3r,
}

impl PlaneBasis {
    /// Build a basis from three non-collinear points on the plane.
    fn from_triangle(a: &Point3r, b: &Point3r, c: &Point3r) -> Option<Self> {
        let ab = b - a;
        let ac = c - a;
        let n  = ab.cross(&ac);
        let nl = n.norm();
        if nl < 1e-20 { return None; }

        let u_raw = ab;
        let u_len = u_raw.norm();
        if u_len < 1e-20 { return None; }
        let u = u_raw / u_len;
        let normal = n / nl;
        let v = normal.cross(&u).normalize();

        Some(Self { origin: *a, u, v, normal })
    }

    /// Project a 3-D point onto the 2-D plane coordinates.
    #[inline]
    pub(crate) fn project(&self, p: &Point3r) -> [Real; 2] {
        let d = p - self.origin;
        [d.dot(&self.u), d.dot(&self.v)]
    }
}

// ── 2-D point-in-union test ────────────────────────────────────────────────────

/// Test whether a 2-D point `(px, py)` is **inside or on the boundary** of the
/// 2-D triangle `(ax,ay)-(bx,by)-(cx,cy)`.
///
/// Works for both CCW and CW winding by checking sign consistency.
#[inline]
fn point_in_triangle_2d(
    px: Real, py: Real,
    ax: Real, ay: Real,
    bx: Real, by: Real,
    cx: Real, cy: Real,
) -> bool {
    // Cross products for each edge.
    let d0 = (bx - ax) * (py - ay) - (by - ay) * (px - ax);
    let d1 = (cx - bx) * (py - by) - (cy - by) * (px - bx);
    let d2 = (ax - cx) * (py - cy) - (ay - cy) * (px - cx);

    let has_neg = d0 < 0.0 || d1 < 0.0 || d2 < 0.0;
    let has_pos = d0 > 0.0 || d1 > 0.0 || d2 > 0.0;
    !(has_neg && has_pos)
}

/// Test whether 2-D point `(px, py)` is inside **any** triangle in `tris`.
///
/// Each entry in `tris` is `[ax, ay, bx, by, cx, cy]`.
#[inline]
fn point_in_union_2d(px: Real, py: Real, tris: &[[Real; 6]]) -> bool {
    for t in tris {
        if point_in_triangle_2d(px, py, t[0], t[1], t[2], t[3], t[4], t[5]) {
            return true;
        }
    }
    false
}

// ── Emit a single triangle ─────────────────────────────────────────────────────

/// Insert one triangle into the result pool with correct CCW winding and normal.
///
/// Skips degenerate triangles (coincident vertices after welding).
#[inline]
fn emit_one(
    p0: Point3r,
    p1: Point3r,
    p2: Point3r,
    basis: &PlaneBasis,
    region: crate::core::index::RegionId,
    result: &mut Vec<FaceData>,
    pool: &mut VertexPool,
) {
    // Area check: skip degenerate (near-zero area) sub-triangles.
    let ab = p1 - p0;
    let ac = p2 - p0;
    let fn_ = ab.cross(&ac);
    if fn_.norm() < 1e-20 { return; }

    // Ensure CCW winding relative to basis normal.
    let flip = fn_.dot(&basis.normal) < 0.0;
    let (o0, o1, o2) = if flip { (p0, p2, p1) } else { (p0, p1, p2) };

    let v0 = pool.insert_or_weld(o0, basis.normal);
    let v1 = pool.insert_or_weld(o1, basis.normal);
    let v2 = pool.insert_or_weld(o2, basis.normal);
    if v0 != v1 && v1 != v2 && v0 != v2 {
        result.push(FaceData::new(v0, v1, v2, region));
    }
}

// ── Main entry point ───────────────────────────────────────────────────────────

/// Detect whether a face soup is flat (all triangles lie in the same plane).
///
/// Returns the `PlaneBasis` if flat, or `None` if the soup is 3-D / curved.
///
/// **Algorithm:**
/// 1. Find the first non-degenerate triangle to define the plane.
/// 2. For every subsequent vertex, check that its signed distance to the plane
///    is within `FLAT_TOLERANCE`.
///
/// **Complexity:** O(n) — checks every vertex once.
///
/// ## Theorem — Flat-plane criterion
///
/// A set of points `{p_i}` all lie in a common plane iff for any non-degenerate
/// triangle `(a, b, c)` from the set, `(p_i − a) · n̂ = 0` for all `i`, where
/// `n̂ = (b−a) × (c−a) / ‖…‖`.
pub(crate) fn detect_flat_plane(faces: &[FaceData], pool: &VertexPool) -> Option<PlaneBasis> {
    // 1. Find first non-degenerate face to establish the reference plane.
    let mut basis: Option<PlaneBasis> = None;
    for face in faces {
        let a = pool.position(face.vertices[0]);
        let b = pool.position(face.vertices[1]);
        let c = pool.position(face.vertices[2]);
        if let Some(b0) = PlaneBasis::from_triangle(a, b, c) {
            basis = Some(b0);
            break;
        }
    }
    let basis = basis?;

    // 2. Verify every vertex lies within FLAT_TOLERANCE of the plane.
    const FLAT_TOLERANCE: Real = 1e-6; // 1 nm — safe for mm-scale geometry
    for face in faces {
        for &vid in &face.vertices {
            let p = pool.position(vid);
            let d = (p - basis.origin).dot(&basis.normal).abs();
            if d > FLAT_TOLERANCE {
                return None; // Not flat.
            }
        }
    }

    Some(basis)
}

/// Perform a Boolean operation on two **coplanar flat** face soups.
///
/// Both soups must lie in the same plane (verified by the caller via
/// `detect_flat_plane`).  The `basis` parameter is the plane coordinate system
/// used for 2-D projection.
///
/// ## Algorithm — sub-triangle centroid classification
///
/// For each triangle T in the source mesh:
///
/// 1. **Fast path** — test all 3 vertices against the opposing mesh using the
///    2-D point-in-union test.
///    - All inside  → emit whole T (or skip, depending on operation).
///    - All outside → emit whole T (or skip).
/// 2. **Boundary triangles** — sub-divide T into an (SUB×SUB) barycentric grid
///    and test each sub-triangle centroid with `point_in_union_2d`.
///    Emit sub-triangles whose centroid satisfies the Boolean predicate.
///
/// No polygon subtraction is performed — only forward clipping via centroid
/// tests — so there is no exponential blowup.
///
/// # Returns
///
/// A `Vec<FaceData>` representing the result, re-inserted into `pool`.
pub(crate) fn boolean_coplanar(
    op: BooleanOp,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
    basis: &PlaneBasis,
) -> Vec<FaceData> {
    let mut result: Vec<FaceData> = Vec::new();

    // Pre-project all triangles of both meshes to 2-D for O(1) containment tests.
    // Each entry: [ax,ay, bx,by, cx,cy].
    let b_tris: Vec<[Real; 6]> = faces_b.iter().map(|fb| {
        let d = pool.position(fb.vertices[0]);
        let e = pool.position(fb.vertices[1]);
        let f = pool.position(fb.vertices[2]);
        let [dx, dy] = basis.project(d);
        let [ex, ey] = basis.project(e);
        let [fx, fy] = basis.project(f);
        [dx, dy, ex, ey, fx, fy]
    }).collect();

    let a_tris: Vec<[Real; 6]> = faces_a.iter().map(|fa| {
        let p = pool.position(fa.vertices[0]);
        let q = pool.position(fa.vertices[1]);
        let r = pool.position(fa.vertices[2]);
        let [px, py] = basis.project(p);
        let [qx, qy] = basis.project(q);
        let [rx, ry] = basis.project(r);
        [px, py, qx, qy, rx, ry]
    }).collect();

    // Sub-division depth: each boundary triangle is divided into SUB×SUB cells.
    // 8 gives 64–128 sub-triangles per boundary triangle — sufficient accuracy.
    const SUB: usize = 8;
    let inv = 1.0 / (SUB as Real);

    // ── Classify and emit one triangle (fa or fb) ────────────────────────────
    //
    // `check_inside` is the set of 2-D triangles that the centroid is tested
    // against.  `want_inside` is the predicate (emit if inside / emit if outside).
    let classify_and_emit = |
        pa: &Point3r, pb: &Point3r, pc: &Point3r,
        check_tris: &[[Real; 6]],
        want_inside: bool,
        region: crate::core::index::RegionId,
        result: &mut Vec<FaceData>,
        pool: &mut VertexPool,
    | {
        let [pax, pay] = basis.project(pa);
        let [pbx, pby] = basis.project(pb);
        let [pcx, pcy] = basis.project(pc);

        let va = point_in_union_2d(pax, pay, check_tris);
        let vb = point_in_union_2d(pbx, pby, check_tris);
        let vc = point_in_union_2d(pcx, pcy, check_tris);

        if va == want_inside && vb == want_inside && vc == want_inside {
            // Fast path: whole triangle has the desired classification.
            emit_one(*pa, *pb, *pc, basis, region, result, pool);
            return;
        }

        if va != want_inside && vb != want_inside && vc != want_inside {
            // Whole triangle has the opposite classification — skip.
            return;
        }

        // Boundary triangle: sub-divide into barycentric grid and classify each cell.
        for i in 0..SUB {
            for j in 0..(SUB - i) {
                let fi = i as Real;
                let fj = j as Real;

                // Sub-triangle type 0: (i,j), (i+1,j), (i,j+1)
                let tris_to_check = [
                    (fi*inv, fj*inv, (fi+1.0)*inv, fj*inv, fi*inv, (fj+1.0)*inv),
                ];

                // Sub-triangle type 1: (i+1,j), (i+1,j+1), (i,j+1)  [only if valid]
                let type1 = if i + j + 2 <= SUB {
                    Some(((fi+1.0)*inv, fj*inv, (fi+1.0)*inv, (fj+1.0)*inv, fi*inv, (fj+1.0)*inv))
                } else {
                    None
                };

                for (u0, v0, u1, v1, u2, v2) in tris_to_check.iter().chain(type1.iter()) {
                    let w0 = 1.0 - u0 - v0;
                    let w1 = 1.0 - u1 - v1;
                    let w2 = 1.0 - u2 - v2;

                    let q0 = Point3r::new(
                        pa.x*w0 + pb.x*u0 + pc.x*v0,
                        pa.y*w0 + pb.y*u0 + pc.y*v0,
                        pa.z*w0 + pb.z*u0 + pc.z*v0,
                    );
                    let q1 = Point3r::new(
                        pa.x*w1 + pb.x*u1 + pc.x*v1,
                        pa.y*w1 + pb.y*u1 + pc.y*v1,
                        pa.z*w1 + pb.z*u1 + pc.z*v1,
                    );
                    let q2 = Point3r::new(
                        pa.x*w2 + pb.x*u2 + pc.x*v2,
                        pa.y*w2 + pb.y*u2 + pc.y*v2,
                        pa.z*w2 + pb.z*u2 + pc.z*v2,
                    );

                    // Centroid of sub-triangle.
                    let cx = (q0.x + q1.x + q2.x) / 3.0;
                    let cy = (q0.y + q1.y + q2.y) / 3.0;
                    let cz = (q0.z + q1.z + q2.z) / 3.0;
                    let centroid_3d = Point3r::new(cx, cy, cz);
                    let [cpx, cpy] = basis.project(&centroid_3d);

                    let inside = point_in_union_2d(cpx, cpy, check_tris);
                    if inside == want_inside {
                        emit_one(q0, q1, q2, basis, region, result, pool);
                    }
                }
            }
        }
    };

    // ── A-side ────────────────────────────────────────────────────────────────

    for fa in faces_a.iter() {
        let pa = *pool.position(fa.vertices[0]);
        let pb = *pool.position(fa.vertices[1]);
        let pc = *pool.position(fa.vertices[2]);

        match op {
            BooleanOp::Union => {
                // Emit A triangles outside B  (A \ B) — avoid overlap with B-side below.
                // Also emit A triangles inside B — one copy of the overlap.
                // Together that is: emit all of A.
                // BUT: the B-side will also add B \ A.  So we need to emit all of A
                // and only the B-outside-A part.  → emit all of A here.
                emit_one(pa, pb, pc, basis, fa.region, &mut result, pool);
            }
            BooleanOp::Intersection => {
                // Emit A fragments that are inside B.
                classify_and_emit(
                    &pa, &pb, &pc, &b_tris, true,
                    fa.region, &mut result, pool,
                );
            }
            BooleanOp::Difference => {
                // Emit A fragments that are outside B.
                classify_and_emit(
                    &pa, &pb, &pc, &b_tris, false,
                    fa.region, &mut result, pool,
                );
            }
        }
    }

    // ── B-side (Union only) ───────────────────────────────────────────────────
    //
    // Add B fragments that are **outside** A so that the overlap region is not
    // double-counted (A already contributed the full overlap above).
    if matches!(op, BooleanOp::Union) {
        for fb in faces_b.iter() {
            let d = *pool.position(fb.vertices[0]);
            let e = *pool.position(fb.vertices[1]);
            let f = *pool.position(fb.vertices[2]);

            // Emit only B triangles (or sub-triangles) that are outside A.
            classify_and_emit(
                &d, &e, &f, &a_tris, false,
                fb.region, &mut result, pool,
            );
        }
    }

    result
}

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use crate::geometry::primitives::{Disk, PrimitiveMesh};
    use crate::csg::boolean::{BooleanOp, csg_boolean_indexed};
    use crate::mesh::IndexedMesh;

    fn signed_area_2d(mesh: &IndexedMesh) -> f64 {
        mesh.faces.iter().map(|f| {
            let a = mesh.vertices.position(f.vertices[0]);
            let b = mesh.vertices.position(f.vertices[1]);
            let c = mesh.vertices.position(f.vertices[2]);
            let ab = b - a;
            let ac = c - a;
            ab.cross(&ac).norm() * 0.5
        }).sum()
    }

    fn make_disk(cx: f64, r: f64, n: usize) -> IndexedMesh {
        use crate::core::scalar::Point3r;
        Disk {
            center:   Point3r::new(cx, 0.0, 0.0),
            radius:   r,
            segments: n,
        }.build().expect("Disk build failed")
    }

    /// Two identical disks — union should equal one disk.
    ///
    /// # Theorem
    ///
    /// `A ∪ A = A`  — union of a set with itself is the set.
    #[test]
    fn identical_disks_union_equals_single() {
        let a = make_disk(0.0, 1.0, 64);
        let b = make_disk(0.0, 1.0, 64);
        let expected = std::f64::consts::PI;

        let union = csg_boolean_indexed(BooleanOp::Union, &a, &b)
            .expect("union should succeed");
        let area = signed_area_2d(&union);
        let err = (area - expected).abs() / expected;
        assert!(err < 0.02,
            "identical disk union area error {:.1}% > 2% (got {area:.4}, expected {expected:.4})",
            err * 100.0);
    }

    /// Two disks offset by d = r.
    ///
    /// # Theorem — Circular segment area
    ///
    /// For two circles of radius r with centres distance d apart (d < 2r):
    /// `θ = arccos(d / 2r)`,  `A_∩ = 2r²(θ − sin θ cos θ)`.
    #[test]
    fn offset_disks_intersection_area() {
        let r = 1.0_f64;
        let d = r;
        let a = make_disk(-d / 2.0, r, 128);
        let b = make_disk( d / 2.0, r, 128);

        let inter = csg_boolean_indexed(BooleanOp::Intersection, &a, &b)
            .expect("intersection should succeed");

        let theta = (d / (2.0 * r)).acos();
        let expected = 2.0 * r * r * (theta - theta.sin() * theta.cos());
        let area = signed_area_2d(&inter);
        let err = (area - expected).abs() / expected;
        assert!(err < 0.02,
            "disk intersection area error {:.1}% > 2% (got {area:.4}, expected {expected:.4})",
            err * 100.0);
    }

    /// Inclusion-exclusion: `area(A) + area(B) ≈ area(A∪B) + area(A∩B)`.
    ///
    /// # Theorem — Inclusion-Exclusion Principle
    ///
    /// For measurable sets A, B: `μ(A) + μ(B) = μ(A∪B) + μ(A∩B)`.
    #[test]
    fn disk_inclusion_exclusion() {
        let r = 1.0_f64;
        let d = r;
        let a = make_disk(-d / 2.0, r, 128);
        let b = make_disk( d / 2.0, r, 128);

        let union = csg_boolean_indexed(BooleanOp::Union, &a, &b)
            .expect("union failed");
        let inter = csg_boolean_indexed(BooleanOp::Intersection, &a, &b)
            .expect("intersection failed");

        let area_a = signed_area_2d(&a);
        let area_b = signed_area_2d(&b);
        let area_u = signed_area_2d(&union);
        let area_i = signed_area_2d(&inter);

        let lhs = area_a + area_b;
        let rhs = area_u + area_i;
        let err = (lhs - rhs).abs() / lhs;
        assert!(err < 0.02,
            "inclusion-exclusion error {:.1}% > 2%  (A+B={lhs:.4}, U+I={rhs:.4})",
            err * 100.0);
    }

    /// Difference `A\B`: `area(A\B) = area(A) − area(A∩B)`.
    #[test]
    fn disk_difference_area() {
        let r = 1.0_f64;
        let d = r;
        let a = make_disk(-d / 2.0, r, 128);
        let b = make_disk( d / 2.0, r, 128);

        let diff = csg_boolean_indexed(BooleanOp::Difference, &a, &b)
            .expect("difference failed");

        let theta   = (d / (2.0 * r)).acos();
        let a_inter = 2.0 * r * r * (theta - theta.sin() * theta.cos());
        let expected = std::f64::consts::PI * r * r - a_inter;
        let area = signed_area_2d(&diff);
        let err = (area - expected).abs() / expected;
        assert!(err < 0.02,
            "disk difference area error {:.1}% > 2% (got {area:.4}, expected {expected:.4})",
            err * 100.0);
    }
}
