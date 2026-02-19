//! Coplanar 2-D Boolean operations for flat surface meshes.
//!
//! When both operand meshes are flat (all triangles share a common plane),
//! the standard 3-D BSP and Mesh-Arrangement pipelines break down because
//! `orient_3d` degenerates for all coplanar query points.
//!
//! ## Algorithm — 2-D Sutherland-Hodgman clipping
//!
//! For each source triangle T:
//!
//! 1. **Project** vertices to the shared 2-D plane via `PlaneBasis::project`.
//! 2. **Fast-path — disjoint**: if T's AABB has no overlap with any opposing
//!    AABB, T is fully disjoint from the opposing mesh → emit (outside) or
//!    skip (inside).
//! 3. **Fast-path — fully inside/outside**: test centroid + 3 vertices against
//!    the AABB-candidate triangles using 2-D point-in-triangle.  If unanimous
//!    inside → emit (inside) or skip (outside); if unanimous outside → skip
//!    (inside) or emit (outside).
//! 4. **Exact boundary clipping**: for each AABB-overlapping opposing triangle
//!    B_i, compute `clip2d_inside_tri(T, B_i)` using the 2-D Sutherland-Hodgman
//!    algorithm (one half-plane clip per B-edge, in 2-D).
//!    - **Intersection**: emit each `T ∩ B_i` fragment (disjoint since B-tris
//!      tile the disk without interior overlap).
//!    - **Outside**: progressive subtraction — start with `remaining = {T}`;
//!      for each B_i remove `inside(B_i)` pieces from every remaining poly.
//!
//! ## Why 2-D clipping is necessary
//!
//! `clip_polygon_to_halfplane` from `csg::clip` uses `orient_3d` for the
//! inside/outside decision.  For points in the shared flat plane, every query
//! returns `Degenerate` (signed tet volume = 0), so all points appear "inside"
//! every half-space — making the 3-D clipper useless for this case.
//!
//! The 2-D Sutherland-Hodgman implementation here uses the 2-D cross product
//! sign to classify query points, which is exact for planar geometry.
//!
//! ## Output quality
//!
//! Boundary vertices are computed as exact 2-D edge-edge intersections and
//! lifted back to 3-D, then registered via `VertexPool::insert_or_weld` →
//! shared seam vertices across adjacent A-triangles produce a manifold,
//! smooth boundary curve with no staircase artefacts.

#![allow(missing_docs)]

use crate::core::scalar::{Point3r, Real, Vector3r};
use crate::csg::boolean::BooleanOp;
use crate::csg::clip::fan_triangulate;
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;

// ── PlaneBasis ─────────────────────────────────────────────────────────────────

pub(crate) struct PlaneBasis {
    pub(crate) origin: Point3r,
    pub(crate) u: Vector3r,
    pub(crate) v: Vector3r,
    pub(crate) normal: Vector3r,
}

impl PlaneBasis {
    fn from_triangle(a: &Point3r, b: &Point3r, c: &Point3r) -> Option<Self> {
        let ab = b - a;
        let ac = c - a;
        let n = ab.cross(&ac);
        let nl = n.norm();
        if nl < 1e-20 { return None; }
        let ul = ab.norm();
        if ul < 1e-20 { return None; }
        let u = ab / ul;
        let normal = n / nl;
        let v = normal.cross(&u).normalize();
        Some(Self { origin: *a, u, v, normal })
    }

    #[inline]
    pub(crate) fn project(&self, p: &Point3r) -> [Real; 2] {
        let d = p - self.origin;
        [d.dot(&self.u), d.dot(&self.v)]
    }

    /// Lift a 2-D point (u,v) back to 3-D.
    #[inline]
    fn lift(&self, u: Real, v: Real) -> Point3r {
        self.origin + self.u * u + self.v * v
    }
}

// ── 2-D point and AABB helpers ─────────────────────────────────────────────────

/// 2-D signed cross product of vectors (b-a) and (p-a).
/// Positive → p is to the left of a→b (CCW); Negative → right (CW).
#[inline]
fn cross2(ax: Real, ay: Real, bx: Real, by: Real, px: Real, py: Real) -> Real {
    (bx - ax) * (py - ay) - (by - ay) * (px - ax)
}

/// Test whether 2-D point (px,py) lies inside or on the boundary of the
/// CCW-wound triangle (ax,ay)→(bx,by)→(cx,cy).
#[inline]
fn point_in_tri_2d(px: Real, py: Real,
                   ax: Real, ay: Real,
                   bx: Real, by: Real,
                   cx: Real, cy: Real) -> bool {
    let d0 = cross2(ax,ay,bx,by,px,py);
    let d1 = cross2(bx,by,cx,cy,px,py);
    let d2 = cross2(cx,cy,ax,ay,px,py);
    // Point is inside iff all cross products have the same sign (or zero).
    let neg = d0<0.0 || d1<0.0 || d2<0.0;
    let pos = d0>0.0 || d1>0.0 || d2>0.0;
    !(neg && pos)
}

/// Test whether 2-D point is inside the union of a set of triangles
/// (given as flat `[ax,ay,bx,by,cx,cy]` arrays).
#[inline]
fn point_in_union_2d(px: Real, py: Real, tris: &[[Real; 6]]) -> bool {
    tris.iter().any(|t| point_in_tri_2d(px, py, t[0], t[1], t[2], t[3], t[4], t[5]))
}

/// 2-D AABB of a triangle: `[min_u, min_v, max_u, max_v]`.
#[inline]
fn aabb2(ax: Real, ay: Real, bx: Real, by: Real, cx: Real, cy: Real) -> [Real; 4] {
    [ax.min(bx).min(cx), ay.min(by).min(cy),
     ax.max(bx).max(cx), ay.max(by).max(cy)]
}

/// True if two 2-D AABBs intersect (inclusive boundary).
#[inline]
fn aabb_overlaps(a: &[Real; 4], b: &[Real; 4]) -> bool {
    a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3]
}

// ── 2-D Sutherland-Hodgman clipping ───────────────────────────────────────────
//
// All clipping is performed in the 2-D projected coordinate system.
// The polygon is represented as a `Vec<[Real; 2]>` of 2-D points.
// After clipping, results are lifted back to 3-D via PlaneBasis::lift.

/// Clip a 2-D polygon against the left half-plane of the directed edge (ax,ay)→(bx,by).
///
/// "Left" means cross product ≥ 0 (i.e. points where `cross2(a,b,p) >= 0`).
/// Sutherland-Hodgman single-edge pass.
fn clip2d_halfplane(poly: &[[Real; 2]], ax: Real, ay: Real, bx: Real, by: Real)
    -> Vec<[Real; 2]>
{
    if poly.len() < 2 { return Vec::new(); }
    let n = poly.len();
    let mut out = Vec::with_capacity(n + 1);

    for i in 0..n {
        let s = poly[i];
        let e = poly[(i + 1) % n];
        let s_cross = cross2(ax, ay, bx, by, s[0], s[1]);
        let e_cross = cross2(ax, ay, bx, by, e[0], e[1]);
        let s_in = s_cross >= 0.0;
        let e_in = e_cross >= 0.0;

        match (s_in, e_in) {
            (true, true) => { out.push(e); }
            (true, false) => {
                // Compute intersection of segment s→e with the clip line.
                let denom = s_cross - e_cross;
                if denom.abs() > 1e-30 {
                    let t = s_cross / denom;
                    out.push([s[0] + (e[0]-s[0])*t, s[1] + (e[1]-s[1])*t]);
                }
                // else parallel — emit nothing (degenerate)
            }
            (false, true) => {
                let denom = s_cross - e_cross;
                if denom.abs() > 1e-30 {
                    let t = s_cross / denom;
                    out.push([s[0] + (e[0]-s[0])*t, s[1] + (e[1]-s[1])*t]);
                }
                out.push(e);
            }
            (false, false) => {}
        }
    }
    out
}

/// Clip a 2-D polygon to the inside of a **CCW** triangle (d,e,f) given by 2-D coords.
///
/// Returns the intersection polygon, or empty if fully outside.
/// Uses three half-plane clips (one per edge) via `clip2d_halfplane`.
fn clip2d_inside_tri(poly: &[[Real; 2]],
                     dx: Real, dy: Real,
                     ex: Real, ey: Real,
                     fx: Real, fy: Real) -> Vec<[Real; 2]>
{
    // Detect CW triangle and swap winding so we always clip to CCW inside.
    let area2 = cross2(dx, dy, ex, ey, fx, fy);
    let (dx, dy, ex, ey, fx, fy) = if area2 < 0.0 {
        (dx, dy, fx, fy, ex, ey) // swap e↔f to make CCW
    } else {
        (dx, dy, ex, ey, fx, fy)
    };
    // If area is exactly 0 (degenerate triangle), return empty.
    if area2 == 0.0 { return Vec::new(); }

    let p = clip2d_halfplane(poly, dx, dy, ex, ey);
    if p.len() < 3 { return Vec::new(); }
    let p = clip2d_halfplane(&p, ex, ey, fx, fy);
    if p.len() < 3 { return Vec::new(); }
    clip2d_halfplane(&p, fx, fy, dx, dy)
}

/// Decompose the 2-D `poly \ triangle(d,e,f)` into at most 3 disjoint convex pieces.
///
/// Uses the complement decomposition (triangle's complement = 3 half-planes):
/// ```text
///   piece_0 = poly ∩ right(d→e)        [outside edge d→e]
///   piece_1 = poly ∩ left(d→e) ∩ right(e→f)
///   piece_2 = poly ∩ left(d→e) ∩ left(e→f) ∩ right(f→d)
/// ```
/// For a **CCW** triangle (d,e,f), these three regions partition the complement.
fn split2d_outside_tri(poly: &[[Real; 2]],
                       dx: Real, dy: Real,
                       ex: Real, ey: Real,
                       fx: Real, fy: Real) -> Vec<Vec<[Real; 2]>>
{
    // Ensure CCW orientation.
    let area2 = cross2(dx, dy, ex, ey, fx, fy);
    let (dx, dy, ex, ey, fx, fy) = if area2 < 0.0 {
        (dx, dy, fx, fy, ex, ey)
    } else {
        (dx, dy, ex, ey, fx, fy)
    };
    if area2 == 0.0 { return vec![poly.to_vec()]; }

    let mut out = Vec::with_capacity(3);

    // piece_0: right of d→e (= outside edge d→e for CCW triangle)
    let p0 = clip2d_halfplane(poly, ex, ey, dx, dy); // reversed = right side
    if p0.len() >= 3 { out.push(p0); }

    // piece_1 & 2: first restrict to left of d→e (inside that edge)
    let in_de = clip2d_halfplane(poly, dx, dy, ex, ey);
    if in_de.len() >= 3 {
        // piece_1: left(d→e) ∩ right(e→f)  → right of e→f = left of f→e
        let p1 = clip2d_halfplane(&in_de, fx, fy, ex, ey);
        if p1.len() >= 3 { out.push(p1); }

        // piece_2: left(d→e) ∩ left(e→f) ∩ right(f→d)
        let in_ef = clip2d_halfplane(&in_de, ex, ey, fx, fy);
        if in_ef.len() >= 3 {
            // right of f→d = left of d→f
            let p2 = clip2d_halfplane(&in_ef, dx, dy, fx, fy);
            if p2.len() >= 3 { out.push(p2); }
        }
    }
    out
}

// ── Lift 2-D polygon to 3-D ───────────────────────────────────────────────────

fn lift_poly(poly2d: &[[Real; 2]], basis: &PlaneBasis) -> Vec<Point3r> {
    poly2d.iter().map(|&[u,v]| basis.lift(u, v)).collect()
}

// ── Emit helpers ───────────────────────────────────────────────────────────────

fn emit_one(p0: Point3r, p1: Point3r, p2: Point3r,
            basis: &PlaneBasis, region: crate::core::index::RegionId,
            result: &mut Vec<FaceData>, pool: &mut VertexPool) {
    let ab = p1 - p0; let ac = p2 - p0;
    let fn_ = ab.cross(&ac);
    if fn_.norm() < 1e-20 { return; }
    let flip = fn_.dot(&basis.normal) < 0.0;
    let (o0,o1,o2) = if flip { (p0,p2,p1) } else { (p0,p1,p2) };
    let v0 = pool.insert_or_weld(o0, basis.normal);
    let v1 = pool.insert_or_weld(o1, basis.normal);
    let v2 = pool.insert_or_weld(o2, basis.normal);
    if v0!=v1 && v1!=v2 && v0!=v2 {
        result.push(FaceData::new(v0, v1, v2, region));
    }
}

fn emit_poly2d(poly2d: &[[Real; 2]], basis: &PlaneBasis,
               region: crate::core::index::RegionId,
               result: &mut Vec<FaceData>, pool: &mut VertexPool) {
    let poly3d = lift_poly(poly2d, basis);
    for [t0,t1,t2] in fan_triangulate(&poly3d) {
        emit_one(t0, t1, t2, basis, region, result, pool);
    }
}

// ── detect_flat_plane ──────────────────────────────────────────────────────────

pub(crate) fn detect_flat_plane(faces: &[FaceData], pool: &VertexPool) -> Option<PlaneBasis> {
    let mut basis: Option<PlaneBasis> = None;
    for face in faces {
        let a = pool.position(face.vertices[0]);
        let b = pool.position(face.vertices[1]);
        let c = pool.position(face.vertices[2]);
        if let Some(b0) = PlaneBasis::from_triangle(a, b, c) { basis = Some(b0); break; }
    }
    let basis = basis?;
    const TOL: Real = 1e-6;
    for face in faces {
        for &vid in &face.vertices {
            if (pool.position(vid) - basis.origin).dot(&basis.normal).abs() > TOL {
                return None;
            }
        }
    }
    Some(basis)
}

// ── Pre-computed triangle data ─────────────────────────────────────────────────

struct TriData {
    coords2d: [Real; 6],    // [ax,ay, bx,by, cx,cy] for point-in-union and clipping
    aabb2d:   [Real; 4],    // [min_u, min_v, max_u, max_v]
    verts3d:  [Point3r; 3], // 3-D positions (needed to emit original triangles)
}

fn build_tri_data(faces: &[FaceData], pool: &VertexPool, basis: &PlaneBasis) -> Vec<TriData> {
    faces.iter().map(|f| {
        let p = *pool.position(f.vertices[0]);
        let q = *pool.position(f.vertices[1]);
        let r = *pool.position(f.vertices[2]);
        let [px,py] = basis.project(&p);
        let [qx,qy] = basis.project(&q);
        let [rx,ry] = basis.project(&r);
        TriData {
            coords2d: [px,py,qx,qy,rx,ry],
            aabb2d:   aabb2(px,py,qx,qy,rx,ry),
            verts3d:  [p,q,r],
        }
    }).collect()
}

// ── Core: process one source triangle against opposing triangles ───────────────

/// Process one source triangle (in 2-D) against opposing triangles.
///
/// `want_inside = true`  → emit src ∩ (∪ opp)   (Intersection)
/// `want_inside = false` → emit src \ (∪ opp)   (Difference / Union B\A)
fn process_triangle(
    src: &[Real; 6],         // [ax,ay,bx,by,cx,cy] of source in 2-D
    src_3d: &[Point3r; 3],   // 3-D positions for fast-path emit
    src_aabb: &[Real; 4],
    opp: &[TriData],
    opp_tris: &[[Real; 6]],  // 2-D coords of ALL opposing triangles
    opp_aabbs: &[[Real; 4]],
    want_inside: bool,
    basis: &PlaneBasis,
    region: crate::core::index::RegionId,
    result: &mut Vec<FaceData>,
    pool: &mut VertexPool,
) {
    // ── Fast path 1: AABB disjoint from all opposing triangles ────────────────
    let candidates: Vec<usize> = opp_aabbs.iter().enumerate()
        .filter(|(_,ob)| aabb_overlaps(src_aabb, ob))
        .map(|(i,_)| i)
        .collect();

    if candidates.is_empty() {
        if !want_inside {
            emit_one(src_3d[0], src_3d[1], src_3d[2], basis, region, result, pool);
        }
        return;
    }

    // ── Fast path 2 & 3: vertex + centroid classification ─────────────────────
    // Build a local view of only the candidate opposing triangles for speed.
    let cand_tris: Vec<[Real; 6]> = candidates.iter().map(|&i| opp_tris[i]).collect();

    let [ax,ay,bx,by,cx,cy] = *src;
    let mx = (ax+bx+cx) / 3.0;
    let my = (ay+by+cy) / 3.0;

    let va = point_in_union_2d(ax, ay, &cand_tris);
    let vb = point_in_union_2d(bx, by, &cand_tris);
    let vc = point_in_union_2d(cx, cy, &cand_tris);
    let vm = point_in_union_2d(mx, my, &cand_tris);

    let all_in  = va && vb && vc && vm;
    let all_out = !va && !vb && !vc && !vm;

    if all_in {
        if want_inside {
            emit_one(src_3d[0], src_3d[1], src_3d[2], basis, region, result, pool);
        }
        return;
    }
    if all_out {
        if !want_inside {
            emit_one(src_3d[0], src_3d[1], src_3d[2], basis, region, result, pool);
        }
        return;
    }

    // ── Boundary: exact 2-D Sutherland-Hodgman clipping ───────────────────────
    let src_poly: Vec<[Real; 2]> = vec![[ax,ay],[bx,by],[cx,cy]];

    if want_inside {
        // Intersection: emit A ∩ B_i for each candidate.
        // B-triangles tile a disk without interior overlap → fragments are disjoint.
        for &ci in &candidates {
            let [dx,dy,ex,ey,fx,fy] = opp[ci].coords2d;
            let inside = clip2d_inside_tri(&src_poly, dx, dy, ex, ey, fx, fy);
            if inside.len() >= 3 {
                emit_poly2d(&inside, basis, region, result, pool);
            }
        }
    } else {
        // Outside: A \ (∪B) via progressive subtraction.
        let mut remaining: Vec<Vec<[Real; 2]>> = vec![src_poly];

        for &ci in &candidates {
            let [dx,dy,ex,ey,fx,fy] = opp[ci].coords2d;
            let mut new_rem: Vec<Vec<[Real; 2]>> = Vec::new();

            for poly in &remaining {
                let inside = clip2d_inside_tri(poly, dx, dy, ex, ey, fx, fy);
                if inside.len() < 3 {
                    new_rem.push(poly.clone());
                    continue;
                }
                // poly straddles B_i: split each fan sub-triangle into outside pieces.
                let poly3d = lift_poly(poly, basis);
                for [t0,t1,t2] in fan_triangulate(&poly3d) {
                    let t2d = [basis.project(&t0), basis.project(&t1), basis.project(&t2)];
                    let t_inside = clip2d_inside_tri(&t2d, dx, dy, ex, ey, fx, fy);
                    if t_inside.len() < 3 {
                        new_rem.push(t2d.to_vec());
                    } else {
                        for piece in split2d_outside_tri(&t2d, dx, dy, ex, ey, fx, fy) {
                            new_rem.push(piece);
                        }
                    }
                }
            }
            remaining = new_rem;
        }

        for poly in &remaining {
            emit_poly2d(poly, basis, region, result, pool);
        }
    }
}

// ── boolean_coplanar ───────────────────────────────────────────────────────────

pub(crate) fn boolean_coplanar(
    op: BooleanOp,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
    basis: &PlaneBasis,
) -> Vec<FaceData> {
    let mut result: Vec<FaceData> = Vec::new();

    let b_data = build_tri_data(faces_b, pool, basis);
    let a_data = build_tri_data(faces_a, pool, basis);

    let b_tris: Vec<[Real; 6]> = b_data.iter().map(|b| b.coords2d).collect();
    let a_tris: Vec<[Real; 6]> = a_data.iter().map(|a| a.coords2d).collect();
    let b_aabbs: Vec<[Real; 4]> = b_data.iter().map(|b| b.aabb2d).collect();
    let a_aabbs: Vec<[Real; 4]> = a_data.iter().map(|a| a.aabb2d).collect();

    // ── A-side ────────────────────────────────────────────────────────────────

    for (ai, fa) in faces_a.iter().enumerate() {
        let src     = &a_tris[ai];
        let src_3d  = &a_data[ai].verts3d;
        let aabb    = &a_aabbs[ai];

        match op {
            BooleanOp::Union => {
                // Emit all of A; B\A is handled in the B-side pass.
                emit_one(src_3d[0], src_3d[1], src_3d[2], basis, fa.region, &mut result, pool);
            }
            BooleanOp::Intersection => {
                process_triangle(
                    src, src_3d, aabb,
                    &b_data, &b_tris, &b_aabbs,
                    true, basis, fa.region, &mut result, pool,
                );
            }
            BooleanOp::Difference => {
                process_triangle(
                    src, src_3d, aabb,
                    &b_data, &b_tris, &b_aabbs,
                    false, basis, fa.region, &mut result, pool,
                );
            }
        }
    }

    // ── B-side (Union only): add B \ A ───────────────────────────────────────

    if matches!(op, BooleanOp::Union) {
        for (bi, fb) in faces_b.iter().enumerate() {
            process_triangle(
                &b_tris[bi], &b_data[bi].verts3d, &b_aabbs[bi],
                &a_data, &a_tris, &a_aabbs,
                false, basis, fb.region, &mut result, pool,
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

    fn area(mesh: &IndexedMesh) -> f64 {
        mesh.faces.iter().map(|f| {
            let a = mesh.vertices.position(f.vertices[0]);
            let b = mesh.vertices.position(f.vertices[1]);
            let c = mesh.vertices.position(f.vertices[2]);
            (b-a).cross(&(c-a)).norm() * 0.5
        }).sum()
    }

    fn disk(cx: f64, r: f64, n: usize) -> IndexedMesh {
        use crate::core::scalar::Point3r;
        Disk { center: Point3r::new(cx,0.,0.), radius: r, segments: n }
            .build().unwrap()
    }

    #[test]
    fn identical_disks_union_equals_single() {
        let a = disk(0., 1., 64); let b = disk(0., 1., 64);
        let u = csg_boolean_indexed(BooleanOp::Union, &a, &b).unwrap();
        let err = (area(&u) - std::f64::consts::PI).abs() / std::f64::consts::PI;
        assert!(err < 0.02, "union err {:.1}%", err*100.);
    }

    #[test]
    fn offset_disks_intersection_area() {
        let (r,d) = (1.0_f64, 1.0_f64);
        let a = disk(-d/2., r, 128); let b = disk(d/2., r, 128);
        let inter = csg_boolean_indexed(BooleanOp::Intersection, &a, &b).unwrap();
        let th = (d/(2.*r)).acos();
        let exp = 2.*r*r*(th - th.sin()*th.cos());
        let err = (area(&inter)-exp).abs()/exp;
        assert!(err < 0.02, "inter err {:.1}%", err*100.);
    }

    #[test]
    fn disk_inclusion_exclusion() {
        let (r,d) = (1.0_f64, 1.0_f64);
        let a = disk(-d/2., r, 128); let b = disk(d/2., r, 128);
        let u = csg_boolean_indexed(BooleanOp::Union, &a, &b).unwrap();
        let i = csg_boolean_indexed(BooleanOp::Intersection, &a, &b).unwrap();
        let lhs = area(&a)+area(&b);
        let rhs = area(&u)+area(&i);
        let err = (lhs-rhs).abs()/lhs;
        assert!(err < 0.02, "incl-excl err {:.1}%", err*100.);
    }

    #[test]
    fn disk_difference_area() {
        let (r,d) = (1.0_f64, 1.0_f64);
        let a = disk(-d/2., r, 128); let b = disk(d/2., r, 128);
        let diff = csg_boolean_indexed(BooleanOp::Difference, &a, &b).unwrap();
        let th = (d/(2.*r)).acos();
        let exp = std::f64::consts::PI*r*r - 2.*r*r*(th-th.sin()*th.cos());
        let err = (area(&diff)-exp).abs()/exp;
        assert!(err < 0.02, "diff err {:.1}%", err*100.);
    }
}
