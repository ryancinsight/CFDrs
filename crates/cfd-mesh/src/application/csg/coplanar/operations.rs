//! Coplanar Boolean Operations Kernel

use crate::domain::core::scalar::{Point3r, Real};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;
use crate::application::csg::boolean::BooleanOp;
use super::basis::PlaneBasis;
use super::geometry2d::{aabb2, aabb_overlaps, point_in_tri_2d_exact, point_in_union_2d_exact};
use crate::application::csg::clip::{clip_polygon_to_triangle, fan_triangulate, split_polygon_outside_triangle};

fn emit_one(
    p0: Point3r,
    p1: Point3r,
    p2: Point3r,
    basis: &PlaneBasis,
    region: crate::domain::core::index::RegionId,
    result: &mut Vec<FaceData>,
    pool: &mut VertexPool,
) {
    let ab = p1 - p0;
    let ac = p2 - p0;
    let fn_ = ab.cross(&ac);
    if fn_.norm() < 1e-20 {
        return;
    }
    let flip = fn_.dot(&basis.normal) < 0.0;
    let (o0, o1, o2) = if flip { (p0, p2, p1) } else { (p0, p1, p2) };
    let v0 = pool.insert_or_weld(o0, basis.normal);
    let v1 = pool.insert_or_weld(o1, basis.normal);
    let v2 = pool.insert_or_weld(o2, basis.normal);
    if v0 != v1 && v1 != v2 && v0 != v2 {
        result.push(FaceData::new(v0, v1, v2, region));
    }
}

fn emit_poly2d(
    poly2d: &[[Real; 2]],
    basis: &PlaneBasis,
    region: crate::domain::core::index::RegionId,
    result: &mut Vec<FaceData>,
    pool: &mut VertexPool,
) {
    let poly3d: Vec<Point3r> = poly2d.iter().map(|&[u, v]| basis.lift(u, v)).collect();
    for [t0, t1, t2] in fan_triangulate(&poly3d) {
        emit_one(t0, t1, t2, basis, region, result, pool);
    }
}

// ── Pre-computed triangle data ─────────────────────────────────────────────────

struct TriData {
    coords2d: [Real; 6],   // [ax,ay, bx,by, cx,cy] for point-in-union and clipping
    aabb2d: [Real; 4],     // [min_u, min_v, max_u, max_v]
    verts3d: [Point3r; 3], // 3-D positions (needed to emit original triangles)
}

fn build_tri_data(faces: &[FaceData], pool: &VertexPool, basis: &PlaneBasis) -> Vec<TriData> {
    faces
        .iter()
        .map(|f| {
            let p = *pool.position(f.vertices[0]);
            let q = *pool.position(f.vertices[1]);
            let r = *pool.position(f.vertices[2]);
            let [px, py] = basis.project(&p);
            let [qx, qy] = basis.project(&q);
            let [rx, ry] = basis.project(&r);
            TriData {
                coords2d: [px, py, qx, qy, rx, ry],
                aabb2d: aabb2(px, py, qx, qy, rx, ry),
                verts3d: [p, q, r],
            }
        })
        .collect()
}

// ── Core: process one source triangle against opposing triangles ───────────────

/// Process one source triangle (in 2-D) against opposing triangles.
///
/// `want_inside = true`  → emit src ∩ (∪ opp)   (Intersection)
/// `want_inside = false` → emit src \ (∪ opp)   (Difference / Union B\A)
fn process_triangle(
    src: &[Real; 6],       // [ax,ay,bx,by,cx,cy] of source in 2-D
    src_3d: &[Point3r; 3], // 3-D positions for fast-path emit
    src_aabb: &[Real; 4],
    opp: &[TriData],
    opp_tris: &[[Real; 6]], // 2-D coords of ALL opposing triangles
    opp_aabbs: &[[Real; 4]],
    want_inside: bool,
    basis: &PlaneBasis,
    region: crate::domain::core::index::RegionId,
    result: &mut Vec<FaceData>,
    pool: &mut VertexPool,
) {
    let candidates: Vec<usize> = opp_aabbs
        .iter()
        .enumerate()
        .filter(|(_, ob)| aabb_overlaps(src_aabb, ob))
        .map(|(i, _)| i)
        .collect();

    if candidates.is_empty() {
        if !want_inside {
            emit_one(src_3d[0], src_3d[1], src_3d[2], basis, region, result, pool);
        }
        return;
    }

    let cand_tris: Vec<[Real; 6]> = candidates.iter().map(|&i| opp_tris[i]).collect();

    let [ax, ay, bx, by, cx, cy] = *src;

    let va = point_in_union_2d_exact(ax, ay, &cand_tris);
    let vb = point_in_union_2d_exact(bx, by, &cand_tris);
    let vc = point_in_union_2d_exact(cx, cy, &cand_tris);

    let all_in = va && vb && vc;
    let all_out_vertices = !va && !vb && !vc;
    let all_out = all_out_vertices
        && cand_tris.iter().all(|tri| {
            let [dx, dy, ex, ey, fx, fy] = *tri;
            !point_in_tri_2d_exact(dx, dy, ax, ay, bx, by, cx, cy)
                && !point_in_tri_2d_exact(ex, ey, ax, ay, bx, by, cx, cy)
                && !point_in_tri_2d_exact(fx, fy, ax, ay, bx, by, cx, cy)
        });

    if all_in {
        if want_inside {
            emit_one(src_3d[0], src_3d[1], src_3d[2], basis, region, result, pool);
        }
        return;
    }
    if all_out && want_inside {
        return;
    }

    let src_poly: Vec<[Real; 2]> = vec![[ax, ay], [bx, by], [cx, cy]];

    if want_inside {
        for &ci in &candidates {
            let [dx, dy, ex, ey, fx, fy] = opp[ci].coords2d;
            let inside = clip_polygon_to_triangle(&src_poly, dx, dy, ex, ey, fx, fy);
            if inside.len() >= 3 {
                emit_poly2d(&inside, basis, region, result, pool);
            }
        }
    } else {
        let mut remaining: Vec<Vec<[Real; 2]>> = vec![src_poly];

        for &ci in &candidates {
            let [dx, dy, ex, ey, fx, fy] = opp[ci].coords2d;
            let mut new_rem: Vec<Vec<[Real; 2]>> = Vec::new();

            for poly in &remaining {
                let inside = clip_polygon_to_triangle(poly, dx, dy, ex, ey, fx, fy);
                if inside.len() < 3 {
                    new_rem.push(poly.clone());
                    continue;
                }
                for piece in split_polygon_outside_triangle(poly, dx, dy, ex, ey, fx, fy) {
                    new_rem.push(piece);
                }
            }
            remaining = new_rem;
        }

        for poly in &remaining {
            emit_poly2d(poly, basis, region, result, pool);
        }
    }
}

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

    for (ai, fa) in faces_a.iter().enumerate() {
        let src = &a_tris[ai];
        let src_3d = &a_data[ai].verts3d;
        let aabb = &a_aabbs[ai];

        match op {
            BooleanOp::Union => {
                process_triangle(src, src_3d, aabb, &b_data, &b_tris, &b_aabbs, false, basis, fa.region, &mut result, pool);
                process_triangle(src, src_3d, aabb, &b_data, &b_tris, &b_aabbs, true, basis, fa.region, &mut result, pool);
            }
            BooleanOp::Intersection => {
                process_triangle(src, src_3d, aabb, &b_data, &b_tris, &b_aabbs, true, basis, fa.region, &mut result, pool);
            }
            BooleanOp::Difference => {
                process_triangle(src, src_3d, aabb, &b_data, &b_tris, &b_aabbs, false, basis, fa.region, &mut result, pool);
            }
        }
    }

    if matches!(op, BooleanOp::Union) {
        for (bi, fb) in faces_b.iter().enumerate() {
            process_triangle(
                &b_tris[bi],
                &b_data[bi].verts3d,
                &b_aabbs[bi],
                &a_data,
                &a_tris,
                &a_aabbs,
                false,
                basis,
                fb.region,
                &mut result,
                pool,
            );
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use crate::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
    use crate::domain::geometry::primitives::{Disk, PrimitiveMesh};
    use crate::domain::mesh::IndexedMesh;
    use super::super::geometry2d::polygon_area_2d;

    fn area(mesh: &IndexedMesh) -> f64 {
        mesh.faces
            .iter()
            .map(|f| {
                let a = mesh.vertices.position(f.vertices[0]);
                let b = mesh.vertices.position(f.vertices[1]);
                let c = mesh.vertices.position(f.vertices[2]);
                (b - a).cross(&(c - a)).norm() * 0.5
            })
            .sum()
    }

    fn disk(cx: f64, r: f64, n: usize) -> IndexedMesh {
        use crate::domain::core::scalar::Point3r;
        Disk {
            center: Point3r::new(cx, 0., 0.),
            radius: r,
            segments: n,
        }
        .build()
        .unwrap()
    }

    #[test]
    fn split_outside_partition_identity() {
        use crate::application::csg::clip::{
            clip_polygon_to_triangle, split_polygon_outside_triangle,
        };

        let poly = vec![[0.0, 0.0], [2.0, 0.0], [1.0, 2.0]];
        let (dx, dy, ex, ey, fx, fy) = (0.5, 0.0, 1.5, 0.0, 1.0, 1.0);
        let a_poly = polygon_area_2d(&poly);
        let a_inside = polygon_area_2d(&clip_polygon_to_triangle(&poly, dx, dy, ex, ey, fx, fy));
        let a_outside: f64 = split_polygon_outside_triangle(&poly, dx, dy, ex, ey, fx, fy)
            .iter()
            .map(|p| polygon_area_2d(p))
            .sum();
        let err1 = ((a_inside + a_outside) - a_poly).abs();
        assert!(err1 < 1e-12, "case 1: partition err={err1:.2e}");
    }

    #[test]
    fn identical_disks_union_equals_single() {
        let a = disk(0., 1., 64);
        let b = disk(0., 1., 64);
        let u = csg_boolean_indexed(BooleanOp::Union, &a, &b).unwrap();
        let err = (area(&u) - std::f64::consts::PI).abs() / std::f64::consts::PI;
        assert!(err < 0.02, "union err {:.1}%", err * 100.);
    }

    #[test]
    fn offset_disks_intersection_area() {
        let (r, d) = (1.0_f64, 1.0_f64);
        let a = disk(-d / 2., r, 128);
        let b = disk(d / 2., r, 128);
        let inter = csg_boolean_indexed(BooleanOp::Intersection, &a, &b).unwrap();
        let th = (d / (2. * r)).acos();
        let exp = 2. * r * r * (th - th.sin() * th.cos());
        let err = (area(&inter) - exp).abs() / exp;
        assert!(err < 0.02, "inter err {:.1}%", err * 100.);
    }

    #[test]
    fn disk_inclusion_exclusion() {
        let (r, d) = (1.0_f64, 1.0_f64);
        let a = disk(-d / 2., r, 128);
        let b = disk(d / 2., r, 128);
        let u = csg_boolean_indexed(BooleanOp::Union, &a, &b).unwrap();
        let i = csg_boolean_indexed(BooleanOp::Intersection, &a, &b).unwrap();
        let lhs = area(&a) + area(&b);
        let rhs = area(&u) + area(&i);
        let err = (lhs - rhs).abs() / lhs;
        assert!(err < 0.02, "incl-excl err {:.1}%", err * 100.);
    }

    #[test]
    fn disk_difference_area() {
        let (r, d) = (1.0_f64, 1.0_f64);
        let a = disk(-d / 2., r, 128);
        let b = disk(d / 2., r, 128);
        let diff = csg_boolean_indexed(BooleanOp::Difference, &a, &b).unwrap();
        let inter = csg_boolean_indexed(BooleanOp::Intersection, &a, &b).unwrap();
        let area_a = area(&a);
        let lhs = area(&diff) + area(&inter);
        let err = (lhs - area_a).abs() / area_a;
        assert!(err < 0.05, "inclusion-exclusion err {:.1}%", err * 100.);
    }
}
