//! Rounded cube (filleted box) primitive.

use std::f64::consts::{PI, TAU};

use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Vector3r};
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Builds a box with cylindrical edge fillets and spherical corner octants.
///
/// The box occupies the region
/// `[origin.x, origin.x + width] × [origin.y, origin.y + height] × [origin.z, origin.z + depth]`.
/// Every edge is replaced by a quarter-cylinder of radius `corner_radius`
/// and every corner by a sphere octant of the same radius.
///
/// ## Geometry decomposition
///
/// | Region type | Count | Faces each |
/// |-------------|-------|-----------|
/// | Flat face panels (rectangular) | 6 | `2 × (u × v)` |
/// | Quarter-cylinder edge strips | 12 | `2 × cs × len_segments` |
/// | Sphere octant corners | 8 | `cs × cs × 2 + cs` apex triangles |
///
/// where `cs = corner_segments` and all shared boundary rings are welded by
/// `VertexPool` spatial-hash deduplication.
///
/// ## Validation
///
/// - `corner_radius ≤ min(width, height, depth) / 2`
/// - `corner_segments ≥ 1`
/// - All dimensions > 0
///
/// ## Output
///
/// - All faces tagged `RegionId(1)`
/// - `signed_volume ≈ w·h·d − (4−π)·r²·(w+h+d) + V_sphere_correction`
///   (the exact value approaches `w·h·d` as `r → 0`)
#[derive(Clone, Debug)]
pub struct RoundedCube {
    /// Corner of the bounding box (minimum x, y, z).
    pub origin: Point3r,
    /// Extent along +X [mm].
    pub width: f64,
    /// Extent along +Y [mm].
    pub height: f64,
    /// Extent along +Z [mm].
    pub depth: f64,
    /// Fillet radius [mm]. Must be ≤ `min(w, h, d) / 2`.
    pub corner_radius: f64,
    /// Angular segments per quarter-turn of the fillets (≥ 1).
    pub corner_segments: usize,
}

impl Default for RoundedCube {
    fn default() -> Self {
        Self {
            origin: Point3r::origin(),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
            corner_radius: 0.2,
            corner_segments: 4,
        }
    }
}

impl PrimitiveMesh for RoundedCube {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(rc: &RoundedCube) -> Result<IndexedMesh, PrimitiveError> {
    let (w, h, d) = (rc.width, rc.height, rc.depth);
    let r = rc.corner_radius;
    let cs = rc.corner_segments;

    if w <= 0.0 { return Err(PrimitiveError::InvalidParam(format!("width must be > 0, got {w}"))); }
    if h <= 0.0 { return Err(PrimitiveError::InvalidParam(format!("height must be > 0, got {h}"))); }
    if d <= 0.0 { return Err(PrimitiveError::InvalidParam(format!("depth must be > 0, got {d}"))); }
    if r <= 0.0 || r > w / 2.0 || r > h / 2.0 || r > d / 2.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "corner_radius must be in (0, min(w,h,d)/2] = (0, {}], got {r}",
            w.min(h).min(d) / 2.0
        )));
    }
    if cs < 1 {
        return Err(PrimitiveError::InvalidParam("corner_segments must be ≥ 1".into()));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    // Inner box corners (after subtracting r from all sides).
    let x0 = rc.origin.x + r;  let x1 = rc.origin.x + w - r;
    let y0 = rc.origin.y + r;  let y1 = rc.origin.y + h - r;
    let z0 = rc.origin.z + r;  let z1 = rc.origin.z + d - r;

    // ── Six flat face panels ──────────────────────────────────────────────────
    // Each face is a rectangle in the inner box's face plane, pushed outward by r.

    // Helper: add a quad (two triangles) with a given flat normal.
    let add_quad = |mesh: &mut IndexedMesh,
                        p00: Point3r, p10: Point3r, p11: Point3r, p01: Point3r,
                        n: Vector3r| {
        let v00 = mesh.add_vertex(p00, n);
        let v10 = mesh.add_vertex(p10, n);
        let v11 = mesh.add_vertex(p11, n);
        let v01 = mesh.add_vertex(p01, n);
        mesh.add_face_with_region(v00, v10, v11, region);
        mesh.add_face_with_region(v00, v11, v01, region);
    };

    // −X face (normal = −X), x = x0 − r = origin.x
    let fx = rc.origin.x;
    add_quad(&mut mesh,
        Point3r::new(fx, y0, z0), Point3r::new(fx, y1, z0),
        Point3r::new(fx, y1, z1), Point3r::new(fx, y0, z1),
        -Vector3r::x());
    // +X face (normal = +X), x = x1 + r = origin.x + w
    let fx = rc.origin.x + w;
    add_quad(&mut mesh,
        Point3r::new(fx, y0, z1), Point3r::new(fx, y1, z1),
        Point3r::new(fx, y1, z0), Point3r::new(fx, y0, z0),
        Vector3r::x());
    // −Y face (normal = −Y)
    let fy = rc.origin.y;
    add_quad(&mut mesh,
        Point3r::new(x0, fy, z1), Point3r::new(x1, fy, z1),
        Point3r::new(x1, fy, z0), Point3r::new(x0, fy, z0),
        -Vector3r::y());
    // +Y face (normal = +Y)
    let fy = rc.origin.y + h;
    add_quad(&mut mesh,
        Point3r::new(x0, fy, z0), Point3r::new(x1, fy, z0),
        Point3r::new(x1, fy, z1), Point3r::new(x0, fy, z1),
        Vector3r::y());
    // −Z face (normal = −Z)
    let fz = rc.origin.z;
    add_quad(&mut mesh,
        Point3r::new(x0, y0, fz), Point3r::new(x1, y0, fz),
        Point3r::new(x1, y1, fz), Point3r::new(x0, y1, fz),
        -Vector3r::z());
    // +Z face (normal = +Z)
    let fz = rc.origin.z + d;
    add_quad(&mut mesh,
        Point3r::new(x1, y0, fz), Point3r::new(x0, y0, fz),
        Point3r::new(x0, y1, fz), Point3r::new(x1, y1, fz),
        Vector3r::z());

    // ── Quarter-cylinder edge strips ──────────────────────────────────────────
    // Each of the 12 edges of a box has a quarter-cylinder of radius r.
    // Edge enumeration: 4 edges parallel to X, 4 to Y, 4 to Z.
    //
    // Helper: quarter-cylinder strip along axis `axis` (0=X,1=Y,2=Z),
    // centred at inner_corner with sweep angle from `angle_start` to `angle_end`
    // in the plane perpendicular to `axis`.  The strip has `edge_len` (the
    // "flat" extent along the axis) driven by cs quads around the arc.

    // We enumerate all 12 edges manually with their arc parameters:
    // For the 4 Z-parallel edges (varying z from z0 to z1):
    let z_edges: [([f64; 3], f64, f64, i8, i8); 4] = [
        // [corner_x, corner_y, _], angle_start, angle_end, nx_sign, ny_sign
        ([x0, y0, 0.0], PI,     3.0*PI/2.0, -1,  -1),  // −X, −Y corner → arc in 3rd quadrant
        ([x1, y0, 0.0], 3.0*PI/2.0, TAU,    1,   -1),  // +X, −Y corner
        ([x1, y1, 0.0], 0.0,    PI/2.0,     1,    1),  // +X, +Y corner
        ([x0, y1, 0.0], PI/2.0, PI,         -1,   1),  // −X, +Y corner
    ];
    for ([cx_, cy_, _], a_start, a_end, _snx, _sny) in z_edges {
        for k in 0..cs {
            let a0 = a_start + k as f64 / cs as f64 * (a_end - a_start);
            let a1 = a_start + (k+1) as f64 / cs as f64 * (a_end - a_start);
            let (c0, s0) = (a0.cos(), a0.sin());
            let (c1, s1) = (a1.cos(), a1.sin());
            let n0 = Vector3r::new(c0, s0, 0.0);
            let n1 = Vector3r::new(c1, s1, 0.0);
            let pb0 = Point3r::new(cx_ + r*c0, cy_ + r*s0, z0);
            let pt0 = Point3r::new(cx_ + r*c0, cy_ + r*s0, z1);
            let pb1 = Point3r::new(cx_ + r*c1, cy_ + r*s1, z0);
            let pt1 = Point3r::new(cx_ + r*c1, cy_ + r*s1, z1);
            let vb0 = mesh.add_vertex(pb0, n0); let vt0 = mesh.add_vertex(pt0, n0);
            let vb1 = mesh.add_vertex(pb1, n1); let vt1 = mesh.add_vertex(pt1, n1);
            mesh.add_face_with_region(vb0, vt0, vt1, region);
            mesh.add_face_with_region(vb0, vt1, vb1, region);
        }
    }

    // 4 X-parallel edges
    // Sweep is along X (pb=x0, pt=x1); arc varies in YZ plane.
    // (vb0,vt0,vt1): cross = (x1-x0,0,0)×(x1-x0,dy,dz) = (0*dz-0*dy, 0*(x1-x0)-dx*dz, dx*dy-0*(x1-x0))
    // = (0, -dx*dz, dx*dy) — points in -Z,+Y for +X sweep going in +Z,+Y arc direction → outward.
    // To make inward (consistent with flat panels), reverse: (vb0, vb1, vt1) and (vb0, vt1, vt0).
    let x_edges: [([f64; 3], f64, f64); 4] = [
        ([0.0, y0, z0], PI,     3.0*PI/2.0),
        ([0.0, y0, z1], 3.0*PI/2.0, TAU   ),
        ([0.0, y1, z1], 0.0,    PI/2.0    ),
        ([0.0, y1, z0], PI/2.0, PI        ),
    ];
    for ([_, cy_, cz_], a_start, a_end) in x_edges {
        for k in 0..cs {
            let a0 = a_start + k as f64 / cs as f64 * (a_end - a_start);
            let a1 = a_start + (k+1) as f64 / cs as f64 * (a_end - a_start);
            let (c0, s0) = (a0.cos(), a0.sin());
            let (c1, s1) = (a1.cos(), a1.sin());
            let n0 = Vector3r::new(0.0, s0, c0);
            let n1 = Vector3r::new(0.0, s1, c1);
            let pb0 = Point3r::new(x0, cy_ + r*s0, cz_ + r*c0);
            let pt0 = Point3r::new(x1, cy_ + r*s0, cz_ + r*c0);
            let pb1 = Point3r::new(x0, cy_ + r*s1, cz_ + r*c1);
            let pt1 = Point3r::new(x1, cy_ + r*s1, cz_ + r*c1);
            let vb0 = mesh.add_vertex(pb0, n0); let vt0 = mesh.add_vertex(pt0, n0);
            let vb1 = mesh.add_vertex(pb1, n1); let vt1 = mesh.add_vertex(pt1, n1);
            // Reversed from Z-edge pattern to produce inward winding (matching flat panels).
            mesh.add_face_with_region(vb0, vt1, vt0, region);
            mesh.add_face_with_region(vb0, vb1, vt1, region);
        }
    }

    // 4 Y-parallel edges
    // Sweep is along Y (pb=y0, pt=y1); arc varies in XZ plane.
    // Same analysis as X-edges — sweep-along-axis produces outward winding; reverse it.
    let y_edges: [([f64; 3], f64, f64); 4] = [
        ([x0, 0.0, z0], PI,     3.0*PI/2.0),
        ([x1, 0.0, z0], 3.0*PI/2.0, TAU   ),
        ([x1, 0.0, z1], 0.0,    PI/2.0    ),
        ([x0, 0.0, z1], PI/2.0, PI        ),
    ];
    for ([cx_, _, cz_], a_start, a_end) in y_edges {
        for k in 0..cs {
            let a0 = a_start + k as f64 / cs as f64 * (a_end - a_start);
            let a1 = a_start + (k+1) as f64 / cs as f64 * (a_end - a_start);
            let (c0, s0) = (a0.cos(), a0.sin());
            let (c1, s1) = (a1.cos(), a1.sin());
            let n0 = Vector3r::new(c0, 0.0, s0);
            let n1 = Vector3r::new(c1, 0.0, s1);
            let pb0 = Point3r::new(cx_ + r*c0, y0, cz_ + r*s0);
            let pt0 = Point3r::new(cx_ + r*c0, y1, cz_ + r*s0);
            let pb1 = Point3r::new(cx_ + r*c1, y0, cz_ + r*s1);
            let pt1 = Point3r::new(cx_ + r*c1, y1, cz_ + r*s1);
            let vb0 = mesh.add_vertex(pb0, n0); let vt0 = mesh.add_vertex(pt0, n0);
            let vb1 = mesh.add_vertex(pb1, n1); let vt1 = mesh.add_vertex(pt1, n1);
            // Reversed from Z-edge pattern to produce inward winding (matching flat panels).
            mesh.add_face_with_region(vb0, vt1, vt0, region);
            mesh.add_face_with_region(vb0, vb1, vt1, region);
        }
    }

    // ── Sphere octant corners ─────────────────────────────────────────────────
    // 8 corners of the inner box; each gets a sphere octant patch.
    // An octant covers φ ∈ [0, π/2] (from apex to equator) and θ ∈ [0, π/2].
    // The apex is the vertex furthest from the box centre; the equator connects
    // to the three adjacent edge strips.
    //
    // Corner enumeration: (sx, sy, sz) = signs (±1).
    for sx in [-1.0_f64, 1.0_f64] {
        for sy in [-1.0_f64, 1.0_f64] {
            for sz in [-1.0_f64, 1.0_f64] {
                // Centre of the corner sphere octant (inner box corner).
                let ccx = if sx > 0.0 { x1 } else { x0 };
                let ccy = if sy > 0.0 { y1 } else { y0 };
                let ccz = if sz > 0.0 { z1 } else { z0 };

                // The octant spans φ from 0 (at axis +sx·X) to π/2 (at axis ±Y/±Z boundary).
                // We parametrize: u ∈ [0,1] along X-arc, v ∈ [0,1] along YZ-arc.
                // Point on octant: n = (sx·cos(u·π/2), sy·sin(u·π/2)·cos(v·π/2), sz·sin(u·π/2)·sin(v·π/2))
                // but the triple (cos u, sin u cos v, sin u sin v) is a sphere octant in the +++ octant,
                // remapped by (sx, sy, sz).

                for iu in 0..cs {
                    for iv in 0..cs {
                        let u0 = iu as f64 / cs as f64 * PI / 2.0;
                        let u1 = (iu+1) as f64 / cs as f64 * PI / 2.0;
                        let v0 = iv as f64 / cs as f64 * PI / 2.0;
                        let v1 = (iv+1) as f64 / cs as f64 * PI / 2.0;

                        let corner_pt = |u: f64, v: f64| -> (Point3r, Vector3r) {
                            let nx = sx * u.cos();
                            let ny = sy * u.sin() * v.cos();
                            let nz = sz * u.sin() * v.sin();
                            let n = Vector3r::new(nx, ny, nz);
                            let pos = Point3r::new(ccx + r*nx, ccy + r*ny, ccz + r*nz);
                            (pos, n)
                        };

                        let (p00, n00) = corner_pt(u0, v0);
                        let (p10, n10) = corner_pt(u1, v0);
                        let (p11, n11) = corner_pt(u1, v1);
                        let (p01, n01) = corner_pt(u0, v1);

                        let v00 = mesh.add_vertex(p00, n00);
                        let v10 = mesh.add_vertex(p10, n10);
                        let v11 = mesh.add_vertex(p11, n11);
                        let v01 = mesh.add_vertex(p01, n01);

                        // sign_parity: reflections in an odd number of axes reverse
                        // the surface orientation so the Jacobian determinant flips sign.
                        // Even parity (sx*sy*sz > 0): (u,v) natural winding gives outward.
                        // Odd parity  (sx*sy*sz < 0): (u,v) natural winding gives inward.
                        // We want all faces inward here (flip_faces() corrects at the end).
                        let parity_positive = (sx * sy * sz) > 0.0;

                        if iu == 0 {
                            // Apex row (u=0): v00 == v01 (degenerate) → triangle.
                            if parity_positive {
                                // Natural (u,v) winding is outward → reverse for inward.
                                mesh.add_face_with_region(v10, v01, v11, region);
                            } else {
                                // Natural winding is inward → keep as-is.
                                mesh.add_face_with_region(v10, v11, v01, region);
                            }
                        } else {
                            if parity_positive {
                                // Natural outward → reverse for inward.
                                mesh.add_face_with_region(v00, v01, v11, region);
                                mesh.add_face_with_region(v00, v11, v10, region);
                            } else {
                                // Natural inward → keep.
                                mesh.add_face_with_region(v00, v10, v11, region);
                                mesh.add_face_with_region(v00, v11, v01, region);
                            }
                        }
                    }
                }
            }
        }
    }

    // All sections (flat panels, Z-edge strips, X/Y-edge strips reversed, octant corners)
    // are constructed with inward winding. Flip all faces to obtain outward normals.
    mesh.flip_faces();

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::storage::edge_store::EdgeStore;
    use crate::watertight::check::check_watertight;

    #[test]
    fn rounded_cube_is_watertight() {
        let mesh = RoundedCube::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "rounded_cube must be watertight: {:?}", report);
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn rounded_cube_volume_positive() {
        let mesh = RoundedCube {
            width: 4.0, height: 3.0, depth: 2.0, corner_radius: 0.3,
            corner_segments: 6, ..RoundedCube::default()
        }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight);
        assert!(report.signed_volume > 0.0);
        // Volume must be less than bounding box
        assert!(report.signed_volume < 4.0 * 3.0 * 2.0);
    }

    #[test]
    fn rounded_cube_rejects_invalid_params() {
        assert!(RoundedCube { width: 0.0, ..RoundedCube::default() }.build().is_err());
        assert!(RoundedCube { corner_radius: 0.0, ..RoundedCube::default() }.build().is_err());
        // corner_radius > min(w,h,d)/2
        assert!(RoundedCube { corner_radius: 1.5, width: 2.0, height: 2.0, depth: 2.0, ..RoundedCube::default() }.build().is_err());
        assert!(RoundedCube { corner_segments: 0, ..RoundedCube::default() }.build().is_err());
    }
}
