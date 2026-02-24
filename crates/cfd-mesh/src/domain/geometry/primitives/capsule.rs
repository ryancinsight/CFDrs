//! Capsule primitive (cylinder + hemispherical end caps).

use std::f64::consts::{PI, TAU};

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a capsule: a closed right cylinder capped with two hemispheres.
///
/// The capsule is aligned with the +Y axis and centred at `center`.
/// The total height is `cylinder_height + 2 × radius`.
///
/// ## Mesh structure
///
/// | Section | Faces |
/// |---------|-------|
/// | Top hemisphere (`hemisphere_stacks` stacks from +Y apex to equator) | `segments × (hemisphere_stacks − 1) × 2 + segments` |
/// | Lateral cylinder | `2 × segments` |
/// | Bottom hemisphere (`hemisphere_stacks` stacks from equator to −Y apex) | `segments × (hemisphere_stacks − 1) × 2 + segments` |
///
/// Shared equatorial rings at `y = center.y ± cylinder_height/2` are
/// automatically welded by `VertexPool` spatial-hash deduplication.
///
/// ## Uses
///
/// Bacteria (E. coli, B. subtilis), pharmaceutical drug-delivery capsules,
/// elongated droplets, flexible particles in Lagrangian–Eulerian blood flow.
///
/// ## Output
///
/// - `signed_volume ≈ π r² (cylinder_height + 4r/3)`
/// - All faces tagged `RegionId(1)`
#[derive(Clone, Debug)]
pub struct Capsule {
    /// Hemisphere and cylinder radius [mm].
    pub radius: f64,
    /// Length of the cylindrical midsection [mm]. May be 0 (sphere).
    pub cylinder_height: f64,
    /// Centre of the capsule.
    pub center: Point3r,
    /// Angular subdivisions around the axis (≥ 3).
    pub segments: usize,
    /// Latitude subdivisions per hemisphere (≥ 1).
    pub hemisphere_stacks: usize,
}

impl Default for Capsule {
    fn default() -> Self {
        Self {
            radius: 0.5,
            cylinder_height: 1.0,
            center: Point3r::origin(),
            segments: 32,
            hemisphere_stacks: 8,
        }
    }
}

impl PrimitiveMesh for Capsule {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(cap: &Capsule) -> Result<IndexedMesh, PrimitiveError> {
    if cap.radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {}",
            cap.radius
        )));
    }
    if cap.cylinder_height < 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "cylinder_height must be ≥ 0, got {}",
            cap.cylinder_height
        )));
    }
    if cap.segments < 3 {
        return Err(PrimitiveError::TooFewSegments(cap.segments));
    }
    if cap.hemisphere_stacks < 1 {
        return Err(PrimitiveError::InvalidParam(
            "hemisphere_stacks must be ≥ 1".into(),
        ));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let r = cap.radius;
    let hl = cap.cylinder_height;
    let cx = cap.center.x;
    let cy = cap.center.y;
    let cz = cap.center.z;
    let ns = cap.segments;
    let hs = cap.hemisphere_stacks;

    // Y-offsets for the two hemisphere centres (= cylinder cap positions).
    let top_cy = cy + hl / 2.0;
    let bot_cy = cy - hl / 2.0;

    // ── Top hemisphere (φ: 0 → π/2, centre at y = top_cy) ───────────────────
    // φ = 0 → north pole (y = top_cy + r)
    // φ = π/2 → equator ring (y = top_cy, r_xy = r) → shared with cylinder top ring
    for i in 0..ns {
        let t0 = i as f64 / ns as f64 * TAU;
        let t1 = (i + 1) as f64 / ns as f64 * TAU;
        for j in 0..hs {
            let phi0 = j as f64 / hs as f64 * PI / 2.0;
            let phi1 = (j + 1) as f64 / hs as f64 * PI / 2.0;

            let vat = |theta: f64, phi: f64| -> (Point3r, Vector3r) {
                let sp = phi.sin();
                let cp = phi.cos();
                let ct = theta.cos();
                let st = theta.sin();
                let n = Vector3r::new(sp * ct, cp, sp * st);
                let p = Point3r::new(cx + r * sp * ct, top_cy + r * cp, cz + r * sp * st);
                (p, n)
            };

            let (p00, n00) = vat(t0, phi0);
            let (p10, n10) = vat(t1, phi0);
            let (p11, n11) = vat(t1, phi1);
            let (p01, n01) = vat(t0, phi1);

            let v00 = mesh.add_vertex(p00, n00);
            let v10 = mesh.add_vertex(p10, n10);
            let v11 = mesh.add_vertex(p11, n11);
            let v01 = mesh.add_vertex(p01, n01);

            if j == 0 {
                // North-pole row: v00 == v10 (φ=0) → triangle
                mesh.add_face_with_region(v10, v11, v01, region);
            } else {
                mesh.add_face_with_region(v00, v10, v11, region);
                mesh.add_face_with_region(v00, v11, v01, region);
            }
        }
    }

    // ── Cylinder lateral (only when cylinder_height > 0) ────────────────────
    if hl > 0.0 {
        for i in 0..ns {
            let t0 = i as f64 / ns as f64 * TAU;
            let t1 = (i + 1) as f64 / ns as f64 * TAU;

            let (c0, s0) = (t0.cos(), t0.sin());
            let (c1, s1) = (t1.cos(), t1.sin());
            let n0 = Vector3r::new(c0, 0.0, s0);
            let n1 = Vector3r::new(c1, 0.0, s1);

            let vb0 = mesh.add_vertex(Point3r::new(cx + r * c0, bot_cy, cz + r * s0), n0);
            let vb1 = mesh.add_vertex(Point3r::new(cx + r * c1, bot_cy, cz + r * s1), n1);
            let vt0 = mesh.add_vertex(Point3r::new(cx + r * c0, top_cy, cz + r * s0), n0);
            let vt1 = mesh.add_vertex(Point3r::new(cx + r * c1, top_cy, cz + r * s1), n1);

            mesh.add_face_with_region(vb0, vt0, vt1, region);
            mesh.add_face_with_region(vb0, vt1, vb1, region);
        }
    }

    // ── Bottom hemisphere (φ: π/2 → π, centre at y = bot_cy) ────────────────
    // φ = π/2 → equator ring (y = bot_cy) → shared with cylinder bottom ring
    // φ = π  → south pole (y = bot_cy − r)
    for i in 0..ns {
        let t0 = i as f64 / ns as f64 * TAU;
        let t1 = (i + 1) as f64 / ns as f64 * TAU;
        for j in 0..hs {
            let phi0 = PI / 2.0 + j as f64 / hs as f64 * PI / 2.0;
            let phi1 = PI / 2.0 + (j + 1) as f64 / hs as f64 * PI / 2.0;

            let vat = |theta: f64, phi: f64| -> (Point3r, Vector3r) {
                let sp = phi.sin();
                let cp = phi.cos();
                let ct = theta.cos();
                let st = theta.sin();
                let n = Vector3r::new(sp * ct, cp, sp * st);
                let p = Point3r::new(cx + r * sp * ct, bot_cy + r * cp, cz + r * sp * st);
                (p, n)
            };

            let (p00, n00) = vat(t0, phi0);
            let (p10, n10) = vat(t1, phi0);
            let (p11, n11) = vat(t1, phi1);
            let (p01, n01) = vat(t0, phi1);

            let v00 = mesh.add_vertex(p00, n00);
            let v10 = mesh.add_vertex(p10, n10);
            let v11 = mesh.add_vertex(p11, n11);
            let v01 = mesh.add_vertex(p01, n01);

            if j == hs - 1 {
                // South-pole row: v11 == v01 (φ=π) → triangle
                mesh.add_face_with_region(v00, v10, v01, region);
            } else {
                mesh.add_face_with_region(v00, v10, v11, region);
                mesh.add_face_with_region(v00, v11, v01, region);
            }
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;
    use std::f64::consts::PI;

    #[test]
    fn capsule_is_watertight() {
        let mesh = Capsule::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "capsule must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn capsule_volume_positive_and_approximately_correct() {
        let r = 0.5_f64;
        let h = 2.0_f64;
        let mesh = Capsule {
            radius: r,
            cylinder_height: h,
            segments: 64,
            hemisphere_stacks: 16,
            ..Capsule::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        // V = π r² (h + 4r/3)
        let expected = PI * r * r * (h + 4.0 * r / 3.0);
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.005,
            "volume error {:.4}% should be < 0.5%",
            error * 100.0
        );
    }

    #[test]
    fn capsule_zero_height_is_sphere() {
        // cylinder_height = 0 → capsule becomes a sphere
        let r = 1.0_f64;
        let mesh = Capsule {
            radius: r,
            cylinder_height: 0.0,
            segments: 64,
            hemisphere_stacks: 16,
            ..Capsule::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight);
        let expected = 4.0 / 3.0 * PI * r * r * r;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.005,
            "sphere-degenerate error {:.4}%",
            error * 100.0
        );
    }

    #[test]
    fn capsule_rejects_invalid_params() {
        assert!(Capsule {
            radius: 0.0,
            ..Capsule::default()
        }
        .build()
        .is_err());
        assert!(Capsule {
            radius: -1.0,
            ..Capsule::default()
        }
        .build()
        .is_err());
        assert!(Capsule {
            cylinder_height: -0.1,
            ..Capsule::default()
        }
        .build()
        .is_err());
        assert!(Capsule {
            segments: 2,
            ..Capsule::default()
        }
        .build()
        .is_err());
        assert!(Capsule {
            hemisphere_stacks: 0,
            ..Capsule::default()
        }
        .build()
        .is_err());
    }
}
