//! Right pyramid primitive (n-gon base).

use std::f64::consts::TAU;

use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Vector3r};
use crate::geometry::normal::triangle_normal;
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Builds a right pyramid with a regular n-gon base.
///
/// The pyramid is aligned with the +Y axis: the base is centred at
/// `base_center` and the apex is at `base_center + (0, height, 0)`.
///
/// ## Geometry
///
/// - Base: regular polygon with `sides` vertices inscribed in circle of `base_radius`.
/// - Apex: single vertex at height `height` above the base centre.
/// - Lateral faces: `sides` equilateral (for `sides = 4`) or isosceles triangles.
///
/// `sides = 4` gives a square pyramid (half an octahedron when `height = √2 r`).
///
/// ## Topology
///
/// - V = sides + 1, E = 2·sides, F = sides + 1
/// - `V − E + F = 2`  (χ = 2, genus 0)
///
/// ## Output
///
/// - `RegionId(1)` on all faces
/// - `signed_volume = (1/3)·base_area·height = (sides·r²·sin(2π/sides)/6)·height`
#[derive(Clone, Debug)]
pub struct Pyramid {
    /// Centre of the base polygon.
    pub base_center: Point3r,
    /// Circumradius of the base polygon [mm].
    pub base_radius: f64,
    /// Pyramid height [mm].
    pub height: f64,
    /// Number of sides of the base polygon (≥ 3).
    pub sides: usize,
}

impl Default for Pyramid {
    fn default() -> Self {
        Self {
            base_center: Point3r::origin(),
            base_radius: 1.0,
            height: 1.5,
            sides: 4,
        }
    }
}

impl PrimitiveMesh for Pyramid {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(p: &Pyramid) -> Result<IndexedMesh, PrimitiveError> {
    if p.base_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "base_radius must be > 0, got {}", p.base_radius
        )));
    }
    if p.height <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "height must be > 0, got {}", p.height
        )));
    }
    if p.sides < 3 {
        return Err(PrimitiveError::TooFewSegments(p.sides));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    let r  = p.base_radius;
    let bx = p.base_center.x;
    let by = p.base_center.y;
    let bz = p.base_center.z;
    let ns = p.sides;

    // Apex
    let apex = Point3r::new(bx, by + p.height, bz);
    // Base polygon vertices (CCW from above = CW from below = correct for outward normals)
    let base: Vec<Point3r> = (0..ns).map(|i| {
        let angle = i as f64 / ns as f64 * TAU;
        Point3r::new(bx + r * angle.cos(), by, bz + r * angle.sin())
    }).collect();

    // ── Lateral faces ─────────────────────────────────────────────────────────
    // Triangle: apex, base[i], base[i+1] — normal points outward.
    // Verify CCW winding by checking sign with triangle_normal.
    for i in 0..ns {
        let j = (i + 1) % ns;
        let n = triangle_normal(&apex, &base[i], &base[j])
            .unwrap_or(Vector3r::zeros());
        // Ensure normal points outward (away from axis, positive horizontal component).
        let va = mesh.add_vertex(apex,    n);
        let v0 = mesh.add_vertex(base[i], n);
        let v1 = mesh.add_vertex(base[j], n);
        mesh.add_face_with_region(va, v0, v1, region);
    }

    // ── Base cap (normal = −Y) ─────────────────────────────────────────────────
    // Fan from base centre, CCW from below.
    {
        let n_down = -Vector3r::y();
        let vc = mesh.add_vertex(Point3r::new(bx, by, bz), n_down);
        for i in 0..ns {
            let j = (i + 1) % ns;
            let v0 = mesh.add_vertex(base[i], n_down);
            let v1 = mesh.add_vertex(base[j], n_down);
            // CCW from below: vc → vj → vi
            mesh.add_face_with_region(vc, v1, v0, region);
        }
    }

    // The fan-triangulated lateral faces and base cap produce consistent
    // inward winding due to the CCW base order and apex placement.
    // Flip all faces so outward normals point away from the interior.
    mesh.flip_faces();

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::storage::edge_store::EdgeStore;
    use crate::watertight::check::check_watertight;
    use std::f64::consts::PI;

    #[test]
    fn pyramid_is_watertight() {
        let mesh = Pyramid::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "pyramid must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn pyramid_volume_positive_and_correct() {
        let r  = 1.0_f64;
        let h  = 2.0_f64;
        let ns = 4_usize;
        let mesh = Pyramid { base_radius: r, height: h, sides: ns, ..Pyramid::default() }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0, "pyramid signed_volume must be positive");
        // V = (1/3) * base_area * h
        // base_area = n * r^2 * sin(2π/n) / 2
        let base_area = ns as f64 * r * r * (TAU / ns as f64).sin() / 2.0;
        let expected = base_area * h / 3.0;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(error < 1e-9, "volume error {:.2e} should be machine eps", error);
    }

    #[test]
    fn pyramid_rejects_invalid_params() {
        assert!(Pyramid { base_radius: 0.0, ..Pyramid::default() }.build().is_err());
        assert!(Pyramid { height: 0.0, ..Pyramid::default() }.build().is_err());
        assert!(Pyramid { sides: 2, ..Pyramid::default() }.build().is_err());
    }
}
