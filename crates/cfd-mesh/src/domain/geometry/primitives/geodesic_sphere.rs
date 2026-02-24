//! Geodesic sphere primitive — icosahedron with subdivided faces.

use super::icosahedron::{icosahedron_vertices, ICOSAHEDRON_FACES};
use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a geodesic sphere by subdividing each icosahedron face `frequency`
/// times and projecting vertices onto a sphere of the given `radius`.
///
/// ## Algorithm
///
/// 1. Build the 12 icosahedron vertices at the given radius (no centre offset yet).
/// 2. For each of the 20 faces, generate a `(frequency+1) × (frequency+2) / 2`
///    triangular grid of barycentric sub-vertices, then project them radially
///    onto the sphere.
/// 3. `VertexPool` spatial-hash deduplication automatically welds the shared
///    edge vertices between adjacent faces.
///
/// ## Why use this instead of [`UvSphere`]?
///
/// - No polar crowding (all triangles have similar solid angles)
/// - Suitable for immersed-boundary particle surface tracking
/// - Better CSG robustness for near-spherical geometry
///
/// ## Topology
///
/// - F = 20·f²,  E = 30·f²,  V = 10·f² + 2  →  χ = 2  (genus 0)
///
/// ## Output
///
/// - All faces tagged `RegionId(1)`
/// - Volume error < 1 % relative to `(4/3)·π·R³` for `frequency ≥ 8`
///
/// [`UvSphere`]: super::sphere::UvSphere
#[derive(Clone, Debug)]
pub struct GeodesicSphere {
    /// Circumsphere radius [mm].
    pub radius: f64,
    /// Centre position.
    pub center: Point3r,
    /// Subdivision frequency (≥ 1). `frequency = 1` yields the base icosahedron
    /// (20 faces); `frequency = 3` gives 180 faces with < 1 % volume error.
    pub frequency: usize,
}

impl Default for GeodesicSphere {
    fn default() -> Self {
        Self {
            radius: 1.0,
            center: Point3r::origin(),
            frequency: 3,
        }
    }
}

impl PrimitiveMesh for GeodesicSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(g: &GeodesicSphere) -> Result<IndexedMesh, PrimitiveError> {
    if g.radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {}",
            g.radius
        )));
    }
    if g.frequency < 1 {
        return Err(PrimitiveError::InvalidParam("frequency must be ≥ 1".into()));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    let r = g.radius;
    let cx = g.center.x;
    let cy = g.center.y;
    let cz = g.center.z;
    let f = g.frequency;

    // Raw icosahedron vertices at radius r (centred at origin; we add center offset below).
    let raw = icosahedron_vertices(r);
    let base_verts: [[f64; 3]; 12] = raw;

    // For each icosahedron face, generate the sub-grid.
    for [ia, ib, ic] in ICOSAHEDRON_FACES {
        let [ax, ay, az] = base_verts[ia];
        let [bx, by, bz] = base_verts[ib];
        let [cx_, cy_, cz_] = base_verts[ic];

        // Barycentric sub-vertices: p(i,j) for i+j <= f.
        // Barycentric coordinates: (u, v, w) = ((f-i-j)/f, i/f, j/f)
        // Position = u*A + v*B + w*C (then projected onto sphere).
        let sub_pt = |si: usize, sj: usize| -> (Point3r, Vector3r) {
            let u = (f - si - sj) as f64 / f as f64;
            let v = si as f64 / f as f64;
            let w = sj as f64 / f as f64;
            let x = u * ax + v * bx + w * cx_;
            let y = u * ay + v * by + w * cy_;
            let z = u * az + v * bz + w * cz_;
            // Project onto sphere: normalize and scale to r.
            let len = (x * x + y * y + z * z).sqrt();
            let (px, py, pz) = if len > 1e-15 {
                (x / len * r, y / len * r, z / len * r)
            } else {
                (0.0, 0.0, r)
            };
            // Outward normal = unit position vector after projection.
            let n = Vector3r::new(px / r, py / r, pz / r);
            let pos = Point3r::new(cx + px, cy + py, cz + pz);
            (pos, n)
        };

        // Iterate over each small triangle in the sub-grid.
        // Two types of triangles per unit rhombus:
        //   "Up"   triangle: (i,j), (i+1,j), (i,j+1)   — i+j+1 <= f
        //   "Down" triangle: (i+1,j), (i+1,j+1), (i,j+1) — i+1+j+1 <= f
        for si in 0..f {
            for sj in 0..f - si {
                // Up triangle
                let (p0, n0) = sub_pt(si, sj);
                let (p1, n1) = sub_pt(si + 1, sj);
                let (p2, n2) = sub_pt(si, sj + 1);
                let v0 = mesh.add_vertex(p0, n0);
                let v1 = mesh.add_vertex(p1, n1);
                let v2 = mesh.add_vertex(p2, n2);
                mesh.add_face_with_region(v0, v1, v2, region);

                // Down triangle (only when there is room)
                if si + 1 + sj + 1 <= f {
                    let (p3, n3) = sub_pt(si + 1, sj + 1);
                    let v3 = mesh.add_vertex(p3, n3);
                    // Down triangle: (i+1,j) → (i+1,j+1) → (i,j+1) — same CCW sense
                    mesh.add_face_with_region(v1, v3, v2, region);
                }
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
    fn geodesic_sphere_f1_is_watertight() {
        // Frequency 1 = icosahedron itself.
        let mesh = GeodesicSphere {
            frequency: 1,
            ..GeodesicSphere::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "geodesic f=1 must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn geodesic_sphere_f3_is_watertight() {
        let mesh = GeodesicSphere::default().build().unwrap(); // f=3
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "geodesic f=3 must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn geodesic_sphere_f8_volume_within_1pct() {
        // f=3 gives only ~5.9% volume error (180 triangles, coarse approximation).
        // Volume error scales as O(1/f²); f=8 (1280 triangles) gives < 1%.
        let r = 1.0_f64;
        let mesh = GeodesicSphere {
            radius: r,
            frequency: 8,
            ..GeodesicSphere::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        let expected = 4.0 / 3.0 * PI * r * r * r;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.01,
            "volume error {:.4}% should be < 1% for f=8",
            error * 100.0
        );
    }

    #[test]
    fn geodesic_sphere_face_count() {
        // F = 20 * f^2
        for f in 1..=4 {
            let mesh = GeodesicSphere {
                frequency: f,
                ..GeodesicSphere::default()
            }
            .build()
            .unwrap();
            assert_eq!(
                mesh.faces.len(),
                20 * f * f,
                "f={}: expected {} faces, got {}",
                f,
                20 * f * f,
                mesh.faces.len()
            );
        }
    }

    #[test]
    fn geodesic_sphere_rejects_invalid_params() {
        assert!(GeodesicSphere {
            radius: 0.0,
            ..GeodesicSphere::default()
        }
        .build()
        .is_err());
        assert!(GeodesicSphere {
            frequency: 0,
            ..GeodesicSphere::default()
        }
        .build()
        .is_err());
    }
}
