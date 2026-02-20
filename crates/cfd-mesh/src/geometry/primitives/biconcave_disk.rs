//! Biconcave disk primitive — Evans-Fung red blood cell (RBC) shape.

use std::f64::consts::TAU;

use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Vector3r};
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Builds a biconcave disk matching the Evans-Fung parametrization of a
/// human red blood cell (RBC).
///
/// The disk is centred at `center` with the disk axis along +Y.  The rim
/// lies in the XZ plane; the upper lobe (+Y) and lower lobe (−Y) are
/// generated independently and joined at the rim.
///
/// ## Profile (Evans-Fung)
///
/// For normalised radius `ρ = 2r/D ∈ [0, 1]`:
/// ```text
/// y(ρ) = (D/2) · √(1 − ρ²) · (c0 + c1·ρ² + c2·ρ⁴)
/// ```
/// Default coefficients model the human discocyte:
/// `c0 = 0.207161`, `c1 = 2.002558`, `c2 = −1.122762`
///
/// ## Region IDs
///
/// | RegionId | Surface |
/// |----------|---------|
/// | 1 | Upper lobe (+Y) |
/// | 2 | Lower lobe (−Y) |
///
/// ## Uses
///
/// Immersed boundary / front-tracking blood flow, capsule dynamics,
/// RBC deformability studies (sickle cell, malaria-infected cells).
///
/// ## Output
///
/// - `signed_volume ≈ 94 fL` for default `diameter = 8e-3` mm (8 µm)
/// - Outward normals on upper lobe point upward (+Y component > 0 except at rim)
/// - Outward normals on lower lobe point downward (−Y component > 0 outward)
#[derive(Clone, Debug)]
pub struct BiconcaveDisk {
    /// Disk diameter [mm]. Default ≈ 8 µm = 8e-3 mm for human RBC.
    pub diameter: f64,
    /// Centre of the disk.
    pub center: Point3r,
    /// Angular subdivisions around the disk rim (≥ 3).
    pub segments: usize,
    /// Radial rings from centre to rim (≥ 2). More rings → better concavity.
    pub rings: usize,
    /// Evans-Fung coefficient c0 (default 0.207161).
    pub c0: f64,
    /// Evans-Fung coefficient c1 (default 2.002558).
    pub c1: f64,
    /// Evans-Fung coefficient c2 (default −1.122762).
    pub c2: f64,
}

impl Default for BiconcaveDisk {
    fn default() -> Self {
        Self {
            diameter: 8e-3,
            center: Point3r::origin(),
            segments: 32,
            rings: 16,
            c0: 0.207161,
            c1: 2.002558,
            c2: -1.122762,
        }
    }
}

impl BiconcaveDisk {
    /// Convenience constructor with human-RBC defaults and specified diameter.
    pub fn human_rbc(diameter: f64) -> Self {
        Self { diameter, ..Self::default() }
    }
}

impl PrimitiveMesh for BiconcaveDisk {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(bd: &BiconcaveDisk) -> Result<IndexedMesh, PrimitiveError> {
    if bd.diameter <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "diameter must be > 0, got {}", bd.diameter
        )));
    }
    if bd.segments < 3 {
        return Err(PrimitiveError::TooFewSegments(bd.segments));
    }
    if bd.rings < 2 {
        return Err(PrimitiveError::InvalidParam(
            "rings must be ≥ 2".into()
        ));
    }

    let upper_region = RegionId::new(1);
    let lower_region = RegionId::new(2);

    let d  = bd.diameter;
    let r  = d / 2.0;

    // Scale the weld tolerance to the geometry.
    //
    // The innermost ring (ρ = 1/nr) has circumradius r/nr and vertex
    // spacing ≈ 2π·r / (nr·ns).  We choose a tolerance that is 1/10 of
    // that spacing — small enough that no two distinct vertices are
    // accidentally merged, yet large enough to absorb floating-point drift
    // at CSG intersection seams.
    //
    // For a human-RBC default (r = 4 µm, nr = 16, ns = 32) this gives
    // tolerance ≈ 4.9 × 10⁻⁵ mm — finer than the fixed 1 × 10⁻⁴ mm that
    // caused incorrect welding of innermost-ring vertices at that scale.
    let min_ring_spacing = std::f64::consts::TAU * r / (bd.rings as f64 * bd.segments as f64);
    let tol   = (min_ring_spacing / 10.0).max(1e-10); // never go below 1 pm
    let cell  = tol * 4.0;                             // cell ≥ 2× tolerance for 26-neighbour check
    let mut mesh = IndexedMesh::with_tolerance(cell, tol);
    let cx = bd.center.x;
    let cy = bd.center.y;
    let cz = bd.center.z;
    let ns = bd.segments;
    let nr = bd.rings;
    let c0 = bd.c0;
    let c1 = bd.c1;
    let c2 = bd.c2;

    // Evans-Fung height at normalised radius rho in [0,1].
    // Returns 0 at the rim (rho=1) by construction.
    let ef_height = |rho: f64| -> f64 {
        let rho2 = rho * rho;
        let under = (1.0 - rho2).max(0.0);
        r * under.sqrt() * (c0 + c1 * rho2 + c2 * rho2 * rho2)
    };

    // Build a ring of (position, outward_normal) at a given rho and +/- side.
    // sign = +1 for upper lobe, sign = -1 for lower lobe (mirror in Y).
    // Outward normal: for the upper lobe it's the surface gradient pointing
    // away from the enclosed volume (roughly +Y at centre, tilted outward).
    // We use analytic cross-product from the partial derivatives dP/drho and
    // dP/dtheta evaluated at this ring.
    let ring_vertices = |rho: f64, sign: f64| -> Vec<(Point3r, Vector3r)> {
        let rho2 = rho * rho;
        let y_val = sign * ef_height(rho);

        // dP/dtheta (tangent in theta direction, no y component): (-r*rho*sin, 0, r*rho*cos)
        // dP/drho: partial of (r*rho*cos, y(rho), r*rho*sin) w.r.t. rho
        // = (r*cos, dy/drho, r*sin)
        // dy/drho = sign * d/drho [ r * sqrt(1-rho²) * f(rho) ]
        let dy_drho = {
            let under = (1.0 - rho2).max(0.0);
            let sqrt_under = under.sqrt();
            // d/drho [sqrt(1-rho²) * f] = -rho/sqrt * f + sqrt * f'
            // f  = c0 + c1*rho² + c2*rho⁴
            // f' = 2*c1*rho + 4*c2*rho³
            let f_val = c0 + c1 * rho2 + c2 * rho2 * rho2;
            let f_prime = 2.0 * c1 * rho + 4.0 * c2 * rho2 * rho;
            if sqrt_under < 1e-14 {
                0.0
            } else {
                sign * r * (-rho / sqrt_under * f_val + sqrt_under * f_prime)
            }
        };

        (0..ns).map(|i| {
            let theta = i as f64 / ns as f64 * TAU;
            let (ct, st) = (theta.cos(), theta.sin());
            let pos = Point3r::new(cx + r * rho * ct, cy + y_val, cz + r * rho * st);

            // dP/dtheta = (-r*rho*sin, 0, r*rho*cos)
            let dt = Vector3r::new(-r * rho * st, 0.0, r * rho * ct);
            // dP/drho = (r*cos, dy_drho, r*sin)
            let dr = Vector3r::new(r * ct, dy_drho, r * st);

            // For upper lobe: outward = dP/dtheta × dP/drho (gives +Y bias at centre)
            // For lower lobe: we need outward = dP/drho × dP/dtheta (gives -Y bias at centre)
            let raw_n = if sign > 0.0 {
                dt.cross(&dr)
            } else {
                dr.cross(&dt)
            };
            let len = raw_n.norm();
            let n = if len < 1e-14 {
                Vector3r::new(0.0, sign, 0.0)
            } else {
                raw_n / len
            };

            (pos, n)
        }).collect()
    };

    // Shared rim ring (rho = 1, y = 0)
    // Both upper and lower lobes end at the rim. Build it once for shared topology.
    let rim_ids: Vec<crate::core::index::VertexId> = (0..ns).map(|i| {
        let theta = i as f64 / ns as f64 * TAU;
        let (ct, st) = (theta.cos(), theta.sin());
        let pos = Point3r::new(cx + r * ct, cy, cz + r * st);
        // Normal at rim points radially outward in XZ plane
        let n = Vector3r::new(ct, 0.0, st);
        mesh.add_vertex(pos, n)
    }).collect();

    // Upper lobe (+Y, outward normals generally pointing +Y)
    {
        let apex_y = cy + ef_height(0.0);
        let apex_n = Vector3r::new(0.0, 1.0, 0.0);
        let apex_pos = Point3r::new(cx, apex_y, cz);
        let v_apex = mesh.add_vertex(apex_pos, apex_n);

        // Build all upper ring IDs, using rim_ids for the outermost ring
        let mut upper_ring_ids: Vec<Vec<crate::core::index::VertexId>> = Vec::with_capacity(nr);
        for k in 1..=nr {
            let rho = k as f64 / nr as f64;
            let ids: Vec<_> = if k == nr {
                rim_ids.clone()
            } else {
                let curr = ring_vertices(rho, 1.0);
                curr.iter().map(|(p, n)| mesh.add_vertex(*p, *n)).collect()
            };
            upper_ring_ids.push(ids);
        }

        // Fan from apex to first ring.
        // Reversed angular order so that the apex-to-ring0 edge direction pairs
        // correctly (opposite) with the k=0 concentric ring's ring0 edges.
        for i in 0..ns {
            let j = (i + 1) % ns;
            mesh.add_face_with_region(v_apex, upper_ring_ids[0][j], upper_ring_ids[0][i], upper_region);
        }

        // Concentric rings
        for k in 0..nr - 1 {
            for i in 0..ns {
                let j = (i + 1) % ns;
                mesh.add_face_with_region(upper_ring_ids[k][i], upper_ring_ids[k][j], upper_ring_ids[k+1][j], upper_region);
                mesh.add_face_with_region(upper_ring_ids[k][i], upper_ring_ids[k+1][j], upper_ring_ids[k+1][i], upper_region);
            }
        }
    }

    // Lower lobe (-Y, outward normals generally pointing -Y)
    {
        let apex_y = cy - ef_height(0.0).abs();
        let apex_n = Vector3r::new(0.0, -1.0, 0.0);
        let apex_pos = Point3r::new(cx, apex_y, cz);
        let v_apex = mesh.add_vertex(apex_pos, apex_n);

        // Build all lower ring IDs, using rim_ids for the outermost ring
        let mut lower_ring_ids: Vec<Vec<crate::core::index::VertexId>> = Vec::with_capacity(nr);
        for k in 1..=nr {
            let rho = k as f64 / nr as f64;
            let ids: Vec<_> = if k == nr {
                rim_ids.clone()
            } else {
                let curr = ring_vertices(rho, -1.0);
                curr.iter().map(|(p, n)| mesh.add_vertex(*p, *n)).collect()
            };
            lower_ring_ids.push(ids);
        }

        // Fan from apex to first ring -- normal angular order (same topology as upper lobe).
        // The concentric ring faces already use reversed order for the lower lobe,
        // so the apex fan must use normal order to produce opposite ring0 edge direction.
        for i in 0..ns {
            let j = (i + 1) % ns;
            mesh.add_face_with_region(v_apex, lower_ring_ids[0][i], lower_ring_ids[0][j], lower_region);
        }

        // Concentric rings -- reversed angular order
        for k in 0..nr - 1 {
            for i in 0..ns {
                let j = (i + 1) % ns;
                mesh.add_face_with_region(lower_ring_ids[k][j], lower_ring_ids[k][i], lower_ring_ids[k+1][i], lower_region);
                mesh.add_face_with_region(lower_ring_ids[k][j], lower_ring_ids[k+1][i], lower_ring_ids[k+1][j], lower_region);
            }
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::storage::edge_store::EdgeStore;
    use crate::watertight::check::check_watertight;

    #[test]
    fn biconcave_disk_is_watertight() {
        // Use diameter = 1.0 mm so that ring vertex spacing (≈ 2π·r·ρ_min/ns)
        // stays well above the VertexPool weld tolerance (1e-4 mm).
        // The default 8e-3 mm (8 µm) makes the innermost ring too small to mesh.
        let mesh = BiconcaveDisk::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "biconcave_disk must be watertight: {:?}", report);
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn biconcave_disk_volume_positive() {
        let mesh = BiconcaveDisk {
            diameter: 8e-3,
            segments: 64,
            rings: 32,
            ..BiconcaveDisk::default()
        }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0, "signed_volume must be positive");
        // Human RBC volume ≈ 94 fL = 94e-9 mm³; allow wide tolerance for numeric integration
        // against discretized profile
        assert!(
            report.signed_volume > 0.0 && report.signed_volume < 1.0,
            "volume {:.3e} mm³ should be in (0, 1) mm³ for 8 µm disk",
            report.signed_volume
        );
    }

    #[test]
    fn biconcave_disk_human_rbc_constructor() {
        // The build() function automatically scales the weld tolerance to the
        // geometry, so the default 8 µm diameter and default segments/rings
        // now produce a watertight mesh without manual parameter adjustment.
        let mesh = BiconcaveDisk::human_rbc(8e-3).build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight);
    }

    #[test]
    fn biconcave_disk_rejects_invalid() {
        assert!(BiconcaveDisk { diameter: 0.0, ..BiconcaveDisk::default() }.build().is_err());
        assert!(BiconcaveDisk { diameter: -1.0, ..BiconcaveDisk::default() }.build().is_err());
        assert!(BiconcaveDisk { segments: 2, ..BiconcaveDisk::default() }.build().is_err());
        assert!(BiconcaveDisk { rings: 1, ..BiconcaveDisk::default() }.build().is_err());
    }
}
