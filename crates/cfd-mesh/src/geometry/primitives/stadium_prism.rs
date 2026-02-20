//! Stadium prism (rounded-rectangle cross-section) primitive.

use std::f64::consts::PI;

use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Vector3r};
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Builds a prism with a stadium (rounded-rectangle) cross-section.
///
/// A *stadium* is a rectangle with semicircular ends (also called an
/// *oblong* or *discorectangle*).  This is the dominant cross-section in
/// soft-lithography microfluidic channels fabricated from SU-8 moulds.
///
/// The prism extends along +Y from `base_center`.
///
/// ## Geometry
///
/// The cross-section in the XZ plane consists of:
/// - Two straight edges of length `flat_length = width − 2·corner_radius`
///   parallel to X, at Z = ±corner_radius.
/// - Two semicircles of radius `corner_radius` centred at
///   `(0, ±corner_radius)` (in the Z direction).
///
/// Total width (X extent) = `2 · corner_radius`.
/// Total depth (Z extent) = `flat_length + 2 · corner_radius = width`.
///
/// ## Validation
///
/// - `width > 0`, `height > 0`
/// - `corner_radius ≤ width / 2`
/// - `corner_segments ≥ 1`
/// - `segments ≥ 3` (total angular subdivisions around the profile)
///
/// ## Topology
///
/// χ = 2 (genus 0).  The flat lateral faces and curved lateral faces share
/// edges that are automatically welded by `VertexPool` deduplication.
///
/// ## Output
///
/// - All faces tagged `RegionId(1)`
/// - `signed_volume = (π·r² + 2·r·flat_length) · height`
///   where `r = corner_radius` and `flat_length = width − 2·r`
#[derive(Clone, Debug)]
pub struct StadiumPrism {
    /// Centre of the base face.
    pub base_center: Point3r,
    /// Total width of the stadium cross-section [mm] (= 2 · corner_radius).
    pub width: f64,
    /// Extrusion height along +Y [mm].
    pub height: f64,
    /// Radius of the semicircular ends [mm]. Must be ≤ width / 2.
    pub corner_radius: f64,
    /// Angular segments per *semicircle* end (≥ 1).
    pub corner_segments: usize,
}

impl Default for StadiumPrism {
    fn default() -> Self {
        Self {
            base_center: Point3r::origin(),
            width: 2.0,
            height: 2.0,
            corner_radius: 0.5,
            corner_segments: 8,
        }
    }
}

impl PrimitiveMesh for StadiumPrism {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(sp: &StadiumPrism) -> Result<IndexedMesh, PrimitiveError> {
    if sp.width <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "width must be > 0, got {}", sp.width
        )));
    }
    if sp.height <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "height must be > 0, got {}", sp.height
        )));
    }
    if sp.corner_radius <= 0.0 || sp.corner_radius > sp.width / 2.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "corner_radius must be in (0, width/2] = (0, {}], got {}",
            sp.width / 2.0, sp.corner_radius
        )));
    }
    if sp.corner_segments < 1 {
        return Err(PrimitiveError::InvalidParam(
            "corner_segments must be ≥ 1".into()
        ));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    let r   = sp.corner_radius;
    let h   = sp.height;
    let bx  = sp.base_center.x;
    let by  = sp.base_center.y;
    let bz  = sp.base_center.z;
    let cs  = sp.corner_segments;

    // flat_length along X between the two semicircle centres
    let flat = sp.width - 2.0 * r;

    // Build the closed 2D profile points in the XZ plane (CCW when viewed from
    // +Y — this gives outward normals for lateral faces with +Y convention).
    //
    // Profile traversal (CCW from above):
    //   Right semicircle (centred at +X=flat/2, Z=0), angle -π/2 → +π/2
    //   Left  semicircle (centred at  X=-flat/2, Z=0), angle +π/2 → 3π/2
    //
    // Each semicircle has cs+1 sample points (cs arcs).
    // The two straight edges connecting the semicircles are implicit in the
    // endpoint sharing between the two arcs.

    let mut profile: Vec<[f64; 2]> = Vec::new(); // [x, z]

    // Right semicircle: centre (+flat/2, 0), angles -π/2 → +π/2 (CW = CCW in XZ)
    // Going from bottom-right (+flat/2, -r) around the right end to top-right (+flat/2, +r)
    // in CCW direction from above:
    for i in 0..=cs {
        let angle = -PI / 2.0 + i as f64 / cs as f64 * PI;
        let x = flat / 2.0 + r * angle.cos();
        let z = r * angle.sin();
        profile.push([x, z]);
    }
    // Left semicircle: centre (-flat/2, 0), angles +π/2 → 3π/2
    // Going from top-left (-flat/2, +r) around the left end to bottom-left (-flat/2, -r)
    for i in 0..=cs {
        let angle = PI / 2.0 + i as f64 / cs as f64 * PI;
        let x = -flat / 2.0 + r * angle.cos();
        let z = r * angle.sin();
        profile.push([x, z]);
    }
    // When flat > 0, the seam points (at z=±r, x=±flat/2) differ between the two
    // semicircles (right ends at +flat/2, left starts at -flat/2), so the straight
    // edges are encoded implicitly between consecutive points.
    // When flat = 0, the seam points coincide (both semicircles end/start at x=0).
    // In that case we must remove the duplicates to avoid degenerate zero-length edges
    // and non-manifold topology in the mesh.
    // Deduplicate consecutive near-coincident profile points (tolerance matches VertexPool).
    let profile: Vec<[f64; 2]> = {
        let mut dedup: Vec<[f64; 2]> = Vec::with_capacity(profile.len());
        let tol_sq = 1e-8_f64;
        for p in &profile {
            if let Some(last) = dedup.last() {
                let dx = p[0] - last[0];
                let dz = p[1] - last[1];
                if dx*dx + dz*dz < tol_sq {
                    continue; // skip duplicate
                }
            }
            dedup.push(*p);
        }
        // Also check wrap-around: last point must not equal first
        if dedup.len() > 1 {
            let first = dedup[0];
            let last = *dedup.last().unwrap();
            let dx = last[0] - first[0];
            let dz = last[1] - first[1];
            if dx*dx + dz*dz < tol_sq {
                dedup.pop();
            }
        }
        dedup
    };

    let np = profile.len();

    // Pre-build bottom and top ring vertex ID arrays for shared topology.
    let mut bot_vids: Vec<crate::core::index::VertexId> = Vec::with_capacity(np);
    let mut top_vids: Vec<crate::core::index::VertexId> = Vec::with_capacity(np);
    for i in 0..np {
        let [x0, z0] = profile[i];
        let j = (i + 1) % np;
        let [x1, z1] = profile[j];
        let dx = x1 - x0;
        let dz = z1 - z0;
        let len = (dx * dx + dz * dz).sqrt();
        let (nx, nz) = if len > 1e-14 { (dz / len, -dx / len) } else { (0.0, 0.0) };
        let n = Vector3r::new(nx, 0.0, nz);
        let pb = Point3r::new(bx + x0, by,     bz + z0);
        let pt = Point3r::new(bx + x0, by + h, bz + z0);
        bot_vids.push(mesh.add_vertex(pb, n));
        top_vids.push(mesh.add_vertex(pt, n));
    }

    // Lateral surface
    for i in 0..np {
        let j = (i + 1) % np;
        let vb0 = bot_vids[i];
        let vb1 = bot_vids[j];
        let vt0 = top_vids[i];
        let vt1 = top_vids[j];
        // CCW from outside: bot0->top0->top1, bot0->top1->bot1
        mesh.add_face_with_region(vb0, vt0, vt1, region);
        mesh.add_face_with_region(vb0, vt1, vb1, region);
    }

    // Bottom cap (y = by, normal -Y)
    // Lateral face 2 (vb0, vt1, vb1) creates bottom edge bj->bi.
    // Cap must provide opposite direction bi->bj to complete the manifold edge.
    // Winding (vc, bi, bj) = (vc, v0, v1) gives cross product with -Y component. 
    {
        let n_down = -Vector3r::y();
        let center_bottom = Point3r::new(bx, by, bz);
        let vc = mesh.add_vertex(center_bottom, n_down);
        for i in 0..np {
            let j = (i + 1) % np;
            let v0 = bot_vids[i];
            let v1 = bot_vids[j];
            // Edge at bottom from lateral: bj->bi; cap provides bi->bj.
            // CCW from below (-Y): vc -> bi -> bj
            mesh.add_face_with_region(vc, v0, v1, region);
        }
    }

    // Top cap (y = by + h, normal +Y)
    // Lateral face 1 (vb0, vt0, vt1) creates top edge ti->tj.
    // Cap must provide opposite direction tj->ti to complete the manifold edge.
    // Winding (vc, tj, ti) = (vc, v1, v0) gives cross product with +Y component.
    {
        let n_up = Vector3r::y();
        let center_top = Point3r::new(bx, by + h, bz);
        let vc = mesh.add_vertex(center_top, n_up);
        for i in 0..np {
            let j = (i + 1) % np;
            let v0 = top_vids[i];
            let v1 = top_vids[j];
            // Edge at top from lateral: ti->tj; cap provides tj->ti.
            // CCW from above (+Y): vc -> tj -> ti
            mesh.add_face_with_region(vc, v1, v0, region);
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::storage::edge_store::EdgeStore;
    use crate::watertight::check::check_watertight;
    use std::f64::consts::PI;

    #[test]
    fn stadium_prism_is_watertight() {
        let mesh = StadiumPrism::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "stadium_prism must be watertight: {:?}", report);
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn stadium_prism_volume_positive_and_approximately_correct() {
        let r  = 0.5_f64;
        let w  = 2.0_f64;  // width = 2*r here (semicircle-only, flat=1.0)
        let h  = 2.0_f64;
        let flat = w - 2.0 * r;
        let mesh = StadiumPrism {
            width: w, height: h, corner_radius: r, corner_segments: 16,
            ..StadiumPrism::default()
        }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        // V = (π·r² + 2·r·flat_length) · h
        let expected = (PI * r * r + 2.0 * r * flat) * h;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(error < 0.01, "volume error {:.4}% should be < 1%", error * 100.0);
    }

    #[test]
    fn stadium_prism_full_circle_matches_cylinder() {
        // When flat_length = 0 (corner_radius = width/2), cross-section is a circle.
        let r = 1.0_f64;
        let h = 2.0_f64;
        let mesh = StadiumPrism {
            width: 2.0 * r, height: h, corner_radius: r, corner_segments: 16,
            ..StadiumPrism::default()
        }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight);
        let expected = PI * r * r * h;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(error < 0.02, "cylinder degenerate error {:.4}%", error * 100.0);
    }

    #[test]
    fn stadium_prism_rejects_invalid_params() {
        assert!(StadiumPrism { width: 0.0, ..StadiumPrism::default() }.build().is_err());
        assert!(StadiumPrism { height: 0.0, ..StadiumPrism::default() }.build().is_err());
        assert!(StadiumPrism { corner_radius: 0.0, ..StadiumPrism::default() }.build().is_err());
        assert!(StadiumPrism { corner_radius: 1.5, width: 2.0, ..StadiumPrism::default() }.build().is_err());
        assert!(StadiumPrism { corner_segments: 0, ..StadiumPrism::default() }.build().is_err());
    }
}
