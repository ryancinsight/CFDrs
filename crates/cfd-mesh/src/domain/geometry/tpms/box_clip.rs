//! AABB-clipped TPMS extraction — fills a rectangular volume with a TPMS mesh.
//!
//! This module mirrors [`build_tpms_sphere`](super::build_tpms_sphere) but clips
//! the marching-cubes extraction to an axis-aligned bounding box instead of a
//! sphere.  This is the geometry kernel used when filling a shell cuboid cavity
//! with a TPMS lattice network.
//!
//! ## Theorem — Grid Convergence
//!
//! The same `O(h)` Hausdorff convergence rate (Lorensen & Cline 1987) applies:
//! as `resolution → ∞`, the extracted surface converges to the true `{F = c}`
//! level-set inside the box.

use crate::domain::geometry::primitives::PrimitiveError;
use crate::domain::geometry::tpms::marching_cubes;
use crate::domain::geometry::tpms::Tpms;
use crate::domain::geometry::tpms::Vector3r;
use crate::domain::mesh::IndexedMesh;

use std::collections::HashMap;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Point3r;

// ── Parameters ────────────────────────────────────────────────────────────────

/// Parameters for AABB-clipped TPMS extraction.
///
/// The TPMS surface is extracted inside the box
/// `[x_min, x_max] × [y_min, y_max] × [z_min, z_max]`.
#[derive(Clone, Debug)]
pub struct TpmsBoxParams {
    /// AABB bounds: `[x_min, y_min, z_min, x_max, y_max, z_max]` in mm.
    pub bounds: [f64; 6],
    /// TPMS unit-cell period [mm].  `k = 2π / period`.
    pub period: f64,
    /// Voxels per axis.  Higher → denser, more accurate.
    pub resolution: usize,
    /// Level-set iso-value.  `0.0` = exact minimal surface mid-sheet.
    pub iso_value: f64,
}

impl TpmsBoxParams {
    /// Validate all parameters.
    ///
    /// # Errors
    ///
    /// Returns [`PrimitiveError::InvalidParam`] for degenerate or non-finite
    /// bounds, non-positive period, or resolution < 4.
    pub fn validate(&self) -> Result<(), PrimitiveError> {
        let [x0, y0, z0, x1, y1, z1] = self.bounds;
        if !x0.is_finite() || !y0.is_finite() || !z0.is_finite()
            || !x1.is_finite() || !y1.is_finite() || !z1.is_finite()
        {
            return Err(PrimitiveError::InvalidParam(
                "all AABB bounds must be finite".to_string(),
            ));
        }
        if x1 <= x0 || y1 <= y0 || z1 <= z0 {
            return Err(PrimitiveError::InvalidParam(format!(
                "AABB must have positive extents: ({x0},{y0},{z0})→({x1},{y1},{z1})"
            )));
        }
        if self.period <= 0.0 {
            return Err(PrimitiveError::InvalidParam(format!(
                "period must be > 0, got {}",
                self.period
            )));
        }
        if self.resolution < 4 {
            return Err(PrimitiveError::InvalidParam(format!(
                "resolution must be >= 4, got {}",
                self.resolution
            )));
        }
        Ok(())
    }
}

// ── Builder ───────────────────────────────────────────────────────────────────

/// Extract a TPMS mesh clipped to an axis-aligned bounding box.
///
/// This is the rectangular counterpart of [`build_tpms_sphere`](super::build_tpms_sphere).
/// The marching-cubes grid spans `[x_min, x_max] × [y_min, y_max] × [z_min, z_max]`
/// with `resolution` voxels per axis.  All extracted triangles lie inside the box
/// (no centroid distance culling).
///
/// # Errors
///
/// Returns [`PrimitiveError::InvalidParam`] on parameter validation failure.
pub fn build_tpms_box<S: Tpms>(
    surface: &S,
    params: &TpmsBoxParams,
) -> Result<IndexedMesh, PrimitiveError> {
    params.validate()?;

    let [x0, y0, z0, x1, y1, z1] = params.bounds;
    let k = std::f64::consts::TAU / params.period;
    let n = params.resolution;
    let iso = params.iso_value;

    let dx = (x1 - x0) / n as f64;
    let dy = (y1 - y0) / n as f64;
    let dz = (z1 - z0) / n as f64;
    // Pad by 1 voxel on each side so the marching cubes bounds enclose the box
    let gs = n + 3;

    let cx = (x0 + x1) * 0.5;
    let cy = (y0 + y1) * 0.5;
    let cz = (z0 + z1) * 0.5;
    let hx = (x1 - x0) * 0.5;
    let hy = (y1 - y0) * 0.5;
    let hz = (z1 - z0) * 0.5;

    // Pre-sample field on padded grid.
    let mut field = vec![0.0_f64; gs * gs * gs];
    let idx = |ix: usize, iy: usize, iz: usize| iz * gs * gs + iy * gs + ix;
    for iz in 0..gs {
        for iy in 0..gs {
            for ix in 0..gs {
                let wx = x0 + (ix as isize - 1) as f64 * dx;
                let wy = y0 + (iy as isize - 1) as f64 * dy;
                let wz = z0 + (iz as isize - 1) as f64 * dz;
                
                let tpms_val = surface.field(wx, wy, wz, k) - iso;
                
                // Box SDF
                let qx = (wx - cx).abs() - hx;
                let qy = (wy - cy).abs() - hy;
                let qz = (wz - cz).abs() - hz;
                let box_sdf = qx.max(0.0).hypot(qy.max(0.0)).hypot(qz.max(0.0))
                            + qx.max(qy).max(qz).min(0.0);
                
                // Solid intersection: max(tpms, box_sdf).
                // Negative = inside fluid, Positive = outside (wall or blocked).
                field[idx(ix, iy, iz)] = tpms_val.max(box_sdf);
            }
        }
    }

    let mut mesh = IndexedMesh::new();
    let mut cache: HashMap<(usize, usize, usize, usize), VertexId> = HashMap::new();

    for iz in 0..(gs - 1) {
        for iy in 0..(gs - 1) {
            for ix in 0..(gs - 1) {
                // Corner field values and sign configuration.
                let mut cube_vals = [0.0_f64; 8];
                let mut cube_cfg: usize = 0;
                for (ci, &(cdx, cdy, cdz)) in marching_cubes::CORNERS.iter().enumerate() {
                    let v = field[idx(
                        ix + cdx as usize,
                        iy + cdy as usize,
                        iz + cdz as usize,
                    )];
                    cube_vals[ci] = v;
                    if v < 0.0 {
                        cube_cfg |= 1 << ci;
                    }
                }

                let emask = marching_cubes::EDGE_TABLE[cube_cfg];
                if emask == 0 {
                    continue;
                }

                // Resolve or create vertex for each intersected edge.
                let mut edge_vids: [Option<VertexId>; 12] = [None; 12];
                for (ei, &[ca, cb]) in marching_cubes::EDGES.iter().enumerate() {
                    if emask & (1 << ei) == 0 {
                        continue;
                    }
                    let vid = *cache.entry((ix, iy, iz, ei)).or_insert_with(|| {
                        let (ax, ay, az) = (
                            ix + marching_cubes::CORNERS[ca].0 as usize,
                            iy + marching_cubes::CORNERS[ca].1 as usize,
                            iz + marching_cubes::CORNERS[ca].2 as usize,
                        );
                        let (bx, by, bz) = (
                            ix + marching_cubes::CORNERS[cb].0 as usize,
                            iy + marching_cubes::CORNERS[cb].1 as usize,
                            iz + marching_cubes::CORNERS[cb].2 as usize,
                        );
                        let va = cube_vals[ca];
                        let vb = cube_vals[cb];
                        let t = if (vb - va).abs() > 1e-15 {
                            (-va / (vb - va)).clamp(0.0, 1.0)
                        } else {
                            0.5
                        };
                        
                        let wx = x0 + (ax as f64 * (1.0 - t) + bx as f64 * t - 1.0) * dx;
                        let wy = y0 + (ay as f64 * (1.0 - t) + by as f64 * t - 1.0) * dy;
                        let wz = z0 + (az as f64 * (1.0 - t) + bz as f64 * t - 1.0) * dz;
                        
                        // Because the boundary is defined by the Box SDF max intersection,
                        // surface normals at the exact box boundary should point inwards (from wall).
                        // If it's a TPMS body surface, use TPMS gradient.
                        // We use a simple numeric SDF check:
                        let qx = (wx - cx).abs() - hx;
                        let qy = (wy - cy).abs() - hy;
                        let qz = (wz - cz).abs() - hz;
                        let box_sdf = qx.max(0.0).hypot(qy.max(0.0)).hypot(qz.max(0.0))
                                    + qx.max(qy).max(qz).min(0.0);
                        
                        let normal = if box_sdf.abs() < 1e-5 {
                            // On box boundary, normal points outward of the fluid (into the box wall)
                            // We construct the normal by seeing which face we're on
                            let mut nx = 0.0;
                            let mut ny = 0.0;
                            let mut nz = 0.0;
                            if (wx - x0).abs() < 1e-5 { nx = -1.0; }
                            else if (wx - x1).abs() < 1e-5 { nx = 1.0; }
                            if (wy - y0).abs() < 1e-5 { ny = -1.0; }
                            else if (wy - y1).abs() < 1e-5 { ny = 1.0; }
                            if (wz - z0).abs() < 1e-5 { nz = -1.0; }
                            else if (wz - z1).abs() < 1e-5 { nz = 1.0; }
                            let v = Vector3r::new(nx, ny, nz);
                            if v.norm_squared() > 1e-6 { v.normalize() } else { surface.gradient(wx, wy, wz, k) }
                        } else {
                            surface.gradient(wx, wy, wz, k)
                        };
                        mesh.add_vertex(Point3r::new(wx, wy, wz), normal)
                    });
                    edge_vids[ei] = Some(vid);
                }

                // Emit triangles
                let tri_row = &marching_cubes::TRI_TABLE[cube_cfg];
                let mut ti = 0;
                while ti + 2 < 16 && tri_row[ti] >= 0 {
                    let e0 = tri_row[ti] as usize;
                    let e1 = tri_row[ti + 1] as usize;
                    let e2 = tri_row[ti + 2] as usize;
                    if let (Some(v0), Some(v1), Some(v2)) =
                        (edge_vids[e0], edge_vids[e1], edge_vids[e2])
                    {
                        mesh.add_face(v0, v1, v2);
                    }
                    ti += 3;
                }
            }
        }
    }

    Ok(mesh)
}

// ── Graded builder ────────────────────────────────────────────────────────────

/// Validate AABB bounds and resolution (shared by `build_tpms_box` and the
/// graded variant).
fn validate_box_bounds(bounds: &[f64; 6], resolution: usize) -> Result<(), PrimitiveError> {
    let [x0, y0, z0, x1, y1, z1] = *bounds;
    if !x0.is_finite() || !y0.is_finite() || !z0.is_finite()
        || !x1.is_finite() || !y1.is_finite() || !z1.is_finite()
    {
        return Err(PrimitiveError::InvalidParam(
            "all AABB bounds must be finite".to_string(),
        ));
    }
    if x1 <= x0 || y1 <= y0 || z1 <= z0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "AABB must have positive extents: ({x0},{y0},{z0})→({x1},{y1},{z1})"
        )));
    }
    if resolution < 4 {
        return Err(PrimitiveError::InvalidParam(format!(
            "resolution must be >= 4, got {resolution}"
        )));
    }
    Ok(())
}

/// Extract a TPMS mesh with **spatially-varying period** inside an AABB.
///
/// This is the adaptive counterpart of [`build_tpms_box`].  Instead of a
/// single global period, the caller provides a closure `period_fn(x, y, z)`
/// that returns the local period [mm] at each world coordinate.
///
/// This enables graded pore structures for size-based cell separation:
/// fine period at the periphery (small pores → block RBCs) and coarse
/// period at the center (large pores → pass WBCs/CTCs).
///
/// # Arguments
///
/// * `surface` — the TPMS implicit surface to extract.
/// * `bounds` — AABB `[x_min, y_min, z_min, x_max, y_max, z_max]`.
/// * `resolution` — voxels per axis (≥ 4).
/// * `iso_value` — level-set threshold to subtract from the field.
/// * `period_fn` — closure `(x, y, z) → period_mm` for spatially-varying
///   period.  Must return a positive, finite value for all inputs.
///
/// # Errors
///
/// Returns [`PrimitiveError::InvalidParam`] on degenerate bounds or low
/// resolution.
pub fn build_tpms_box_graded<S: Tpms>(
    surface: &S,
    bounds: [f64; 6],
    resolution: usize,
    iso_value: f64,
    period_fn: impl Fn(f64, f64, f64) -> f64,
) -> Result<IndexedMesh, PrimitiveError> {
    validate_box_bounds(&bounds, resolution)?;

    let [x0, y0, z0, x1, y1, z1] = bounds;
    let n = resolution;
    let iso = iso_value;

    let dx = (x1 - x0) / n as f64;
    let dy = (y1 - y0) / n as f64;
    let dz = (z1 - z0) / n as f64;
    let gs = n + 3;

    let cx = (x0 + x1) * 0.5;
    let cy = (y0 + y1) * 0.5;
    let cz = (z0 + z1) * 0.5;
    let hx = (x1 - x0) * 0.5;
    let hy = (y1 - y0) * 0.5;
    let hz = (z1 - z0) * 0.5;

    // Pre-sample field on (n+1)³ grid with spatially-varying k.
    let mut field = vec![0.0_f64; gs * gs * gs];
    let idx = |ix: usize, iy: usize, iz: usize| iz * gs * gs + iy * gs + ix;
    for iz in 0..gs {
        for iy in 0..gs {
            for ix in 0..gs {
                let wx = x0 + (ix as isize - 1) as f64 * dx;
                let wy = y0 + (iy as isize - 1) as f64 * dy;
                let wz = z0 + (iz as isize - 1) as f64 * dz;
                let local_period = period_fn(wx, wy, wz).max(1e-12);
                let local_k = std::f64::consts::TAU / local_period;
                let tpms_val = surface.field(wx, wy, wz, local_k) - iso;

                let qx = (wx - cx).abs() - hx;
                let qy = (wy - cy).abs() - hy;
                let qz = (wz - cz).abs() - hz;
                let box_sdf = qx.max(0.0).hypot(qy.max(0.0)).hypot(qz.max(0.0))
                            + qx.max(qy).max(qz).min(0.0);
                
                field[idx(ix, iy, iz)] = tpms_val.max(box_sdf);
            }
        }
    }

    let mut mesh = IndexedMesh::new();
    let mut cache: HashMap<(usize, usize, usize, usize), VertexId> = HashMap::new();

    for iz in 0..(gs - 1) {
        for iy in 0..(gs - 1) {
            for ix in 0..(gs - 1) {
                let mut cube_vals = [0.0_f64; 8];
                let mut cube_cfg: usize = 0;
                for (ci, &(cdx, cdy, cdz)) in marching_cubes::CORNERS.iter().enumerate() {
                    let v = field[idx(
                        ix + cdx as usize,
                        iy + cdy as usize,
                        iz + cdz as usize,
                    )];
                    cube_vals[ci] = v;
                    if v < 0.0 {
                        cube_cfg |= 1 << ci;
                    }
                }

                let emask = marching_cubes::EDGE_TABLE[cube_cfg];
                if emask == 0 {
                    continue;
                }

                let mut edge_vids: [Option<VertexId>; 12] = [None; 12];
                for (ei, &[ca, cb]) in marching_cubes::EDGES.iter().enumerate() {
                    if emask & (1 << ei) == 0 {
                        continue;
                    }
                    let vid = *cache.entry((ix, iy, iz, ei)).or_insert_with(|| {
                        let (ax, ay, az) = (
                            ix + marching_cubes::CORNERS[ca].0 as usize,
                            iy + marching_cubes::CORNERS[ca].1 as usize,
                            iz + marching_cubes::CORNERS[ca].2 as usize,
                        );
                        let (bx, by, bz) = (
                            ix + marching_cubes::CORNERS[cb].0 as usize,
                            iy + marching_cubes::CORNERS[cb].1 as usize,
                            iz + marching_cubes::CORNERS[cb].2 as usize,
                        );
                        let va = cube_vals[ca];
                        let vb = cube_vals[cb];
                        let t = if (vb - va).abs() > 1e-15 {
                            (-va / (vb - va)).clamp(0.0, 1.0)
                        } else {
                            0.5
                        };
                        let wx = x0 + (ax as f64 * (1.0 - t) + bx as f64 * t - 1.0) * dx;
                        let wy = y0 + (ay as f64 * (1.0 - t) + by as f64 * t - 1.0) * dy;
                        let wz = z0 + (az as f64 * (1.0 - t) + bz as f64 * t - 1.0) * dz;
                        
                        let qx = (wx - cx).abs() - hx;
                        let qy = (wy - cy).abs() - hy;
                        let qz = (wz - cz).abs() - hz;
                        let box_sdf = qx.max(0.0).hypot(qy.max(0.0)).hypot(qz.max(0.0))
                                    + qx.max(qy).max(qz).min(0.0);
                        
                        let normal = if box_sdf.abs() < 1e-5 {
                            let mut nx = 0.0;
                            let mut ny = 0.0;
                            let mut nz = 0.0;
                            if (wx - x0).abs() < 1e-5 { nx = -1.0; }
                            else if (wx - x1).abs() < 1e-5 { nx = 1.0; }
                            if (wy - y0).abs() < 1e-5 { ny = -1.0; }
                            else if (wy - y1).abs() < 1e-5 { ny = 1.0; }
                            if (wz - z0).abs() < 1e-5 { nz = -1.0; }
                            else if (wz - z1).abs() < 1e-5 { nz = 1.0; }
                            let v = Vector3r::new(nx, ny, nz);
                            if v.norm_squared() > 1e-6 { v.normalize() } else { 
                                let local_period = period_fn(wx, wy, wz).max(1e-12);
                                let local_k = std::f64::consts::TAU / local_period;
                                surface.gradient(wx, wy, wz, local_k) 
                            }
                        } else {
                            // Gradient uses local k at the interpolated position.
                            let local_period = period_fn(wx, wy, wz).max(1e-12);
                            let local_k = std::f64::consts::TAU / local_period;
                            surface.gradient(wx, wy, wz, local_k)
                        };
                        mesh.add_vertex(Point3r::new(wx, wy, wz), normal)
                    });
                    edge_vids[ei] = Some(vid);
                }

                let tri_row = &marching_cubes::TRI_TABLE[cube_cfg];
                let mut ti = 0;
                while ti + 2 < 16 && tri_row[ti] >= 0 {
                    let e0 = tri_row[ti] as usize;
                    let e1 = tri_row[ti + 1] as usize;
                    let e2 = tri_row[ti + 2] as usize;
                    if let (Some(v0), Some(v1), Some(v2)) =
                        (edge_vids[e0], edge_vids[e1], edge_vids[e2])
                    {
                        mesh.add_face(v0, v1, v2);
                    }
                    ti += 3;
                }
            }
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::geometry::tpms::Gyroid;

    #[test]
    fn box_clip_produces_nonempty_mesh() {
        let params = TpmsBoxParams {
            bounds: [-5.0, -5.0, -5.0, 5.0, 5.0, 5.0],
            period: 2.5,
            resolution: 16,
            iso_value: 0.0,
        };
        let mesh = build_tpms_box(&Gyroid, &params).expect("should succeed");
        assert!(
            mesh.face_count() > 0,
            "gyroid-in-box must produce at least one face"
        );
        assert!(
            mesh.vertex_count() > 0,
            "gyroid-in-box must produce at least one vertex"
        );
    }

    #[test]
    fn box_clip_validates_degenerate_bounds() {
        let params = TpmsBoxParams {
            bounds: [5.0, 0.0, 0.0, 5.0, 10.0, 10.0], // x_max == x_min
            period: 2.5,
            resolution: 16,
            iso_value: 0.0,
        };
        assert!(build_tpms_box(&Gyroid, &params).is_err());
    }

    #[test]
    fn box_clip_validates_low_resolution() {
        let params = TpmsBoxParams {
            bounds: [0.0, 0.0, 0.0, 10.0, 10.0, 10.0],
            period: 2.5,
            resolution: 2,
            iso_value: 0.0,
        };
        assert!(build_tpms_box(&Gyroid, &params).is_err());
    }

    #[test]
    fn box_clip_all_vertices_within_bounds() {
        let params = TpmsBoxParams {
            bounds: [-3.0, -2.0, -1.0, 4.0, 5.0, 6.0],
            period: 3.0,
            resolution: 16,
            iso_value: 0.0,
        };
        let mesh = build_tpms_box(&Gyroid, &params).expect("should succeed");
        let eps = params.period / params.resolution as f64; // one voxel tolerance
        for vid in 0..mesh.vertex_count() {
            let p = mesh.vertices.position(VertexId(vid as u32));
            assert!(
                p.x >= -3.0 - eps && p.x <= 4.0 + eps,
                "vertex x={} out of bounds",
                p.x,
            );
            assert!(
                p.y >= -2.0 - eps && p.y <= 5.0 + eps,
                "vertex y={} out of bounds",
                p.y,
            );
            assert!(
                p.z >= -1.0 - eps && p.z <= 6.0 + eps,
                "vertex z={} out of bounds",
                p.z,
            );
        }
    }

    // ── Graded builder tests ──────────────────────────────────────────────

    #[test]
    fn graded_uniform_matches_box_clip() {
        // A graded builder with constant period should produce the same mesh
        // topology as build_tpms_box (same face count ± small tolerance from
        // floating point differences).
        let bounds = [-5.0, -5.0, -5.0, 5.0, 5.0, 5.0];
        let period = 2.5;
        let params = TpmsBoxParams {
            bounds,
            period,
            resolution: 16,
            iso_value: 0.0,
        };
        let uniform = build_tpms_box(&Gyroid, &params).unwrap();
        let graded = build_tpms_box_graded(
            &Gyroid,
            bounds,
            16,
            0.0,
            |_x, _y, _z| period,
        )
        .unwrap();
        assert_eq!(
            uniform.face_count(),
            graded.face_count(),
            "constant-period graded must produce same face count as uniform"
        );
    }

    #[test]
    fn graded_mesh_nonempty() {
        // A graded mesh with period varying from 1.5 (walls) to 5.0 (center)
        // should produce a non-empty mesh.
        let bounds = [0.0, 0.0, 0.0, 10.0, 10.0, 5.0];
        let mesh = build_tpms_box_graded(
            &Gyroid,
            bounds,
            20,
            0.0,
            |_x, y, _z| {
                // Y ranges [0, 10]: center at 5.0
                let y_frac = y / 10.0;
                let wall_dist = (2.0 * (y_frac - 0.5)).abs();
                5.0 * (1.0 - wall_dist) + 1.5 * wall_dist
            },
        )
        .unwrap();
        assert!(
            mesh.face_count() > 0,
            "graded gyroid must produce faces"
        );
    }

    #[test]
    fn graded_rejects_degenerate_bounds() {
        assert!(build_tpms_box_graded(
            &Gyroid,
            [0.0, 0.0, 0.0, 0.0, 10.0, 10.0],
            16,
            0.0,
            |_, _, _| 3.0,
        )
        .is_err());
    }
}
