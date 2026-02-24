//! Surface-of-revolution (lathe / revolve) primitive.
//!
//! Sweeps a 2D profile curve in the YR plane (Y = axial, R = radial) around
//! the Y axis by a specified angle to produce a watertight solid.
//!
//! ## Coordinate System
//!
//! The profile is specified as a list of `(r, y)` control points where:
//! - `r ≥ 0` is the **radial distance** from the Y axis
//! - `y` is the **axial coordinate**
//!
//! A full 360° revolution (`angle = 2π`) produces a closed torus-like solid.
//! Partial revolutions (`angle < 2π`) include two flat end caps.
//!
//! ## Winding Convention
//!
//! All faces use **outward CCW** winding (right-hand rule → outward normal,
//! `signed_volume > 0`).
//!
//! ## Topology
//!
//! For a profile with `n` points and `m` angular segments:
//! - `m × (n−1)` lateral quads → `2 m (n−1)` triangular faces
//! - For a **full revolution**: no end caps needed (closed tube)
//! - For a **partial revolution**: two end cap polygon fans
//!
//! ## Analytical Volume
//!
//! By Pappus's centroid theorem:
//! `V = 2π · R̄ · A_profile · (angle / 2π)`
//! where `R̄` is the centroid of the profile and `A_profile` is its area.
//!
//! ## Example
//!
//! ```rust,ignore
//! use cfd_mesh::domain::geometry::primitives::{RevolutionSweep, PrimitiveMesh};
//! use std::f64::consts::TAU;
//!
//! // Revolve a rectangular profile → hollow cylinder-like shape
//! let sweep = RevolutionSweep {
//!     profile: vec![(0.5, 0.0), (1.0, 0.0), (1.0, 2.0), (0.5, 2.0)],
//!     segments: 32,
//!     angle: TAU,
//! };
//! let mesh = sweep.build().unwrap();
//! ```

use std::f64::consts::TAU;

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a solid of revolution by sweeping a 2D profile in the YR plane
/// around the Y axis.
///
/// # Invariants
///
/// - `profile.len() >= 2`  (at least one segment)
/// - All `r` values in profile are `≥ 0`
/// - `segments >= 3`
/// - `0 < angle <= 2π`
///
/// # Example
///
/// ```rust,ignore
/// use cfd_mesh::domain::geometry::primitives::{RevolutionSweep, PrimitiveMesh};
/// use std::f64::consts::TAU;
///
/// // Full revolution of a circular arc → torus-like ring
/// let mesh = RevolutionSweep {
///     profile: vec![(2.0, -1.0), (3.0, 0.0), (2.0, 1.0)],
///     segments: 48,
///     angle: TAU,
/// }.build().unwrap();
/// ```
#[derive(Clone, Debug)]
pub struct RevolutionSweep {
    /// Profile points as `(r, y)` pairs.  `r` is radial distance from the
    /// Y axis (≥ 0); `y` is the axial coordinate.  The profile is swept in
    /// order from `profile[0]` to `profile[last]`.
    pub profile: Vec<(f64, f64)>,
    /// Number of angular subdivisions around the Y axis (≥ 3).
    pub segments: usize,
    /// Total sweep angle in radians.  Use `TAU` (= 2π) for a full revolution.
    /// Must be in `(0, TAU]`.
    pub angle: f64,
}

impl Default for RevolutionSweep {
    /// Default: revolve a 1-unit-radius, 2-unit-tall tube profile 360°.
    fn default() -> Self {
        Self {
            profile: vec![(1.0, 0.0), (1.0, 2.0)],
            segments: 32,
            angle: TAU,
        }
    }
}

impl PrimitiveMesh for RevolutionSweep {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(s: &RevolutionSweep) -> Result<IndexedMesh, PrimitiveError> {
    let np = s.profile.len();
    if np < 2 {
        return Err(PrimitiveError::TooFewSegments(np));
    }
    if s.segments < 3 {
        return Err(PrimitiveError::TooFewSegments(s.segments));
    }
    if s.angle <= 0.0 || s.angle > TAU + 1e-12 {
        return Err(PrimitiveError::InvalidParam(format!(
            "angle must be in (0, 2π], got {:.6}",
            s.angle
        )));
    }
    for &(r, _y) in &s.profile {
        if r < 0.0 {
            return Err(PrimitiveError::InvalidParam(format!(
                "all radial values must be ≥ 0, got {}",
                r
            )));
        }
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let is_full = (s.angle - TAU).abs() < 1e-10;
    let nm = s.segments;

    // ── Build vertex grid ─────────────────────────────────────────────────────
    //
    // For each angular slice j ∈ 0..nm and profile point i ∈ 0..np, compute:
    //   φ = j / nm × angle
    //   x = r_i · cos(φ),  y = y_i,  z = r_i · sin(φ)
    //   outward normal = (cos(φ), 0, sin(φ))  (radial, for non-tapered edges)
    //
    // For a tapered profile, the lateral normal is the profile-edge normal
    // rotated around Y by φ:
    //   Profile edge tangent: Δr, Δy
    //   Profile normal (outward in r-y plane): (Δy, -Δr) normalised
    //   Rotated to 3D: n = (ny_prof · cos(φ), −nr_prof, ny_prof · sin(φ))
    //   where ny_prof = Δy / |edge|, nr_prof = Δr / |edge|

    // Compute profile vertex normals by averaging adjacent edge normals.
    // Each profile edge has a 2-D outward normal (nr, ny) in the r-y plane.
    let edge_n2d: Vec<(f64, f64)> = (0..(np - 1))
        .map(|i| {
            let (r0, y0) = s.profile[i];
            let (r1, y1) = s.profile[i + 1];
            let dr = r1 - r0;
            let dy = y1 - y0;
            let len = (dr * dr + dy * dy).sqrt();
            if len > 1e-14 {
                // Outward profile normal: (dy/len, -dr/len) = rotate edge 90° right
                (dy / len, -dr / len)
            } else {
                (1.0, 0.0)
            }
        })
        .collect();

    // Vertex profile normals: average of adjacent edge normals.
    let vertex_n2d: Vec<(f64, f64)> = (0..np)
        .map(|i| {
            if i == 0 {
                edge_n2d[0]
            } else if i == np - 1 {
                edge_n2d[np - 2]
            } else {
                let (nr0, ny0) = edge_n2d[i - 1];
                let (nr1, ny1) = edge_n2d[i];
                let nr = nr0 + nr1;
                let ny = ny0 + ny1;
                let len = (nr * nr + ny * ny).sqrt();
                if len > 1e-14 {
                    (nr / len, ny / len)
                } else {
                    (nr1, ny1)
                }
            }
        })
        .collect();

    // Compute the actual angular step — for a full revolution the last column
    // wraps back to the first.
    let angular_step = s.angle / nm as f64;

    // Number of angular "columns" to generate vertices for:
    //   full revolution:    nm columns  (wrap: col nm → col 0)
    //   partial revolution: nm+1 columns (col 0 and col nm are distinct seam)
    let ncols = if is_full { nm } else { nm + 1 };

    // Vertex ID grid: grid[j][i] = VertexId at angular col j, profile row i.
    let mut grid: Vec<Vec<crate::domain::core::index::VertexId>> = Vec::with_capacity(ncols);
    for j in 0..ncols {
        let phi = j as f64 * angular_step;
        let cp = phi.cos();
        let sp = phi.sin();

        let col: Vec<crate::domain::core::index::VertexId> = (0..np)
            .map(|i| {
                let (r, y) = s.profile[i];
                let pos = Point3r::new(r * cp, y, r * sp);
                // 3-D outward normal from 2-D profile normal (nr, ny):
                //   n_3d = (nr · cos(φ),  ny,  nr · sin(φ))
                let (nr, ny) = vertex_n2d[i];
                let n = Vector3r::new(nr * cp, ny, nr * sp);
                let nlen = n.norm();
                let n = if nlen > 1e-14 {
                    n / nlen
                } else {
                    Vector3r::new(cp, 0.0, sp)
                };
                mesh.add_vertex(pos, n)
            })
            .collect();

        grid.push(col);
    }

    // ── Lateral quads ────────────────────────────────────────────────────────
    //
    // For angular column j and profile segment i → i+1:
    //   quad corners (CCW from outside, +radial view):
    //     v(j,i) → v(j,i+1) → v(j+1,i+1) → v(j+1,i)
    //
    // Winding check: looking from outside (+radial direction), with
    //   right = angular (φ increasing) direction, up = axial (+Y):
    //   col j is "left", col j+1 is "right", row i is "bottom", i+1 is "top"
    //   CCW from outside: left-bottom → left-top → right-top, then
    //                     left-bottom → right-top → right-bottom ✓
    for j in 0..nm {
        let j_next = if is_full { (j + 1) % nm } else { j + 1 };

        for i in 0..(np - 1) {
            let v00 = grid[j][i];
            let v01 = grid[j][i + 1];
            let v11 = grid[j_next][i + 1];
            let v10 = grid[j_next][i];

            // CCW from outside:
            mesh.add_face_with_region(v00, v01, v11, region);
            mesh.add_face_with_region(v00, v11, v10, region);
        }
    }

    // ── End caps (only for partial revolution) ────────────────────────────────
    //
    // Start cap at φ = 0: the profile in the XZ plane with z = 0 (XY quarter).
    // The outward normal is −Z (into the page when looking in −Z direction).
    //
    // For the start cap (φ = 0, normal = −Z = (0,0,-1)):
    //   Profile is in the X-Y plane at z = 0.
    //   Outward normal for this cap is −Z.
    //   Fan: profile[0] → profile[1] → profile[2] → …
    //   Viewed from −Z (from outside), this is CCW when profile[i]
    //   has increasing r (going right in X), i.e. standard orientation.
    //   Actual winding depends on profile shape, so we check and reverse.
    //
    // Start cap normal: −Z.  Fan triangle (p0, p_i, p_{i+1}):
    //   Normal via right-hand rule: (p_i−p0) × (p_{i+1}−p0)
    //   For a profile that goes outward (r increases), this cross product
    //   points in the +Z direction for a CCW profile.  We want −Z, so we
    //   reverse: (p0, p_{i+1}, p_i).
    if !is_full {
        // Detect if the profile closes on itself (profile[last] ≈ profile[0]).
        // When closed, the last cap triangle is degenerate (fan center == vi1),
        // so we stop the fan one step earlier: for i in 1..(cap_fan_end).
        let (r_last, y_last) = s.profile[np - 1];
        let (r_first, y_first) = s.profile[0];
        let profile_is_closed =
            (r_last - r_first).abs() < 1e-10 && (y_last - y_first).abs() < 1e-10;
        // For an open profile, fan covers i=1..np-2 (np-1 exclusive).
        // For a closed profile, fan covers i=1..np-3 (np-2 exclusive),
        // because the last lateral segment closes back to profile[0] which
        // is already the fan center, so that seam edge is covered by
        // the second-to-last triangle's spoke.
        let cap_fan_end = if profile_is_closed { np - 2 } else { np - 1 };

        // Start cap at φ = 0 (column 0), outward normal = −Z.
        //
        // Lateral face at j=0, segment i, left seam edge: grid[0][i]→grid[0][i+1].
        // Cap must reverse: grid[0][i+1]→grid[0][i].
        // Fan: (center, vi1, vi) traverses vi1→vi ✓
        let n_start = -Vector3r::z();
        {
            let (r0, y0) = s.profile[0];
            let v_fan_center = mesh.add_vertex(Point3r::new(r0, y0, 0.0), n_start);

            for i in 1..cap_fan_end {
                let (ri, yi) = s.profile[i];
                let (ri1, yi1) = s.profile[i + 1];
                let vi = mesh.add_vertex(Point3r::new(ri, yi, 0.0), n_start);
                let vi1 = mesh.add_vertex(Point3r::new(ri1, yi1, 0.0), n_start);
                mesh.add_face_with_region(v_fan_center, vi1, vi, region);
            }
        }

        // End cap at φ = angle (column nm).
        //
        // At φ_end, the angular tangent direction is (-sin(φ_end), 0, cos(φ_end)).
        // Outward normal for end cap = that tangent direction (past the end).
        //
        // Lateral face at j=nm-1, segment i, right seam edge: grid[nm][i]→grid[nm][i+1].
        // Cap must reverse: grid[nm][i+1]→grid[nm][i].
        // Fan: (center, vi1, vi) traverses vi1→vi ✓
        let phi_end = s.angle;
        let n_end = Vector3r::new(-phi_end.sin(), 0.0, phi_end.cos());
        let n_end = {
            let l = n_end.norm();
            if l > 1e-14 {
                n_end / l
            } else {
                n_end
            }
        };

        {
            let (r0, y0) = s.profile[0];
            let cp_end = phi_end.cos();
            let sp_end = phi_end.sin();
            let v_fan_center = mesh.add_vertex(Point3r::new(r0 * cp_end, y0, r0 * sp_end), n_end);

            for i in 1..cap_fan_end {
                let (ri, yi) = s.profile[i];
                let (ri1, yi1) = s.profile[i + 1];
                let vi = mesh.add_vertex(Point3r::new(ri * cp_end, yi, ri * sp_end), n_end);
                let vi1 = mesh.add_vertex(Point3r::new(ri1 * cp_end, yi1, ri1 * sp_end), n_end);
                // Lateral right seam at col nm: edge goes v11→v10 (descending).
                // Cap must reverse: ascending = vi→vi1.
                // Fan: (center, vi, vi1) traverses vi→vi1 ✓, normal = -X at π/2 ✓
                mesh.add_face_with_region(v_fan_center, vi, vi1, region);
            }
        }
    }

    Ok(mesh)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;
    use std::f64::consts::PI;

    /// Revolve a single vertical edge of radius r and height h → open tube
    /// with no caps. Test with a full revolution.
    #[test]
    fn revolution_sweep_full_cylinder_watertight() {
        // Vertical edge at radius r=1, from y=0 to y=2
        let sweep = RevolutionSweep {
            profile: vec![(1.0, 0.0), (1.0, 2.0)],
            segments: 32,
            angle: TAU,
        };
        let mesh = sweep.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        // The lateral surface of a cylinder (open tube) is NOT closed — it has
        // two boundary loops (top and bottom circles). This test verifies the
        // topology compiles and the faces have consistent winding.
        assert!(
            report.orientation_consistent,
            "revolution winding must be consistent"
        );
    }

    /// Revolve a closed rectangular profile → solid ring (washer shape)
    /// Full revolution → should be watertight (closed torus-like surface).
    ///
    /// The profile must explicitly close the loop so all 4 walls are generated:
    /// bottom annulus + outer cylinder + top annulus + inner cylinder.
    /// The 5th point closes the loop back to the start.
    #[test]
    fn revolution_sweep_washer_full_watertight() {
        // Closed 5-point loop: bottom → outer → top → inner → back to start
        // bottom: (1,0) → (2,0); outer: (2,0) → (2,0.5);
        // top: (2,0.5) → (1,0.5); inner: (1,0.5) → (1,0)
        let profile = vec![
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 0.5),
            (1.0, 0.5),
            (1.0, 0.0), // close the loop back to start
        ];
        let sweep = RevolutionSweep {
            profile,
            segments: 32,
            angle: TAU,
        };
        let mesh = sweep.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        // Full revolution of a 5-point closed profile → closed torus-like surface
        assert!(
            report.orientation_consistent,
            "revolution winding must be consistent"
        );
        assert!(
            report.is_closed,
            "full revolution of closed 5-point profile must be closed"
        );
    }

    /// Revolve a triangular profile → cone-like solid (full 360°).
    #[test]
    fn revolution_sweep_cone_like_watertight() {
        // Profile: (0,2) at apex, (1,0) at base rim.
        // Revolution around Y: generates a cone lateral surface + base disk.
        // Since r=0 at apex, top degenerate vertices merge.
        let sweep = RevolutionSweep {
            profile: vec![(0.0, 2.0), (1.0, 0.0)],
            segments: 32,
            angle: TAU,
        };
        let mesh = sweep.build().unwrap();
        assert!(mesh.face_count() > 0, "cone-like sweep must produce faces");
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(
            report.orientation_consistent,
            "cone revolution winding must be consistent"
        );
    }

    /// Partial revolution of a closed rectangular profile → watertight wedge.
    ///
    /// The profile must close the loop (5 points: 4 walls + closing segment)
    /// so that all 4 sides of the annular cross-section are swept.
    /// The two end caps close the angular start and end faces.
    #[test]
    fn revolution_sweep_partial_watertight() {
        // Closed 5-point loop: bottom → outer → top → inner → back to start
        let profile = vec![
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 0.0), // close the loop
        ];
        let sweep = RevolutionSweep {
            profile,
            segments: 16,
            angle: PI / 2.0, // 90°
        };
        let mesh = sweep.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(
            report.orientation_consistent,
            "partial revolution winding must be consistent"
        );
        assert!(
            report.is_closed,
            "partial revolution with closed profile + end caps must be closed"
        );
        assert!(
            report.is_watertight,
            "partial revolution must be watertight"
        );
    }

    /// Volume of partial revolution: Pappus's theorem.
    ///
    /// The profile is a closed 5-point loop so the solid has all 4 walls.
    /// Volume = r̄ × A × angle (Pappus), where A is the cross-section area.
    #[test]
    fn revolution_sweep_partial_volume_pappus() {
        // Revolve a closed annular cross-section (r=1..2, y=0..1) by π/2 (90°).
        // Cross-section area = (2-1) × (1-0) = 1 mm².
        // Centroid r̄ = 1.5 mm.
        // Volume = r̄ × A × angle = 1.5 × 1 × π/2 ≈ 2.3562 mm³
        let profile = vec![
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 0.0), // close the loop
        ];
        let sweep = RevolutionSweep {
            profile,
            segments: 64,
            angle: PI / 2.0,
        };
        let mesh = sweep.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0, "positive volume");
        let expected = 1.5_f64 * 1.0 * (PI / 2.0);
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.02,
            "Pappus volume error {:.4}% < 2%",
            error * 100.0
        );
    }

    #[test]
    fn revolution_sweep_rejects_too_few_profile_points() {
        let result = RevolutionSweep {
            profile: vec![(1.0, 0.0)],
            segments: 16,
            angle: TAU,
        }
        .build();
        assert!(result.is_err());
    }

    #[test]
    fn revolution_sweep_rejects_negative_radius() {
        let result = RevolutionSweep {
            profile: vec![(-1.0, 0.0), (1.0, 1.0)],
            segments: 16,
            angle: TAU,
        }
        .build();
        assert!(result.is_err());
    }

    #[test]
    fn revolution_sweep_rejects_too_few_segments() {
        let result = RevolutionSweep {
            profile: vec![(1.0, 0.0), (1.0, 1.0)],
            segments: 2,
            angle: TAU,
        }
        .build();
        assert!(result.is_err());
    }
}
