//! Linear-extrusion sweep primitive.
//!
//! Sweeps a closed 2D polygon profile along the +Y axis to produce a
//! watertight prismatic solid.
//!
//! ## Parametrisation
//!
//! Given a profile of `n` vertices `p_0 … p_{n-1}` in the XZ plane (y = 0),
//! the sweep extrudes them by `height` along +Y.  The polygon must be
//! **convex or simple** (no self-intersections) and vertices must be wound
//! **counter-clockwise when viewed from −Y** so that the bottom cap face
//! normals point downward (outward).
//!
//! ## Topology
//!
//! ```text
//! Faces = 2·(n−2)   (bottom cap — triangle fan)
//!       + 2·(n−2)   (top cap    — triangle fan)
//!       + 2·n       (lateral quads → 2 triangles each)
//! Total = 4n − 4    faces
//! ```
//!
//! All faces carry `RegionId(1)`.
//!
//! ## Winding Convention
//!
//! Profile vertices are supplied in the XZ plane (y = 0).  They must be
//! ordered **CCW when viewed from −Y** (i.e. standard mathematical
//! counter-clockwise in the XZ plane with X right, Z up, looking in the +Y
//! direction, which is the same as CW if you think of X right, Z down).
//!
//! Concretely: the signed area of the profile computed as
//! `Σ (x_i · z_{i+1} − x_{i+1} · z_i) / 2` must be **positive** for a
//! valid (outward bottom normal) winding.  The builder checks this and
//! returns [`PrimitiveError::InvalidParam`] otherwise.
//!
//! ## Analytical Volume
//!
//! `V = A_profile × height`
//! where `A_profile` is the signed area of the polygon.
//!
//! ## Example
//!
//! ```rust,ignore
//! use cfd_mesh::domain::geometry::primitives::{LinearSweep, PrimitiveMesh};
//! use cfd_mesh::domain::core::scalar::Point2r;
//!
//! // Square cross-section 2×2, swept 3 mm tall
//! let profile = vec![
//!     Point2r::new(-1.0, -1.0),
//!     Point2r::new( 1.0, -1.0),
//!     Point2r::new( 1.0,  1.0),
//!     Point2r::new(-1.0,  1.0),
//! ];
//! let mesh = LinearSweep { profile, height: 3.0 }.build().unwrap();
//! ```

use std::f64::consts::TAU;

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

// Re-use nalgebra's 2-D point via f64 scalars.
/// 2-D point in the XZ plane.
pub type Point2 = nalgebra::Point2<f64>;

/// Builds a watertight prismatic solid by linearly sweeping a 2D polygon
/// profile along the +Y axis.
///
/// # Invariants
///
/// - `profile.len() >= 3`
/// - Profile polygon has positive signed area (CCW when viewed from −Y)
/// - `height > 0`
///
/// # Example
///
/// ```rust,ignore
/// use cfd_mesh::domain::geometry::primitives::{LinearSweep, PrimitiveMesh};
///
/// let triangle = LinearSweep {
///     profile: vec![
///         LinearSweep::pt(-1.0, 0.0),
///         LinearSweep::pt( 1.0, 0.0),
///         LinearSweep::pt( 0.0, 1.732),
///     ],
///     height: 2.0,
/// };
/// let mesh = triangle.build().unwrap();
/// ```
#[derive(Clone, Debug)]
pub struct LinearSweep {
    /// Profile vertices in the XZ plane (y = 0), ordered CCW when viewed
    /// from −Y.  The polygon is automatically closed (last vertex connects
    /// back to the first).
    pub profile: Vec<Point2>,
    /// Extrusion length along +Y [mm].
    pub height: f64,
}

impl LinearSweep {
    /// Convenience constructor for a profile point (x, z).
    pub fn pt(x: f64, z: f64) -> Point2 {
        Point2::new(x, z)
    }

    /// Build a regular n-gon (convex polygon) profile of the given radius.
    /// The first vertex is at angle 0 (i.e. `(radius, 0)` in XZ).
    pub fn regular_polygon(n: usize, radius: f64) -> Vec<Point2> {
        (0..n)
            .map(|i| {
                let a = i as f64 / n as f64 * TAU;
                Point2::new(radius * a.cos(), radius * a.sin())
            })
            .collect()
    }

    /// Compute the signed area of the profile polygon.
    ///
    /// Positive → CCW in XZ plane (correct winding for outward bottom normal).
    pub fn signed_area(profile: &[Point2]) -> f64 {
        let n = profile.len();
        let mut area = 0.0_f64;
        for i in 0..n {
            let j = (i + 1) % n;
            area += profile[i].x * profile[j].y - profile[j].x * profile[i].y;
        }
        area * 0.5
    }
}

impl Default for LinearSweep {
    /// Default: unit square (1×1) profile, 1 mm tall.
    fn default() -> Self {
        Self {
            profile: vec![
                Point2::new(-0.5, -0.5),
                Point2::new(0.5, -0.5),
                Point2::new(0.5, 0.5),
                Point2::new(-0.5, 0.5),
            ],
            height: 1.0,
        }
    }
}

impl PrimitiveMesh for LinearSweep {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(s: &LinearSweep) -> Result<IndexedMesh, PrimitiveError> {
    let n = s.profile.len();
    if n < 3 {
        return Err(PrimitiveError::TooFewSegments(n));
    }
    if s.height <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "height must be > 0, got {}",
            s.height
        )));
    }

    let area = LinearSweep::signed_area(&s.profile);
    if area <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "profile signed area must be > 0 (CCW winding), got {:.6}",
            area
        )));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let h = s.height;

    // ── Edge normals for the lateral surface ─────────────────────────────────
    //
    // For each edge (p_i → p_{i+1}) in the XZ plane, the outward lateral
    // normal is perpendicular to the edge and pointing away from the interior.
    //
    // Edge vector: d = (dx, dz) = (p_{i+1}.x - p_i.x, p_{i+1}.z - p_i.z)
    // Outward perpendicular (CCW polygon): n = (dz, -dx) normalised.
    let edge_normals: Vec<Vector3r> = (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            let dx = s.profile[j].x - s.profile[i].x;
            let dz = s.profile[j].y - s.profile[i].y;
            // 2-D outward normal: (dz, -dx) — points right of the CCW edge direction
            let nx = dz;
            let nz = -dx;
            let len = (nx * nx + nz * nz).sqrt();
            if len > 1e-14 {
                Vector3r::new(nx / len, 0.0, nz / len)
            } else {
                Vector3r::new(1.0, 0.0, 0.0) // degenerate edge fallback
            }
        })
        .collect();

    // ── Lateral surface ───────────────────────────────────────────────────────
    //
    // For edge i → j, the lateral quad is:
    //   p_i at y=0, p_j at y=0, p_j at y=h, p_i at y=h
    //
    // Outward CCW from outside (normal = edge_normals[i]):
    //   tri 1: bot_i → bot_j → top_j
    //   tri 2: bot_i → top_j → top_i
    //
    // Verify: looking from outside (from the normal direction), with
    //   X-right = tangent along edge, Y-up = up axis:
    //   bot_i (left-bottom) → bot_j (right-bottom) → top_j (right-top) is CCW ✓
    for i in 0..n {
        let j = (i + 1) % n;
        let en = edge_normals[i];

        // Use the average of the two vertex normals for the edge midpoint.
        // For smooth lighting each vertex carries its own edge normal.
        let ni = {
            // vertex i is shared between edges (i-1, i) and (i, i+1);
            // for a flat-shaded prism, edge normal suffices.
            let prev_i = (i + n - 1) % n;
            let en_prev = edge_normals[prev_i];
            let avg = en + en_prev;
            let len = avg.norm();
            if len > 1e-14 {
                avg / len
            } else {
                en
            }
        };
        let nj = {
            let next_j = (j + 1) % n;
            let en_next = edge_normals[j];
            let avg = en + en_next;
            let _ = next_j; // suppress unused warning
            let len = avg.norm();
            if len > 1e-14 {
                avg / len
            } else {
                en
            }
        };

        let pi = &s.profile[i];
        let pj = &s.profile[j];

        let p_bot_i = Point3r::new(pi.x, 0.0, pi.y);
        let p_bot_j = Point3r::new(pj.x, 0.0, pj.y);
        let p_top_i = Point3r::new(pi.x, h, pi.y);
        let p_top_j = Point3r::new(pj.x, h, pj.y);

        let v_bot_i = mesh.add_vertex(p_bot_i, ni);
        let v_bot_j = mesh.add_vertex(p_bot_j, nj);
        let v_top_i = mesh.add_vertex(p_top_i, ni);
        let v_top_j = mesh.add_vertex(p_top_j, nj);

        // Outward CCW from outside:
        // tri1: (bot_i, top_j, bot_j) — outward normal = right-perpendicular of edge ✓
        // tri2: (bot_i, top_i, top_j)
        mesh.add_face_with_region(v_bot_i, v_top_j, v_bot_j, region);
        mesh.add_face_with_region(v_bot_i, v_top_i, v_top_j, region);
    }

    // ── Bottom cap (y = 0, normal −Y) ─────────────────────────────────────────
    //
    // Lateral face bottom seam edge (from first tri) goes  bot_j → bot_i.
    // The cap must traverse that same edge in reverse: bot_i → bot_j.
    //
    // Fan triangle for profile segment i (i=1..n-2):
    //   (p_0, p_i, p_{i+1})   ← forward fan
    //   edge covered: p_i → p_{i+1}  (= bot_i → bot_j) ✓
    //
    // Normal check: (p_i-p_0) × (p_{i+1}-p_0) in the XZ plane.
    //   profile is CCW when viewed from -Y (from below), so this cross product
    //   points in -Y → outward bottom normal ✓
    {
        let n_down = -Vector3r::y();
        let p0 = &s.profile[0];
        let v0 = mesh.add_vertex(Point3r::new(p0.x, 0.0, p0.y), n_down);

        for i in 1..(n - 1) {
            let pi = &s.profile[i];
            let pi1 = &s.profile[i + 1];
            let vi = mesh.add_vertex(Point3r::new(pi.x, 0.0, pi.y), n_down);
            let vi1 = mesh.add_vertex(Point3r::new(pi1.x, 0.0, pi1.y), n_down);
            // Forward fan: (v0, v_i, v_{i+1}) → covers edge bot_i→bot_j ✓, normal -Y ✓
            mesh.add_face_with_region(v0, vi, vi1, region);
        }
    }

    // ── Top cap (y = h, normal +Y) ────────────────────────────────────────────
    //
    // Lateral face top seam edge (from second tri) goes top_i → top_j.
    // The cap must traverse that same edge in reverse: top_j → top_i.
    //
    // Fan triangle for profile segment i (i=1..n-2):
    //   (p_0, p_{i+1}, p_i)   ← reversed fan
    //   edge covered: p_{i+1} → p_i  (= top_j → top_i) ✓
    //
    // Normal check: (p_{i+1}-p_0) × (p_i-p_0) → points +Y → outward top cap ✓
    {
        let n_up = Vector3r::y();
        let p0 = &s.profile[0];
        let v0 = mesh.add_vertex(Point3r::new(p0.x, h, p0.y), n_up);

        for i in 1..(n - 1) {
            let pi = &s.profile[i];
            let pi1 = &s.profile[i + 1];
            let vi = mesh.add_vertex(Point3r::new(pi.x, h, pi.y), n_up);
            let vi1 = mesh.add_vertex(Point3r::new(pi1.x, h, pi1.y), n_up);
            // Reversed fan: (v0, v_{i+1}, v_i) → covers edge top_j→top_i ✓, normal +Y ✓
            mesh.add_face_with_region(v0, vi1, vi, region);
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

    #[test]
    fn linear_sweep_square_is_watertight() {
        // Square profile: 2×2, swept 3 mm tall → should be identical to a box
        let profile = vec![
            Point2::new(-1.0, -1.0),
            Point2::new(1.0, -1.0),
            Point2::new(1.0, 1.0),
            Point2::new(-1.0, 1.0),
        ];
        let mesh = LinearSweep {
            profile,
            height: 3.0,
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_closed, "sweep must be closed");
        assert!(
            report.orientation_consistent,
            "sweep must be consistently oriented"
        );
        assert!(report.is_watertight, "sweep must be watertight");
        assert_eq!(report.euler_characteristic, Some(2), "genus-0 solid: χ=2");
    }

    #[test]
    fn linear_sweep_square_volume_correct() {
        // 2×2 square × 3 tall → volume = 4 × 3 = 12 mm³
        let profile = vec![
            Point2::new(-1.0, -1.0),
            Point2::new(1.0, -1.0),
            Point2::new(1.0, 1.0),
            Point2::new(-1.0, 1.0),
        ];
        let mesh = LinearSweep {
            profile,
            height: 3.0,
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(
            report.signed_volume > 0.0,
            "outward normals → positive volume"
        );
        let expected = 12.0_f64; // exact for a box
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 1e-10,
            "square sweep volume must be exact, error={:.2e}",
            error
        );
    }

    #[test]
    fn linear_sweep_triangle_is_watertight() {
        // Equilateral-ish triangle
        let profile = vec![
            Point2::new(-1.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.732),
        ];
        let mesh = LinearSweep {
            profile,
            height: 2.0,
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "triangular prism must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn linear_sweep_hexagonal_prism_watertight_and_volume_correct() {
        // Regular hexagon r=1, swept h=2 → A = 3√3/2 ≈ 2.5981, V ≈ 5.1962
        let n = 6_usize;
        let r = 1.0_f64;
        let h = 2.0_f64;
        let profile = LinearSweep::regular_polygon(n, r);
        let mesh = LinearSweep { profile, height: h }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "hexagonal prism must be watertight");

        // Exact hexagon area: 3√3/2 * r²
        let expected_area = 3.0 * (3.0_f64).sqrt() / 2.0 * r * r;
        let expected_vol = expected_area * h;
        let error = (report.signed_volume - expected_vol).abs() / expected_vol;
        assert!(report.signed_volume > 0.0, "positive volume");
        assert!(
            error < 0.01,
            "hexagonal prism volume error {:.4}% < 1%",
            error * 100.0
        );
    }

    #[test]
    fn linear_sweep_rejects_cw_profile() {
        // CW profile (negative area) should be rejected
        let profile = vec![
            Point2::new(-1.0, 1.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, -1.0),
            Point2::new(-1.0, -1.0),
        ];
        let result = LinearSweep {
            profile,
            height: 1.0,
        }
        .build();
        assert!(result.is_err(), "CW profile should return InvalidParam");
    }

    #[test]
    fn linear_sweep_rejects_zero_height() {
        let profile = vec![
            Point2::new(-1.0, -1.0),
            Point2::new(1.0, -1.0),
            Point2::new(1.0, 1.0),
        ];
        let result = LinearSweep {
            profile,
            height: 0.0,
        }
        .build();
        assert!(result.is_err());
    }

    #[test]
    fn linear_sweep_rejects_too_few_points() {
        let profile = vec![Point2::new(-1.0, -1.0), Point2::new(1.0, -1.0)];
        let result = LinearSweep {
            profile,
            height: 1.0,
        }
        .build();
        assert!(result.is_err());
    }
}
