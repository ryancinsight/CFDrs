//! Flat circular disk — a 2-D surface mesh in the XZ plane.
//!
//! A `Disk` produces a triangulated flat fan of `segments` triangles, all
//! lying in the plane `y = center.y`.  Unlike a `Cylinder` (which is a closed
//! 3-D solid), a disk is an **open surface**: it has a boundary (its rim) and
//! zero signed volume.  It is the natural 2-D operand for coplanar Boolean
//! operations.
//!
//! ## Winding and normals
//!
//! All triangles are wound counter-clockwise when viewed from `+Y` (the
//! default upward normal direction).  Every vertex carries a stored normal of
//! `(0, +1, 0)`.  If you need the disk to face `−Y`, call `mesh.flip_faces()`
//! after building.
//!
//! ## Geometry
//!
//! ```text
//!  centre = (cx, y, cz)        radius = r        segments = N
//!
//!  vertex i  (i = 0..N−1):
//!    x = cx + r · cos(2π·i / N)
//!    y = y
//!    z = cz + r · sin(2π·i / N)
//!
//!  face i: (centre, rim[i], rim[(i+1) % N])   — CCW from +Y
//! ```
//!
//! ## Area
//!
//! `A = π r²`  (exact in the limit N → ∞; discretisation error ~ O(1/N²))
//!
//! At `N = 64` the area error is < 0.04%.
//!
//! ## Theorem — Convergence of Polygon Area to π r²
//!
//! The area of a regular N-gon inscribed in a circle of radius r is
//! `A_N = (N/2) r² sin(2π/N)`.  As N → ∞:
//! `lim A_N = (N/2) r² · (2π/N) = π r²`  (Taylor expansion of sin).

use std::f64::consts::TAU;

use crate::core::scalar::{Point3r, Vector3r};
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Flat circular disk in the plane `y = base_center.y`.
///
/// # Invariants
///
/// - `radius > 0`
/// - `segments >= 3`
/// - All output triangles lie in the plane `y = base_center.y`
/// - All vertex normals are `(0, +1, 0)`
/// - Triangle winding is CCW as viewed from `+Y`
#[derive(Debug, Clone)]
pub struct Disk {
    /// Centre of the disk (the `y` component sets the plane height).
    pub center: Point3r,
    /// Disk radius (> 0).
    pub radius: f64,
    /// Number of rim vertices / triangles (>= 3).
    pub segments: usize,
}

impl Default for Disk {
    fn default() -> Self {
        Self {
            center:   Point3r::origin(),
            radius:   1.0,
            segments: 32,
        }
    }
}

impl PrimitiveMesh for Disk {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        if self.radius <= 0.0 {
            return Err(PrimitiveError::InvalidParam(
                format!("Disk radius must be > 0, got {}", self.radius),
            ));
        }
        if self.segments < 3 {
            return Err(PrimitiveError::TooFewSegments(self.segments));
        }

        let mut mesh = IndexedMesh::new();

        // All vertices face +Y.
        let up = Vector3r::new(0.0, 1.0, 0.0);

        let cx = self.center.x;
        let cy = self.center.y;
        let cz = self.center.z;

        // Centre vertex.
        let centre_id = mesh.add_vertex(self.center, up);

        // Rim vertices: CCW from +Y means angle advances positively around +Y axis.
        // Convention: x = cos(θ), z = −sin(θ) gives CCW winding when viewed from +Y
        // (right-hand rule: n = (rim[i]−c) × (rim[i+1]−c) points in +Y direction).
        let rim: Vec<_> = (0..self.segments)
            .map(|i| {
                let theta = TAU * (i as f64) / (self.segments as f64);
                let pos = Point3r::new(
                    cx + self.radius * theta.cos(),
                    cy,
                    cz - self.radius * theta.sin(),
                );
                mesh.add_vertex(pos, up)
            })
            .collect();

        // Triangulate as a fan from centre.
        for i in 0..self.segments {
            let j = (i + 1) % self.segments;
            mesh.add_face(centre_id, rim[i], rim[j]);
        }

        Ok(mesh)
    }
}

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn build_disk(r: f64, n: usize) -> IndexedMesh {
        Disk { center: Point3r::origin(), radius: r, segments: n }
            .build()
            .expect("Disk::build failed")
    }

    /// # Theorem
    ///
    /// A regular N-gon fan has exactly N triangles (one per rim edge).
    #[test]
    fn disk_face_count_equals_segments() {
        let mesh = build_disk(1.0, 32);
        assert_eq!(mesh.face_count(), 32);
    }

    /// Vertex count = 1 (centre) + segments (rim).
    #[test]
    fn disk_vertex_count() {
        let n = 48usize;
        let mesh = build_disk(1.0, n);
        assert_eq!(mesh.vertices.len(), n + 1);
    }

    /// All vertices lie in the plane y = 0.
    #[test]
    fn disk_all_vertices_coplanar_y0() {
        let mesh = build_disk(2.0, 64);
        for (_, v) in mesh.vertices.iter() {
            assert_relative_eq!(v.position.y, 0.0, epsilon = 1e-12);
        }
    }

    /// Area converges to π r² from below (inscribed polygon).
    ///
    /// # Theorem
    ///
    /// `A_N = (N/2) r² sin(2π/N) → π r²`  as `N → ∞`
    #[test]
    fn disk_area_converges_to_pi_r_squared() {
        let r = 2.5_f64;
        let mesh = build_disk(r, 256);
        let area: f64 = mesh.faces.iter().map(|f| {
            let a = mesh.vertices.position(f.vertices[0]);
            let b = mesh.vertices.position(f.vertices[1]);
            let c = mesh.vertices.position(f.vertices[2]);
            let ab = b - a;
            let ac = c - a;
            ab.cross(&ac).norm() * 0.5
        }).sum();
        let expected = std::f64::consts::PI * r * r;
        let err = (area - expected).abs() / expected;
        assert!(err < 0.001, "area error {:.4}% > 0.1%", err * 100.0);
    }

    /// All vertex normals point in +Y.
    #[test]
    fn disk_normals_point_up() {
        let mesh = build_disk(1.0, 32);
        for (_, v) in mesh.vertices.iter() {
            assert_relative_eq!(v.normal.x, 0.0, epsilon = 1e-12);
            assert_relative_eq!(v.normal.y, 1.0, epsilon = 1e-12);
            assert_relative_eq!(v.normal.z, 0.0, epsilon = 1e-12);
        }
    }

    /// Disk rejects non-positive radius.
    #[test]
    fn disk_rejects_zero_radius() {
        let err = Disk { center: Point3r::origin(), radius: 0.0, segments: 32 }.build();
        assert!(err.is_err());
    }

    /// Disk rejects too few segments.
    #[test]
    fn disk_rejects_too_few_segments() {
        let err = Disk { center: Point3r::origin(), radius: 1.0, segments: 2 }.build();
        assert!(err.is_err());
    }

    /// Disk centred at non-origin lies in plane y = cy.
    #[test]
    fn disk_offset_center() {
        let cy = 3.7_f64;
        let disk = Disk {
            center: Point3r::new(1.0, cy, -2.0),
            radius: 0.5,
            segments: 16,
        }.build().expect("build failed");
        for (_, v) in disk.vertices.iter() {
            assert_relative_eq!(v.position.y, cy, epsilon = 1e-12);
        }
    }
}
