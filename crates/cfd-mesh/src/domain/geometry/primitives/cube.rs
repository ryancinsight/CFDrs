//! Axis-aligned box (cuboid) primitive.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds an axis-aligned box (cuboid) with the given dimensions.
///
/// The box is positioned with one corner at `origin` and extends
/// `(width, height, depth)` along the `+X`, `+Y`, `+Z` axes respectively.
///
/// ## Output
///
/// - 8 unique vertices, 12 faces (6 quads × 2 triangles)
/// - `RegionId(1)` on all faces
/// - `signed_volume = width × height × depth > 0`
///
/// ## Example
///
/// ```rust,ignore
/// // Unit cube at origin
/// let mesh = Cube::unit().build().unwrap();
/// ```
#[derive(Clone, Debug)]
pub struct Cube {
    /// Corner closest to (−∞, −∞, −∞).
    pub origin: Point3r,
    /// Extent along +X [mm].
    pub width: f64,
    /// Extent along +Y [mm].
    pub height: f64,
    /// Extent along +Z [mm].
    pub depth: f64,
}

impl Cube {
    /// Unit cube `[0, 1]³` at the origin.
    pub fn unit() -> Self {
        Self {
            origin: Point3r::origin(),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
    }

    /// Cube of side `s` centred at the origin.
    pub fn centred(s: f64) -> Self {
        let h = s / 2.0;
        Self {
            origin: Point3r::new(-h, -h, -h),
            width: s,
            height: s,
            depth: s,
        }
    }
}

impl Default for Cube {
    fn default() -> Self {
        Self::unit()
    }
}

impl PrimitiveMesh for Cube {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(c: &Cube) -> Result<IndexedMesh, PrimitiveError> {
    let (w, h, d) = (c.width, c.height, c.depth);
    if w <= 0.0 || h <= 0.0 || d <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "all dimensions must be > 0, got ({w}, {h}, {d})"
        )));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    let ox = c.origin.x;
    let oy = c.origin.y;
    let oz = c.origin.z;

    // 8 corner positions
    let p = [
        Point3r::new(ox, oy, oz),             // 0 left-bottom-back
        Point3r::new(ox + w, oy, oz),         // 1 right-bottom-back
        Point3r::new(ox + w, oy + h, oz),     // 2 right-top-back
        Point3r::new(ox, oy + h, oz),         // 3 left-top-back
        Point3r::new(ox, oy, oz + d),         // 4 left-bottom-front
        Point3r::new(ox + w, oy, oz + d),     // 5 right-bottom-front
        Point3r::new(ox + w, oy + h, oz + d), // 6 right-top-front
        Point3r::new(ox, oy + h, oz + d),     // 7 left-top-front
    ];

    // 6 quads: (corner indices [CCW from outside], outward normal)
    // Winding is CCW when viewed from outside the face.
    let quads: &[([usize; 4], Vector3r)] = &[
        ([0, 3, 2, 1], -Vector3r::z()), // −Z back face
        ([4, 5, 6, 7], Vector3r::z()),  // +Z front face
        ([0, 1, 5, 4], -Vector3r::y()), // −Y bottom face
        ([3, 7, 6, 2], Vector3r::y()),  // +Y top face
        ([0, 4, 7, 3], -Vector3r::x()), // −X left face
        ([1, 2, 6, 5], Vector3r::x()),  // +X right face
    ];

    for &(idx, normal) in quads {
        let [i0, i1, i2, i3] = idx;
        let v0 = mesh.add_vertex(p[i0], normal);
        let v1 = mesh.add_vertex(p[i1], normal);
        let v2 = mesh.add_vertex(p[i2], normal);
        let v3 = mesh.add_vertex(p[i3], normal);
        // Split quad into 2 CCW triangles
        mesh.add_face_with_region(v0, v1, v2, region);
        mesh.add_face_with_region(v0, v2, v3, region);
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;
    use approx::assert_relative_eq;

    #[test]
    fn cube_is_watertight() {
        let mesh = Cube::unit().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "cube must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn cube_volume_correct() {
        let mesh = Cube {
            origin: Point3r::origin(),
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert_relative_eq!(report.signed_volume, 24.0, epsilon = 1e-10);
    }

    #[test]
    fn cube_invalid_dimensions() {
        let result = Cube {
            origin: Point3r::origin(),
            width: -1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build();
        assert!(result.is_err());
    }
}
