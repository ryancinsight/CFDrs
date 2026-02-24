//! Substrate (chip body) generation.
//!
//! Creates the outer chip body as a cuboid, with channels subtracted via CSG.
//! Equivalent to blue2mesh's cuboid substrate generation.

use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Builder for millifluidic chip substrates.
pub struct SubstrateBuilder {
    /// Width of the substrate (X dimension).
    pub width: Real,
    /// Depth of the substrate (Y dimension).
    pub depth: Real,
    /// Height of the substrate (Z dimension).
    pub height: Real,
    /// Origin corner (min X, min Y, min Z).
    pub origin: Point3r,
}

impl SubstrateBuilder {
    /// Create a new substrate builder with the given dimensions.
    pub fn new(width: Real, depth: Real, height: Real) -> Self {
        Self {
            width,
            depth,
            height,
            origin: Point3r::origin(),
        }
    }

    /// Set the origin corner.
    pub fn with_origin(mut self, origin: Point3r) -> Self {
        self.origin = origin;
        self
    }

    /// Generate the cuboid substrate mesh.
    ///
    /// Produces 12 triangles (2 per face × 6 faces) with outward normals.
    pub fn build(&self, vertex_pool: &mut VertexPool, region: RegionId) -> Vec<FaceData> {
        let o = self.origin;
        let w = self.width;
        let d = self.depth;
        let h = self.height;

        // 8 corners of the cuboid
        let corners = [
            Point3r::new(o.x, o.y, o.z),             // 0: min
            Point3r::new(o.x + w, o.y, o.z),         // 1
            Point3r::new(o.x + w, o.y + d, o.z),     // 2
            Point3r::new(o.x, o.y + d, o.z),         // 3
            Point3r::new(o.x, o.y, o.z + h),         // 4
            Point3r::new(o.x + w, o.y, o.z + h),     // 5
            Point3r::new(o.x + w, o.y + d, o.z + h), // 6
            Point3r::new(o.x, o.y + d, o.z + h),     // 7
        ];

        let normals = [
            Vector3r::new(0.0, 0.0, -1.0), // bottom
            Vector3r::new(0.0, 0.0, 1.0),  // top
            Vector3r::new(0.0, -1.0, 0.0), // front
            Vector3r::new(0.0, 1.0, 0.0),  // back
            Vector3r::new(-1.0, 0.0, 0.0), // left
            Vector3r::new(1.0, 0.0, 0.0),  // right
        ];

        // Face quads (vertex indices into `corners`), each split into 2 triangles.
        // The quad [a,b,c,d] is fan-triangulated as (a,b,c) + (a,c,d).
        // Winding must be CCW when viewed from outside so both triangles agree.
        let face_quads: [([usize; 4], usize); 6] = [
            ([0, 3, 2, 1], 0), // bottom (-Z): 0→3→2→1
            ([4, 5, 6, 7], 1), // top    (+Z): 4→5→6→7
            ([0, 1, 5, 4], 2), // front  (-Y): 0→1→5→4
            ([2, 3, 7, 6], 3), // back   (+Y): 2→3→7→6
            ([0, 4, 7, 3], 4), // left   (-X): 0→4→7→3
            ([1, 2, 6, 5], 5), // right  (+X): 1→2→6→5
        ];

        let mut faces = Vec::with_capacity(12);

        for (quad, normal_idx) in &face_quads {
            let n = normals[*normal_idx];
            let vids: Vec<_> = quad
                .iter()
                .map(|&ci| vertex_pool.insert_or_weld(corners[ci], n))
                .collect();

            faces.push(FaceData {
                vertices: [vids[0], vids[1], vids[2]],
                region,
            });
            faces.push(FaceData {
                vertices: [vids[0], vids[2], vids[3]],
                region,
            });
        }

        faces
    }
}
