//! Junction geometry for channel intersections.
//!
//! Handles T-junctions, Y-junctions, and cross-junctions common in
//! millifluidic chip designs.

use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Type of junction.
#[derive(Clone, Debug)]
pub enum JunctionType {
    /// T-junction: one channel meets another at 90°.
    Tee {
        /// Radius of the main channel.
        main_radius: Real,
        /// Radius of the branch channel.
        branch_radius: Real,
    },
    /// Y-junction: two channels merge into one.
    Wye {
        /// Radius of the inlet channels.
        inlet_radius: Real,
        /// Radius of the outlet channel.
        outlet_radius: Real,
        /// Angle between the two inlet arms (radians).
        angle: Real,
    },
    /// Cross-junction: four channels meeting at a point.
    Cross {
        /// Radius of all channels.
        radius: Real,
    },
}

impl JunctionType {
    /// Generate junction mesh faces.
    ///
    /// `center`: position of the junction center.
    /// `direction`: the primary flow direction.
    ///
    /// Returns the generated faces. This is a stub that creates a sphere-like
    /// junction body; full implementation will blend channel profiles smoothly.
    pub fn generate(
        &self,
        center: &Point3r,
        _direction: &Vector3r,
        vertex_pool: &mut VertexPool,
        region: RegionId,
    ) -> Vec<FaceData> {
        let radius = match self {
            JunctionType::Tee { main_radius, .. } => *main_radius,
            JunctionType::Wye { outlet_radius, .. } => *outlet_radius,
            JunctionType::Cross { radius } => *radius,
        };

        // Generate a simple sphere approximation at the junction
        generate_icosphere_faces(center, radius, 1, vertex_pool, region)
    }
}

/// Generate an icosphere (subdivision level `depth`) at `center` with `radius`.
fn generate_icosphere_faces(
    center: &Point3r,
    radius: Real,
    _depth: usize,
    pool: &mut VertexPool,
    region: RegionId,
) -> Vec<FaceData> {
    // Icosahedron base vertices
    let phi: Real = (1.0 + (5.0 as Real).sqrt()) / 2.0;
    let base_verts = [
        Vector3r::new(-1.0, phi, 0.0).normalize(),
        Vector3r::new(1.0, phi, 0.0).normalize(),
        Vector3r::new(-1.0, -phi, 0.0).normalize(),
        Vector3r::new(1.0, -phi, 0.0).normalize(),
        Vector3r::new(0.0, -1.0, phi).normalize(),
        Vector3r::new(0.0, 1.0, phi).normalize(),
        Vector3r::new(0.0, -1.0, -phi).normalize(),
        Vector3r::new(0.0, 1.0, -phi).normalize(),
        Vector3r::new(phi, 0.0, -1.0).normalize(),
        Vector3r::new(phi, 0.0, 1.0).normalize(),
        Vector3r::new(-phi, 0.0, -1.0).normalize(),
        Vector3r::new(-phi, 0.0, 1.0).normalize(),
    ];

    #[rustfmt::skip]
    let base_faces: [[usize; 3]; 20] = [
        [0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],
        [1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],
        [3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],
        [4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1],
    ];

    let vids: Vec<_> = base_verts
        .iter()
        .map(|n| {
            let pos = Point3r::from(center.coords + *n * radius);
            pool.insert_or_weld(pos, *n)
        })
        .collect();

    base_faces
        .iter()
        .map(|f| FaceData {
            vertices: [vids[f[0]], vids[f[1]], vids[f[2]]],
            region,
        })
        .collect()
}
