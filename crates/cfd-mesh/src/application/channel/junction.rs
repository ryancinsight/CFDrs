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
    /// Returns the generated faces as an icosphere approximation centered
    /// at the junction.
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
///
/// # Geometry
///
/// Starts from a regular icosahedron (20 faces) and performs `depth` rounds
/// of Loop-like subdivision: each triangle is split into 4 by introducing
/// midpoints on each edge, re-projected onto the unit sphere. After `d`
/// subdivisions the mesh has `20 · 4^d` faces and `10 · 4^d + 2` vertices.
fn generate_icosphere_faces(
    center: &Point3r,
    radius: Real,
    depth: usize,
    pool: &mut VertexPool,
    region: RegionId,
) -> Vec<FaceData> {
    // Icosahedron base vertices
    let phi: Real = f64::midpoint(1.0, (5.0 as Real).sqrt());
    let base_verts: Vec<Vector3r> = [
        Vector3r::new(-1.0, phi, 0.0),
        Vector3r::new(1.0, phi, 0.0),
        Vector3r::new(-1.0, -phi, 0.0),
        Vector3r::new(1.0, -phi, 0.0),
        Vector3r::new(0.0, -1.0, phi),
        Vector3r::new(0.0, 1.0, phi),
        Vector3r::new(0.0, -1.0, -phi),
        Vector3r::new(0.0, 1.0, -phi),
        Vector3r::new(phi, 0.0, -1.0),
        Vector3r::new(phi, 0.0, 1.0),
        Vector3r::new(-phi, 0.0, -1.0),
        Vector3r::new(-phi, 0.0, 1.0),
    ]
    .iter()
    .map(|v| v.normalize())
    .collect();

    #[rustfmt::skip]
    let base_faces: Vec<[usize; 3]> = vec![
        [0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],
        [1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],
        [3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],
        [4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1],
    ];

    let mut normals = base_verts;
    let mut faces = base_faces;

    // Subdivide: each triangle → 4 triangles via edge midpoint splitting
    for _ in 0..depth {
        let mut next_faces = Vec::with_capacity(faces.len() * 4);
        let mut midpoint_cache = std::collections::HashMap::new();

        for tri in &faces {
            let mids: [usize; 3] = [
                get_midpoint(tri[0], tri[1], &mut normals, &mut midpoint_cache),
                get_midpoint(tri[1], tri[2], &mut normals, &mut midpoint_cache),
                get_midpoint(tri[2], tri[0], &mut normals, &mut midpoint_cache),
            ];
            next_faces.push([tri[0], mids[0], mids[2]]);
            next_faces.push([tri[1], mids[1], mids[0]]);
            next_faces.push([tri[2], mids[2], mids[1]]);
            next_faces.push([mids[0], mids[1], mids[2]]);
        }

        faces = next_faces;
    }

    // Map unit-sphere normals to world-space vertices
    let vids: Vec<_> = normals
        .iter()
        .map(|n| {
            let pos = Point3r::from(center.coords + *n * radius);
            pool.insert_or_weld(pos, *n)
        })
        .collect();

    faces
        .iter()
        .map(|f| FaceData {
            vertices: [vids[f[0]], vids[f[1]], vids[f[2]]],
            region,
        })
        .collect()
}

/// Return the index of the midpoint between vertices `a` and `b`,
/// creating and projecting it onto the unit sphere if it doesn't exist.
fn get_midpoint(
    a: usize,
    b: usize,
    normals: &mut Vec<Vector3r>,
    cache: &mut std::collections::HashMap<(usize, usize), usize>,
) -> usize {
    let key = if a < b { (a, b) } else { (b, a) };
    if let Some(&idx) = cache.get(&key) {
        return idx;
    }
    let mid = ((normals[a] + normals[b]) * 0.5).normalize();
    let idx = normals.len();
    normals.push(mid);
    cache.insert(key, idx);
    idx
}
