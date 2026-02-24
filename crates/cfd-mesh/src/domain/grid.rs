//! Structured grid builder.
//!
//! Generates a regular Cartesian grid over the unit cube [0,1]³,
//! subdivided into nx×ny×nz hexahedra (each decomposed to 5 tetrahedra).
//!
//! This module is a **volume/FEM tool** — it intentionally uses `Mesh<T>` for
//! hexahedral cell topology and is exempt from the surface-mesh deprecation.

// Volume tool: Mesh<T> is the correct type here.
#![allow(deprecated)]

use crate::domain::mesh::IndexedMesh;
use crate::domain::topology::Cell;
use crate::domain::core::index::VertexId;
use nalgebra::Point3;

/// Error type for grid building.
#[derive(Debug)]
pub struct GridError(pub String);

impl std::fmt::Display for GridError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "grid error: {}", self.0)
    }
}

impl std::error::Error for GridError {}

/// Builds a structured hexahedral grid over the unit cube.
///
/// `nx`, `ny`, `nz` are the number of *cells* (not nodes) along each axis.
pub struct StructuredGridBuilder {
    nx: usize,
    ny: usize,
    nz: usize,
}

impl StructuredGridBuilder {
    /// Create a builder with `nx × ny × nz` cells.
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self { nx, ny, nz }
    }

    /// Build the mesh.
    pub fn build(self) -> Result<IndexedMesh<f64>, GridError> {
        build_structured_grid(self.nx, self.ny, self.nz)
    }
}

fn build_structured_grid(nx: usize, ny: usize, nz: usize) -> Result<IndexedMesh<f64>, GridError> {
    let nx = nx.max(1);
    let ny = ny.max(1);
    let nz = nz.max(1);

    let vnx = nx + 1;
    let vny = ny + 1;
    let vnz = nz + 1;

    let mut mesh = IndexedMesh::<f64>::new();
    let mut v_ids = Vec::with_capacity(vnx * vny * vnz);
    
    // Hash map to ensure face deduplication (critical for 3D volumetrics)
    let mut face_map: std::collections::HashMap<[crate::domain::core::index::VertexId; 3], crate::domain::core::index::FaceId> = std::collections::HashMap::new();

    let mut add_tri = |v0: crate::domain::core::index::VertexId, v1: crate::domain::core::index::VertexId, v2: crate::domain::core::index::VertexId, mesh: &mut IndexedMesh<f64>| -> crate::domain::core::index::FaceId {
        let mut key = [v0, v1, v2];
        key.sort_unstable_by_key(|vid| vid.as_usize());
        *face_map.entry(key).or_insert_with(|| mesh.add_face(key[0], key[1], key[2]))
    };

    // Create corner vertices on a regular grid.
    for iz in 0..vnz {
        for iy in 0..vny {
            for ix in 0..vnx {
                let x = ix as f64 / nx as f64;
                let y = iy as f64 / ny as f64;
                let z = iz as f64 / nz as f64;
                v_ids.push(mesh.add_vertex_pos(Point3::new(x, y, z)));
            }
        }
    }

    let v_idx = |ix: usize, iy: usize, iz: usize| v_ids[iz * vny * vnx + iy * vnx + ix];

    // Create cells: each hex cell is split into 5 tetrahedra.
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                // 8 corner indices of the hex cell.
                let v: [VertexId; 8] = [
                    v_idx(ix, iy, iz),
                    v_idx(ix + 1, iy, iz),
                    v_idx(ix + 1, iy + 1, iz),
                    v_idx(ix, iy + 1, iz),
                    v_idx(ix, iy, iz + 1),
                    v_idx(ix + 1, iy, iz + 1),
                    v_idx(ix + 1, iy + 1, iz + 1),
                    v_idx(ix, iy + 1, iz + 1),
                ];

                // Alternating 5-tet decomposition to ensure conforming faces.
                if (ix + iy + iz) % 2 == 0 {
                    let tets_a: [[VertexId; 4]; 5] = [
                        [v[0], v[1], v[3], v[4]],
                        [v[1], v[2], v[3], v[6]],
                        [v[4], v[5], v[6], v[1]],
                        [v[4], v[7], v[6], v[3]],
                        [v[1], v[3], v[4], v[6]],
                    ];
                    for tet in &tets_a {
                        let f0 = add_tri(tet[0], tet[1], tet[2], &mut mesh).as_usize();
                        let f1 = add_tri(tet[0], tet[1], tet[3], &mut mesh).as_usize();
                        let f2 = add_tri(tet[0], tet[2], tet[3], &mut mesh).as_usize();
                        let f3 = add_tri(tet[1], tet[2], tet[3], &mut mesh).as_usize();
                        mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
                    }
                } else {
                    let tets_b: [[VertexId; 4]; 5] = [
                        [v[1], v[0], v[2], v[5]],
                        [v[3], v[0], v[2], v[7]],
                        [v[4], v[0], v[5], v[7]],
                        [v[6], v[2], v[5], v[7]],
                        [v[0], v[2], v[5], v[7]],
                    ];
                    for tet in &tets_b {
                        let f0 = add_tri(tet[0], tet[1], tet[2], &mut mesh).as_usize();
                        let f1 = add_tri(tet[0], tet[1], tet[3], &mut mesh).as_usize();
                        let f2 = add_tri(tet[0], tet[2], tet[3], &mut mesh).as_usize();
                        let f3 = add_tri(tet[1], tet[2], tet[3], &mut mesh).as_usize();
                        mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
                    }
                }
            }
        }
    }

    // Label boundary faces.
    let mut boundary_updates = Vec::new();
    for (f_idx, face) in mesh.faces.iter_enumerated() {
        let verts: Vec<_> = face
            .vertices
            .iter()
            .map(|&vi| mesh.vertices.position(vi))
            .collect();
        if verts.is_empty() {
            continue;
        }
        let all_bottom = verts.iter().all(|p| p.z < 1e-9);
        let all_top = verts.iter().all(|p| p.z > 1.0 - 1e-9);
        let all_front = verts.iter().all(|p| p.y < 1e-9);
        let all_back = verts.iter().all(|p| p.y > 1.0 - 1e-9);
        let all_left = verts.iter().all(|p| p.x < 1e-9);
        let all_right = verts.iter().all(|p| p.x > 1.0 - 1e-9);
        if all_bottom {
            boundary_updates.push((f_idx, "inlet"));
        } else if all_top {
            boundary_updates.push((f_idx, "outlet"));
        } else if all_front || all_back || all_left || all_right {
            boundary_updates.push((f_idx, "wall"));
        }
    }
    for (f_idx, label) in boundary_updates {
        mesh.mark_boundary(f_idx, label);
    }

    Ok(mesh)
}
