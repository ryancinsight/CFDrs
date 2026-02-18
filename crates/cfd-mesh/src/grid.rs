//! Structured grid builder.
//!
//! Generates a regular Cartesian grid over the unit cube [0,1]³,
//! subdivided into nx×ny×nz hexahedra (each decomposed to 5 tetrahedra).

use crate::mesh::Mesh;
use crate::topology::{Cell, Face, Vertex};
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
    pub fn build(self) -> Result<Mesh<f64>, GridError> {
        build_structured_grid(self.nx, self.ny, self.nz)
    }
}

fn build_structured_grid(nx: usize, ny: usize, nz: usize) -> Result<Mesh<f64>, GridError> {
    let nx = nx.max(1);
    let ny = ny.max(1);
    let nz = nz.max(1);

    let vnx = nx + 1;
    let vny = ny + 1;
    let vnz = nz + 1;

    let mut mesh = Mesh::<f64>::new();

    // Create corner vertices on a regular grid.
    for iz in 0..vnz {
        for iy in 0..vny {
            for ix in 0..vnx {
                let x = ix as f64 / nx as f64;
                let y = iy as f64 / ny as f64;
                let z = iz as f64 / nz as f64;
                mesh.add_vertex(Vertex::new(Point3::new(x, y, z)));
            }
        }
    }

    let v_idx = |ix: usize, iy: usize, iz: usize| iz * vny * vnx + iy * vnx + ix;

    // Create cells: each hex cell is split into 5 tetrahedra.
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                // 8 corner indices of the hex cell.
                let v: [usize; 8] = [
                    v_idx(ix,   iy,   iz  ),
                    v_idx(ix+1, iy,   iz  ),
                    v_idx(ix+1, iy+1, iz  ),
                    v_idx(ix,   iy+1, iz  ),
                    v_idx(ix,   iy,   iz+1),
                    v_idx(ix+1, iy,   iz+1),
                    v_idx(ix+1, iy+1, iz+1),
                    v_idx(ix,   iy+1, iz+1),
                ];

                // Standard 5-tet decomposition of a hex.
                let tets: [[usize; 4]; 5] = [
                    [v[0], v[1], v[3], v[4]],
                    [v[1], v[4], v[5], v[6]],
                    [v[1], v[3], v[6], v[7]],  // corrected
                    [v[3], v[4], v[6], v[7]],
                    [v[1], v[3], v[4], v[6]],
                ];

                for tet in &tets {
                    let f0 = mesh.add_face(Face::triangle(tet[0], tet[1], tet[2]));
                    let f1 = mesh.add_face(Face::triangle(tet[0], tet[1], tet[3]));
                    let f2 = mesh.add_face(Face::triangle(tet[0], tet[2], tet[3]));
                    let f3 = mesh.add_face(Face::triangle(tet[1], tet[2], tet[3]));
                    mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
                }
            }
        }
    }

    // Label boundary faces.
    let n_faces_before_labeling = mesh.face_count();
    for f_idx in 0..n_faces_before_labeling {
        if let Some(face) = mesh.face(f_idx) {
            let verts: Vec<_> = face.vertices.iter()
                .filter_map(|&vi| mesh.vertex(vi))
                .map(|v| v.position)
                .collect();
            if verts.is_empty() { continue; }
            let all_bottom = verts.iter().all(|p| p.z < 1e-9);
            let all_top = verts.iter().all(|p| p.z > 1.0 - 1e-9);
            let all_front = verts.iter().all(|p| p.y < 1e-9);
            let all_back = verts.iter().all(|p| p.y > 1.0 - 1e-9);
            let all_left = verts.iter().all(|p| p.x < 1e-9);
            let all_right = verts.iter().all(|p| p.x > 1.0 - 1e-9);
            if all_bottom { mesh.mark_boundary(f_idx, "inlet".to_string()); }
            else if all_top { mesh.mark_boundary(f_idx, "outlet".to_string()); }
            else if all_front || all_back || all_left || all_right {
                mesh.mark_boundary(f_idx, "wall".to_string());
            }
        }
    }

    Ok(mesh)
}
