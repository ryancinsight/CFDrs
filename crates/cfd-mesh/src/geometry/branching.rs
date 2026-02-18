//! Branching (bifurcation / trifurcation) mesh builder.
//!
//! Builds a structured mesh for a Y-shaped or T-shaped branching passage.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

use crate::mesh::Mesh;
use crate::topology::{Cell, Face, Vertex};
use crate::geometry::venturi::BuildError;
use nalgebra::Point3;

/// Builds a branching (bifurcation) flow passage mesh.
///
/// All length/geometry parameters are in metres.
#[derive(Clone, Debug)]
pub struct BranchingMeshBuilder<T: Copy + RealField> {
    d_parent: T,
    l_parent: T,
    d_daughter: T,
    l_daughter: T,
    branching_angle: T,
    resolution: usize,
    n_daughters: usize,
}

impl<T: Copy + RealField + Float + FromPrimitive> BranchingMeshBuilder<T> {
    /// Create a symmetric bifurcation (1 parent, 2 daughters).
    pub fn bifurcation(
        d_parent: T,
        l_parent: T,
        d_daughter: T,
        l_daughter: T,
        branching_angle: T,
        resolution: usize,
    ) -> Self {
        Self {
            d_parent,
            l_parent,
            d_daughter,
            l_daughter,
            branching_angle,
            resolution,
            n_daughters: 2,
        }
    }

    /// Create a symmetric trifurcation (1 parent, 3 daughters).
    pub fn trifurcation(
        d_parent: T,
        l_parent: T,
        d_daughter: T,
        l_daughter: T,
        branching_angle: T,
        resolution: usize,
    ) -> Self {
        Self {
            d_parent,
            l_parent,
            d_daughter,
            l_daughter,
            branching_angle,
            resolution,
            n_daughters: 3,
        }
    }

    /// Build the mesh.
    pub fn build(self) -> Result<Mesh<T>, BuildError> {
        build_branching_mesh(self)
    }
}

/// Decompose a triangular prism defined by two triangles (b0,a0,a1) and (b1,c0,c1)
/// into 3 tetrahedra and add them to the mesh.
fn add_prism_tets<T: Copy + RealField>(
    mesh: &mut Mesh<T>,
    b0: usize, a0: usize, a1: usize,
    b1: usize, c0: usize, c1: usize,
) {
    // Tet 1: [b0, a0, a1, b1]
    let f0 = mesh.add_face(Face::triangle(b0, a0, a1));
    let f1 = mesh.add_face(Face::triangle(b0, a0, b1));
    let f2 = mesh.add_face(Face::triangle(b0, a1, b1));
    let f3 = mesh.add_face(Face::triangle(a0, a1, b1));
    mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
    // Tet 2: [a0, a1, b1, c0]
    let g0 = mesh.add_face(Face::triangle(a0, a1, b1));
    let g1 = mesh.add_face(Face::triangle(a0, a1, c0));
    let g2 = mesh.add_face(Face::triangle(a0, b1, c0));
    let g3 = mesh.add_face(Face::triangle(a1, b1, c0));
    mesh.add_cell(Cell::tetrahedron(g0, g1, g2, g3));
    // Tet 3: [a1, b1, c0, c1]
    let h0 = mesh.add_face(Face::triangle(a1, b1, c0));
    let h1 = mesh.add_face(Face::triangle(a1, b1, c1));
    let h2 = mesh.add_face(Face::triangle(a1, c0, c1));
    let h3 = mesh.add_face(Face::triangle(b1, c0, c1));
    mesh.add_cell(Cell::tetrahedron(h0, h1, h2, h3));
}

fn build_branching_mesh<T: Copy + RealField + Float + FromPrimitive>(
    b: BranchingMeshBuilder<T>,
) -> Result<Mesh<T>, BuildError> {
    let n_ax = b.resolution.max(4);
    let n_ang = 8_usize;
    let two_pi = T::from_f64(std::f64::consts::TAU)
        .ok_or_else(|| BuildError("float conv".into()))?;
    let r_parent = b.d_parent / T::from_f64(2.0).ok_or_else(|| BuildError("float conv".into()))?;
    let r_daughter = b.d_daughter / T::from_f64(2.0).ok_or_else(|| BuildError("float conv".into()))?;

    let mut mesh = Mesh::new();
    let verts_per_ring = n_ang + 1;

    // Build parent tube along +z.
    for iz in 0..n_ax {
        let t = T::from_usize(iz).unwrap() / T::from_usize(n_ax - 1).unwrap();
        let z = b.l_parent * t;
        mesh.add_vertex(Vertex::new(Point3::new(T::zero(), T::zero(), z)));
        for ia in 0..n_ang {
            let theta = two_pi * T::from_usize(ia).unwrap() / T::from_usize(n_ang).unwrap();
            let x = r_parent * Float::cos(theta);
            let y = r_parent * Float::sin(theta);
            mesh.add_vertex(Vertex::new(Point3::new(x, y, z)));
        }
    }

    // Parent cells — decompose each annular prism sector into 3 tetrahedra.
    for iz in 0..(n_ax - 1) {
        let base0 = iz * verts_per_ring;
        let base1 = (iz + 1) * verts_per_ring;
        for ia in 0..n_ang {
            let ia1 = (ia + 1) % n_ang;
            add_prism_tets(&mut mesh, base0, base0+1+ia, base0+1+ia1, base1, base1+1+ia, base1+1+ia1);
        }
    }

    // Inlet face
    let base0 = 0;
    for ia in 0..n_ang {
        let ia1 = (ia + 1) % n_ang;
        let fi = mesh.add_face(Face::triangle(base0+1+ia, base0+1+ia1, base0));
        mesh.mark_boundary(fi, "inlet".to_string());
    }

    // Build daughter tubes
    for d in 0..b.n_daughters {
        let angle_step = if b.n_daughters == 1 {
            T::zero()
        } else {
            b.branching_angle * T::from_f64(
                (d as f64 - (b.n_daughters - 1) as f64 / 2.0)
            ).unwrap()
        };
        let sin_a = Float::sin(angle_step);
        let cos_a = Float::cos(angle_step);

        let vertex_offset = mesh.vertex_count();
        let junction_z = b.l_parent;
        for iz in 0..n_ax {
            let t = T::from_usize(iz).unwrap() / T::from_usize(n_ax - 1).unwrap();
            let dz = b.l_daughter * t * cos_a;
            let dx = b.l_daughter * t * sin_a;
            let z = junction_z + dz;
            let x_off = dx;
            mesh.add_vertex(Vertex::new(Point3::new(x_off, T::zero(), z)));
            for ia in 0..n_ang {
                let theta = two_pi * T::from_usize(ia).unwrap() / T::from_usize(n_ang).unwrap();
                let x = x_off + r_daughter * Float::cos(theta);
                let y = r_daughter * Float::sin(theta);
                mesh.add_vertex(Vertex::new(Point3::new(x, y, z)));
            }
        }

        for iz in 0..(n_ax - 1) {
            let base0 = vertex_offset + iz * verts_per_ring;
            let base1 = vertex_offset + (iz + 1) * verts_per_ring;
            for ia in 0..n_ang {
                let ia1 = (ia + 1) % n_ang;
                add_prism_tets(&mut mesh, base0, base0+1+ia, base0+1+ia1, base1, base1+1+ia, base1+1+ia1);
            }
        }

        // Outlet
        let last_base = vertex_offset + (n_ax - 1) * verts_per_ring;
        for ia in 0..n_ang {
            let ia1 = (ia + 1) % n_ang;
            let fi = mesh.add_face(Face::triangle(last_base+1+ia, last_base+1+ia1, last_base));
            mesh.mark_boundary(fi, format!("outlet_{d}"));
        }
    }

    Ok(mesh)
}
