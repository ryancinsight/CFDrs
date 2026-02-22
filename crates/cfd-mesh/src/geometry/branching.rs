//! Branching (bifurcation / trifurcation) mesh builder.
//!
//! Builds a structured mesh for a Y-shaped or T-shaped branching passage.
//! Use [`BranchingMeshBuilder::build_surface`] for the modern [`IndexedMesh`]
//! boundary-surface output.

use nalgebra::RealField;
use nalgebra::Point3;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::mesh::{Mesh, IndexedMesh};
use crate::topology::{Cell, Face, Vertex};
use crate::geometry::venturi::BuildError;
use crate::core::index::RegionId;
use crate::core::scalar::{Real, Point3r, Vector3r};

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

    /// Build the mesh (volume mesh with tetrahedra).
    ///
    /// This returns a `Mesh<T>` suitable for FEM/FVM simulations.
    pub fn build(&self) -> Result<Mesh<T>, BuildError> {
        build_branching_mesh(self.clone())
    }

    /// Build a watertight surface mesh (parent + daughter walls, inlet, and outlet caps).
    ///
    /// Region IDs:
    /// - `RegionId(0)` — wall (all tube surfaces)
    /// - `RegionId(1)` — inlet cap (parent inlet)
    /// - `RegionId(2+d)` — outlet cap for daughter `d`
    pub fn build_surface(&self) -> Result<IndexedMesh, BuildError> {
        build_branching_surface(self)
    }
}

fn build_branching_surface<T: Copy + RealField + Float + FromPrimitive + ToPrimitive>(
    b: &BranchingMeshBuilder<T>,
) -> Result<IndexedMesh, BuildError> {
    let f = |v: T| v.to_f64().ok_or_else(|| BuildError("float conv".into()));

    let d_parent = f(b.d_parent)?;
    let l_parent = f(b.l_parent)?;
    let d_daughter = f(b.d_daughter)?;
    let l_daughter = f(b.l_daughter)?;
    let branching_angle = f(b.branching_angle)?;

    let r_parent = d_parent / 2.0;
    let r_daughter = d_daughter / 2.0;
    let n_ax = b.resolution.max(4);
    let n_ang: usize = 8;

    let wall_region = RegionId::from_usize(0);
    let inlet_region = RegionId::from_usize(1);

    let mut mesh = IndexedMesh::new();

    // Helper: add a tube section along an axis defined by (origin, direction, radius)
    // and return the rings of VertexIds.
    let build_tube = |mesh: &mut IndexedMesh,
                      origin: (Real, Real, Real),
                      dir: (Real, Real, Real),
                      r: Real,
                      n_steps: usize|
     -> Vec<Vec<crate::core::index::VertexId>> {
        let (ox, oy, oz) = origin;
        let (dx, dy, dz) = dir;
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        let (udx, udy, udz) = (dx / len, dy / len, dz / len);

        // Build a local frame perpendicular to the tube direction.
        // If dir is nearly along z, use x as "up"; otherwise use z.
        let (ex, ey, ez) = if udz.abs() < 0.9 {
            let (lx, ly, lz) = (0.0, 0.0, 1.0);
            let dot = udx * lx + udy * ly + udz * lz;
            let (sx, sy, sz) = (lx - dot * udx, ly - dot * udy, lz - dot * udz);
            let slen = (sx * sx + sy * sy + sz * sz).sqrt();
            (sx / slen, sy / slen, sz / slen)
        } else {
            let (lx, ly, lz) = (1.0, 0.0, 0.0);
            let dot = udx * lx + udy * ly + udz * lz;
            let (sx, sy, sz) = (lx - dot * udx, ly - dot * udy, lz - dot * udz);
            let slen = (sx * sx + sy * sy + sz * sz).sqrt();
            (sx / slen, sy / slen, sz / slen)
        };
        // Third frame axis = dir × e
        let (fx, fy, fz) = (
            udy * ez - udz * ey,
            udz * ex - udx * ez,
            udx * ey - udy * ex,
        );

        let mut rings = Vec::with_capacity(n_steps);
        for i in 0..n_steps {
            let t = i as Real / (n_steps - 1) as Real;
            let cx = ox + dx * t;
            let cy = oy + dy * t;
            let cz = oz + dz * t;
            let mut ring = Vec::with_capacity(n_ang);
            for ia in 0..n_ang {
                let theta = std::f64::consts::TAU * ia as Real / n_ang as Real;
                let (sin_t, cos_t) = theta.sin_cos();
                let nx_v = cos_t * ex + sin_t * fx;
                let ny_v = cos_t * ey + sin_t * fy;
                let nz_v = cos_t * ez + sin_t * fz;
                let vid = mesh.add_vertex(
                    Point3r::new(cx + r * nx_v, cy + r * ny_v, cz + r * nz_v),
                    Vector3r::new(nx_v, ny_v, nz_v),
                );
                ring.push(vid);
            }
            rings.push(ring);
        }
        rings
    };

    // Parent tube along +z.
    let parent_rings = build_tube(
        &mut mesh,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, l_parent),
        r_parent,
        n_ax,
    );

    // Wall faces for parent.
    for iz in 0..(n_ax - 1) {
        for ia in 0..n_ang {
            let ia1 = (ia + 1) % n_ang;
            let v00 = parent_rings[iz][ia];
            let v01 = parent_rings[iz][ia1];
            let v10 = parent_rings[iz + 1][ia];
            let v11 = parent_rings[iz + 1][ia1];
            mesh.add_face_with_region(v00, v10, v01, wall_region);
            mesh.add_face_with_region(v01, v10, v11, wall_region);
        }
    }

    // Inlet cap at z = 0.
    let ic = mesh.add_vertex(Point3r::new(0.0, 0.0, 0.0), Vector3r::new(0.0, 0.0, -1.0));
    for ia in 0..n_ang {
        let ia1 = (ia + 1) % n_ang;
        mesh.add_face_with_region(ic, parent_rings[0][ia1], parent_rings[0][ia], inlet_region);
    }

    // Daughter tubes.
    for d in 0..b.n_daughters {
        let angle_step = if b.n_daughters == 1 {
            0.0
        } else {
            branching_angle * (d as Real - (b.n_daughters - 1) as Real / 2.0)
        };
        let sin_a = angle_step.sin();
        let cos_a = angle_step.cos();

        let daughter_rings = build_tube(
            &mut mesh,
            (0.0, 0.0, l_parent),
            (l_daughter * sin_a, 0.0, l_daughter * cos_a),
            r_daughter,
            n_ax,
        );

        for iz in 0..(n_ax - 1) {
            for ia in 0..n_ang {
                let ia1 = (ia + 1) % n_ang;
                let v00 = daughter_rings[iz][ia];
                let v01 = daughter_rings[iz][ia1];
                let v10 = daughter_rings[iz + 1][ia];
                let v11 = daughter_rings[iz + 1][ia1];
                mesh.add_face_with_region(v00, v10, v01, wall_region);
                mesh.add_face_with_region(v01, v10, v11, wall_region);
            }
        }

        // Outlet cap for this daughter.
        let outlet_region = RegionId::from_usize(2 + d);
        let last = n_ax - 1;
        let end_x = l_daughter * sin_a;
        let end_z = l_parent + l_daughter * cos_a;
        let oc = mesh.add_vertex(
            Point3r::new(end_x, 0.0, end_z),
            Vector3r::new(sin_a, 0.0, cos_a),
        );
        for ia in 0..n_ang {
            let ia1 = (ia + 1) % n_ang;
            mesh.add_face_with_region(oc, daughter_rings[last][ia], daughter_rings[last][ia1], outlet_region);
        }
    }

    Ok(mesh)
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
                d as f64 - (b.n_daughters - 1) as f64 / 2.0
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
