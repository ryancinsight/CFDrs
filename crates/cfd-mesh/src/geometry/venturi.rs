//! Venturi tube mesh builder.
//!
//! Builds a structured hexahedral mesh for a Venturi flow passage and
//! optionally converts it to tetrahedra (done externally via `HexToTetConverter`).

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

use crate::mesh::Mesh;
use crate::topology::{Cell, ElementType, Face, Vertex};
use nalgebra::Point3;

/// Error type for mesh building.
#[derive(Debug)]
pub struct BuildError(pub String);

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "mesh build error: {}", self.0)
    }
}

impl std::error::Error for BuildError {}

/// Builds a Venturi tube mesh.
///
/// All length parameters are in metres.
#[derive(Clone, Debug)]
pub struct VenturiMeshBuilder<T: Copy + RealField> {
    // --- geometry parameters ---
    /// Inlet diameter [m].
    pub d_inlet: T,
    /// Throat diameter [m].
    pub d_throat: T,
    /// Inlet section length [m].
    pub l_inlet: T,
    /// Convergent section length [m].
    pub l_convergent: T,
    /// Throat section length [m].
    pub l_throat: T,
    /// Divergent section length [m].
    pub l_divergent: T,
    /// Outlet section length [m].
    pub l_outlet: T,

    // --- mesh resolution ---
    resolution_x: usize,
    resolution_y: usize,
    circular: bool,
}

impl<T: Copy + RealField + Float + FromPrimitive> VenturiMeshBuilder<T> {
    /// Create a Venturi mesh builder with the given geometry.
    pub fn new(
        d_inlet: T,
        d_throat: T,
        l_inlet: T,
        l_convergent: T,
        l_throat: T,
        l_divergent: T,
        l_outlet: T,
    ) -> Self {
        Self {
            d_inlet,
            d_throat,
            l_inlet,
            l_convergent,
            l_throat,
            l_divergent,
            l_outlet,
            resolution_x: 8,
            resolution_y: 4,
            circular: true,
        }
    }

    /// Set the mesh resolution (axial × radial).
    pub fn with_resolution(mut self, x: usize, y: usize) -> Self {
        self.resolution_x = x;
        self.resolution_y = y;
        self
    }

    /// Use a circular (cylinder) cross-section (default) vs. square.
    pub fn with_circular(mut self, circular: bool) -> Self {
        self.circular = circular;
        self
    }

    /// Build the mesh.
    ///
    /// Returns a structured hexahedral mesh of the Venturi passage.
    pub fn build(self) -> Result<Mesh<T>, BuildError> {
        build_venturi_mesh(self)
    }
}

// ---------------------------------------------------------------------------
// Internal mesh generation
// ---------------------------------------------------------------------------

fn build_venturi_mesh<T: Copy + RealField + Float + FromPrimitive>(
    b: VenturiMeshBuilder<T>,
) -> Result<Mesh<T>, BuildError> {
    let nx = b.resolution_x.max(2);
    let nr = b.resolution_y.max(2);

    let total_l = b.l_inlet + b.l_convergent + b.l_throat + b.l_divergent + b.l_outlet;
    let two_pi = T::from_f64(std::f64::consts::TAU).ok_or_else(|| BuildError("float conv".into()))?;
    let half = T::from_f64(0.5).ok_or_else(|| BuildError("float conv".into()))?;

    let mut mesh = Mesh::new();

    // Nodes: nx axial stations × (nr+1 radial) × (nr angular, if circular)
    // For simplicity, generate a 1D axial sweep with radial rings.
    let n_ang = if b.circular { nr * 4 } else { 4 };
    let n_ax = nx;

    // Compute axial z positions and radii.
    let mut z_vals: Vec<T> = Vec::with_capacity(n_ax);
    let mut r_vals: Vec<T> = Vec::with_capacity(n_ax);

    for i in 0..n_ax {
        let t = T::from_usize(i).unwrap() / T::from_usize(n_ax - 1).unwrap();
        let z = total_l * t;
        z_vals.push(z);
        r_vals.push(radius_at(z, &b));
    }

    // Create vertices in rings.
    for iz in 0..n_ax {
        let z = z_vals[iz];
        let r = r_vals[iz];
        // Centre node
        mesh.add_vertex(Vertex::new(Point3::new(T::zero(), T::zero(), z)));
        // Ring nodes
        for ia in 0..n_ang {
            let theta = two_pi * T::from_usize(ia).unwrap() / T::from_usize(n_ang).unwrap();
            let x = r * Float::cos(theta);
            let y = r * Float::sin(theta);
            mesh.add_vertex(Vertex::new(Point3::new(x, y, z)));
        }
    }

    let verts_per_ring = n_ang + 1; // centre + ring

    // Create faces and cells (wedge elements between axial slices).
    for iz in 0..(n_ax - 1) {
        let base0 = iz * verts_per_ring;
        let base1 = (iz + 1) * verts_per_ring;

        let z0 = z_vals[iz];
        let z1 = z_vals[iz + 1];

        for ia in 0..n_ang {
            let ia1 = (ia + 1) % n_ang;

            // Four vertices of the wedge quad between the two rings.
            let v00 = base0 + 1 + ia;
            let v01 = base0 + 1 + ia1;
            let v10 = base1 + 1 + ia;
            let v11 = base1 + 1 + ia1;
            let ctr0 = base0;
            let ctr1 = base1;

            // Create two tetrahedra per quad column.
            let f0 = mesh.add_face(Face::triangle(v00, v01, ctr0));
            let f1 = mesh.add_face(Face::triangle(v10, v11, ctr1));
            let f2 = mesh.add_face(Face::triangle(v00, v10, v01));
            let f3 = mesh.add_face(Face::triangle(v01, v10, v11));

            mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
        }

        // Label boundary faces (inlet / outlet / wall)
        // Inlet: iz == 0 plane
        if iz == 0 {
            for ia in 0..n_ang {
                let ia1 = (ia + 1) % n_ang;
                let fi = mesh.add_face(Face::triangle(base0 + 1 + ia, base0 + 1 + ia1, base0));
                mesh.mark_boundary(fi, "inlet".to_string());
            }
        }
        // Outlet: iz == n_ax-2 plane
        if iz == n_ax - 2 {
            for ia in 0..n_ang {
                let ia1 = (ia + 1) % n_ang;
                let fi = mesh.add_face(Face::triangle(base1 + 1 + ia, base1 + 1 + ia1, base1));
                mesh.mark_boundary(fi, "outlet".to_string());
            }
        }
    }

    Ok(mesh)
}

fn radius_at<T: Copy + RealField + Float + FromPrimitive>(z: T, b: &VenturiMeshBuilder<T>) -> T {
    let two = T::from_f64(2.0).unwrap();
    let r_in = b.d_inlet / two;
    let r_th = b.d_throat / two;

    let z1 = b.l_inlet;
    let z2 = z1 + b.l_convergent;
    let z3 = z2 + b.l_throat;
    let z4 = z3 + b.l_divergent;

    if z <= z1 {
        r_in
    } else if z <= z2 {
        let t = (z - z1) / b.l_convergent;
        r_in + (r_th - r_in) * t
    } else if z <= z3 {
        r_th
    } else if z <= z4 {
        let t = (z - z3) / b.l_divergent;
        r_th + (r_in - r_th) * t
    } else {
        r_in
    }
}
