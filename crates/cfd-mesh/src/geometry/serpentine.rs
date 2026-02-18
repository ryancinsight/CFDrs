//! Serpentine channel mesh builder.
//!
//! Builds a structured mesh for a sinuous (serpentine) microchannel.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

use crate::mesh::Mesh;
use crate::topology::{Cell, Face, Vertex};
use crate::geometry::venturi::BuildError;
use nalgebra::Point3;

/// Builds a serpentine channel mesh.
///
/// All length/geometry parameters are in metres.
#[derive(Clone, Debug)]
pub struct SerpentineMeshBuilder<T: Copy + RealField> {
    /// Channel diameter (or hydraulic diameter) [m].
    pub diameter: T,
    /// Amplitude of the sinusoidal path [m].
    pub amplitude: T,
    /// Wavelength (period) of the channel [m].
    pub wavelength: T,
    /// Number of full periods.
    pub num_periods: usize,
    resolution_x: usize,
    resolution_y: usize,
    circular: bool,
}

impl<T: Copy + RealField + Float + FromPrimitive> SerpentineMeshBuilder<T> {
    /// Create a serpentine builder with given diameter, amplitude, wavelength.
    pub fn new(diameter: T, amplitude: T, wavelength: T) -> Self {
        Self {
            diameter,
            amplitude,
            wavelength,
            num_periods: 3,
            resolution_x: 32,
            resolution_y: 4,
            circular: true,
        }
    }

    /// Set the number of full sinusoidal periods.
    pub fn with_periods(mut self, periods: usize) -> Self {
        self.num_periods = periods;
        self
    }

    /// Set the mesh resolution (axial × radial).
    pub fn with_resolution(mut self, x: usize, y: usize) -> Self {
        self.resolution_x = x;
        self.resolution_y = y;
        self
    }

    /// Use a circular cross-section.
    pub fn with_circular(mut self, circular: bool) -> Self {
        self.circular = circular;
        self
    }

    /// Build the mesh.
    pub fn build(self) -> Result<Mesh<T>, BuildError> {
        build_serpentine_mesh(self)
    }
}

fn build_serpentine_mesh<T: Copy + RealField + Float + FromPrimitive>(
    b: SerpentineMeshBuilder<T>,
) -> Result<Mesh<T>, BuildError> {
    let n_ax = b.resolution_x.max(4) * b.num_periods;
    let n_ang = (b.resolution_y.max(2) * 4).max(4);
    let two_pi = T::from_f64(std::f64::consts::TAU)
        .ok_or_else(|| BuildError("float conv".into()))?;
    let r = b.diameter / T::from_f64(2.0).ok_or_else(|| BuildError("float conv".into()))?;
    let total_len = b.wavelength * T::from_usize(b.num_periods)
        .ok_or_else(|| BuildError("float conv".into()))?;

    let mut mesh = Mesh::new();
    let verts_per_ring = n_ang + 1;

    // Spine path: y(t) = amplitude * sin(2π * t / wavelength)
    let mut centres: Vec<(T, T, T)> = Vec::with_capacity(n_ax);
    for i in 0..n_ax {
        let t = T::from_usize(i).unwrap() / T::from_usize(n_ax - 1).unwrap();
        let z = total_len * t;
        let y = b.amplitude * Float::sin(two_pi * z / b.wavelength);
        centres.push((T::zero(), y, z));
    }

    // Vertices
    for (cx, cy, cz) in &centres {
        mesh.add_vertex(Vertex::new(Point3::new(*cx, *cy, *cz)));
        for ia in 0..n_ang {
            let theta = two_pi * T::from_usize(ia).unwrap() / T::from_usize(n_ang).unwrap();
            let x = *cx + r * Float::cos(theta);
            let y = *cy + r * Float::sin(theta);
            mesh.add_vertex(Vertex::new(Point3::new(x, y, *cz)));
        }
    }

    // Cells
    for iz in 0..(n_ax - 1) {
        let base0 = iz * verts_per_ring;
        let base1 = (iz + 1) * verts_per_ring;
        for ia in 0..n_ang {
            let ia1 = (ia + 1) % n_ang;
            let v00 = base0 + 1 + ia;
            let v01 = base0 + 1 + ia1;
            let v10 = base1 + 1 + ia;
            let v11 = base1 + 1 + ia1;
            let c0 = base0;
            let c1 = base1;
            let f0 = mesh.add_face(Face::triangle(v00, v01, c0));
            let f1 = mesh.add_face(Face::triangle(v10, v11, c1));
            let f2 = mesh.add_face(Face::triangle(v00, v10, v01));
            let f3 = mesh.add_face(Face::triangle(v01, v10, v11));
            mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
        }
        if iz == 0 {
            for ia in 0..n_ang {
                let ia1 = (ia + 1) % n_ang;
                let fi = mesh.add_face(Face::triangle(base0 + 1 + ia, base0 + 1 + ia1, base0));
                mesh.mark_boundary(fi, "inlet".to_string());
            }
        }
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
