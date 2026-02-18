//! Serpentine channel mesh builder.
//!
//! Builds a structured mesh for a sinuous (serpentine) microchannel.
//! Use [`SerpentineMeshBuilder::build_surface`] for the modern [`IndexedMesh`]
//! boundary-surface output.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::mesh::{Mesh, IndexedMesh};
use crate::topology::{Cell, Face, Vertex};
use crate::geometry::venturi::BuildError;
use crate::core::index::RegionId;
use crate::core::scalar::{Real, Point3r, Vector3r};
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
    #[deprecated(note = "use `build_surface()` which returns an `IndexedMesh` boundary surface")]
    pub fn build(self) -> Result<Mesh<T>, BuildError> {
        build_serpentine_mesh(self)
    }

    /// Build a watertight surface mesh (wall + inlet + outlet caps).
    ///
    /// Returns an [`IndexedMesh`] with three named regions:
    /// - `RegionId(0)` — outer wall
    /// - `RegionId(1)` — inlet cap
    /// - `RegionId(2)` — outlet cap
    pub fn build_surface(&self) -> Result<IndexedMesh, BuildError> {
        build_serpentine_surface(self)
    }
}

fn build_serpentine_surface<T: Copy + RealField + Float + FromPrimitive + ToPrimitive>(
    b: &SerpentineMeshBuilder<T>,
) -> Result<IndexedMesh, BuildError> {
    let f = |v: T| v.to_f64().ok_or_else(|| BuildError("float conv".into()));

    let diameter = f(b.diameter)?;
    let amplitude = f(b.amplitude)?;
    let wavelength = f(b.wavelength)?;

    let r = diameter / 2.0;
    let n_ax = b.resolution_x.max(4) * b.num_periods;
    let n_ang: usize = (b.resolution_y.max(2) * 4).max(4);
    let total_len = wavelength * b.num_periods as Real;

    let wall_region = RegionId::from_usize(0);
    let inlet_region = RegionId::from_usize(1);
    let outlet_region = RegionId::from_usize(2);

    let mut mesh = IndexedMesh::new();

    // Spine: y(z) = amplitude * sin(2π * z / wavelength).
    let spine: Vec<(Real, Real, Real)> = (0..n_ax)
        .map(|i| {
            let t = i as Real / (n_ax - 1) as Real;
            let z = total_len * t;
            let y = amplitude * (std::f64::consts::TAU * z / wavelength).sin();
            (0.0, y, z)
        })
        .collect();

    // Build rings of vertices along spine.
    let mut rings: Vec<Vec<crate::core::index::VertexId>> = Vec::with_capacity(n_ax);
    for &(cx, cy, cz) in &spine {
        let mut ring = Vec::with_capacity(n_ang);
        for ia in 0..n_ang {
            let theta = std::f64::consts::TAU * ia as Real / n_ang as Real;
            let (sin_t, cos_t) = theta.sin_cos();
            let vid = mesh.add_vertex(
                Point3r::new(cx + r * cos_t, cy + r * sin_t, cz),
                Vector3r::new(cos_t, sin_t, 0.0),
            );
            ring.push(vid);
        }
        rings.push(ring);
    }

    // Wall: quad strip between adjacent rings.
    for iz in 0..(n_ax - 1) {
        for ia in 0..n_ang {
            let ia1 = (ia + 1) % n_ang;
            let v00 = rings[iz][ia];
            let v01 = rings[iz][ia1];
            let v10 = rings[iz + 1][ia];
            let v11 = rings[iz + 1][ia1];
            mesh.add_face_with_region(v00, v10, v01, wall_region);
            mesh.add_face_with_region(v01, v10, v11, wall_region);
        }
    }

    // Inlet cap at first spine point.
    let (icx, icy, icz) = spine[0];
    let ic = mesh.add_vertex(Point3r::new(icx, icy, icz), Vector3r::new(0.0, 0.0, -1.0));
    for ia in 0..n_ang {
        let ia1 = (ia + 1) % n_ang;
        mesh.add_face_with_region(ic, rings[0][ia1], rings[0][ia], inlet_region);
    }

    // Outlet cap at last spine point.
    let last = n_ax - 1;
    let (ocx, ocy, ocz) = spine[last];
    let oc = mesh.add_vertex(Point3r::new(ocx, ocy, ocz), Vector3r::new(0.0, 0.0, 1.0));
    for ia in 0..n_ang {
        let ia1 = (ia + 1) % n_ang;
        mesh.add_face_with_region(oc, rings[last][ia], rings[last][ia1], outlet_region);
    }

    Ok(mesh)
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
