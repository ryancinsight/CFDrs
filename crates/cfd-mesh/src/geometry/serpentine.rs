//! Serpentine channel mesh builder.
//!
//! Builds a structured mesh for a sinuous (serpentine) microchannel.
//! Use [`SerpentineMeshBuilder::build_surface`] for the modern [`IndexedMesh`]
//! boundary-surface output.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::mesh::{IndexedMesh, Mesh};
use crate::geometry::venturi::BuildError;
use crate::core::index::RegionId;
use crate::core::scalar::{Real, Point3r, Vector3r};

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

    /// Build a watertight surface mesh (wall + inlet + outlet caps).
    ///
    /// Returns an [`IndexedMesh`] with three named regions:
    /// - `RegionId(0)` — outer wall
    /// - `RegionId(1)` — inlet cap
    /// - `RegionId(2)` — outlet cap
    pub fn build_surface(&self) -> Result<IndexedMesh, BuildError> {
        build_serpentine_surface(self)
    }

    /// Legacy build method for volumetric structured mesh (placeholder).
    /// Returns an error as this functionality has been removed in favor of `build_surface`.
    pub fn build(&self) -> Result<Mesh<T>, BuildError> {
        Err(BuildError("Legacy build() removed. Use build_surface()".into()))
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
