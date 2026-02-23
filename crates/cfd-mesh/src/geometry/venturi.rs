//! Venturi tube mesh builder.
//!
//! Builds a structured mesh for a Venturi flow passage.
//! Use [`VenturiMeshBuilder::build_surface`] for the [`IndexedMesh`]
//! boundary-surface output.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::mesh::IndexedMesh;
use crate::core::index::RegionId;
use crate::core::scalar::{Real, Point3r, Vector3r};

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

    /// Build a watertight surface mesh (wall + inlet + outlet caps).
    ///
    /// Returns an [`IndexedMesh`] with three named regions:
    /// - `RegionId(0)` — outer wall
    /// - `RegionId(1)` — inlet cap
    /// - `RegionId(2)` — outlet cap
    pub fn build_surface(&self) -> Result<IndexedMesh, BuildError> {
        build_venturi_surface(self)
    }

    /// Alias for [`build_surface`] to support legacy API calls.
    pub fn build(&self) -> Result<IndexedMesh, BuildError> {
        self.build_surface()
    }
}

// ---------------------------------------------------------------------------
// Internal mesh generation
// ---------------------------------------------------------------------------

fn build_venturi_surface<T: Copy + RealField + Float + FromPrimitive + ToPrimitive>(
    b: &VenturiMeshBuilder<T>,
) -> Result<IndexedMesh, BuildError> {
    let f = |v: T| v.to_f64().ok_or_else(|| BuildError("float conv".into()));

    let d_inlet = f(b.d_inlet)?;
    let d_throat = f(b.d_throat)?;
    let l_inlet = f(b.l_inlet)?;
    let l_convergent = f(b.l_convergent)?;
    let l_throat = f(b.l_throat)?;
    let l_divergent = f(b.l_divergent)?;
    let l_outlet = f(b.l_outlet)?;

    let nx = b.resolution_x.max(2);
    let n_ang: usize = if b.circular { b.resolution_y.max(2) * 4 } else { 4 };
    let total_l = l_inlet + l_convergent + l_throat + l_divergent + l_outlet;

    let wall_region = RegionId::from_usize(0);
    let inlet_region = RegionId::from_usize(1);
    let outlet_region = RegionId::from_usize(2);

    // Radius at axial position z (all in f64).
    let radius_at_f64 = |z: Real| -> Real {
        let r_in = d_inlet / 2.0;
        let r_th = d_throat / 2.0;
        let z1 = l_inlet;
        let z2 = z1 + l_convergent;
        let z3 = z2 + l_throat;
        let z4 = z3 + l_divergent;
        if z <= z1 { r_in }
        else if z <= z2 { let t = (z - z1) / l_convergent; r_in + (r_th - r_in) * t }
        else if z <= z3 { r_th }
        else if z <= z4 { let t = (z - z3) / l_divergent; r_th + (r_in - r_th) * t }
        else { r_in }
    };

    let mut mesh = IndexedMesh::new();

    // Build rings of vertices (no center node needed for surface mesh).
    let mut rings: Vec<Vec<crate::core::index::VertexId>> = Vec::with_capacity(nx);
    for i in 0..nx {
        let t = i as Real / (nx - 1) as Real;
        let z = total_l * t;
        let r = radius_at_f64(z);
        let mut ring = Vec::with_capacity(n_ang);
        for ia in 0..n_ang {
            let theta = std::f64::consts::TAU * ia as Real / n_ang as Real;
            let (sin_t, cos_t) = theta.sin_cos();
            let vid = mesh.add_vertex(
                Point3r::new(r * cos_t, r * sin_t, z),
                Vector3r::new(cos_t, sin_t, 0.0),
            );
            ring.push(vid);
        }
        rings.push(ring);
    }

    // Wall: quad strip between adjacent rings → 2 outward-facing triangles per quad.
    //
    // Ring vertices lie in the XY plane (x = r·cosθ, y = r·sinθ, z = axial).
    // Viewed from outside (+radial), the outward CCW quad is:
    //   v00 → v01 → v11 → v10   (angle increases first, then axially back)
    // giving cross product (angular) × (−axial) = outward radial. ✓
    for iz in 0..(nx - 1) {
        for ia in 0..n_ang {
            let ia1 = (ia + 1) % n_ang;
            let v00 = rings[iz][ia];       // (θ, z_iz)
            let v01 = rings[iz][ia1];      // (θ+1, z_iz)
            let v10 = rings[iz + 1][ia];   // (θ, z_{iz+1})
            let v11 = rings[iz + 1][ia1];  // (θ+1, z_{iz+1})
            // CCW from outside: v00 → v01 → v11, then v00 → v11 → v10
            mesh.add_face_with_region(v00, v01, v11, wall_region);
            mesh.add_face_with_region(v00, v11, v10, wall_region);
        }
    }

    // Inlet cap at z = 0 (outward normal = −z).
    // Wall bottom edge runs rings[0][ia] → rings[0][ia1] (positive-θ).
    // Cap must reverse that shared edge: rings[0][ia1] → rings[0][ia].
    // Fan: ic → rings[0][ia1] → rings[0][ia]  →  normal = −z ✓
    let ic = mesh.add_vertex(Point3r::new(0.0, 0.0, 0.0), Vector3r::new(0.0, 0.0, -1.0));
    for ia in 0..n_ang {
        let ia1 = (ia + 1) % n_ang;
        mesh.add_face_with_region(ic, rings[0][ia1], rings[0][ia], inlet_region);
    }

    // Outlet cap at z = total_l (outward normal = +z).
    // Wall top edge at iz=nx-2 runs rings[last][ia1] → rings[last][ia] (negative-θ).
    // Cap must reverse that shared edge: rings[last][ia] → rings[last][ia1].
    // Fan: oc → rings[last][ia] → rings[last][ia1]  →  normal = +z ✓
    let oc = mesh.add_vertex(Point3r::new(0.0, 0.0, total_l), Vector3r::new(0.0, 0.0, 1.0));
    let last = nx - 1;
    for ia in 0..n_ang {
        let ia1 = (ia + 1) % n_ang;
        mesh.add_face_with_region(oc, rings[last][ia], rings[last][ia1], outlet_region);
    }

    Ok(mesh)
}

