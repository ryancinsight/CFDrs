//! Branching (bifurcation / trifurcation) mesh builder.
//!
//! Builds a structured mesh for a Y-shaped or T-shaped branching passage.
//! Use [`BranchingMeshBuilder::build_surface`] for the modern [`IndexedMesh`]
//! boundary-surface output.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::mesh::{IndexedMesh, Mesh};
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

    /// Build a watertight surface mesh (parent + daughter walls, inlet, and outlet caps).
    ///
    /// Region IDs:
    /// - `RegionId(0)` — wall (all tube surfaces)
    /// - `RegionId(1)` — inlet cap (parent inlet)
    /// - `RegionId(2+d)` — outlet cap for daughter `d`
    pub fn build_surface(&self) -> Result<IndexedMesh, BuildError> {
        build_branching_surface(self)
    }

    /// Alias for [`build_surface`](Self::build_surface).
    ///
    /// # Deprecated
    /// This method converts the [`IndexedMesh`] to the legacy [`Mesh<T>`] type for
    /// compatibility with `cfd-3d` volume meshing tools.
    #[allow(deprecated)]
    pub fn build(&self) -> Result<Mesh<T>, BuildError>
    where
        T: 'static,
    {
        self.build_surface().map(|m| {
            let m_f64: Mesh<Real> = m.into();
            m_f64.cast::<T>()
        })
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

