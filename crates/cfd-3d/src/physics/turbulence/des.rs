//! Detached-Eddy Simulation (DES) hybrid RANS-LES model.
//!
//! # Theorem — DES Length Scale (Spalart et al. 1997)
//!
//! DES modifies the RANS destruction length scale by replacing the wall
//! distance d with:
//!
//! ```text
//! l_DES = min(d, C_DES · Δ_max)
//! ```
//!
//! where Δ_max = max(dx, dy, dz) is the largest cell dimension and
//! C_DES = 0.65 (calibrated from homogeneous turbulence).
//!
//! **Physical interpretation.** Near walls (d < C_DES · Δ), l_DES = d and the
//! model reduces to RANS (Spalart-Allmaras).  Away from walls
//! (d > C_DES · Δ), l_DES = C_DES · Δ and the model operates as LES.
//! The transition is governed by Δ_max / d, providing automatic mode switching.
//!
//! **Proof sketch.** For d ≫ C_DES Δ, the SA production–destruction balance
//! reduces to ε ~ ν̃ · |S| / (C_DES Δ)², giving νₜ ~ (C_DES Δ)² |S|,
//! identical to the Smagorinsky SGS model.  (Shur et al. 1999, §2.)
//!
//! ## Algorithm
//!
//! ```text
//! For each grid point:
//!   1. Evaluate wall distance d_i (provided at construction via wall_distances).
//!   2. Compute Δ_max_i = max(dx, dy, dz) (provided at construction).
//!   3. l_DES_i = min(d_i, C_DES · Δ_max_i)
//!   4. Compute RANS νₜ using SA model with modified destruction length l_DES.
//!      (Approximated here by Smagorinsky with Δ = l_DES for tractability.)
//! ```
//!
//! ## References
//!
//! - Spalart, P.R., Jou, W.H., Strelets, M. & Allmaras, S.R. (1997).
//!   "Comments on the feasibility of LES for wings, and on a hybrid RANS/LES
//!   approach." *Advances in DNS/LES* 1:4–8.
//! - Shur, M., Spalart, P.R., Strelets, M. & Travin, A. (1999). "Detached-
//!   eddy simulation of an airfoil at high angle of attack." *IUTAM Symp.*

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::{DES_C_DES, SMAGORINSKY_CS_DEFAULT};

/// DES hybrid RANS-LES model (Spalart et al. 1997).
///
/// Blends RANS behaviour near walls with LES behaviour away from walls
/// by selecting the minimum of the RANS wall-distance and LES filter scale.
#[derive(Debug, Clone)]
pub struct DESModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> {
    /// DES constant C_DES = 0.65.
    pub c_des: T,
    /// Wall distance for each grid point (order: x-innermost, then y, then z).
    ///
    /// Must have length == nx·ny·nz.  Set via [`DESModel::with_wall_distances`].
    /// A uniform large value (e.g. 1.0 m) produces pure LES everywhere.
    pub wall_distances: Vec<T>,
    /// Maximum cell dimension Δ_max = max(dx, dy, dz) for each grid point.
    ///
    /// For uniform Cartesian grids, this is a constant.  For non-uniform grids
    /// supply the per-cell value.
    pub delta_max: Vec<T>,
    /// Smagorinsky constant C_s used for the LES sub-model away from walls.
    pub cs: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> DESModel<T> {
    /// Create a DES model with uniform wall distance and filter width.
    ///
    /// # Arguments
    /// * `n_points` — number of grid points (nx·ny·nz)
    /// * `uniform_wall_distance` — wall distance [m] applied uniformly (use a
    ///   large value for pure LES, 0 for pure RANS)
    /// * `dx`, `dy`, `dz` — uniform cell dimensions [m]
    pub fn new(n_points: usize, uniform_wall_distance: T, dx: T, dy: T, dz: T) -> Self {
        let c_des = <T as FromPrimitive>::from_f64(DES_C_DES).unwrap_or_else(T::one);
        let delta_max = num_traits::Float::max(dx, num_traits::Float::max(dy, dz));
        Self {
            c_des,
            wall_distances: vec![uniform_wall_distance; n_points],
            delta_max: vec![delta_max; n_points],
            cs: <T as FromPrimitive>::from_f64(SMAGORINSKY_CS_DEFAULT).unwrap_or_else(T::one),
        }
    }

    /// Create a DES model with per-point wall distances and delta_max values.
    pub fn with_wall_distances(wall_distances: Vec<T>, delta_max: Vec<T>, cs: T) -> Self {
        let c_des = <T as FromPrimitive>::from_f64(DES_C_DES).unwrap_or_else(T::one);
        Self {
            c_des,
            wall_distances,
            delta_max,
            cs,
        }
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for DESModel<T>
{
    /// Compute DES eddy viscosity.
    ///
    /// At each grid point, uses l_DES = min(d, C_DES·Δ_max) as the effective
    /// length scale for the Smagorinsky LES sub-model.
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let n = nx * ny * nz;
        let mut viscosity = Vec::with_capacity(n);
        let two = <T as FromPrimitive>::from_f64(2.0).unwrap_or_else(T::one);

        for idx in 0..n {
            let d = if idx < self.wall_distances.len() {
                self.wall_distances[idx]
            } else {
                T::one()
            };
            let delta = if idx < self.delta_max.len() {
                self.delta_max[idx]
            } else {
                T::one()
            };

            // l_DES = min(d, C_DES · Δ_max)
            let l_des = num_traits::Float::min(d, self.c_des * delta);

            let i = idx % nx;
            let j = (idx / nx) % ny;
            let k = idx / (nx * ny);

            // Compute strain rate magnitude at (i,j,k) using central differences with spacing l_des.
            let mut s_sq = T::zero();
            if i > 0 && i < nx - 1 {
                if let (Some(vp), Some(vm)) = (
                    flow_field.velocity.get(i + 1, j, k),
                    flow_field.velocity.get(i - 1, j, k),
                ) {
                    let s11 = (vp.x - vm.x) / (two * l_des);
                    s_sq = s_sq + s11 * s11;
                }
            }
            if j > 0 && j < ny - 1 {
                if let (Some(vp), Some(vm)) = (
                    flow_field.velocity.get(i, j + 1, k),
                    flow_field.velocity.get(i, j - 1, k),
                ) {
                    let s22 = (vp.y - vm.y) / (two * l_des);
                    s_sq = s_sq + s22 * s22;
                }
            }
            if k > 0 && k < nz - 1 {
                if let (Some(vp), Some(vm)) = (
                    flow_field.velocity.get(i, j, k + 1),
                    flow_field.velocity.get(i, j, k - 1),
                ) {
                    let s33 = (vp.z - vm.z) / (two * l_des);
                    s_sq = s_sq + s33 * s33;
                }
            }
            let s_mag = num_traits::Float::sqrt(two * s_sq);
            viscosity.push(self.cs * self.cs * l_des * l_des * s_mag);
        }
        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        self.turbulent_viscosity(flow_field)
    }

    fn name(&self) -> &'static str {
        "DES"
    }
}
