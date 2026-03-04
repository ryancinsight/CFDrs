//! Vreman subgrid-scale model for LES.
//!
//! # Theorem — Vreman SGS Model (Vreman 2004)
//!
//! The eddy viscosity is given by:
//!
//! ```text
//! νₜ = C_V √(B_β / α_ij α_ij)
//! ```
//!
//! where:
//! - α_ij = ∂u_i / ∂x_j  (velocity gradient tensor)
//! - B_β = β₁₁β₂₂ − β₁₂² + β₁₁β₃₃ − β₁₃² + β₂₂β₃₃ − β₂₃²  (second invariant of β)
//! - β_mn = Δ_m² α_im α_in  (filtered gradient metric)
//! - C_V ≈ 2.5 C_s² for C_s ≈ 0.1 → C_V ≈ 0.025  (Vreman 2004 Table 1)
//!
//! **Properties** (Vreman 2004, §2):
//! - Exactly zero for pure shear flow (no SGS dissipation in laminar-like shear)
//! - Exactly zero for solid-body rotation
//! - Avoids the need for wall damping functions
//! - Computationally cheaper than the Dynamic Smagorinsky (no test filtering)
//!
//! ## Algorithm
//!
//! ```text
//! For each grid point:
//!   1. Compute α_ij = ∂u_i/∂x_j via central differences (spacing Δ)
//!   2. Compute α_ij α_ij = sum over all i,j of α_ij²
//!   3. Compute β_mn = Δ² · sum_i(α_im · α_in)  (isotropic Δ simplification)
//!   4. Compute B_β from upper-triangular 3×3 metric B_β = I₂(β)
//!   5. νₜ = C_V · sqrt(max(B_β, 0) / max(α_ij α_ij, ε))
//! ```
//!
//! ## References
//!
//! - Vreman, A.W. (2004). "An eddy-viscosity subgrid-scale model for
//!   turbulent shear flow." *Phys. Fluids* 16(10):3670–3681.

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::{DEARDORFF_ONE_THIRD, VREMAN_CV};

/// Vreman SGS model for LES (Vreman 2004).
///
/// Automatically vanishes for pure shear and solid-body rotation without
/// ad hoc damping functions.
#[derive(Debug, Clone)]
pub struct VremanModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Vreman constant C_V ≈ 0.025 (Vreman 2004 Table 1).
    pub c_v: T,
    /// Physical LES filter width Δ = (dx·dy·dz)^(1/3) [m].
    pub filter_width: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> VremanModel<T> {
    /// Create a Vreman model with default constant and unit filter width.
    pub fn new() -> Self {
        Self {
            c_v: <T as FromPrimitive>::from_f64(VREMAN_CV).unwrap_or_else(T::one),
            filter_width: T::one(),
        }
    }

    /// Create a Vreman model with physically correct filter width Δ = (dx·dy·dz)^(1/3).
    pub fn with_filter_width(dx: T, dy: T, dz: T) -> Self {
        let one_third = <T as FromPrimitive>::from_f64(DEARDORFF_ONE_THIRD).unwrap_or_else(T::one);
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self {
            c_v: <T as FromPrimitive>::from_f64(VREMAN_CV).unwrap_or_else(T::one),
            filter_width,
        }
    }

    /// Compute Vreman eddy viscosity at grid point (i,j,k).
    fn vreman_viscosity_at(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize) -> T {
        let (nx, ny, nz) = flow.velocity.dimensions;
        let delta = self.filter_width;
        let two = <T as FromPrimitive>::from_f64(2.0).unwrap_or_else(T::one);
        let eps = <T as FromPrimitive>::from_f64(1e-30).unwrap_or_else(T::zero);

        // α_ij = ∂u_i / ∂x_j  (central differences where interior, one-sided at boundaries)
        let mut alpha = [[T::zero(); 3]; 3]; // alpha[i][j] = ∂u_i/∂x_j

        // x-derivatives (j=0 column): ∂u/∂x, ∂v/∂x, ∂w/∂x
        if i > 0 && i < nx - 1 {
            if let (Some(vp), Some(vm)) = (
                flow.velocity.get(i + 1, j, k),
                flow.velocity.get(i - 1, j, k),
            ) {
                alpha[0][0] = (vp.x - vm.x) / (two * delta);
                alpha[1][0] = (vp.y - vm.y) / (two * delta);
                alpha[2][0] = (vp.z - vm.z) / (two * delta);
            }
        }
        // y-derivatives (j=1 column)
        if j > 0 && j < ny - 1 {
            if let (Some(vp), Some(vm)) = (
                flow.velocity.get(i, j + 1, k),
                flow.velocity.get(i, j - 1, k),
            ) {
                alpha[0][1] = (vp.x - vm.x) / (two * delta);
                alpha[1][1] = (vp.y - vm.y) / (two * delta);
                alpha[2][1] = (vp.z - vm.z) / (two * delta);
            }
        }
        // z-derivatives (j=2 column)
        if k > 0 && k < nz - 1 {
            if let (Some(vp), Some(vm)) = (
                flow.velocity.get(i, j, k + 1),
                flow.velocity.get(i, j, k - 1),
            ) {
                alpha[0][2] = (vp.x - vm.x) / (two * delta);
                alpha[1][2] = (vp.y - vm.y) / (two * delta);
                alpha[2][2] = (vp.z - vm.z) / (two * delta);
            }
        }

        // Compute |α|² = Σ α_ij²
        let mut alpha_sq = T::zero();
        for ai in &alpha {
            for &a in ai {
                alpha_sq = alpha_sq + a * a;
            }
        }

        if alpha_sq <= eps {
            return T::zero();
        }

        // β_mn = Δ² Σ_i α_im α_in  (isotropic Δ, i.e. Δ_m = Δ_n = Δ)
        let delta_sq = delta * delta;
        let mut beta = [[T::zero(); 3]; 3];
        for m in 0..3 {
            for n in 0..3 {
                let mut s = T::zero();
                for ii in 0..3 {
                    s = s + alpha[ii][m] * alpha[ii][n];
                }
                beta[m][n] = delta_sq * s;
            }
        }

        // B_β = β₁₁β₂₂ − β₁₂² + β₁₁β₃₃ − β₁₃² + β₂₂β₃₃ − β₂₃²
        let b_beta = beta[0][0] * beta[1][1] - beta[0][1] * beta[1][0] + beta[0][0] * beta[2][2]
            - beta[0][2] * beta[2][0]
            + beta[1][1] * beta[2][2]
            - beta[1][2] * beta[2][1];

        if b_beta <= T::zero() {
            return T::zero();
        }

        self.c_v * num_traits::Float::sqrt(b_beta / alpha_sq)
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> Default
    for VremanModel<T>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for VremanModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let mut viscosity = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    viscosity.push(self.vreman_viscosity_at(flow_field, i, j, k));
                }
            }
        }
        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        self.turbulent_viscosity(flow_field)
    }

    fn name(&self) -> &'static str {
        "Vreman"
    }
}
