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
//!   3. Compute β = A diag(Δx², Δy², Δz²) Aᵀ  (anisotropic Cartesian grid)
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

use super::constants::VREMAN_CV;
use super::field_ops::velocity_gradient_tensor;

/// Vreman SGS model for LES (Vreman 2004).
///
/// Automatically vanishes for pure shear and solid-body rotation without
/// ad hoc damping functions.
#[derive(Debug, Clone)]
pub struct VremanModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Vreman constant C_V ≈ 0.025 (Vreman 2004 Table 1).
    pub c_v: T,
    /// Grid spacing in the x direction [m].
    pub dx: T,
    /// Grid spacing in the y direction [m].
    pub dy: T,
    /// Grid spacing in the z direction [m].
    pub dz: T,
    /// Physical LES filter width Δ = (dx·dy·dz)^(1/3) [m].
    pub filter_width: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> VremanModel<T> {
    /// Create a Vreman model with default constant and unit filter width.
    pub fn new() -> Self {
        Self {
            c_v: <T as FromPrimitive>::from_f64(VREMAN_CV)
                .expect("VREMAN_CV is an IEEE 754 representable f64 constant"),
            dx: T::one(),
            dy: T::one(),
            dz: T::one(),
            filter_width: T::one(),
        }
    }

    /// Create a Vreman model with physically correct filter width Δ = (dx·dy·dz)^(1/3).
    pub fn with_filter_width(dx: T, dy: T, dz: T) -> Self {
        let one_third = T::one() / (T::one() + T::one() + T::one());
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self {
            c_v: <T as FromPrimitive>::from_f64(VREMAN_CV)
                .expect("VREMAN_CV is an IEEE 754 representable f64 constant"),
            dx,
            dy,
            dz,
            filter_width,
        }
    }

    /// Compute Vreman eddy viscosity at grid point (i,j,k).
    fn vreman_viscosity_at(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize) -> T {
        let (nx, ny, nz) = flow.velocity.dimensions;
        let gradient = velocity_gradient_tensor(
            &flow.velocity.components,
            nx,
            ny,
            nz,
            i,
            j,
            k,
            self.dx,
            self.dy,
            self.dz,
        );
        let eps = <T as FromPrimitive>::from_f64(1e-30)
            .expect("1e-30 is an IEEE 754 representable f64 constant");

        let alpha = gradient;

        // Compute |α|² = Σ α_ij²
        let mut alpha_sq = T::zero();
        for ai in &alpha {
            for &a in ai {
                alpha_sq += a * a;
            }
        }

        if alpha_sq <= eps {
            return T::zero();
        }

        // β = A diag(Δx², Δy², Δz²) Aᵀ for the row-major velocity gradient A.
        let spacing_sq = [self.dx * self.dx, self.dy * self.dy, self.dz * self.dz];
        let mut beta = [[T::zero(); 3]; 3];
        for m in 0..3 {
            for n in 0..3 {
                let mut s = T::zero();
                for p in 0..3 {
                    s += spacing_sq[p] * alpha[m][p] * alpha[n][p];
                }
                beta[m][n] = s;
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    fn fill_velocity_field<F>(flow: &mut FlowField<f64>, mut generator: F)
    where
        F: FnMut(f64, f64, f64) -> Vector3<f64>,
    {
        let (nx, ny, nz) = flow.velocity.dimensions;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = k * nx * ny + j * nx + i;
                    flow.velocity.components[idx] = generator(i as f64, j as f64, k as f64);
                }
            }
        }
    }

    #[test]
    fn pure_shear_has_zero_viscosity() {
        let mut flow = FlowField::<f64>::new(4, 4, 4);
        fill_velocity_field(&mut flow, |_, y, _z| Vector3::new(y, 0.0, 0.0));

        let model = VremanModel::<f64>::with_filter_width(1.0, 1.0, 1.0);
        let viscosity = model.turbulent_viscosity(&flow);
        let center = 2 * 4 * 4 + 2 * 4 + 2;
        assert_relative_eq!(viscosity[center], 0.0, epsilon = 1e-12);
    }
}
