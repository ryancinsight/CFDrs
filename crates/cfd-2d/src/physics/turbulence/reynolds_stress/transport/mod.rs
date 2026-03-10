//! Full RSTM time-advancement (transport equations).
//!
//! This file owns:
//! - `update_reynolds_stresses` — primary public entry point (optimised)
//! - `update_reynolds_stresses_optimized` — block-cached, SIMD-friendly loop
//! - Wall boundary condition application
//! - Velocity/stress gradient helpers
//! - Scalar Laplacian helper
//! - [`TurbulenceModel`] trait implementation for framework compatibility
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \overline{u_i^\prime u_j^\prime}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\overline{u_i^\prime u_i^\prime} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

mod solver;

use super::diffusion::{dissipation_tensor, turbulent_transport};
use super::model::ReynoldsStressModel;
use super::production::production_term as prod_term;
use super::tensor::ReynoldsStressTensor;
use cfd_core::error::Result;
use nalgebra::{DMatrix, RealField, Vector2};
use num_traits::{FromPrimitive, ToPrimitive};

// Numerical stability floor for ε
const EPSILON_MIN: f64 = 1e-12;

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> ReynoldsStressModel<T> {
    fn c(v: f64) -> T {
        T::from_f64(v).expect("transport constant must be representable")
    }

    /// Dissolve Reynolds stresses to zero at wall cells.
    pub fn apply_wall_boundary_conditions(
        &self,
        xx: &mut DMatrix<T>,
        xy: &mut DMatrix<T>,
        yy: &mut DMatrix<T>,
        k: &mut DMatrix<T>,
        epsilon: &mut DMatrix<T>,
    ) {
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            xx[(i, 0)] = T::zero();
            xy[(i, 0)] = T::zero();
            yy[(i, 0)] = T::zero();
            k[(i, 0)] = T::zero();
            epsilon[(i, 0)] = T::zero();
            xx[(i, ny - 1)] = T::zero();
            xy[(i, ny - 1)] = T::zero();
            yy[(i, ny - 1)] = T::zero();
            k[(i, ny - 1)] = T::zero();
            epsilon[(i, ny - 1)] = T::zero();
        }
    }

    // --- grid helpers ---

    /// Compute the Reynolds stress production tensor component $P_{xy}$.
    pub fn production_term(
        &self,
        rs: &ReynoldsStressTensor<T>,
        velocity_gradient: &[[T; 2]; 2],
        i: usize,
        j: usize,
        x: usize,
        y: usize,
    ) -> T {
        prod_term(rs, velocity_gradient, i, j, x, y)
    }

    /// Compute the isotropic dissipation tensor component $\varepsilon_{xy}$.
    pub fn dissipation_tensor(
        &self,
        rs: &ReynoldsStressTensor<T>,
        i: usize,
        j: usize,
        x: usize,
        y: usize,
    ) -> T {
        dissipation_tensor(rs, i, j, x, y)
    }

    /// Compute the turbulent transport (Daly-Harlow) for a stress component.
    pub fn turbulent_transport(
        &self,
        _rs: &ReynoldsStressTensor<T>,
        k: T,
        epsilon: T,
        stress_gradient: &[[T; 2]; 2],
        i: usize,
        j: usize,
    ) -> T {
        turbulent_transport(k, epsilon, stress_gradient, i, j)
    }

    /// Compute the 2D velocity gradient tensor via central differences.
    pub fn calculate_velocity_gradients(
        &self,
        velocity: &[DMatrix<T>; 2],
        i: usize,
        j: usize,
        dx: T,
        dy: T,
    ) -> [[T; 2]; 2] {
        let half = Self::c(0.5);
        let dx_inv = T::one() / dx;
        let dy_inv = T::one() / dy;
        [
            [
                dx_inv * (velocity[0][(i + 1, j)] - velocity[0][(i - 1, j)]) * half,
                dy_inv * (velocity[0][(i, j + 1)] - velocity[0][(i, j - 1)]) * half,
            ],
            [
                dx_inv * (velocity[1][(i + 1, j)] - velocity[1][(i - 1, j)]) * half,
                dy_inv * (velocity[1][(i, j + 1)] - velocity[1][(i, j - 1)]) * half,
            ],
        ]
    }

    /// Decompose the velocity gradient into strain-rate and rotation-rate tensors.
    pub fn calculate_strain_rotation_rates(
        &self,
        velocity_gradient: &[[T; 2]; 2],
    ) -> ([[T; 2]; 2], [[T; 2]; 2]) {
        let half = Self::c(0.5);
        let du_dx = velocity_gradient[0][0];
        let du_dy = velocity_gradient[0][1];
        let dv_dx = velocity_gradient[1][0];
        let dv_dy = velocity_gradient[1][1];
        let s12 = half * (du_dy + dv_dx);
        let w12 = half * (du_dy - dv_dx);
        (
            [[du_dx, s12], [s12, dv_dy]],
            [[T::zero(), w12], [-w12, T::zero()]],
        )
    }

    /// Compute the spatial gradient of a Reynolds stress component.
    pub fn calculate_stress_gradients(
        &self,
        rs: &ReynoldsStressTensor<T>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
    ) -> [[T; 2]; 2] {
        let half = Self::c(0.5);
        let dx_inv = T::one() / dx;
        let dy_inv = T::one() / dy;
        [
            [
                dx_inv * (rs.xx[(i + 1, j)] - rs.xx[(i - 1, j)]) * half,
                dy_inv * (rs.xx[(i, j + 1)] - rs.xx[(i, j - 1)]) * half,
            ],
            [
                dx_inv * (rs.xy[(i + 1, j)] - rs.xy[(i - 1, j)]) * half,
                dy_inv * (rs.xy[(i, j + 1)] - rs.xy[(i, j - 1)]) * half,
            ],
        ]
    }

    /// Compute the 5-point Laplacian $\nabla^2 \phi$ for a scalar field.
    pub fn calculate_scalar_laplacian(
        &self,
        scalar: &DMatrix<T>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
    ) -> T {
        let two = Self::c(2.0);
        (scalar[(i + 1, j)] - two * scalar[(i, j)] + scalar[(i - 1, j)]) / (dx * dx)
            + (scalar[(i, j + 1)] - two * scalar[(i, j)] + scalar[(i, j - 1)]) / (dy * dy)
    }
}

// TurbulenceModel trait compatibility implementation
use super::super::traits::TurbulenceModel;

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> TurbulenceModel<T>
    for ReynoldsStressModel<T>
{
    fn turbulent_viscosity(&self, k: T, epsilon: T, _density: T) -> T {
        self.c_mu * k * k / epsilon
    }

    fn production_term(
        &self,
        velocity_gradient: &[[T; 2]; 2],
        turbulent_viscosity: T,
        _turbulence_variable: T,
        _wall_distance: T,
        _molecular_viscosity: T,
    ) -> T {
        let s11 = velocity_gradient[0][0];
        let s12 = velocity_gradient[0][1];
        let s22 = velocity_gradient[1][1];
        let strain_sq = Self::c(2.0) * (s11 * s11 + Self::c(2.0) * s12 * s12 + s22 * s22);
        turbulent_viscosity * strain_sq
    }

    fn dissipation_term(&self, _k: T, epsilon: T) -> T {
        epsilon
    }

    fn update(
        &mut self,
        _k: &mut [T],
        _epsilon_or_omega: &mut [T],
        _velocity: &[Vector2<T>],
        _density: T,
        _molecular_viscosity: T,
        _dt: T,
        _dx: T,
        _dy: T,
    ) -> Result<()> {
        Err(cfd_core::error::Error::InvalidConfiguration(
            "Reynolds stress model requires full tensor update via update_reynolds_stresses()"
                .to_string(),
        ))
    }

    fn name(&self) -> &'static str {
        "Reynolds Stress Transport Model (RSTM)"
    }

    fn is_valid_for_reynolds(&self, reynolds: T) -> bool {
        reynolds >= Self::c(1000.0)
    }
}
