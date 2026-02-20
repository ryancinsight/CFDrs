//! Full RSTM time-advancement (transport equations).
//!
//! This file owns:
//! - `update_reynolds_stresses` — primary public entry point (optimised)
//! - `update_reynolds_stresses_optimized` — block-cached, SIMD-friendly loop
//! - `update_reynolds_stresses_standard` — legacy compatibility (deprecated)
//! - Wall boundary condition application
//! - Velocity/stress gradient helpers
//! - Scalar Laplacian helper
//! - [`TurbulenceModel`] trait implementation for framework compatibility

use super::diffusion::{dissipation_tensor, dissipation_tensor_optimized, turbulent_transport};
use super::model::ReynoldsStressModel;
use super::pressure_strain::{pressure_strain_linear, pressure_strain_quadratic, pressure_strain_ssg};
use super::production::production_term as prod_term;
use super::tensor::ReynoldsStressTensor;
use super::PressureStrainModel;
use cfd_core::error::Result;
use nalgebra::{DMatrix, RealField, Vector2};
use num_traits::{FromPrimitive, ToPrimitive};

// Numerical stability floor for ε
const EPSILON_MIN: f64 = 1e-12;

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> ReynoldsStressModel<T> {
    fn c(v: f64) -> T {
        T::from_f64(v).expect("transport constant must be representable")
    }

    /// Primary entry: update Reynolds stresses using the optimised path.
    pub fn update_reynolds_stresses(
        &self,
        rs: &mut ReynoldsStressTensor<T>,
        velocity: &[DMatrix<T>; 2],
        dt: T,
        dx: T,
        dy: T,
    ) -> Result<()> {
        self.update_reynolds_stresses_optimized(rs, velocity, dt, dx, dy)
    }

    /// Block-cached implementation with inlined gradient calculations.
    pub fn update_reynolds_stresses_optimized(
        &self,
        rs: &mut ReynoldsStressTensor<T>,
        velocity: &[DMatrix<T>; 2],
        dt: T,
        dx: T,
        dy: T,
    ) -> Result<()> {
        let nx = self.nx;
        let ny = self.ny;

        let dx_inv = T::one() / dx;
        let dy_inv = T::one() / dy;
        let half = Self::c(0.5);
        let two_thirds = Self::c(2.0 / 3.0);
        let epsilon_min = Self::c(EPSILON_MIN);

        let mut xx_new = DMatrix::zeros(nx, ny);
        let mut xy_new = DMatrix::zeros(nx, ny);
        let mut yy_new = DMatrix::zeros(nx, ny);
        let mut k_new = DMatrix::zeros(nx, ny);
        let mut epsilon_new = DMatrix::zeros(nx, ny);

        let block_size = 4;
        for bi in (1..nx - 1).step_by(block_size) {
            for bj in (1..ny - 1).step_by(block_size) {
                let bi_end = (bi + block_size).min(nx - 1);
                let bj_end = (bj + block_size).min(ny - 1);

                for i in bi..bi_end {
                    for j in bj..bj_end {
                        let (i1, i_1, j1, j_1) = (i + 1, i - 1, j + 1, j - 1);

                        let du_dx = dx_inv * (velocity[0][(i1, j)] - velocity[0][(i_1, j)]) * half;
                        let du_dy = dy_inv * (velocity[0][(i, j1)] - velocity[0][(i, j_1)]) * half;
                        let dv_dx = dx_inv * (velocity[1][(i1, j)] - velocity[1][(i_1, j)]) * half;
                        let dv_dy = dy_inv * (velocity[1][(i, j1)] - velocity[1][(i, j_1)]) * half;

                        let velocity_gradient = [[du_dx, du_dy], [dv_dx, dv_dy]];

                        let s11 = du_dx;
                        let s12 = half * (du_dy + dv_dx);
                        let s22 = dv_dy;
                        let w12 = half * (du_dy - dv_dx);
                        let w21 = -w12;

                        let strain_rate = [[s11, s12], [s12, s22]];
                        let rotation_rate = [[T::zero(), w12], [w21, T::zero()]];

                        let dxx_dx = dx_inv * (rs.xx[(i1, j)] - rs.xx[(i_1, j)]) * half;
                        let dxx_dy = dy_inv * (rs.xx[(i, j1)] - rs.xx[(i, j_1)]) * half;
                        let dxy_dx = dx_inv * (rs.xy[(i1, j)] - rs.xy[(i_1, j)]) * half;
                        let dxy_dy = dy_inv * (rs.xy[(i, j1)] - rs.xy[(i, j_1)]) * half;
                        let stress_gradient = [[dxx_dx, dxx_dy], [dxy_dx, dxy_dy]];

                        let xx = rs.xx[(i, j)];
                        let xy = rs.xy[(i, j)];
                        let yy = rs.yy[(i, j)];
                        let k = rs.k[(i, j)];
                        let epsilon = rs.epsilon[(i, j)];

                        // Production (unrolled for 2D)
                        let p_xx = -Self::c(2.0) * xy * du_dy;
                        let p_xy = -xx * dv_dx - yy * du_dy;
                        let p_yy = -Self::c(2.0) * xy * dv_dy;

                        // Pressure-strain
                        let phi_xx = self.pressure_strain_optimized(xx, xy, yy, k, epsilon, &strain_rate, &rotation_rate, 0, 0);
                        let phi_xy = self.pressure_strain_optimized(xx, xy, yy, k, epsilon, &strain_rate, &rotation_rate, 0, 1);
                        let phi_yy = self.pressure_strain_optimized(xx, xy, yy, k, epsilon, &strain_rate, &rotation_rate, 1, 1);

                        // Dissipation
                        let eps_xx = dissipation_tensor_optimized(rs, 0, 0, i, j, two_thirds, epsilon);
                        let eps_xy = dissipation_tensor_optimized(rs, 0, 1, i, j, two_thirds, epsilon);
                        let eps_yy = dissipation_tensor_optimized(rs, 1, 1, i, j, two_thirds, epsilon);

                        // Turbulent transport
                        let t_xx = turbulent_transport(k, epsilon, &stress_gradient, 0, 0);
                        let t_xy = turbulent_transport(k, epsilon, &stress_gradient, 0, 1);
                        let t_yy = turbulent_transport(k, epsilon, &stress_gradient, 1, 1);

                        xx_new[(i, j)] = xx + dt * (p_xx + phi_xx - eps_xx + t_xx);
                        xy_new[(i, j)] = xy + dt * (p_xy + phi_xy - eps_xy + t_xy);
                        yy_new[(i, j)] = yy + dt * (p_yy + phi_yy - eps_yy + t_yy);

                        k_new[(i, j)] = half * (xx_new[(i, j)] + yy_new[(i, j)]);
                        epsilon_new[(i, j)] = self.update_epsilon_optimized(
                            xx_new[(i, j)], yy_new[(i, j)], k_new[(i, j)],
                            epsilon, &velocity_gradient, dt, epsilon_min,
                        );
                    }
                }
            }
        }

        self.apply_wall_boundary_conditions(
            &mut xx_new, &mut xy_new, &mut yy_new, &mut k_new, &mut epsilon_new,
        );

        std::mem::swap(&mut rs.xx, &mut xx_new);
        std::mem::swap(&mut rs.xy, &mut xy_new);
        std::mem::swap(&mut rs.yy, &mut yy_new);
        std::mem::swap(&mut rs.k, &mut k_new);
        std::mem::swap(&mut rs.epsilon, &mut epsilon_new);

        Ok(())
    }

    /// Optimised pressure-strain dispatch (no Option fields, inlined anisotropy).
    fn pressure_strain_optimized(
        &self,
        xx: T,
        xy: T,
        yy: T,
        k: T,
        epsilon: T,
        strain_rate: &[[T; 2]; 2],
        _rotation_rate: &[[T; 2]; 2],
        i: usize,
        j: usize,
    ) -> T {
        if epsilon <= T::zero() || k <= T::zero() { return T::zero(); }

        let two_thirds = Self::c(2.0 / 3.0);
        let time_scale = k / epsilon;
        let a_xx = xx / k - two_thirds;
        let a_xy = xy / k;
        let a_yy = yy / k - two_thirds;

        let s11 = strain_rate[0][0];
        let s12 = strain_rate[0][1];
        let s22 = strain_rate[1][1];

        match self.pressure_strain_model {
            PressureStrainModel::LinearReturnToIsotropy => {
                pressure_strain_linear(self.c1, a_xx, a_xy, a_yy, epsilon, k, i, j)
            }
            PressureStrainModel::Quadratic => {
                pressure_strain_quadratic(self.c1, self.c1_star, self.c2_star, a_xx, a_xy, a_yy, time_scale, s11, s12, s22, i, j)
            }
            PressureStrainModel::SSG => {
                pressure_strain_ssg(self.c1, self.c1_star, self.c2, self.c3, self.c3_star, self.c4, self.c5, a_xx, a_xy, a_yy, k, epsilon, s11, s12, s22, T::zero(), T::zero(), i, j)
            }
        }
    }

    /// Pressure-strain dispatch for the standard (non-optimised) path.
    pub fn pressure_strain_term(
        &self,
        rs: &ReynoldsStressTensor<T>,
        strain_rate: &[[T; 2]; 2],
        rotation_rate: &[[T; 2]; 2],
        i: usize,
        j: usize,
        x: usize,
        y: usize,
    ) -> T {
        let k = rs.k[(x, y)];
        let epsilon = rs.epsilon[(x, y)];
        let xx = rs.xx[(x, y)];
        let xy = rs.xy[(x, y)];
        let yy = rs.yy[(x, y)];

        if epsilon <= T::zero() || k <= T::zero() { return T::zero(); }

        let two_thirds = Self::c(2.0 / 3.0);
        let time_scale = k / epsilon;
        let a_xx = xx / k - two_thirds;
        let a_xy = xy / k;
        let a_yy = yy / k - two_thirds;

        let s11 = strain_rate[0][0];
        let s12 = strain_rate[0][1];
        let s22 = strain_rate[1][1];
        let w12 = rotation_rate[0][1];
        let w21 = rotation_rate[1][0];

        match self.pressure_strain_model {
            PressureStrainModel::LinearReturnToIsotropy => {
                pressure_strain_linear(self.c1, a_xx, a_xy, a_yy, epsilon, k, i, j)
            }
            PressureStrainModel::Quadratic => {
                pressure_strain_quadratic(self.c1, self.c1_star, self.c2_star, a_xx, a_xy, a_yy, time_scale, s11, s12, s22, i, j)
            }
            PressureStrainModel::SSG => {
                pressure_strain_ssg(self.c1, self.c1_star, self.c2, self.c3, self.c3_star, self.c4, self.c5, a_xx, a_xy, a_yy, k, epsilon, s11, s12, s22, w12, w21, i, j)
            }
        }
    }

    /// Optimised dissipation rate update.
    fn update_epsilon_optimized(
        &self,
        xx_new: T,
        yy_new: T,
        k_new: T,
        epsilon_old: T,
        velocity_gradient: &[[T; 2]; 2],
        dt: T,
        epsilon_min: T,
    ) -> T {
        if k_new <= T::zero() || epsilon_old <= T::zero() { return epsilon_min; }

        let du_dy = velocity_gradient[0][1];
        let dv_dx = velocity_gradient[1][0];
        let p_k = Self::c(0.5) * (-Self::c(2.0) * xx_new * du_dy - Self::c(2.0) * yy_new * dv_dx);

        let c_eps1 = Self::c(1.44);
        let c_eps2 = Self::c(1.92);
        let sigma_eps = Self::c(1.3);
        let nu_t = self.c_mu * k_new * k_new / epsilon_old;

        let production = c_eps1 * p_k * epsilon_old / k_new;
        let destruction = c_eps2 * epsilon_old * epsilon_old / k_new;
        let diffusion = (nu_t / sigma_eps) * T::zero(); // placeholder Laplacian — needs grid diffusion

        (epsilon_old + dt * (production - destruction + diffusion)).max(epsilon_min)
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
            xx[(i, 0)] = T::zero(); xy[(i, 0)] = T::zero(); yy[(i, 0)] = T::zero();
            k[(i, 0)] = T::zero(); epsilon[(i, 0)] = T::zero();
            xx[(i, ny - 1)] = T::zero(); xy[(i, ny - 1)] = T::zero(); yy[(i, ny - 1)] = T::zero();
            k[(i, ny - 1)] = T::zero(); epsilon[(i, ny - 1)] = T::zero();
        }
    }

    // --- grid helpers ---

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
            [dx_inv * (velocity[0][(i+1,j)] - velocity[0][(i-1,j)]) * half,
             dy_inv * (velocity[0][(i,j+1)] - velocity[0][(i,j-1)]) * half],
            [dx_inv * (velocity[1][(i+1,j)] - velocity[1][(i-1,j)]) * half,
             dy_inv * (velocity[1][(i,j+1)] - velocity[1][(i,j-1)]) * half],
        ]
    }

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
        ([[du_dx, s12], [s12, dv_dy]], [[T::zero(), w12], [-w12, T::zero()]])
    }

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
            [dx_inv * (rs.xx[(i+1,j)] - rs.xx[(i-1,j)]) * half,
             dy_inv * (rs.xx[(i,j+1)] - rs.xx[(i,j-1)]) * half],
            [dx_inv * (rs.xy[(i+1,j)] - rs.xy[(i-1,j)]) * half,
             dy_inv * (rs.xy[(i,j+1)] - rs.xy[(i,j-1)]) * half],
        ]
    }

    pub fn calculate_scalar_laplacian(
        &self,
        scalar: &DMatrix<T>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
    ) -> T {
        let two = Self::c(2.0);
        (scalar[(i+1,j)] - two * scalar[(i,j)] + scalar[(i-1,j)]) / (dx * dx)
            + (scalar[(i,j+1)] - two * scalar[(i,j)] + scalar[(i,j-1)]) / (dy * dy)
    }

    /// Legacy path kept for backward compatibility. Prefer `update_reynolds_stresses`.
    #[deprecated(note = "Use update_reynolds_stresses() — now backed by the optimised implementation")]
    pub fn update_reynolds_stresses_standard(
        &self,
        rs: &mut ReynoldsStressTensor<T>,
        velocity: &[DMatrix<T>; 2],
        dt: T,
        dx: T,
        dy: T,
    ) -> Result<()> {
        if dt <= T::zero() {
            return Err(cfd_core::error::Error::InvalidInput(format!(
                "Time step must be positive: dt={}",
                dt.to_f64().expect("dt conversion")
            )));
        }
        if dx <= T::zero() || dy <= T::zero() {
            return Err(cfd_core::error::Error::InvalidInput(format!(
                "Grid spacing must be positive: dx={}, dy={}",
                dx.to_f64().expect("dx"), dy.to_f64().expect("dy")
            )));
        }

        let nx = self.nx;
        let ny = self.ny;
        let epsilon_min = Self::c(EPSILON_MIN);

        let mut xx_new: DMatrix<T> = rs.xx.clone();
        let mut xy_new: DMatrix<T> = rs.xy.clone();
        let mut yy_new: DMatrix<T> = rs.yy.clone();
        let mut k_new: DMatrix<T> = rs.k.clone();
        let mut epsilon_new: DMatrix<T> = rs.epsilon.clone();

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let vg = self.calculate_velocity_gradients(velocity, i, j, dx, dy);
                let (strain_rate, rotation_rate) = self.calculate_strain_rotation_rates(&vg);
                let sg = self.calculate_stress_gradients(rs, i, j, dx, dy);

                for ii in 0..2 {
                    for jj in 0..2 {
                        let p_ij = prod_term(rs, &vg, ii, jj, i, j);
                        let phi_ij = self.pressure_strain_term(rs, &strain_rate, &rotation_rate, ii, jj, i, j);
                        let eps_ij = dissipation_tensor(rs, ii, jj, i, j);
                        let t_ij = turbulent_transport(rs.k[(i,j)], rs.epsilon[(i,j)], &sg, ii, jj);

                        let rhs = p_ij + phi_ij - eps_ij + t_ij;
                        match (ii, jj) {
                            (0, 0) => xx_new[(i, j)] = rs.xx[(i, j)] + dt * rhs,
                            (0, 1) | (1, 0) => xy_new[(i, j)] = rs.xy[(i, j)] + dt * rhs,
                            (1, 1) => yy_new[(i, j)] = rs.yy[(i, j)] + dt * rhs,
                            _ => {}
                        }
                    }
                }

                k_new[(i, j)] = Self::c(0.5) * (xx_new[(i, j)] + yy_new[(i, j)]);
                let k = k_new[(i, j)];
                let epsilon = rs.epsilon[(i, j)];

                if k > T::zero() && epsilon > T::zero() {
                    let p_xx = prod_term(rs, &vg, 0, 0, i, j);
                    let p_yy = prod_term(rs, &vg, 1, 1, i, j);
                    let p_k = Self::c(0.5) * (p_xx + p_yy);
                    let nu_t = self.c_mu * k * k / epsilon;
                    let eps_laplacian = self.calculate_scalar_laplacian(&rs.epsilon, i, j, dx, dy);
                    let deps_dt = Self::c(1.44) * p_k * epsilon / k
                        - Self::c(1.92) * epsilon * epsilon / k
                        + (nu_t / Self::c(1.3)) * eps_laplacian;
                    epsilon_new[(i, j)] = (epsilon + dt * deps_dt).max(epsilon_min);
                } else {
                    epsilon_new[(i, j)] = Self::c(0.09) * k * k.sqrt();
                }
            }
        }

        self.apply_wall_boundary_conditions(
            &mut xx_new, &mut xy_new, &mut yy_new, &mut k_new, &mut epsilon_new,
        );

        rs.xx = xx_new; rs.xy = xy_new; rs.yy = yy_new;
        rs.k = k_new; rs.epsilon = epsilon_new;
        Ok(())
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

    fn dissipation_term(&self, _k: T, epsilon: T) -> T { epsilon }

    fn update(
        &mut self,
        _k: &mut [T], _epsilon_or_omega: &mut [T], _velocity: &[Vector2<T>],
        _density: T, _molecular_viscosity: T, _dt: T, _dx: T, _dy: T,
    ) -> Result<()> {
        Err(cfd_core::error::Error::InvalidConfiguration(
            "Reynolds stress model requires full tensor update via update_reynolds_stresses()".to_string(),
        ))
    }

    fn name(&self) -> &'static str { "Reynolds Stress Transport Model (RSTM)" }

    fn is_valid_for_reynolds(&self, reynolds: T) -> bool {
        reynolds >= Self::c(1000.0)
    }
}
