//! Core Reynolds stress transport solver — optimised block-cached implementation.

use super::super::diffusion::{dissipation_tensor_optimized, turbulent_transport};
use super::super::model::ReynoldsStressModel;
use super::super::pressure_strain::{
    pressure_strain_linear, pressure_strain_quadratic, pressure_strain_ssg,
};
use super::super::tensor::ReynoldsStressTensor;
use super::super::PressureStrainModel;
use super::EPSILON_MIN;
use cfd_core::error::Result;
use nalgebra::{DMatrix, RealField};
use num_traits::{FromPrimitive, ToPrimitive};

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> ReynoldsStressModel<T> {
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
    ///
    /// # Theorem (Realizability of Reynolds Stress Tensor)
    ///
    /// The Reynolds stress tensor $R_{ij} = \langle u_i' u_j' \rangle$ must satisfy:
    ///
    /// 1. **Non-negativity of diagonal**: $R_{ii} \ge 0$ (variance is non-negative)
    /// 2. **Schwarz inequality**: $|R_{ij}|^2 \le R_{ii} R_{jj}$ (Cauchy–Schwarz)
    /// 3. **Non-negative TKE**: $k = \tfrac{1}{2}\mathrm{tr}(R) \ge 0$
    ///
    /// **Proof sketch**: Since $R_{ii} = \langle u_i'^2 \rangle \ge 0$ by definition
    /// of variance, and $|\langle u_i' u_j' \rangle|^2 \le \langle u_i'^2 \rangle
    /// \langle u_j'^2 \rangle$ by the Cauchy–Schwarz inequality for $L^2$ inner products,
    /// the three constraints follow from the positive semi-definiteness of the
    /// covariance matrix. The explicit Euler update can violate these for large $\Delta t$;
    /// we enforce them by clamping diagonal components to zero and projecting the
    /// off-diagonal onto the Schwarz-admissible interval $[-\sqrt{R_{ii}R_{jj}},\,
    /// \sqrt{R_{ii}R_{jj}}]$ after each time step.
    ///
    /// **Reference**: Schumann, U. (1977). "Realizability of Reynolds-stress
    /// turbulence models." *Physics of Fluids*, 20(5), 721–725.
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
                        let phi_xx = self.pressure_strain_optimized(
                            xx,
                            xy,
                            yy,
                            k,
                            epsilon,
                            &strain_rate,
                            &rotation_rate,
                            0,
                            0,
                        );
                        let phi_xy = self.pressure_strain_optimized(
                            xx,
                            xy,
                            yy,
                            k,
                            epsilon,
                            &strain_rate,
                            &rotation_rate,
                            0,
                            1,
                        );
                        let phi_yy = self.pressure_strain_optimized(
                            xx,
                            xy,
                            yy,
                            k,
                            epsilon,
                            &strain_rate,
                            &rotation_rate,
                            1,
                            1,
                        );

                        // Dissipation
                        let eps_xx =
                            dissipation_tensor_optimized(rs, 0, 0, i, j, two_thirds, epsilon);
                        let eps_xy =
                            dissipation_tensor_optimized(rs, 0, 1, i, j, two_thirds, epsilon);
                        let eps_yy =
                            dissipation_tensor_optimized(rs, 1, 1, i, j, two_thirds, epsilon);

                        // Turbulent transport
                        let t_xx = turbulent_transport(k, epsilon, &stress_gradient, 0, 0);
                        let t_xy = turbulent_transport(k, epsilon, &stress_gradient, 0, 1);
                        let t_yy = turbulent_transport(k, epsilon, &stress_gradient, 1, 1);

                        xx_new[(i, j)] = (xx + dt * (p_xx + phi_xx - eps_xx + t_xx))
                            .max(T::zero());
                        xy_new[(i, j)] = xy + dt * (p_xy + phi_xy - eps_xy + t_xy);
                        yy_new[(i, j)] = (yy + dt * (p_yy + phi_yy - eps_yy + t_yy))
                            .max(T::zero());

                        // Schwarz inequality: |⟨u'v'⟩|² ≤ ⟨u'u'⟩⟨v'v'⟩
                        let xy_max = (xx_new[(i, j)] * yy_new[(i, j)]).sqrt();
                        xy_new[(i, j)] = xy_new[(i, j)].clamp(-xy_max, xy_max);

                        k_new[(i, j)] = half * (xx_new[(i, j)] + yy_new[(i, j)]);
                        epsilon_new[(i, j)] = self.update_epsilon_optimized(
                            xx_new[(i, j)],
                            yy_new[(i, j)],
                            k_new[(i, j)],
                            epsilon,
                            &velocity_gradient,
                            &rs.epsilon,
                            i,
                            j,
                            dx,
                            dy,
                            dt,
                            epsilon_min,
                        );
                    }
                }
            }
        }

        self.apply_wall_boundary_conditions(
            &mut xx_new,
            &mut xy_new,
            &mut yy_new,
            &mut k_new,
            &mut epsilon_new,
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
        if epsilon <= T::zero() || k <= T::zero() {
            return T::zero();
        }

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
            PressureStrainModel::Quadratic => pressure_strain_quadratic(
                self.c1,
                self.c1_star,
                self.c2_star,
                a_xx,
                a_xy,
                a_yy,
                time_scale,
                s11,
                s12,
                s22,
                i,
                j,
            ),
            PressureStrainModel::SSG => pressure_strain_ssg(
                self.c1,
                self.c1_star,
                self.c2,
                self.c3,
                self.c3_star,
                self.c4,
                self.c5,
                a_xx,
                a_xy,
                a_yy,
                k,
                epsilon,
                s11,
                s12,
                s22,
                T::zero(),
                T::zero(),
                i,
                j,
            ),
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

        if epsilon <= T::zero() || k <= T::zero() {
            return T::zero();
        }

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
            PressureStrainModel::Quadratic => pressure_strain_quadratic(
                self.c1,
                self.c1_star,
                self.c2_star,
                a_xx,
                a_xy,
                a_yy,
                time_scale,
                s11,
                s12,
                s22,
                i,
                j,
            ),
            PressureStrainModel::SSG => pressure_strain_ssg(
                self.c1,
                self.c1_star,
                self.c2,
                self.c3,
                self.c3_star,
                self.c4,
                self.c5,
                a_xx,
                a_xy,
                a_yy,
                k,
                epsilon,
                s11,
                s12,
                s22,
                w12,
                w21,
                i,
                j,
            ),
        }
    }

    /// Optimised dissipation rate update.
    ///
    /// Computes the ε transport equation:
    ///
    /// $$\frac{\partial \varepsilon}{\partial t} = C_{\varepsilon 1} \frac{P_k \varepsilon}{k} - C_{\varepsilon 2} \frac{\varepsilon^2}{k} + \frac{\partial}{\partial x_j}\!\left(\frac{\nu_t}{\sigma_\varepsilon} \frac{\partial \varepsilon}{\partial x_j}\right)$$
    ///
    /// # Theorem
    ///
    /// For $k > 0$ and $\varepsilon > 0$, the Daly–Harlow gradient-diffusion model
    /// with the 5-point Laplacian stencil is second-order accurate and satisfies
    /// discrete maximum-principle consistency provided the CFL condition holds.
    ///
    /// **Proof sketch**: The central-difference Laplacian
    /// $\nabla^2 \varepsilon \approx (\varepsilon_{i+1,j} - 2\varepsilon_{i,j} +
    /// \varepsilon_{i-1,j})/\Delta x^2 + (\varepsilon_{i,j+1} - 2\varepsilon_{i,j} +
    /// \varepsilon_{i,j-1})/\Delta y^2$ has a non-negative stencil coefficient structure
    /// when $\Delta t \le \frac{1}{2(\nu_t/\sigma_\varepsilon)(1/\Delta x^2 + 1/\Delta y^2)}$,
    /// ensuring no spurious extrema are introduced.
    fn update_epsilon_optimized(
        &self,
        xx_new: T,
        yy_new: T,
        k_new: T,
        epsilon_old: T,
        velocity_gradient: &[[T; 2]; 2],
        epsilon_field: &DMatrix<T>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
        dt: T,
        epsilon_min: T,
    ) -> T {
        if k_new <= T::zero() || epsilon_old <= T::zero() {
            return epsilon_min;
        }

        let du_dy = velocity_gradient[0][1];
        let dv_dx = velocity_gradient[1][0];
        let p_k = Self::c(0.5) * (-Self::c(2.0) * xx_new * du_dy - Self::c(2.0) * yy_new * dv_dx);

        let c_eps1 = Self::c(1.44);
        let c_eps2 = Self::c(1.92);
        let sigma_eps = Self::c(1.3);
        let nu_t = self.c_mu * k_new * k_new / epsilon_old;

        let production = c_eps1 * p_k * epsilon_old / k_new;
        let destruction = c_eps2 * epsilon_old * epsilon_old / k_new;

        // 5-point Laplacian stencil: ∇²ε = (ε_{i+1,j} - 2ε_{i,j} + ε_{i-1,j})/dx²
        //                                 + (ε_{i,j+1} - 2ε_{i,j} + ε_{i,j-1})/dy²
        let two = Self::c(2.0);
        let eps_laplacian = (epsilon_field[(i + 1, j)] - two * epsilon_field[(i, j)]
            + epsilon_field[(i - 1, j)])
            / (dx * dx)
            + (epsilon_field[(i, j + 1)] - two * epsilon_field[(i, j)]
                + epsilon_field[(i, j - 1)])
                / (dy * dy);
        let diffusion = (nu_t / sigma_eps) * eps_laplacian;

        (epsilon_old + dt * (production - destruction + diffusion)).max(epsilon_min)
    }
}
