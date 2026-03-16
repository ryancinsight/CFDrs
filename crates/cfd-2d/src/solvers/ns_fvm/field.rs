//! `FlowField2D` — velocity, pressure, viscosity, and mask field storage.
//!
//! Stores fields on the staggered grid:
//! - `u[i][j]`  at vertical faces  (u ∈ [0..nx+1] × [0..ny])
//! - `v[i][j]`  at horizontal faces (v ∈ [0..nx] × [0..ny+1])
//! - `p[i][j]`  at cell centers     (p ∈ [0..nx] × [0..ny])
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

use crate::solvers::ns_fvm::array2d::{Array2D, Mask2D};
use super::boundary::BloodModel;
use super::grid::StaggeredGrid2D;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// 2D velocity and pressure fields on a staggered grid.
#[derive(Debug, Clone)]
pub struct FlowField2D<T: RealField + Copy> {
    /// u-velocity at vertical faces [nx+1][ny]
    pub u: Array2D<T>,
    /// v-velocity at horizontal faces [nx][ny+1]
    pub v: Array2D<T>,
    /// Pressure at cell centers [nx][ny]
    pub p: Array2D<T>,
    /// Viscosity at cell centers [nx][ny]
    pub mu: Array2D<T>,
    /// Shear rate at cell centers [nx][ny]
    pub gamma_dot: Array2D<T>,
    /// Fluid mask (true = fluid, false = solid)
    pub mask: Mask2D,
}

impl<T: RealField + Copy + FromPrimitive + Float> FlowField2D<T> {
    /// Create a new zero-initialised flow field for an `nx × ny` grid.
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            u: Array2D::new(nx + 1, ny, T::zero()),
            v: Array2D::new(nx, ny + 1, T::zero()),
            p: Array2D::new(nx, ny, T::zero()),
            mu: Array2D::new(nx, ny, T::zero()),
            gamma_dot: Array2D::new(nx, ny, T::zero()),
            mask: Array2D::new(nx, ny, true),
        }
    }

    /// Compute shear rate magnitude `γ̇ = √(2 D:D)` at cell `(i,j)`.
    ///
    /// For 2D incompressible flow:
    /// `γ̇ = √(2[D₁₁² + D₂₂² + 2D₁₂²])`, where `D = (∇u + ∇uᵀ)/2`.
    pub fn compute_shear_rate(&self, i: usize, j: usize, dx: T, dy: T) -> T {
        let two = T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);

        let dudx = if i < self.u.rows() - 1 {
            (self.u[(i + 1, j)] - self.u[(i, j)]) / dx
        } else {
            T::zero()
        };
        let dvdy = if j < self.v.cols() - 1 {
            (self.v[(i, j + 1)] - self.v[(i, j)]) / dy
        } else {
            T::zero()
        };
        let dudy = if j > 0 && j < self.u.cols() - 1 && i < self.u.rows() - 1 {
            let u_up = (self.u[(i, j + 1)] + self.u[(i + 1, j + 1)]) / two;
            let u_dn = (self.u[(i, j - 1)] + self.u[(i + 1, j - 1)]) / two;
            (u_up - u_dn) / (two * dy)
        } else {
            T::zero()
        };
        let dvdx = if i > 0 && i < self.v.rows() - 1 && j < self.v.cols() - 1 {
            let v_r = (self.v[(i + 1, j)] + self.v[(i + 1, j + 1)]) / two;
            let v_l = (self.v[(i - 1, j)] + self.v[(i - 1, j + 1)]) / two;
            (v_r - v_l) / (two * dx)
        } else {
            T::zero()
        };

        let d12 = (dudy + dvdx) / two;
        let gamma_sq = two * (dudx * dudx + dvdy * dvdy + two * d12 * d12);
        Float::sqrt(Float::max(gamma_sq, T::zero()))
    }

    /// Update viscosity field from shear rate using the active `BloodModel`.
    ///
    /// `alpha_mu` under-relaxes the update: μ_new = (1−α)·μ_old + α·μ_computed.
    /// For non-Newtonian models (Casson, Carreau-Yasuda) this prevents viscosity
    /// oscillation between SIMPLE iterations (Patankar 1980, §6.7).
    /// For Newtonian blood the factor has no effect (μ_computed is constant).
    pub fn update_viscosity(
        &mut self,
        grid: &StaggeredGrid2D<T>,
        blood: &BloodModel<T>,
        alpha_mu: T,
    ) {
        let one = T::one();
        for i in 0..grid.nx {
            for j in 0..grid.ny {
                let gamma = self.compute_shear_rate(i, j, grid.dx, grid.dy);
                self.gamma_dot[(i, j)] = gamma;
                let mu_computed = match blood {
                    BloodModel::Casson(m) => m.apparent_viscosity(gamma),
                    BloodModel::CarreauYasuda(m) => m.apparent_viscosity(gamma),
                    BloodModel::Newtonian(mu) => *mu,
                };
                self.mu[(i, j)] = self.mu[(i, j)] * (one - alpha_mu) + mu_computed * alpha_mu;
            }
        }
    }

    /// CFL-based maximum velocity magnitude.
    pub fn max_velocity(&self) -> T {
        let max_u = self
            .u
            .iter()
            .map(|&v| Float::abs(v))
            .max_by(|a: &T, b: &T| a.partial_cmp(b).unwrap())
            .unwrap_or(T::zero());
        let max_v = self
            .v
            .iter()
            .map(|&v| Float::abs(v))
            .max_by(|a: &T, b: &T| a.partial_cmp(b).unwrap())
            .unwrap_or(T::zero());
        Float::max(max_u, max_v)
    }
}
