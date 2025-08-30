//! Velocity predictor step for PISO algorithm

use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;

// Named constants for numerical operations
const HALF: f64 = 0.5;
const TWO: f64 = 2.0;

/// Velocity predictor for PISO algorithm
pub struct VelocityPredictor<T: RealField + Copy> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Grid spacing
    dx: T,
    dy: T,
    /// Under-relaxation factor
    relaxation_factor: T,
}

impl<T: RealField + Copy + FromPrimitive + Copy> VelocityPredictor<T> {
    /// Create new velocity predictor
    pub fn new(grid: &StructuredGrid2D<T>, relaxation_factor: T) -> Self {
        Self {
            nx: grid.nx,
            ny: grid.ny,
            dx: grid.dx,
            dy: grid.dy,
            relaxation_factor,
        }
    }

    /// Predict velocity field without pressure gradient
    /// This is the first step in PISO algorithm
    pub fn predict(&self, fields: &mut SimulationFields<T>, dt: T) -> Result<()> {
        let mut u_star = Field2D::new(self.nx, self.ny, T::zero());
        let mut v_star = Field2D::new(self.nx, self.ny, T::zero());

        // Solve momentum equations without pressure gradient
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                // Get velocities at cell faces (linear interpolation)
                let half = T::from_f64(HALF).unwrap_or_else(T::one);
                let ue = (fields.u.at(i, j) + fields.u.at(i + 1, j)) * half;
                let uw = (fields.u.at(i - 1, j) + fields.u.at(i, j)) * half;
                let un = (fields.u.at(i, j) + fields.u.at(i, j + 1)) * half;
                let us = (fields.u.at(i, j - 1) + fields.u.at(i, j)) * half;

                let ve = (fields.v.at(i, j) + fields.v.at(i + 1, j)) * half;
                let vw = (fields.v.at(i - 1, j) + fields.v.at(i, j)) * half;
                let vn = (fields.v.at(i, j) + fields.v.at(i, j + 1)) * half;
                let vs = (fields.v.at(i, j - 1) + fields.v.at(i, j)) * half;

                // Convective terms (using upwind scheme)
                let conv_u = self.calculate_convection_u(fields, i, j, ue, uw, un, us);
                let conv_v = self.calculate_convection_v(fields, i, j, ve, vw, vn, vs);

                // Diffusive terms (central difference)
                let diff_u = self.calculate_diffusion_u(fields, i, j);
                let diff_v = self.calculate_diffusion_v(fields, i, j);

                // Time integration (explicit Euler for now)
                let u_current = fields.u.at(i, j);
                let v_current = fields.v.at(i, j);

                if let Some(u) = u_star.at_mut(i, j) {
                    *u = u_current + dt * (diff_u - conv_u) / fields.density.at(i, j);
                }
                if let Some(v) = v_star.at_mut(i, j) {
                    *v = v_current + dt * (diff_v - conv_v) / fields.density.at(i, j);
                }
            }
        }

        // Apply under-relaxation
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let u_relaxed = self.relaxation_factor * u_star.at(i, j)
                    + (T::one() - self.relaxation_factor) * fields.u.at(i, j);
                let v_relaxed = self.relaxation_factor * v_star.at(i, j)
                    + (T::one() - self.relaxation_factor) * fields.v.at(i, j);

                if let Some(u) = fields.u.at_mut(i, j) {
                    *u = u_relaxed;
                }
                if let Some(v) = fields.v.at_mut(i, j) {
                    *v = v_relaxed;
                }
            }
        }

        Ok(())
    }

    /// Calculate convection term for u-velocity
    fn calculate_convection_u(
        &self,
        fields: &SimulationFields<T>,
        i: usize,
        j: usize,
        ue: T,
        uw: T,
        un: T,
        us: T,
    ) -> T {
        let rho = fields.density.at(i, j);

        // Upwind scheme for stability
        let fe = if ue > T::zero() {
            rho * ue * fields.u.at(i, j)
        } else {
            rho * ue * fields.u.at(i + 1, j)
        };

        let fw = if uw > T::zero() {
            rho * uw * fields.u.at(i - 1, j)
        } else {
            rho * uw * fields.u.at(i, j)
        };

        let f_n = if un > T::zero() {
            rho * un * fields.u.at(i, j)
        } else {
            rho * un * fields.u.at(i, j + 1)
        };

        let f_s = if us > T::zero() {
            rho * us * fields.u.at(i, j - 1)
        } else {
            rho * us * fields.u.at(i, j)
        };

        (fe - fw) / self.dx + (f_n - f_s) / self.dy
    }

    /// Calculate convection term for v-velocity
    fn calculate_convection_v(
        &self,
        fields: &SimulationFields<T>,
        i: usize,
        j: usize,
        ve: T,
        vw: T,
        vn: T,
        vs: T,
    ) -> T {
        let rho = fields.density.at(i, j);

        // Upwind scheme
        let fe = if ve > T::zero() {
            rho * ve * fields.v.at(i, j)
        } else {
            rho * ve * fields.v.at(i + 1, j)
        };

        let fw = if vw > T::zero() {
            rho * vw * fields.v.at(i - 1, j)
        } else {
            rho * vw * fields.v.at(i, j)
        };

        let f_n = if vn > T::zero() {
            rho * vn * fields.v.at(i, j)
        } else {
            rho * vn * fields.v.at(i, j + 1)
        };

        let f_s = if vs > T::zero() {
            rho * vs * fields.v.at(i, j - 1)
        } else {
            rho * vs * fields.v.at(i, j)
        };

        (fe - fw) / self.dx + (f_n - f_s) / self.dy
    }

    /// Calculate diffusion term for u-velocity
    fn calculate_diffusion_u(&self, fields: &SimulationFields<T>, i: usize, j: usize) -> T {
        let mu = fields.viscosity.at(i, j);
        let two = T::from_f64(TWO).unwrap_or_else(|| T::one() + T::one());

        let d2u_dx2 = (fields.u.at(i + 1, j) - two * fields.u.at(i, j) + fields.u.at(i - 1, j))
            / (self.dx * self.dx);
        let d2u_dy2 = (fields.u.at(i, j + 1) - two * fields.u.at(i, j) + fields.u.at(i, j - 1))
            / (self.dy * self.dy);

        mu * (d2u_dx2 + d2u_dy2)
    }

    /// Calculate diffusion term for v-velocity
    fn calculate_diffusion_v(&self, fields: &SimulationFields<T>, i: usize, j: usize) -> T {
        let mu = fields.viscosity.at(i, j);
        let two = T::from_f64(TWO).unwrap_or_else(|| T::one() + T::one());
        let d2v_dx2 = (fields.v.at(i + 1, j) - two * fields.v.at(i, j) + fields.v.at(i - 1, j))
            / (self.dx * self.dx);
        let d2v_dy2 = (fields.v.at(i, j + 1) - two * fields.v.at(i, j) + fields.v.at(i, j - 1))
            / (self.dy * self.dy);

        mu * (d2v_dx2 + d2v_dy2)
    }
}
