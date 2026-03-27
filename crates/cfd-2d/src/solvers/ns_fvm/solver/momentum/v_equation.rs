//! v-momentum Gauss-Seidel solver.

use crate::error::Error;
use crate::solvers::ns_fvm::boundary::BoundaryCondition;
use crate::solvers::ns_fvm::solver::NavierStokesSolver2D;
use crate::grid::array2d::Array2D;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

impl<T: RealField + Copy + Float + FromPrimitive> NavierStokesSolver2D<T> {
    /// Solves the v-momentum equation on the staggered grid.
    ///
    /// Applies convection, diffusion, and boundary conditions to compute the
    /// intermediate vertical velocity field (`v*`) prior to pressure correction.
    pub fn solve_v_momentum(
        &mut self,
        bc_south: &BoundaryCondition<T>,
        bc_north: &BoundaryCondition<T>,
    ) -> Result<(), Error> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let rho = self.density;
        let alpha = self.config.alpha_u;
        let one = T::one();
        let half = one / (one + one);
        let zero = T::zero();

        let u_old = self.field.u.clone();
        let v_old = self.field.v.clone();
        let mut a_p_v = Array2D::new(nx, ny + 1, T::one());

        for i in 0..nx {
            for j in 1..ny {
                let dy_jm1 = self.grid.dy_at(j - 1);
                let dy_j = self.grid.dy_at(j);
                let dy_cv = (dy_jm1 + dy_j) * half;

                // Effective viscosity = molecular + turbulent.
                let mu_n = if j < ny {
                    self.field.mu[(i, j)] + rho * self.field.nu_t[(i, j)]
                } else {
                    self.field.mu[(i, ny - 1)] + rho * self.field.nu_t[(i, ny - 1)]
                };
                let mu_s = if j > 0 {
                    self.field.mu[(i, j - 1)] + rho * self.field.nu_t[(i, j - 1)]
                } else {
                    self.field.mu[(i, 0)] + rho * self.field.nu_t[(i, 0)]
                };
                let mu_face = (mu_n + mu_s) / (one + one);

                let d_e = mu_face * dy_cv / dx;
                let d_w = d_e;
                let d_n = mu_face * dx / dy_j;
                let d_s = mu_face * dx / dy_jm1;

                let u_e = if i < nx - 1 {
                    (u_old[(i + 1, j - 1)] + u_old[(i + 1, j)]) * half
                } else {
                    zero
                };
                let u_w = if i > 0 {
                    (u_old[(i, j - 1)] + u_old[(i, j)]) * half
                } else {
                    zero
                };
                let f_e = rho * u_e * dy_cv;
                let f_w = rho * u_w * dy_cv;

                let f_n = if j < ny {
                    rho * (v_old[(i, j)] + v_old[(i, j + 1)]) * half * dx
                } else {
                    rho * v_old[(i, j)] * dx
                };
                let f_s = if j > 1 {
                    rho * (v_old[(i, j - 1)] + v_old[(i, j)]) * half * dx
                } else {
                    rho * v_old[(i, j)] * dx
                };

                // Hybrid convection: east/west remain upwind (cross-stream),
                // north/south use hybrid (streamwise for v-equation).
                let a_e = Float::max(d_e - f_e * half, zero) + Float::max(-f_e, zero);
                let a_w = Float::max(d_w + f_w * half, zero) + Float::max(f_w, zero);
                let two = one + one;
                let pe_n = Float::abs(f_n) / d_n;
                let pe_s = Float::abs(f_s) / d_s;
                let a_n = if pe_n <= two {
                    d_n - f_n * half
                } else {
                    d_n + Float::max(zero, -f_n)
                };
                let a_s = if pe_s <= two {
                    d_s + f_s * half
                } else {
                    d_s + Float::max(zero, f_s)
                };

                let mut a_p = a_e + a_w + a_n + a_s + (f_e - f_w) + (f_n - f_s);
                if a_p < T::from_f64(1e-30).unwrap_or(zero) {
                    a_p = T::from_f64(1e-10).unwrap_or(one);
                }

                let p_bot = if j > 0 && j - 1 < ny {
                    self.field.p[(i, j - 1)]
                } else {
                    self.field.p[(i, 0)]
                };
                let p_top = if j < ny {
                    self.field.p[(i, j)]
                } else {
                    self.field.p[(i, ny - 1)]
                };
                let pressure_source = (p_bot - p_top) * dx;

                let v_e = if i + 1 < nx {
                    self.field.v[(i + 1, j)]
                } else {
                    self.field.v[(i, j)]
                };
                let v_w = if i >= 1 {
                    self.field.v[(i - 1, j)]
                } else {
                    self.field.v[(i, j)]
                };
                let v_n = if j < ny {
                    self.field.v[(i, j + 1)]
                } else {
                    self.field.v[(i, j)]
                };
                let v_s = if j >= 1 {
                    self.field.v[(i, j - 1)]
                } else {
                    self.field.v[(i, j)]
                };

                let is_fluid = if j > 0 && j < ny {
                    self.field.mask[(i, j - 1)] && self.field.mask[(i, j)]
                } else {
                    true
                };

                if is_fluid {
                    let v_star =
                        (a_e * v_e + a_w * v_w + a_n * v_n + a_s * v_s + pressure_source) / a_p;
                    self.field.v[(i, j)] = self.field.v[(i, j)] * (one - alpha) + v_star * alpha;
                    a_p_v[(i, j)] = a_p;
                } else {
                    self.field.v[(i, j)] = zero;
                    a_p_v[(i, j)] = one;
                }
            }
        }

        if let BoundaryCondition::Wall { .. } = bc_south {
            for i in 0..nx {
                self.field.v[(i, 0)] = zero;
            }
        } else if bc_south.is_neumann() {
            if let BoundaryCondition::Symmetry = bc_south {
                for i in 0..nx {
                    self.field.v[(i, 0)] = zero;
                }
            }
        }

        if let BoundaryCondition::Wall { .. } = bc_north {
            for i in 0..nx {
                self.field.v[(i, ny)] = zero;
            }
        } else if let BoundaryCondition::Symmetry = bc_north {
            for i in 0..nx {
                self.field.v[(i, ny)] = zero;
            }
        }

        self.a_p_v = a_p_v;
        Ok(())
    }
}
