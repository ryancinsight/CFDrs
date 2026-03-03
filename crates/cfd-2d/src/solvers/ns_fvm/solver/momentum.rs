//! u-momentum and v-momentum Gauss-Seidel solvers with SIMPLE coupling.

use super::NavierStokesSolver2D;
use crate::error::Error;
use crate::solvers::ns_fvm::boundary::{BCType, BoundaryCondition};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

impl<T: RealField + Copy + Float + FromPrimitive> NavierStokesSolver2D<T> {
    /// Solve u-momentum equation using Gauss-Seidel with SIMPLE pressure coupling.
    ///
    /// Hybrid scheme (Patankar, 1980 §6.3–6.4).
    /// Supports non-uniform y-spacing via `grid.dy_at(j)`.
    pub fn solve_u_momentum(
        &mut self,
        bc_west: &BoundaryCondition<T>,
        bc_east: &BoundaryCondition<T>,
        _u_inlet: T,
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
        let mut a_p_u = vec![vec![T::one(); ny]; nx + 1];

        for i in 1..nx {
            for j in 0..ny {
                let dy_j = self.grid.dy_at(j);

                let mu_e = if i < nx {
                    self.field.mu[i][j]
                } else {
                    self.field.mu[nx - 1][j]
                };
                let mu_w = if i > 0 {
                    self.field.mu[i - 1][j]
                } else {
                    self.field.mu[0][j]
                };
                let mu_face = (mu_e + mu_w) / (one + one);

                // Diffusion: east/west faces have area dy_j; north/south use
                // centre-to-centre distance for the gradient.
                let d_e = mu_face * dy_j / dx;
                let d_w = d_e;
                let d_n = if j + 1 < ny {
                    mu_face * dx / self.grid.dy_face(j)
                } else {
                    // Half-cell distance to the north wall (no-slip).
                    mu_face * dx / (dy_j * half)
                };
                let d_s = if j > 0 {
                    mu_face * dx / self.grid.dy_face(j - 1)
                } else {
                    // Half-cell distance to the south wall (no-slip).
                    mu_face * dx / (dy_j * half)
                };

                let f_e = if i < nx {
                    rho * (u_old[i][j] + u_old[i + 1][j]) * half * dy_j
                } else {
                    rho * u_old[i][j] * dy_j
                };
                let f_w = if i > 1 {
                    rho * (u_old[i - 1][j] + u_old[i][j]) * half * dy_j
                } else {
                    rho * u_old[i][j] * dy_j
                };

                let v_n = if j < ny - 1 {
                    (v_old[i - 1][j + 1] + v_old[i][j + 1]) * half
                } else {
                    T::zero()
                };
                let v_s = if j > 0 {
                    (v_old[i - 1][j] + v_old[i][j]) * half
                } else {
                    T::zero()
                };
                let f_n = rho * v_n * dx;
                let f_s = rho * v_s * dx;

                let a_e = Float::max(d_e - f_e * half, zero) + Float::max(-f_e, zero);
                let a_w = Float::max(d_w + f_w * half, zero) + Float::max(f_w, zero);
                let a_n = Float::max(d_n - f_n * half, zero) + Float::max(-f_n, zero);
                let a_s = Float::max(d_s + f_s * half, zero) + Float::max(f_s, zero);

                let mut a_p = a_e + a_w + a_n + a_s + (f_e - f_w) + (f_n - f_s);
                if a_p < T::from_f64(1e-30).unwrap_or(zero) {
                    a_p = T::from_f64(1e-10).unwrap_or(one);
                }

                let p_left = if i > 0 && i - 1 < nx {
                    self.field.p[i - 1][j]
                } else {
                    self.field.p[0][j]
                };
                let p_right = if i < nx {
                    self.field.p[i][j]
                } else {
                    self.field.p[nx - 1][j]
                };
                let pressure_source = (p_left - p_right) * dy_j;

                let u_e = if i < nx {
                    self.field.u[i + 1][j]
                } else {
                    self.field.u[i][j]
                };
                let u_w = if i >= 1 {
                    self.field.u[i - 1][j]
                } else {
                    self.field.u[i][j]
                };
                let u_n = if j + 1 < ny {
                    self.field.u[i][j + 1]
                } else {
                    // No-slip wall: u = 0 at the north boundary.
                    zero
                };
                let u_s = if j >= 1 {
                    self.field.u[i][j - 1]
                } else {
                    // No-slip wall: u = 0 at the south boundary.
                    zero
                };

                let is_fluid = if i > 0 && i < nx {
                    self.field.mask[i - 1][j] && self.field.mask[i][j]
                } else {
                    true
                };

                if is_fluid {
                    let u_star =
                        (a_e * u_e + a_w * u_w + a_n * u_n + a_s * u_s + pressure_source) / a_p;
                    self.field.u[i][j] = self.field.u[i][j] * (one - alpha) + u_star * alpha;
                    a_p_u[i][j] = a_p;
                } else {
                    self.field.u[i][j] = zero;
                    a_p_u[i][j] = one;
                }
            }
        }

        // Apply West BC
        if bc_west.is_dirichlet() {
            if let BoundaryCondition::Wall { .. } = bc_west {
                for j in 0..ny {
                    self.field.u[0][j] = zero;
                }
            }
        }

        // Apply East BC
        match bc_east.fundamental_type() {
            BCType::Neumann => {
                for j in 0..ny {
                    self.field.u[nx][j] = self.field.u[nx - 1][j];
                }
            }
            BCType::Dirichlet => {
                if let BoundaryCondition::Wall { .. } = bc_east {
                    for j in 0..ny {
                        self.field.u[nx][j] = zero;
                    }
                }
            }
            _ => {}
        }

        self.a_p_u = a_p_u;
        Ok(())
    }

    /// Solve v-momentum equation using Gauss-Seidel with SIMPLE pressure coupling.
    /// Supports non-uniform y-spacing via `grid.dy_at(j)`.
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
        let mut a_p_v = vec![vec![T::one(); ny + 1]; nx];

        for i in 0..nx {
            for j in 1..ny {
                // v-face j sits between pressure cells j-1 and j.
                // The control volume around v-face j has height = dy_face(j-1)
                // = (dy[j-1] + dy[j]) / 2.
                let dy_jm1 = self.grid.dy_at(j - 1);
                let dy_j = self.grid.dy_at(j);
                let dy_cv = (dy_jm1 + dy_j) * half;

                let mu_n = if j < ny {
                    self.field.mu[i][j]
                } else {
                    self.field.mu[i][ny - 1]
                };
                let mu_s = if j > 0 {
                    self.field.mu[i][j - 1]
                } else {
                    self.field.mu[i][0]
                };
                let mu_face = (mu_n + mu_s) / (one + one);

                let d_e = mu_face * dy_cv / dx;
                let d_w = d_e;
                let d_n = mu_face * dx / dy_j;
                let d_s = mu_face * dx / dy_jm1;

                let u_e = if i < nx - 1 {
                    (u_old[i + 1][j - 1] + u_old[i + 1][j]) * half
                } else {
                    zero
                };
                let u_w = if i > 0 {
                    (u_old[i][j - 1] + u_old[i][j]) * half
                } else {
                    zero
                };
                let f_e = rho * u_e * dy_cv;
                let f_w = rho * u_w * dy_cv;

                let f_n = if j < ny {
                    rho * (v_old[i][j] + v_old[i][j + 1]) * half * dx
                } else {
                    rho * v_old[i][j] * dx
                };
                let f_s = if j > 1 {
                    rho * (v_old[i][j - 1] + v_old[i][j]) * half * dx
                } else {
                    rho * v_old[i][j] * dx
                };

                let a_e = Float::max(d_e - f_e * half, zero) + Float::max(-f_e, zero);
                let a_w = Float::max(d_w + f_w * half, zero) + Float::max(f_w, zero);
                let a_n = Float::max(d_n - f_n * half, zero) + Float::max(-f_n, zero);
                let a_s = Float::max(d_s + f_s * half, zero) + Float::max(f_s, zero);

                let mut a_p = a_e + a_w + a_n + a_s + (f_e - f_w) + (f_n - f_s);
                if a_p < T::from_f64(1e-30).unwrap_or(zero) {
                    a_p = T::from_f64(1e-10).unwrap_or(one);
                }

                let p_bot = if j > 0 && j - 1 < ny {
                    self.field.p[i][j - 1]
                } else {
                    self.field.p[i][0]
                };
                let p_top = if j < ny {
                    self.field.p[i][j]
                } else {
                    self.field.p[i][ny - 1]
                };
                let pressure_source = (p_bot - p_top) * dx;

                let v_e = if i + 1 < nx {
                    self.field.v[i + 1][j]
                } else {
                    self.field.v[i][j]
                };
                let v_w = if i >= 1 {
                    self.field.v[i - 1][j]
                } else {
                    self.field.v[i][j]
                };
                let v_n = if j < ny {
                    self.field.v[i][j + 1]
                } else {
                    self.field.v[i][j]
                };
                let v_s = if j >= 1 {
                    self.field.v[i][j - 1]
                } else {
                    self.field.v[i][j]
                };

                let is_fluid = if j > 0 && j < ny {
                    self.field.mask[i][j - 1] && self.field.mask[i][j]
                } else {
                    true
                };

                if is_fluid {
                    let v_star =
                        (a_e * v_e + a_w * v_w + a_n * v_n + a_s * v_s + pressure_source) / a_p;
                    self.field.v[i][j] = self.field.v[i][j] * (one - alpha) + v_star * alpha;
                    a_p_v[i][j] = a_p;
                } else {
                    self.field.v[i][j] = zero;
                    a_p_v[i][j] = one;
                }
            }
        }

        // Apply South BC
        if let BoundaryCondition::Wall { .. } = bc_south {
            for i in 0..nx {
                self.field.v[i][0] = zero;
            }
        } else if bc_south.is_neumann() {
            if let BoundaryCondition::Symmetry = bc_south {
                for i in 0..nx {
                    self.field.v[i][0] = zero;
                }
            }
        }

        // Apply North BC
        if let BoundaryCondition::Wall { .. } = bc_north {
            for i in 0..nx {
                self.field.v[i][ny] = zero;
            }
        } else if let BoundaryCondition::Symmetry = bc_north {
            for i in 0..nx {
                self.field.v[i][ny] = zero;
            }
        }

        self.a_p_v = a_p_v;
        Ok(())
    }
}
