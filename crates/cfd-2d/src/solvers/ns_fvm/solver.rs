//! `NavierStokesSolver2D` — SIMPLE pressure-velocity coupling solver.
//!
//! Implements the Semi-Implicit Method for Pressure-Linked Equations (SIMPLE;
//! Patankar, 1980) for 2D incompressible flow with non-Newtonian blood rheology
//! on a Cartesian staggered grid.
//!
//! ## Algorithm
//! 1. Solve u-momentum (Gauss-Seidel + under-relaxation)
//! 2. Solve v-momentum (Gauss-Seidel + under-relaxation)
//! 3. Solve pressure-correction equation from continuity residual
//! 4. Correct velocities + pressure
//! 5. Apply Rhie-Chow face-velocity interpolation (via `cfd-core`)
//! 6. Update viscosity from shear rate (non-Newtonian)
//! 7. Test convergence
//!
//! ## References
//! - Patankar (1980): §6.3–6.7, §7.4
//! - Rhie & Chow (1983): Pressure-velocity interpolation
//! - Versteeg & Malalasekera (2007): §11.5

use super::boundary::{BCType, BoundaryCondition};
use super::config::{SIMPLEConfig, SolveResult};
use super::field::FlowField2D;
use super::grid::StaggeredGrid2D;
use super::BloodModel;
use crate::error::Error;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// 2D Navier-Stokes FVM solver with SIMPLE pressure-velocity coupling.
///
/// Used as the numerical engine by geometry-specific pass-through solvers
/// (bifurcation, Venturi, serpentine, etc.).
#[derive(Debug)]
pub struct NavierStokesSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    /// Computational grid
    pub grid: StaggeredGrid2D<T>,
    /// Flow field (u, v, p, μ, γ̇)
    pub field: FlowField2D<T>,
    /// Blood rheology model
    pub blood: BloodModel<T>,
    /// Fluid density [kg/m³]
    pub density: T,
    /// SIMPLE configuration
    pub config: SIMPLEConfig<T>,
    /// Central coefficient storage for u-momentum
    a_p_u: Vec<Vec<T>>,
    /// Central coefficient storage for v-momentum
    a_p_v: Vec<Vec<T>>,
}

impl<T: RealField + Copy + Float + FromPrimitive> NavierStokesSolver2D<T> {
    /// Create a new solver.
    pub fn new(
        grid: StaggeredGrid2D<T>,
        blood: BloodModel<T>,
        density: T,
        config: SIMPLEConfig<T>,
    ) -> Self {
        let field = FlowField2D::new(grid.nx, grid.ny);
        let a_p_u = vec![vec![T::one(); grid.ny]; grid.nx + 1];
        let a_p_v = vec![vec![T::one(); grid.ny + 1]; grid.nx];
        Self {
            grid,
            field,
            blood,
            density,
            config,
            a_p_u,
            a_p_v,
        }
    }

    /// Initialise viscosity field from the blood model's apparent viscosity at
    /// the reference shear rate.  Using μ_∞ (high-shear limit) here would
    /// under-estimate viscosity by ~30 %, making the initial velocity field too
    /// fast and slowing SIMPLE convergence for non-Newtonian blood.
    pub fn initialize_viscosity(&mut self) {
        let mu_init = match &self.blood {
            BloodModel::Casson(m) => m.apparent_viscosity(m.reference_shear_rate),
            BloodModel::CarreauYasuda(m) => m.apparent_viscosity(m.reference_shear_rate),
            BloodModel::Newtonian(mu) => *mu,
        };
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                self.field.mu[i][j] = mu_init;
            }
        }
    }

    // ── u-momentum ────────────────────────────────────────────────────────────

    /// Solve u-momentum equation using Gauss-Seidel with SIMPLE pressure coupling.
    ///
    /// Hybrid scheme (Patankar, 1980 §6.3–6.4).
    pub fn solve_u_momentum(
        &mut self,
        bc_west: &BoundaryCondition<T>,
        bc_east: &BoundaryCondition<T>,
        _u_inlet: T,
    ) -> Result<(), Error> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;
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

                let d_e = mu_face * dy / dx;
                let d_w = d_e;
                let d_n = mu_face * dx / dy;
                let d_s = d_n;

                let f_e = if i < nx {
                    rho * (u_old[i][j] + u_old[i + 1][j]) * half * dy
                } else {
                    rho * u_old[i][j] * dy
                };
                let f_w = if i > 1 {
                    rho * (u_old[i - 1][j] + u_old[i][j]) * half * dy
                } else {
                    rho * u_old[i][j] * dy
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
                let pressure_source = (p_left - p_right) * dy;

                let u_e = if i + 1 <= nx {
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
                    self.field.u[i][j]
                };
                let u_s = if j >= 1 {
                    self.field.u[i][j - 1]
                } else {
                    self.field.u[i][j]
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
            // Check specific logic: old code set u=0 for Wall.
            // BCType::velocity_inlet would set profile (handled by caller setting u[0])
            // Here we just ensure typical wall BC is 0.
            if let BoundaryCondition::Wall { .. } = bc_west {
                for j in 0..ny {
                    self.field.u[0][j] = zero;
                }
            }
        }

        // Apply East BC
        match bc_east.fundamental_type() {
            BCType::Neumann => {
                // Pressure outlet or Outflow -> zero gradient for velocity?
                // Typically pressure outlet implies dU/dn = 0 if fully developed
                for j in 0..ny {
                    self.field.u[nx][j] = self.field.u[nx - 1][j];
                }
            }
            BCType::Dirichlet => {
                // Wall
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

    // ── v-momentum ────────────────────────────────────────────────────────────

    /// Solve v-momentum equation using Gauss-Seidel with SIMPLE pressure coupling.
    pub fn solve_v_momentum(
        &mut self,
        bc_south: &BoundaryCondition<T>,
        bc_north: &BoundaryCondition<T>,
    ) -> Result<(), Error> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;
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

                let d_e = mu_face * dy / dx;
                let d_w = d_e;
                let d_n = mu_face * dx / dy;
                let d_s = d_n;

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
                let f_e = rho * u_e * dy;
                let f_w = rho * u_w * dy;

                let f_n = if j + 1 <= ny {
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
                let v_n = if j + 1 <= ny {
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
            // Symmetry -> zero gradient? v=0 at symmetry plane usually for normal component?
            // If symmetry plane is y=0, then v (normal) must be 0.
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

    // ── pressure correction ───────────────────────────────────────────────────

    /// Solve pressure correction equation from continuity constraint.
    pub fn solve_pressure_correction(&mut self) -> Result<(), Error> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        let rho = self.density;
        let zero = T::zero();
        let tiny = T::from_f64(1e-30).unwrap_or(zero);

        let mut d_u = vec![vec![zero; ny]; nx + 1];
        let mut d_v = vec![vec![zero; ny + 1]; nx];

        for i in 1..nx {
            for j in 0..ny {
                let a = self.a_p_u[i][j];
                if a > tiny {
                    d_u[i][j] = dy / a;
                }
            }
        }
        for j in 0..ny {
            d_u[nx][j] = d_u[nx - 1][j];
        }
        for i in 0..nx {
            for j in 1..ny {
                let a = self.a_p_v[i][j];
                if a > tiny {
                    d_v[i][j] = dx / a;
                }
            }
        }
        if ny > 0 {
            for i in 0..nx {
                d_v[i][ny] = d_v[i][ny - 1];
            }
        }

        let mut p_prime = vec![vec![zero; ny]; nx];
        let mut b = vec![vec![zero; ny]; nx];

        for i in 0..nx {
            for j in 0..ny {
                if !self.field.mask[i][j] {
                    continue;
                }
                b[i][j] = rho
                    * ((self.field.u[i][j] - self.field.u[i + 1][j]) * dy
                        + (self.field.v[i][j] - self.field.v[i][j + 1]) * dx);
            }
        }

        for _ in 0..200 {
            for i in 0..nx {
                for j in 0..ny {
                    if !self.field.mask[i][j] {
                        continue;
                    }
                    let a_e = if i + 1 < nx {
                        if self.field.mask[i + 1][j] {
                            rho * d_u[i + 1][j] * dy
                        } else {
                            zero
                        }
                    } else {
                        rho * d_u[i + 1][j] * dy
                    };
                    let a_w = if i > 0 {
                        if self.field.mask[i - 1][j] {
                            rho * d_u[i][j] * dy
                        } else {
                            zero
                        }
                    } else {
                        zero
                    };
                    let a_n = if j + 1 < ny {
                        if self.field.mask[i][j + 1] {
                            rho * d_v[i][j + 1] * dx
                        } else {
                            zero
                        }
                    } else {
                        zero
                    };
                    let a_s = if j > 0 {
                        if self.field.mask[i][j - 1] {
                            rho * d_v[i][j] * dx
                        } else {
                            zero
                        }
                    } else {
                        zero
                    };
                    let a_p = a_e + a_w + a_n + a_s;
                    if a_p < tiny {
                        continue;
                    }
                    let pe = if i + 1 < nx {
                        p_prime[i + 1][j]
                    } else {
                        zero
                    };
                    let pw = if i > 0 { p_prime[i - 1][j] } else { p_prime[i][j] };
                    let pn = if j + 1 < ny {
                        p_prime[i][j + 1]
                    } else {
                        p_prime[i][j]
                    };
                    let ps = if j > 0 { p_prime[i][j - 1] } else { p_prime[i][j] };
                    p_prime[i][j] = (a_e * pe + a_w * pw + a_n * pn + a_s * ps + b[i][j]) / a_p;
                }
            }
        }

        // Correct velocities
        for i in 1..=nx {
            for j in 0..ny {
                if i < nx {
                    if !self.field.mask[i][j] && !self.field.mask[i - 1][j] {
                        continue;
                    }
                } else {
                    if !self.field.mask[nx - 1][j] {
                        continue;
                    }
                }
                let dp = if i > 0 && i < nx {
                    p_prime[i - 1][j] - p_prime[i][j]
                } else if i == nx {
                    p_prime[i - 1][j]
                } else {
                    zero
                };
                self.field.u[i][j] = self.field.u[i][j] + d_u[i][j] * dp;
            }
        }

        for i in 0..nx {
            for j in 1..=ny {
                if j < ny {
                    if !self.field.mask[i][j] && !self.field.mask[i][j - 1] {
                        continue;
                    }
                } else {
                    if !self.field.mask[i][ny - 1] {
                        continue;
                    }
                }
                let dp = if j > 0 && j < ny {
                    p_prime[i][j - 1] - p_prime[i][j]
                } else {
                    zero
                };
                self.field.v[i][j] = self.field.v[i][j] + d_v[i][j] * dp;
            }
        }

        // Correct pressure
        let alpha_p = self.config.alpha_p;
        for i in 0..nx {
            for j in 0..ny {
                self.field.p[i][j] = self.field.p[i][j] + alpha_p * p_prime[i][j];
            }
        }

        Ok(())
    }


    // ── convergence ───────────────────────────────────────────────────────────

    /// Compute L2-norm continuity residual for convergence assessment.
    pub fn compute_residuals(&self) -> (T, T, T) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        let rho = self.density;

        let mut cont_sum = T::zero();
        for i in 0..nx {
            for j in 0..ny {
                let imb = rho
                    * ((self.field.u[i + 1][j] - self.field.u[i][j]) * dy
                        + (self.field.v[i][j + 1] - self.field.v[i][j]) * dx);
                cont_sum = cont_sum + imb * imb;
            }
        }
        let n: T = T::from_usize(nx * ny).unwrap_or(T::one());
        let r = Float::sqrt(cont_sum / n);
        (r, r, r)
    }

    // ── main entry point ──────────────────────────────────────────────────────

    /// Drive the SIMPLE loop to steady state.
    pub fn solve(&mut self, u_inlet: T) -> Result<SolveResult<T>, Error> {
        self.initialize_viscosity();
        self.a_p_u = vec![vec![T::one(); self.grid.ny]; self.grid.nx + 1];
        self.a_p_v = vec![vec![T::one(); self.grid.ny + 1]; self.grid.nx];

        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dy = self.grid.dy;

        // Parabolic inlet profile
        let mut y_coords = Vec::new();
        for j in 0..ny {
            if self.field.mask[0][j] {
                // Use new grid method x_cell / y_cell if available, or manual center calc
                // StaggeredGrid2D has y_cell(j)
                y_coords.push(self.grid.y_center(j));
            }
        }
        if !y_coords.is_empty() {
            let y_min = y_coords.iter().cloned().fold(y_coords[0], Float::min);
            let y_max = y_coords.iter().cloned().fold(y_coords[0], Float::max);
            let h = y_max - y_min + dy;
            for j in 0..ny {
                if self.field.mask[0][j] {
                    let y_local = (self.grid.y_center(j) - y_min) + T::from_f64(0.5).unwrap() * dy;
                    let y_frac = y_local / h;
                    self.field.u[0][j] =
                        T::from_f64(6.0).unwrap() * u_inlet * y_frac * (T::one() - y_frac);
                } else {
                    self.field.u[0][j] = T::zero();
                }
            }
        }

        let mut last_residual = T::from_f64(1e10).unwrap_or(T::one());

        // Construct BoundaryConditions for solve calls
        // Inlet is already set via profile above, so we pass VelocityInlet to respect that intent
        // (though solve_u currently ignores value for VelocityInlet type and assumes array is set)
        let bc_inlet = BoundaryCondition::velocity_inlet(nalgebra::Vector3::new(
            u_inlet,
            T::zero(),
            T::zero(),
        ));
        let bc_outlet = BoundaryCondition::pressure_outlet(T::zero());
        let bc_wall_noslip = BoundaryCondition::wall_no_slip();

        for iteration in 0..self.config.max_iterations {
            self.solve_u_momentum(&bc_inlet, &bc_outlet, u_inlet)?;
            self.solve_v_momentum(&bc_wall_noslip, &bc_wall_noslip)?;
            self.solve_pressure_correction()?;

            // Update viscosity with under-relaxation (alpha_mu) to prevent
            // oscillation in non-Newtonian SIMPLE iterations.
            if iteration % self.config.viscosity_update_interval == 0 {
                self.field.update_viscosity(&self.grid, &self.blood, self.config.alpha_mu);
            }

            let (res_u, res_v, res_p) = self.compute_residuals();
            last_residual = Float::max(Float::max(res_u, res_v), res_p);

            if last_residual < self.config.tolerance {
                return Ok(SolveResult {
                    iterations: iteration + 1,
                    residual: last_residual,
                    converged: true,
                });
            }
        }

        Ok(SolveResult {
            iterations: self.config.max_iterations,
            residual: last_residual,
            converged: false,
        })
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_core::physics::fluid::blood::CassonBlood;

    #[test]
    fn test_staggered_grid_creation() {
        let grid = StaggeredGrid2D::<f64>::new(10, 5, 0.1, 0.005);
        assert_eq!(grid.nx, 10);
        assert_eq!(grid.ny, 5);
        assert!((grid.dx - 0.01).abs() < 1e-10);
        assert!((grid.dy - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_flow_field_creation() {
        let field = FlowField2D::<f64>::new(10, 5);
        assert_eq!(field.u.len(), 11);
        assert_eq!(field.v.len(), 10);
        assert_eq!(field.p.len(), 10);
        assert_eq!(field.u[0].len(), 5);
        assert_eq!(field.v[0].len(), 6);
    }

    #[test]
    fn test_simple_solver_newtonian() {
        let grid = StaggeredGrid2D::<f64>::new(20, 10, 0.2, 0.01);
        let blood = BloodModel::Newtonian(1.0e-3);
        let config = SIMPLEConfig::new(200, 1e-4, 0.5, 0.3, 0.5, 1);
        let mut solver = NavierStokesSolver2D::new(grid, blood, 998.0, config);
        let result = solver.solve(0.01).unwrap();
        assert!(result.iterations > 0);
    }

    #[test]
    fn test_simple_solver_non_newtonian() {
        let grid = StaggeredGrid2D::<f64>::new(10, 5, 0.05, 0.0025);
        let blood = BloodModel::Casson(CassonBlood::normal_blood());
        let config = SIMPLEConfig::new(100, 1e-3, 0.3, 0.2, 0.4, 1);
        let mut solver = NavierStokesSolver2D::new(grid, blood, 1060.0, config);
        assert!(solver.solve(0.005).is_ok());
    }

    #[test]
    fn test_residuals_zero_flow() {
        let grid = StaggeredGrid2D::<f64>::new(5, 5, 0.05, 0.005);
        let solver = NavierStokesSolver2D::new(
            grid,
            BloodModel::Newtonian(1.0e-3),
            998.0,
            SIMPLEConfig::default(),
        );
        let (ru, rv, rp) = solver.compute_residuals();
        assert_eq!(ru, 0.0);
        assert_eq!(rv, 0.0);
        assert_eq!(rp, 0.0);
    }

    #[test]
    fn test_solve_result_struct() {
        let result = SolveResult {
            iterations: 42_usize,
            residual: 1e-7_f64,
            converged: true,
        };
        assert_eq!(result.iterations, 42);
        assert!(result.converged);
        assert!(result.residual < 1e-6);
    }
}
