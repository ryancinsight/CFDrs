use crate::scalar;
use crate::solvers::ns_fvm::boundary::{BloodModel, BoundaryCondition};
use crate::solvers::ns_fvm::config::SolveResult;
use crate::solvers::ns_fvm::solver::NavierStokesSolver2D;
use cfd_core::error::Error;
use cfd_core::CfdScalar;
use eunomia::{FloatElement, NumericElement};
use leto::geometry::Vector2;

impl<T: CfdScalar + Copy + FloatElement> NavierStokesSolver2D<T> {
    /// Drive the SIMPLE loop to steady state.
    pub fn solve(&mut self, u_inlet: T) -> Result<SolveResult<T>, Error> {
        self.initialize_viscosity();

        let ny = self.grid.ny;

        // Parabolic inlet profile (supports non-uniform y-spacing)
        let mut y_coords = Vec::with_capacity(ny);
        let mut dy_cells = Vec::with_capacity(ny);
        for j in 0..ny {
            if self.field.mask[(0, j)] {
                y_coords.push(self.grid.y_center(j));
                dy_cells.push(self.grid.dy_at(j));
            }
        }
        if !y_coords.is_empty() {
            let y_min = y_coords
                .iter()
                .copied()
                .fold(y_coords[0], NumericElement::min_scalar);
            let y_max = y_coords
                .iter()
                .copied()
                .fold(y_coords[0], NumericElement::max_scalar);
            // Channel height = span from first to last fluid centre + half-cells at edges.
            let h = y_max - y_min
                + (dy_cells[0] + dy_cells[dy_cells.len().saturating_sub(1)])
                    * <T as FloatElement>::from_f64(0.5);

            let mut discrete_sum: T = scalar::zero();
            let half: T = <T as FloatElement>::from_f64(0.5);
            let six: T = <T as FloatElement>::from_f64(6.0);
            let one: T = scalar::one();
            let zero: T = scalar::zero();

            for j in 0..ny {
                if self.field.mask[(0, j)] {
                    let dy_j = self.grid.dy_at(j);
                    let y_local = (self.grid.y_center(j) - y_min) + half * dy_j;
                    let y_frac = y_local / h;
                    let u_val = six * u_inlet * y_frac * (one - y_frac);
                    self.field.u[(0, j)] = u_val;
                    discrete_sum += u_val * dy_j;
                } else {
                    self.field.u[(0, j)] = zero;
                }
            }

            // Normalise the discrete profile so numerical mass flux matches continuous theory precisely
            let target_sum = u_inlet * h;
            let tiny: T = <T as FloatElement>::from_f64(1e-30);
            if discrete_sum > tiny {
                let normalize_factor = target_sum / discrete_sum;
                for j in 0..ny {
                    if self.field.mask[(0, j)] {
                        self.field.u[(0, j)] *= normalize_factor;
                    }
                }
            }
        }

        let mut last_residual: T = <T as FloatElement>::from_f64(1e10);
        let zero: T = scalar::zero();

        let bc_inlet =
            BoundaryCondition::velocity_inlet(leto::geometry::Vector3::new(u_inlet, zero, zero));
        let bc_outlet = BoundaryCondition::pressure_outlet(zero);
        let bc_wall_noslip = BoundaryCondition::wall_no_slip();
        let newtonian_blood = matches!(&self.blood, BloodModel::Newtonian(_));

        for iteration in 0..self.config.max_iterations {
            self.solve_u_momentum(&bc_inlet, &bc_outlet, u_inlet)?;
            self.solve_v_momentum(&bc_wall_noslip, &bc_wall_noslip)?;

            // Execute PISO pressure correction loops (n_correctors = 1 for SIMPLE, >1 for PISO)
            for _ in 0..self.config.n_correctors {
                self.solve_pressure_correction()?;
            }

            // Global mass-flux correction: scale outlet face velocities so
            // that Q_outlet = Q_inlet exactly.  Applied after the initial
            // development phase (iteration > 50) to avoid interfering with
            // the pressure field while it's still establishing the flow
            // pattern (Versteeg & Malalasekera 2007, §11.9).
            // if iteration > 50 {
            //     self.apply_mass_flux_correction();
            // }

            // Turbulence model update: solve k and omega transport equations
            // and compute nu_t.  Only runs when turbulence is enabled and
            // at the specified update interval.
            if let Some(ref mut turb) = self.turbulence {
                if iteration % turb.update_interval == 0 && iteration > 10 {
                    // Build velocity vector from staggered u,v fields.
                    let nx = self.grid.nx;
                    let ny = self.grid.ny;
                    let zero: T = scalar::zero();
                    let mut velocity = vec![Vector2::new(zero, zero); nx * ny];
                    let half: T = <T as FloatElement>::from_f64(0.5);
                    for i in 0..nx {
                        for j in 0..ny {
                            let u_cc = (self.field.u[(i, j)] + self.field.u[(i + 1, j)]) * half;
                            let v_cc = (self.field.v[(i, j)] + self.field.v[(i, j + 1)]) * half;
                            velocity[j * nx + i] = Vector2::new(u_cc, v_cc);
                        }
                    }
                    let mu_mol = self.field.mu[(0, 0)]; // reference molecular viscosity
                    let dt_pseudo: T = <T as FloatElement>::from_f64(1e-3);
                    let _ = turb.model.update(
                        &mut turb.k,
                        &mut turb.omega,
                        &velocity,
                        self.density,
                        mu_mol / self.density, // kinematic viscosity
                        dt_pseudo,
                        self.grid.dx,
                        self.grid.dy_at(0),
                    );
                    // Update nu_t field from k and omega.
                    use crate::physics::turbulence::TurbulenceModel;
                    for i in 0..nx {
                        for j in 0..ny {
                            let idx = j * nx + i;
                            let nu_t = turb.model.turbulent_viscosity(
                                turb.k[idx],
                                turb.omega[idx],
                                self.density,
                            );
                            let nu_t_val = nu_t / self.density;
                            self.field.nu_t[(i, j)] = if nu_t_val > zero { nu_t_val } else { zero };
                        }
                    }
                }
            }

            // Update viscosity with under-relaxation (alpha_mu) to prevent
            // oscillation in non-Newtonian SIMPLE iterations.
            if !newtonian_blood && iteration % self.config.viscosity_update_interval == 0 {
                self.field
                    .update_viscosity(&self.grid, &self.blood, self.config.alpha_mu);
            }

            let (res_cont, res_cont_l1, res_max) = self.compute_residuals_inner();
            last_residual = res_cont.max_scalar(res_cont_l1);

            if tracing::enabled!(tracing::Level::DEBUG) {
                let mut max_u: T = scalar::zero();
                let mut max_v: T = scalar::zero();
                let mut max_p: T = scalar::zero();
                for i in 0..=self.grid.nx {
                    for j in 0..self.grid.ny {
                        let val = <T as NumericElement>::abs(self.field.u[(i, j)]);
                        if val > max_u {
                            max_u = val;
                        }
                    }
                }
                for i in 0..self.grid.nx {
                    for j in 0..=self.grid.ny {
                        let val = <T as NumericElement>::abs(self.field.v[(i, j)]);
                        if val > max_v {
                            max_v = val;
                        }
                    }
                }
                for i in 0..self.grid.nx {
                    for j in 0..self.grid.ny {
                        let val = <T as NumericElement>::abs(self.field.p[(i, j)]);
                        if val > max_p {
                            max_p = val;
                        }
                    }
                }

                tracing::debug!(
                    iteration,
                    cont = <T as NumericElement>::to_f64(res_cont),
                    cont_l1 = <T as NumericElement>::to_f64(res_cont_l1),
                    max_pointwise = <T as NumericElement>::to_f64(res_max),
                    max_u = <T as NumericElement>::to_f64(max_u),
                    max_v = <T as NumericElement>::to_f64(max_v),
                    max_p = <T as NumericElement>::to_f64(max_p),
                    "SIMPLE residuals"
                );
            }

            if self.check_divergence() || !<T as NumericElement>::is_finite(last_residual) {
                return Err(Error::Solver(
                    "SIMPLE solver diverged: NaN or Inf detected in velocity field".to_string(),
                ));
            }
            // Residual growth guard: after the startup transient, a very large
            // pointwise residual usually indicates oscillation or divergence.
            // Skip this check during the initial flow establishment phase so
            // high-Re or highly branched cases can settle before being judged.
            if iteration > 50 {
                let growth_limit: T = <T as FloatElement>::from_f64(1e6);
                if res_max > growth_limit {
                    return Err(Error::Solver(format!(
                        "SIMPLE solver diverged: max pointwise residual {} exceeds limit",
                        <T as NumericElement>::to_f64(res_max)
                    )));
                }
            }

            if last_residual < self.config.tolerance {
                if newtonian_blood {
                    self.field.update_shear_rate(&self.grid);
                }
                return Ok(SolveResult {
                    iterations: iteration + 1,
                    residual: last_residual,
                    converged: true,
                });
            }
        }

        if newtonian_blood {
            self.field.update_shear_rate(&self.grid);
        }
        Ok(SolveResult {
            iterations: self.config.max_iterations,
            residual: last_residual,
            converged: false,
        })
    }
}
