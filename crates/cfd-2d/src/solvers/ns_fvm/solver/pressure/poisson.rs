//! Pressure correction Poisson solver.

use crate::scalar;
use crate::scalar::Cfd2dScalar;
use crate::solvers::ns_fvm::BloodModel;
use crate::solvers::ns_fvm::solver::NavierStokesSolver2D;
use cfd_core::error::Error;
use eunomia::FloatElement;

// A conservative fixed SOR factor for the masked pressure-correction grid.
// The classical Poisson optimum approaches 2 on large regular grids; 1.7
// retains stability when the Venturi mask removes neighboring coefficients.
const PRESSURE_SOR_RELAXATION: f64 = 1.7;
const PRESSURE_SOR_MAX_SWEEPS: usize = 32;
const PRESSURE_SOR_REFINED_MAX_SWEEPS: usize = 200;
// Channel profiles use alpha_p <= 0.12; the general SIMPLE profiles use
// alpha_p >= 0.2. The midpoint preserves the deeper correction budget for
// higher pressure-relaxation profiles without penalizing the channel path.
const PRESSURE_SOR_REFINEMENT_ALPHA_P: f64 = 0.15;

impl<T: Cfd2dScalar + eunomia::RealField + Copy + FloatElement> NavierStokesSolver2D<T> {
    /// Solves the pressure-correction Poisson equation.
    ///
    /// Formulates the continuity residual from the intermediate velocity fields
    /// and solves for the pressure correction `p'` iteratively via SOR.
    pub fn solve_pressure_correction(&mut self) -> Result<(), Error> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let rho = self.density;
        let zero: T = scalar::zero();
        let tiny: T = scalar::from_f64(1e-30);

        let a_p_u = &self.a_p_u;
        let a_p_v = &self.a_p_v;
        let d_u = &mut self.pressure_poisson_d_u;
        let d_v = &mut self.pressure_poisson_d_v;
        let p_prime = &mut self.pressure_poisson_p_prime;
        let b = &mut self.pressure_poisson_rhs;
        let a_e_workspace = &mut self.pressure_poisson_a_e;
        let a_w_workspace = &mut self.pressure_poisson_a_w;
        let a_n_workspace = &mut self.pressure_poisson_a_n;
        let a_s_workspace = &mut self.pressure_poisson_a_s;
        let a_p_workspace = &mut self.pressure_poisson_a_p;

        d_u.fill(zero);
        d_v.fill(zero);
        p_prime.fill(zero);
        b.fill(zero);
        a_e_workspace.fill(zero);
        a_w_workspace.fill(zero);
        a_n_workspace.fill(zero);
        a_s_workspace.fill(zero);
        a_p_workspace.fill(zero);

        {
            let field = &self.field;

            for i in 1..nx {
                for j in 0..ny {
                    let a = a_p_u[(i, j)];
                    if a > tiny {
                        d_u[(i, j)] = self.grid.dy_at(j) / a;
                    }
                }
            }
            for j in 0..ny {
                d_u[(nx, j)] = d_u[(nx - 1, j)];
            }
            for i in 0..nx {
                for j in 1..ny {
                    let a = a_p_v[(i, j)];
                    if a > tiny {
                        d_v[(i, j)] = dx / a;
                    }
                }
            }
            if ny > 0 {
                for i in 0..nx {
                    d_v[(i, ny)] = d_v[(i, ny - 1)];
                }
            }

            for i in 0..nx {
                for j in 0..ny {
                    if !field.mask[(i, j)] {
                        continue;
                    }
                    let dy_j = self.grid.dy_at(j);
                    b[(i, j)] = rho
                        * ((field.u[(i, j)] - field.u[(i + 1, j)]) * dy_j
                            + (field.v[(i, j)] - field.v[(i, j + 1)]) * dx);
                }
            }

            for i in 0..nx {
                for j in 0..ny {
                    if !field.mask[(i, j)] {
                        continue;
                    }
                    let dy_j = self.grid.dy_at(j);
                    let a_e = if i + 1 < nx {
                        if field.mask[(i + 1, j)] {
                            rho * d_u[(i + 1, j)] * dy_j
                        } else {
                            zero
                        }
                    } else {
                        rho * d_u[(i + 1, j)] * dy_j
                    };
                    let a_w = if i > 0 {
                        if field.mask[(i - 1, j)] {
                            rho * d_u[(i, j)] * dy_j
                        } else {
                            zero
                        }
                    } else {
                        zero
                    };
                    let a_n = if j + 1 < ny {
                        if field.mask[(i, j + 1)] {
                            rho * d_v[(i, j + 1)] * dx
                        } else {
                            zero
                        }
                    } else {
                        zero
                    };
                    let a_s = if j > 0 {
                        if field.mask[(i, j - 1)] {
                            rho * d_v[(i, j)] * dx
                        } else {
                            zero
                        }
                    } else {
                        zero
                    };
                    a_e_workspace[(i, j)] = a_e;
                    a_w_workspace[(i, j)] = a_w;
                    a_n_workspace[(i, j)] = a_n;
                    a_s_workspace[(i, j)] = a_s;
                    a_p_workspace[(i, j)] = a_e + a_w + a_n + a_s;
                }
            }

            let relaxation = match &self.blood {
                BloodModel::Newtonian(_) => scalar::one(),
                BloodModel::Casson(_) | BloodModel::CarreauYasuda(_) => {
                    scalar::from_f64::<T>(PRESSURE_SOR_RELAXATION)
                }
            };
            // SIMPLE repeats this correction and uses the field continuity
            // residual below as the convergence oracle. The cap limits inner
            // work without treating a partially converged pressure correction
            // as a converged outer solution.
            let max_sweeps =
                if self.config.alpha_p > scalar::from_f64::<T>(PRESSURE_SOR_REFINEMENT_ALPHA_P) {
                    PRESSURE_SOR_REFINED_MAX_SWEEPS
                } else {
                    PRESSURE_SOR_MAX_SWEEPS
                };
            for _ in 0..max_sweeps {
                for i in 0..nx {
                    for j in 0..ny {
                        let a_e = a_e_workspace[(i, j)];
                        let a_w = a_w_workspace[(i, j)];
                        let a_n = a_n_workspace[(i, j)];
                        let a_s = a_s_workspace[(i, j)];
                        let a_p = a_p_workspace[(i, j)];
                        if a_p < tiny {
                            continue;
                        }
                        let pe = if i + 1 < nx {
                            p_prime[(i + 1, j)]
                        } else {
                            zero
                        };
                        let pw = if i > 0 {
                            p_prime[(i - 1, j)]
                        } else {
                            p_prime[(i, j)]
                        };
                        let pn = if j + 1 < ny {
                            p_prime[(i, j + 1)]
                        } else {
                            p_prime[(i, j)]
                        };
                        let ps = if j > 0 {
                            p_prime[(i, j - 1)]
                        } else {
                            p_prime[(i, j)]
                        };
                        let gauss_seidel_value =
                            (a_e * pe + a_w * pw + a_n * pn + a_s * ps + b[(i, j)]) / a_p;
                        let previous = p_prime[(i, j)];
                        p_prime[(i, j)] = previous + relaxation * (gauss_seidel_value - previous);
                    }
                }
            }
        }

        for i in 1..=nx {
            for j in 0..ny {
                if i < nx {
                    if !self.field.mask[(i, j)] && !self.field.mask[(i - 1, j)] {
                        continue;
                    }
                } else if !self.field.mask[(nx - 1, j)] {
                    continue;
                }
                let dp = if i > 0 && i < nx {
                    p_prime[(i - 1, j)] - p_prime[(i, j)]
                } else if i == nx {
                    p_prime[(i - 1, j)]
                } else {
                    zero
                };
                self.field.u[(i, j)] += d_u[(i, j)] * dp;
            }
        }

        for i in 0..nx {
            for j in 1..=ny {
                if j < ny {
                    if !self.field.mask[(i, j)] && !self.field.mask[(i, j - 1)] {
                        continue;
                    }
                } else if !self.field.mask[(i, ny - 1)] {
                    continue;
                }
                let dp = if j > 0 && j < ny {
                    p_prime[(i, j - 1)] - p_prime[(i, j)]
                } else {
                    zero
                };
                self.field.v[(i, j)] += d_v[(i, j)] * dp;
            }
        }

        let alpha_p = self.config.alpha_p;
        let alpha_u = self.config.alpha_u;
        for i in 0..nx {
            for j in 0..ny {
                self.field.p[(i, j)] += (alpha_p / alpha_u) * p_prime[(i, j)];
            }
        }

        Ok(())
    }
}
