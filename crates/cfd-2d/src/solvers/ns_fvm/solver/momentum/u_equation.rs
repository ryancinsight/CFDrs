//! u-momentum Gauss-Seidel solver.

use crate::error::Error;
use crate::solvers::ns_fvm::boundary::{BCType, BoundaryCondition};
use crate::solvers::ns_fvm::solver::NavierStokesSolver2D;
use crate::grid::array2d::Array2D;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

impl<T: RealField + Copy + Float + FromPrimitive> NavierStokesSolver2D<T> {
    /// Solves the u-momentum equation on the staggered grid.
    ///
    /// Calculates upwind convective fluxes and central differenced diffusive
    /// fluxes. Formulates the tridiagonal system and updates the intermediate `u*` field.
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
        let mut a_p_u = Array2D::new(nx + 1, ny, T::one());

        for i in 1..nx {
            for j in 0..ny {
                let dy_j = self.grid.dy_at(j);

                // Effective viscosity = molecular + turbulent (rho * nu_t).
                // For laminar flows nu_t = 0 and this reduces to mu alone.
                let mu_e = if i < nx {
                    self.field.mu[(i, j)] + rho * self.field.nu_t[(i, j)]
                } else {
                    self.field.mu[(nx - 1, j)] + rho * self.field.nu_t[(nx - 1, j)]
                };
                let mu_w = if i > 0 {
                    self.field.mu[(i - 1, j)] + rho * self.field.nu_t[(i - 1, j)]
                } else {
                    self.field.mu[(0, j)] + rho * self.field.nu_t[(0, j)]
                };
                let mu_face = (mu_e + mu_w) / (one + one);

                let d_e = mu_face * dy_j / dx;
                let d_w = d_e;
                let d_n = if j + 1 < ny {
                    mu_face * dx / self.grid.dy_face(j)
                } else {
                    mu_face * dx / (dy_j * half)
                };
                let d_s = if j > 0 {
                    mu_face * dx / self.grid.dy_face(j - 1)
                } else {
                    mu_face * dx / (dy_j * half)
                };

                let f_e = if i < nx {
                    rho * (u_old[(i, j)] + u_old[(i + 1, j)]) * half * dy_j
                } else {
                    rho * u_old[(i, j)] * dy_j
                };
                let f_w = if i > 1 {
                    rho * (u_old[(i - 1, j)] + u_old[(i, j)]) * half * dy_j
                } else {
                    rho * u_old[(i, j)] * dy_j
                };

                let v_n = if j < ny - 1 {
                    (v_old[(i - 1, j + 1)] + v_old[(i, j + 1)]) * half
                } else {
                    T::zero()
                };
                let v_s = if j > 0 {
                    (v_old[(i - 1, j)] + v_old[(i, j)]) * half
                } else {
                    T::zero()
                };
                let f_n = rho * v_n * dx;
                let f_s = rho * v_s * dx;

                // Hybrid convection scheme: central differencing at Pe < 2,
                // upwind at Pe >= 2.  This gives second-order accuracy in
                // low-Pe regions (channel cores) while maintaining stability
                // near walls and in high-velocity throat regions.
                let two = one + one;
                let pe_e = Float::abs(f_e) / d_e;
                let pe_w = Float::abs(f_w) / d_w;
                let a_e = if pe_e <= two {
                    d_e - f_e * half // central
                } else {
                    d_e + Float::max(zero, -f_e) // upwind
                };
                let a_w = if pe_w <= two {
                    d_w + f_w * half
                } else {
                    d_w + Float::max(zero, f_w)
                };
                // North/south faces remain upwind (cross-stream direction
                // typically has low velocity so upwind is stable and adequate).
                let a_n = Float::max(d_n - f_n * half, zero) + Float::max(-f_n, zero);
                let a_s = Float::max(d_s + f_s * half, zero) + Float::max(f_s, zero);

                let mut a_p = a_e + a_w + a_n + a_s + (f_e - f_w) + (f_n - f_s);
                if a_p < T::from_f64(1e-30).unwrap_or(zero) {
                    a_p = T::from_f64(1e-10).unwrap_or(one);
                }

                let p_left = if i > 0 && i - 1 < nx {
                    self.field.p[(i - 1, j)]
                } else {
                    self.field.p[(0, j)]
                };
                let p_right = if i < nx {
                    self.field.p[(i, j)]
                } else {
                    self.field.p[(nx - 1, j)]
                };
                let pressure_source = (p_left - p_right) * dy_j;

                let u_e = if i < nx {
                    self.field.u[(i + 1, j)]
                } else {
                    self.field.u[(i, j)]
                };
                let u_w = if i >= 1 {
                    self.field.u[(i - 1, j)]
                } else {
                    self.field.u[(i, j)]
                };
                let u_n = if j + 1 < ny {
                    self.field.u[(i, j + 1)]
                } else {
                    zero
                };
                let u_s = if j >= 1 {
                    self.field.u[(i, j - 1)]
                } else {
                    zero
                };

                let is_fluid = if i > 0 && i < nx {
                    self.field.mask[(i - 1, j)] && self.field.mask[(i, j)]
                } else {
                    true
                };

                if is_fluid {
                    let u_star =
                        (a_e * u_e + a_w * u_w + a_n * u_n + a_s * u_s + pressure_source) / a_p;
                    self.field.u[(i, j)] = self.field.u[(i, j)] * (one - alpha) + u_star * alpha;
                    a_p_u[(i, j)] = a_p;
                } else {
                    self.field.u[(i, j)] = zero;
                    a_p_u[(i, j)] = one;
                }
            }
        }

        if bc_west.is_dirichlet() {
            match bc_west {
                BoundaryCondition::Wall { .. } => {
                    for j in 0..ny {
                        self.field.u[(0, j)] = zero;
                    }
                }
                BoundaryCondition::VelocityInlet { velocity } => {
                    for j in 0..ny {
                        self.field.u[(0, j)] = velocity[0];
                    }
                }
                BoundaryCondition::PressureInlet {
                    velocity_direction: Some(velocity_direction),
                    ..
                }
                | BoundaryCondition::CharacteristicInlet {
                    velocity: Some(velocity_direction),
                    ..
                } => {
                    for j in 0..ny {
                        self.field.u[(0, j)] = velocity_direction[0];
                    }
                }
                BoundaryCondition::Dirichlet {
                    value,
                    component_values,
                } => {
                    let inlet_value = component_values
                        .as_ref()
                        .and_then(|components| components.first().copied().flatten())
                        .unwrap_or(*value);
                    for j in 0..ny {
                        self.field.u[(0, j)] = inlet_value;
                    }
                }
                _ => {}
            }
        }

        match bc_east {
            BoundaryCondition::PressureOutlet { .. } => {
                for j in 0..ny {
                    self.field.u[(nx, j)] = self.field.u[(nx - 1, j)];
                }
            }
            _ => match bc_east.fundamental_type() {
                BCType::Neumann => {
                    for j in 0..ny {
                        self.field.u[(nx, j)] = self.field.u[(nx - 1, j)];
                    }
                }
                BCType::Dirichlet => {
                    if let BoundaryCondition::Wall { .. } = bc_east {
                        for j in 0..ny {
                            self.field.u[(nx, j)] = zero;
                        }
                    }
                }
                _ => {}
            },
        }

        self.a_p_u = a_p_u;
        Ok(())
    }
}
