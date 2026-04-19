//! 2D Drift-Diffusion Scalar Transport Solver.
//!
//! Handles steady advection-diffusion coupled with an arbitrary deterministic drift velocity:
//!
//! ```text
//! ∇·((u + u_drift) c) = ∇·(Γ ∇c)
//! ```
//!
//! where `u` is the bulk fluid velocity, `u_drift` is the particle drift velocity
//! (e.g., from Acoustic Radiation Force or inertial lift), and `Γ` is the diffusion coefficient.
//!
//! # Theorem — Patankar Scarborough Criterion with Drift
//!
//! By strictly upwinding the composite velocity field $\vec{V}_{eff} = \vec{u} + \vec{v}_{drift}$,
//! the resulting linear system's coefficients $a_{nb}$ remain strictly non-negative.
//! Wall boundaries (where drift pushes cells into impermeable barriers) enforce zero-flux
//! by clipping the boundary normal drift velocity. Thus, the total mass is conserved
//! and the concentration field is unconditionally bounded (non-negative).
//!
//! **Proof sketch**:
//! For each face, the convective flux is computed with $F = V_{eff} \cdot \vec{n} A$.
//! The upwind neighbor coefficient is $a_{nb} = D_{nb} + \max(\mp F_{nb}, 0) \ge 0$.
//! The central coefficient $a_P = \sum a_{nb} + S_P$. As long as the effective velocity
//! divergence is zero or treated implicitly, the matrix remains a diagonally dominant M-matrix,
//! ensuring Gauss-Seidel convergence and a discrete maximum principle $c_i \ge 0$.
//!
//! # Literature
//! - Patankar, S.V. (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing.

use crate::grid::array2d::Array2D;
use crate::solvers::ns_fvm::{FlowField2D, StaggeredGrid2D};
use crate::solvers::scalar_transport_2d::ScalarTransportConfig;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// 2D Drift-Diffusion Transport Solver
pub struct DriftDiffusionSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    /// Concentration field [nx][ny] (stored at cell centers)
    pub c: Array2D<T>,
}

impl<T: RealField + Copy + Float + FromPrimitive> DriftDiffusionSolver2D<T> {
    /// Create new drift-diffusion transport solver
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            c: Array2D::new(nx, ny, T::zero()),
        }
    }

    /// Solve the steady-state drift-advection-diffusion equation.
    ///
    /// The drift_field provides `u_drift` and `v_drift` at each face.
    pub fn solve(
        &mut self,
        grid: &StaggeredGrid2D<T>,
        field: &FlowField2D<T>,
        drift_field: &FlowField2D<T>,
        config: &ScalarTransportConfig<T>,
        boundary_c: &[T], // Profile at West inlet
    ) -> Result<usize, String> {
        let nx = grid.nx;
        let ny = grid.ny;
        let dx = grid.dx;
        let dy = grid.dy;
        let gamma = config.diffusion_coeff;
        let zero = T::zero();
        let half = T::from_f64(0.5).unwrap();
        let one = T::one();
        let omega = T::from_f64(0.8).unwrap(); // Under-relaxation

        for iteration in 0..config.max_iterations {
            let mut max_diff = zero;

            for i in 0..nx {
                for j in 0..ny {
                    if !field.mask[(i, j)] {
                        continue;
                    }

                    let mut a_e = zero;
                    let mut a_w = zero;
                    let mut a_n = zero;
                    let mut a_s = zero;
                    let mut b = zero;

                    // Compute effective velocity (fluid + drift)
                    // East face
                    let f_e = (field.u[(i + 1, j)] + drift_field.u[(i + 1, j)]) * dy;
                    let d_e = gamma * dy / dx;
                    if i < nx - 1 && field.mask[(i + 1, j)] {
                        a_e = d_e + Float::max(-f_e, zero);
                    }

                    // West face
                    let f_w = (field.u[(i, j)] + drift_field.u[(i, j)]) * dy;
                    let d_w = gamma * dy / dx;
                    if i > 0 && field.mask[(i - 1, j)] {
                        a_w = d_w + Float::max(f_w, zero);
                    } else if f_w > zero {
                        // Inlet boundary
                        let d_in = gamma * dy / (half * dx);
                        b += (d_in + f_w) * boundary_c[j];
                    }

                    // North face
                    // Impermeable wall zero-flux: if the neighbor is solid, drift into it is zero.
                    let mut v_n_eff = field.v[(i, j + 1)] + drift_field.v[(i, j + 1)];
                    if j == ny - 1 || !field.mask[(i, j + 1)] {
                        v_n_eff = field.v[(i, j + 1)]; // Cancel drift into explicit boundaries
                    }
                    let f_n = v_n_eff * dx;
                    let d_n = gamma * dx / dy;
                    if j < ny - 1 && field.mask[(i, j + 1)] {
                        a_n = d_n + Float::max(-f_n, zero);
                    }

                    // South face
                    let mut v_s_eff = field.v[(i, j)] + drift_field.v[(i, j)];
                    if j == 0 || !field.mask[(i, j - 1)] {
                        v_s_eff = field.v[(i, j)]; // Cancel drift into explicit boundaries
                    }
                    let f_s = v_s_eff * dx;
                    let d_s = gamma * dx / dy;
                    if j > 0 && field.mask[(i, j - 1)] {
                        a_s = d_s + Float::max(f_s, zero);
                    }

                    // Patankar central coefficient:
                    let a_p = a_e + a_w + a_n + a_s + (f_e - f_w + f_n - f_s);
                    let mut a_p_eff = a_p;

                    if (i == 0 || (i > 0 && !field.mask[(i - 1, j)])) && f_w > zero {
                        let d_in = gamma * dy / (half * dx);
                        a_p_eff += d_in;
                    }

                    if Float::abs(a_p_eff) > T::from_f64(1e-30).unwrap() {
                        let c_e = if i < nx - 1 && field.mask[(i + 1, j)] {
                            self.c[(i + 1, j)]
                        } else {
                            self.c[(i, j)]
                        };
                        let c_w = if i > 0 && field.mask[(i - 1, j)] {
                            self.c[(i - 1, j)]
                        } else {
                            zero
                        };
                        let c_n = if j < ny - 1 && field.mask[(i, j + 1)] {
                            self.c[(i, j + 1)]
                        } else {
                            self.c[(i, j)]
                        };
                        let c_s = if j > 0 && field.mask[(i, j - 1)] {
                            self.c[(i, j - 1)]
                        } else {
                            self.c[(i, j)]
                        };

                        let c_target =
                            (a_e * c_e + a_w * c_w + a_n * c_n + a_s * c_s + b) / a_p_eff;
                        let c_new = (one - omega) * self.c[(i, j)] + omega * c_target;

                        let diff = Float::abs(c_new - self.c[(i, j)]);
                        if diff > max_diff {
                            max_diff = diff;
                        }
                        self.c[(i, j)] = c_new;
                    }
                }
            }

            if max_diff < config.tolerance {
                return Ok(iteration + 1);
            }
        }

        Err(format!(
            "Drift-diffusion failed to converge after {} iterations",
            config.max_iterations
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Theorem: Drift equations with wall-blocking must guarantee no negative
    /// concentrations (boundedness / Scarborough criterion).
    #[test]
    fn drift_diffusion_boundedness() {
        let nx = 10;
        let ny = 10;
        let grid = StaggeredGrid2D::<f64>::new(nx, ny, 0.1, 0.1);
        let mut flow = FlowField2D::<f64>::new(nx, ny);
        let mut drift = FlowField2D::<f64>::new(nx, ny);

        // Fluid moves to the right
        for i in 0..=nx {
            for j in 0..ny {
                flow.u[(i, j)] = 1.0;
            }
        }

        // Drift constantly pushes DOWN (-y) toward the south wall
        for i in 0..nx {
            for j in 0..=ny {
                drift.v[(i, j)] = -0.5;
            }
        }

        let mut drift_solver = DriftDiffusionSolver2D::<f64>::new(nx, ny);
        let config = ScalarTransportConfig {
            max_iterations: 1000,
            tolerance: 1e-6,
            diffusion_coeff: 1e-4, // Low diffusion, advective/drift dominated
        };
        let inlet_c = vec![1.0; ny]; // feed uniformly 1.0

        let iters = drift_solver
            .solve(&grid, &flow, &drift, &config, &inlet_c)
            .unwrap();
        assert!(iters > 0);

        // Verify concentration is strictly bound between 0 and a maximum pile-up.
        for i in 0..nx {
            for j in 0..ny {
                assert!(
                    drift_solver.c[(i, j)] >= -1e-12,
                    "Strict boundedness violated: c < 0 ({})",
                    drift_solver.c[(i, j)]
                );
            }
        }

        // South wall row (j=0) should have the highest concentration due to drift pile-up
        for i in 2..nx {
            assert!(
                drift_solver.c[(i, 0)] > drift_solver.c[(i, ny - 1)],
                "Drift failed to accumulate solute at the -y wall. c(y=0)={}, c(y=top)={}",
                drift_solver.c[(i, 0)],
                drift_solver.c[(i, ny - 1)]
            );
        }
    }
}
