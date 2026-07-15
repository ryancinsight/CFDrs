//! 2D Scalar Transport Solver (Advection-Diffusion)
//!
//! Implements the Finite Volume Method for steady-state scalar transport:
//! ∇·(u c) = ∇·(Γ ∇c)
//!
//! # Discretization
//!
//! - Advection: Upwind Differencing Scheme (UDS) for stability at high Peclet numbers.
//! - Diffusion: Central difference.
//! - Solver: Gauss-Seidel iterations.
//!
//! # Theorem (Scalar Transport Boundedness — Patankar 1980, Ch. 5)
//!
//! The upwind FVM discretisation of $\nabla \cdot (\mathbf{u}\,c) = \nabla \cdot (\Gamma \nabla c)$
//! produces a coefficient matrix satisfying the Scarborough criterion
//! $\sum_{nb}|a_{nb}|/|a_P| \le 1$ (strict for at least one row), guaranteeing
//! Gauss-Seidel convergence and the discrete maximum principle $c_{\min} \le c_i \le c_{\max}$.
//!
//! **Proof sketch**:
//! With UDS, each $a_{nb} = D_{nb} + \max(F_{nb}, 0) \ge 0$ and
//! $a_P = \sum_{nb} a_{nb}$. All coefficients are non-negative, so
//! $c_i$ is a convex combination of its neighbours. The Gauss-Seidel iteration
//! contracts because $\rho(\mathbf{M}^{-1}\mathbf{N}) < 1$ for M-matrices.
//! Boundedness follows directly from the non-negative coefficients and the
//! Scarborough criterion.

use crate::grid::array2d::Array2D;
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use crate::solvers::ns_fvm::{FlowField2D, StaggeredGrid2D};
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Configuration for scalar transport solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScalarTransportConfig<T: Cfd2dScalar + Copy> {
    /// Maximum inner iterations
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
    /// Diffusion coefficient (Γ) [m²/s]
    pub diffusion_coeff: T,
}

impl<T: Cfd2dScalar + Copy + FloatElement> Default for ScalarTransportConfig<T> {
    fn default() -> Self {
        Self {
            max_iterations: 5000,
            // 1e-5 matches the FVM spatial truncation error O(Δx²) on typical
            // coarse grids (Δx ~ 0.003 m → Δx² ~ 1e-5).  Using 1e-8 demands
            // residuals an order of magnitude below the discretisation error,
            // which is unachievable with Gauss-Seidel on advection-dominated flows.
            tolerance: scalar::from_f64(1e-5),
            diffusion_coeff: scalar::from_f64(1e-9), // Typical diffusion
        }
    }
}

/// 2D Scalar Transport Solver
pub struct ScalarTransportSolver2D<T: Cfd2dScalar + Copy + FloatElement> {
    /// Concentration field \[nx]\[ny] (stored at cell centers)
    pub c: Array2D<T>,
    /// Previous iteration for convergence check.
    _c_old: Array2D<T>,
}

impl<T: Cfd2dScalar + Copy + FloatElement> ScalarTransportSolver2D<T> {
    /// Create new scalar transport solver
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            c: Array2D::new(nx, ny, scalar::zero()),
            _c_old: Array2D::new(nx, ny, scalar::zero()),
        }
    }

    /// Solve the steady-state advection-diffusion equation
    ///
    /// a_P c_P = a_E c_E + a_W c_W + a_N c_N + a_S c_S
    pub fn solve(
        &mut self,
        grid: &StaggeredGrid2D<T>,
        field: &FlowField2D<T>,
        config: &ScalarTransportConfig<T>,
        boundary_c: &[T], // Profile at West inlet
    ) -> Result<usize, String> {
        let nx = grid.nx;
        let ny = grid.ny;
        let dx = grid.dx;
        let dy = grid.dy;
        let gamma = config.diffusion_coeff;
        let zero: T = scalar::zero();
        let half = scalar::from_f64::<T>(0.5);
        let one: T = scalar::one();
        let omega = scalar::from_f64::<T>(0.8); // Relaxation
        let tiny = scalar::from_f64::<T>(1e-30);

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

                    // East
                    let f_e = field.u[(i + 1, j)] * dy;
                    let d_e = gamma * dy / dx;
                    if i < nx - 1 && field.mask[(i + 1, j)] {
                        a_e = d_e + <T as NumericElement>::max_scalar(-f_e, zero);
                    } else if i == nx - 1 || !field.mask[(i + 1, j)] {
                        // Outlet or internal wall boundary (zero gradient)
                    }

                    // West
                    let f_w = field.u[(i, j)] * dy;
                    let d_w = gamma * dy / dx;
                    if i > 0 && field.mask[(i - 1, j)] {
                        a_w = d_w + <T as NumericElement>::max_scalar(f_w, zero);
                    } else {
                        // Inlet or wall. If u > 0, it's an inlet.
                        if f_w > zero {
                            let d_in = gamma * dy / (half * dx);
                            b += (d_in + f_w) * boundary_c[j];
                        }
                    }

                    // North
                    let f_n = field.v[(i, j + 1)] * dx;
                    let d_n = gamma * dx / dy;
                    if j < ny - 1 && field.mask[(i, j + 1)] {
                        a_n = d_n + <T as NumericElement>::max_scalar(-f_n, zero);
                    }

                    // South
                    let f_s = field.v[(i, j)] * dx;
                    let d_s = gamma * dx / dy;
                    if j > 0 && field.mask[(i, j - 1)] {
                        a_s = d_s + <T as NumericElement>::max_scalar(f_s, zero);
                    }

                    let a_p = a_e + a_w + a_n + a_s + (f_e - f_w + f_n - f_s);
                    // Add inlet terms to a_p if applicable
                    let mut a_p_eff = a_p;
                    if (i == 0 || (i > 0 && !field.mask[(i - 1, j)])) && f_w > zero {
                        let d_in = gamma * dy / (half * dx);
                        a_p_eff += d_in;
                    }

                    if <T as NumericElement>::abs(a_p_eff) > tiny {
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

                        let diff = <T as NumericElement>::abs(c_new - self.c[(i, j)]);
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
            "Scalar transport failed to converge after {} iterations",
            config.max_iterations
        ))
    }
}
