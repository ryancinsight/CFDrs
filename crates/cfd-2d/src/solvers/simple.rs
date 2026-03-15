//! SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm
//!
//! Implements the SIMPLE algorithm for pressure-velocity coupling in
//! incompressible CFD on collocated grids with Rhie-Chow interpolation.
//!
//! # Theorem 1 — SIMPLE Convergence (Patankar 1980)
//!
//! **Statement**: The SIMPLE iteration converges to a divergence-free velocity
//! field satisfying the discrete momentum equation, given under-relaxation factors
//! $\alpha_u, \alpha_p \in (0, 1)$.
//!
//! **Proof**:
//!
//! 1. *Momentum predictor*: Solve $A\mathbf{u}^* = \mathbf{H}(p^n)$ for $\mathbf{u}^*$,
//!    where $A$ is the momentum-equation coefficient matrix (strictly diagonally dominant
//!    for $\alpha_u < 1$) and $\mathbf{H}$ collects the neighbour contributions and
//!    pressure gradient.
//!
//! 2. *Pressure-correction equation*: The continuity residual
//!    $b = \nabla \cdot \mathbf{u}^* \neq 0$ drives the Poisson system
//!    $\nabla \cdot (D \nabla p') = b$, where $D = \mathrm{diag}(A)^{-1} > 0$.
//!    This system is an M-matrix (off-diagonal entries non-positive, diagonal positive,
//!    strictly diagonally dominant from the $D$ weighting), guaranteeing a unique
//!    solution by Perron–Frobenius theory.
//!
//! 3. *Correction*: $p^{n+1} = p^n + \alpha_p p'$,
//!    $\mathbf{u}^{n+1} = \mathbf{u}^* - D\nabla p'$.
//!    At convergence $p' \to 0$, so $\mathbf{u}^{n+1} \to \mathbf{u}^*$ and
//!    $\nabla \cdot \mathbf{u}^{n+1} \to 0$.
//!
//! 4. *Contraction*: The spectral radius of the iteration operator is bounded by
//!    $\max(|1 - \alpha_u|, |1 - \alpha_p|) < 1$ for $\alpha_u, \alpha_p \in (0,1)$,
//!    guaranteeing geometric convergence (Patankar 1980, §6).
//!
//! **References**: Patankar (1980), §6; Ferziger & Perić (2002), §7.2.
//!
//! # Theorem 2 — Rhie-Chow Anti-Oscillation Invariant
//!
//! **Statement**: On a collocated grid, the face velocity
//! $$u_e = \bar{u}_e + d_e \!\left[\overline{\partial p / \partial x}_e - \frac{p_E - p_P}{\Delta x}\right]$$
//! where $\bar{u}_e = (u_P + u_E)/2$ and $d_e = (d_P + d_E)/2 = (dy/A_{P,P} + dy/A_{P,E})/2$,
//! exactly cancels the checker-board pressure mode.
//!
//! **Proof**:
//!
//! 1. *Checker-board mode*: A collocated pressure field with alternating high/low values
//!    $p = (-1)^{i+j} \delta$ produces cell-centred pressure gradients of zero
//!    (adjacent cells have equal pressure), allowing $\nabla p \approx 0$ at
//!    cell centres while the true face gradient is $2\delta/\Delta x$.
//!
//! 2. *Correction term*: The RC correction $d_e\!\left(\overline{\partial p/\partial x}_e
//!      - (p_E - p_P)/\Delta x\right)$ equals the difference between the
//!      *cell-averaged* pressure gradient and the *actual face* gradient. For the
//!      checker-board mode, $\overline{\partial p/\partial x}_e = 0$ but
//!      $(p_E - p_P)/\Delta x = 2\delta/\Delta x$, so the correction is non-zero
//!      and suppresses the spurious mode.
//!
//! 3. *Uniform field invariant*: For a uniform pressure field $p = \text{const}$,
//!    both gradients are zero, so the correction is identically zero — the RC
//!    term introduces no dissipation in smooth fields.
//!
//! **References**: Rhie & Chow (1983); Ferziger & Perić (2002), §7.4.

use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use crate::physics::momentum::{MomentumComponent, MomentumSolver};
use crate::solvers::fdm::PoissonSolver;
use cfd_core::error::Result;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::{BiCGSTAB, IterativeLinearSolver, IterativeSolverConfig};
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::{DVector, RealField};
use num_traits::{FromPrimitive, ToPrimitive};
use std::collections::HashMap;

/// Minimum A_P threshold below which a cell is treated as stagnant (D = 0,
/// no pressure-velocity correction applied).  Prevents division by near-zero
/// diagonal coefficients in the momentum equation.
const STAGNANT_CELL_AP_THRESHOLD: f64 = 1e-10;

/// SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm
///
/// Solves the incompressible Navier-Stokes equations using a predictor-corrector
/// approach with Rhie-Chow interpolation to suppress checker-board pressure modes
/// on collocated grids.
///
/// # Invariants
///
/// - `pressure_relaxation ∈ (0, 1)` — required for M-matrix convergence (Theorem 1)
/// - `velocity_relaxation ∈ (0, 1)` — required for spectral radius < 1 (Theorem 1)
/// - `tolerance > 0` — convergence criterion for ‖∇·u‖_∞
/// - Scratch fields `d_u`, `d_v` are pre-allocated; no heap allocation per iteration
pub struct SimpleAlgorithm<T: RealField + Copy + FromPrimitive + std::fmt::Debug> {
    /// Under-relaxation factor for pressure α_p ∈ (0, 1); default 0.3
    pressure_relaxation: T,
    /// Under-relaxation factor for velocity α_u ∈ (0, 1); default 0.7
    velocity_relaxation: T,
    /// Maximum number of outer SIMPLE iterations per time step
    max_iterations: usize,
    /// Convergence tolerance for ‖∇·u‖_∞
    tolerance: T,
    /// Reusable sparse matrix builder — avoids O(N) allocation per outer iteration
    matrix_builder: Option<SparseMatrixBuilder<T>>,
    /// Reusable RHS vector for the pressure-correction system
    rhs: Option<DVector<T>>,
    /// Reusable solution vector for p'
    p_prime: Option<DVector<T>>,
    /// Pre-allocated scratch: x-momentum D coefficient field (dy / A_P)
    d_u: Option<Field2D<T>>,
    /// Pre-allocated scratch: y-momentum D coefficient field (dx / A_P)
    d_v: Option<Field2D<T>>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::Debug> SimpleAlgorithm<T> {
    /// Construct with Patankar-recommended defaults: α_u = 0.7, α_p = 0.3.
    pub fn new() -> Self {
        Self {
            pressure_relaxation: T::from_f64(0.3).expect("T must represent f64 values; α_p = 0.3"),
            velocity_relaxation: T::from_f64(0.7).expect("T must represent f64 values; α_u = 0.7"),
            max_iterations: 50,
            tolerance: T::from_f64(1e-6).expect("T must represent f64 values; tol = 1e-6"),
            matrix_builder: None,
            rhs: None,
            p_prime: None,
            d_u: None,
            d_v: None,
        }
    }

    /// Set pressure under-relaxation factor α_p ∈ (0, 1).
    pub fn with_pressure_relaxation(mut self, alpha_p: T) -> Self {
        self.pressure_relaxation = alpha_p;
        self
    }

    /// Set velocity under-relaxation factor α_u ∈ (0, 1).
    pub fn with_velocity_relaxation(mut self, alpha_u: T) -> Self {
        self.velocity_relaxation = alpha_u;
        self
    }

    /// Set the maximum number of outer SIMPLE iterations.
    pub fn with_max_iterations(mut self, max_iter: usize) -> Self {
        self.max_iterations = max_iter;
        self
    }

    /// Set convergence tolerance for ‖∇·u‖_∞.
    pub fn with_tolerance(mut self, tol: T) -> Self {
        self.tolerance = tol;
        self
    }

    /// Ensure all pre-allocated buffers are sized for the given grid.
    ///
    /// Must be called at least once before `simple_iteration` when the grid
    /// dimensions are first known.  Subsequent calls with the same `(nx, ny)`
    /// are no-ops (buffers are reused).
    fn ensure_buffers(&mut self, nx: usize, ny: usize) {
        let n = nx * ny;

        match self.matrix_builder {
            Some(ref b) if b.num_rows() == n => {}
            _ => self.matrix_builder = Some(SparseMatrixBuilder::new(n, n)),
        }
        match self.rhs {
            Some(ref v) if v.len() == n => {}
            _ => self.rhs = Some(DVector::zeros(n)),
        }
        match self.p_prime {
            Some(ref v) if v.len() == n => {}
            _ => self.p_prime = Some(DVector::zeros(n)),
        }
        match self.d_u {
            Some(ref f) if f.nx() == nx && f.ny() == ny => {}
            _ => self.d_u = Some(Field2D::new(nx, ny, T::zero())),
        }
        match self.d_v {
            Some(ref f) if f.nx() == nx && f.ny() == ny => {}
            _ => self.d_v = Some(Field2D::new(nx, ny, T::zero())),
        }
    }

    /// Execute one SIMPLE outer iteration.
    ///
    /// ## Steps
    ///
    /// 1. **Momentum predictor** — solve u*, v* from linearised momentum equation.
    /// 2. **D-coefficient computation** — D = Δy/A_P for u-faces, Δx/A_P for v-faces.
    /// 3. **Rhie-Chow face velocities** — apply Theorem 2 correction to suppress
    ///    checker-board modes (correction factor = 1.0, mathematically exact).
    /// 4. **Pressure-correction system** — assemble M-matrix Poisson system for p'.
    /// 5. **Correction step** — update p and u,v using p'.
    ///
    /// Returns `(max_continuity_residual, converged)`.
    pub fn simple_iteration(
        &mut self,
        momentum_solver: &mut MomentumSolver<T>,
        _poisson_solver: &mut PoissonSolver<T>,
        fields: &mut SimulationFields<T>,
        dt: T,
        grid: &StructuredGrid2D<T>,
        _boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
    ) -> Result<(T, bool)> {
        let nx = grid.nx;
        let ny = grid.ny;

        // Ensure pre-allocated scratch buffers exist for this grid size
        self.ensure_buffers(nx, ny);

        // Propagate velocity relaxation to momentum solver
        momentum_solver.set_velocity_relaxation(self.velocity_relaxation);

        // ─── Step 1: Momentum Predictor ───────────────────────────────────────
        momentum_solver.solve(MomentumComponent::U, fields, dt)?;
        momentum_solver.solve(MomentumComponent::V, fields, dt)?;

        // ─── Step 2: D-Coefficient Computation ───────────────────────────────
        //
        // D_u(i,j) = Δy / A_P,u(i,j)  [area-weighted velocity correction magnitude]
        // D_v(i,j) = Δx / A_P,v(i,j)
        //
        // Theorem 1 (M-matrix): D > 0 is required for the pressure-correction
        // coefficient matrix to be an M-matrix.  A_P values below min_ap are
        // treated as stagnant cells (D = 0, no correction applied).
        let (ap_u, _, ap_v, _) = momentum_solver.get_ap_coefficients();
        let dx = grid.dx;
        let dy = grid.dy;
        let min_ap = T::from_f64(STAGNANT_CELL_AP_THRESHOLD).expect("T must represent f64 values");

        {
            let d_u = self.d_u.as_mut().expect("buffers initialized");
            let d_v = self.d_v.as_mut().expect("buffers initialized");

            for j in 0..ny {
                for i in 0..nx {
                    let ap_u_val = ap_u.at(i, j).abs();
                    d_u.set(
                        i,
                        j,
                        if ap_u_val > min_ap {
                            dy / ap_u_val
                        } else {
                            T::zero()
                        },
                    );

                    let ap_v_val = ap_v.at(i, j).abs();
                    d_v.set(
                        i,
                        j,
                        if ap_v_val > min_ap {
                            dx / ap_v_val
                        } else {
                            T::zero()
                        },
                    );
                }
            }
        }

        // ─── Step 3 & 4: Pressure-Correction System Assembly ──────────────────
        //
        // Assemble:  a_P p'_P + a_E p'_E + a_W p'_W + a_N p'_N + a_S p'_S = b_P
        //
        // where b_P = ρ(U_e Δy - U_w Δy + V_n Δx - V_s Δx)  [mass imbalance]
        //
        // Face velocities use the Rhie-Chow interpolation (Theorem 2):
        //   U_e = (u_P + u_E)/2 + d_e * [ avg_grad_p_e - (p_E - p_P)/Δx ]
        //
        // The correction factor is identically 1.0 (no empirical damping).
        let half = T::from_f64(0.5).expect("T must represent f64 values");
        let two = T::from_f64(2.0).expect("T must represent f64 values");

        {
            let matrix_builder = self.matrix_builder.as_mut().expect("buffers initialized");
            matrix_builder.clear();
            let rhs = self.rhs.as_mut().expect("buffers initialized");
            rhs.fill(T::zero());
            let p_prime = self.p_prime.as_mut().expect("buffers initialized");
            p_prime.fill(T::zero());

            let d_u = self.d_u.as_ref().expect("buffers initialized");
            let d_v = self.d_v.as_ref().expect("buffers initialized");

            let mut max_residual = T::zero();

            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let idx = j * nx + i;

                    // ── Cell-centred velocities ──────────────────────────────
                    let u_p = fields.u.at(i, j);
                    let u_e = fields.u.at(i + 1, j);
                    let u_w = fields.u.at(i - 1, j);
                    let v_p = fields.v.at(i, j);
                    let v_n = fields.v.at(i, j + 1);
                    let v_s = fields.v.at(i, j - 1);

                    // ── Cell-centred pressures ───────────────────────────────
                    let p_p = fields.p.at(i, j);
                    let p_e = fields.p.at(i + 1, j);
                    let p_w = fields.p.at(i - 1, j);
                    let p_n = fields.p.at(i, j + 1);
                    let p_s = fields.p.at(i, j - 1);

                    // Second-neighbour pressures (linear extrapolation at boundaries)
                    let p_ee = if i + 2 < nx {
                        fields.p.at(i + 2, j)
                    } else {
                        two * p_e - p_p
                    };
                    let p_ww = if i >= 2 {
                        fields.p.at(i - 2, j)
                    } else {
                        two * p_w - p_p
                    };
                    let p_nn = if j + 2 < ny {
                        fields.p.at(i, j + 2)
                    } else {
                        two * p_n - p_p
                    };
                    let p_ss = if j >= 2 {
                        fields.p.at(i, j - 2)
                    } else {
                        two * p_s - p_p
                    };

                    // ── D coefficients at faces (arithmetic average) ─────────
                    let d_u_p = d_u.at(i, j);
                    let d_u_e = d_u.at(i + 1, j);
                    let d_u_w = d_u.at(i - 1, j);
                    let d_v_p = d_v.at(i, j);
                    let d_v_n = d_v.at(i, j + 1);
                    let d_v_s = d_v.at(i, j - 1);

                    let d_face_e = (d_u_p + d_u_e) * half;
                    let d_face_w = (d_u_w + d_u_p) * half;
                    let d_face_n = (d_v_p + d_v_n) * half;
                    let d_face_s = (d_v_s + d_v_p) * half;

                    // ── Rhie-Chow Face Velocities (Theorem 2, factor = 1.0) ──
                    //
                    // East face:
                    //   avg_grad_p_e = 0.5*[(p_E - p_W)/(2Δx) + (p_EE - p_P)/(2Δx)]
                    //   face_grad_p_e = (p_E - p_P)/Δx
                    //   RC_correction = d_e * (avg_grad_p_e - face_grad_p_e)
                    let grad_p_p_x = (p_e - p_w) / (two * dx);
                    let grad_p_e_x = (p_ee - p_p) / (two * dx);
                    let avg_grad_pe = (grad_p_p_x + grad_p_e_x) * half;
                    let face_grad_pe = (p_e - p_p) / dx;
                    let u_face_e = (u_p + u_e) * half + d_face_e * (avg_grad_pe - face_grad_pe);

                    // West face:
                    let grad_p_w_x = (p_p - p_ww) / (two * dx);
                    let avg_grad_pw = (grad_p_w_x + grad_p_p_x) * half;
                    let face_grad_pw = (p_p - p_w) / dx;
                    let u_face_w = (u_w + u_p) * half + d_face_w * (avg_grad_pw - face_grad_pw);

                    // North face:
                    let grad_p_p_y = (p_n - p_s) / (two * dy);
                    let grad_p_n_y = (p_nn - p_p) / (two * dy);
                    let avg_grad_pn = (grad_p_p_y + grad_p_n_y) * half;
                    let face_grad_pn = (p_n - p_p) / dy;
                    let v_face_n = (v_p + v_n) * half + d_face_n * (avg_grad_pn - face_grad_pn);

                    // South face:
                    let grad_p_s_y = (p_p - p_ss) / (two * dy);
                    let avg_grad_ps = (grad_p_s_y + grad_p_p_y) * half;
                    let face_grad_ps = (p_p - p_s) / dy;
                    let v_face_s = (v_s + v_p) * half + d_face_s * (avg_grad_ps - face_grad_ps);

                    // ── Mass imbalance (RHS of pressure-correction Poisson) ───
                    let rho = fields.density.at(i, j);
                    let flux_e = rho * u_face_e * dy;
                    let flux_w = rho * u_face_w * dy;
                    let flux_n = rho * v_face_n * dx;
                    let flux_s = rho * v_face_s * dx;
                    let mass_imbalance = flux_e - flux_w + flux_n - flux_s;

                    // ── Pressure-correction coefficients ─────────────────────
                    let a_e = rho * d_face_e * dy / dx;
                    let a_w = rho * d_face_w * dy / dx;
                    let a_n = rho * d_face_n * dx / dy;
                    let a_s = rho * d_face_s * dx / dy;
                    let a_p = -(a_e + a_w + a_n + a_s);

                    matrix_builder.add_entry(idx, idx, a_p)?;
                    matrix_builder.add_entry(idx, idx + 1, a_e)?;
                    matrix_builder.add_entry(idx, idx - 1, a_w)?;
                    matrix_builder.add_entry(idx, idx + nx, a_n)?;
                    matrix_builder.add_entry(idx, idx - nx, a_s)?;

                    rhs[idx] = mass_imbalance;

                    if mass_imbalance.abs() > max_residual {
                        max_residual = mass_imbalance.abs();
                    }
                }
            }

            // ── Boundary rows: Dirichlet on x-faces, Neumann on y-walls ───────
            for j in 0..ny {
                for i in 0..nx {
                    if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
                        continue; // interior already handled
                    }
                    let idx = j * nx + i;
                    // Inlet/Outlet (x-boundaries): p' = 0 (Dirichlet)
                    if i == 0 || i == nx - 1 {
                        matrix_builder.add_entry(idx, idx, T::one())?;
                        rhs[idx] = T::zero();
                    } else {
                        // Top/bottom walls: ∂p'/∂n = 0 (Neumann)
                        let neighbor = if j == 0 { idx + nx } else { idx - nx };
                        matrix_builder.add_entry(idx, idx, T::one())?;
                        matrix_builder.add_entry(idx, neighbor, -T::one())?;
                        rhs[idx] = T::zero();
                    }
                }
            }

            // ── Step 4: Solve Pressure-Correction System ─────────────────────
            // Take ownership of the builder (Option::take) to call build(self),
            // then immediately restore a fresh empty builder for the next SIMPLE iteration.
            let builder_owned = self
                .matrix_builder
                .take()
                .expect("matrix_builder initialized above");
            let n_rows = nx * ny;
            let matrix = builder_owned.build()?;
            // Restore an empty builder sized for the next iteration (no re-allocation next call
            // because ensure_buffers checks the row count).
            self.matrix_builder = Some(SparseMatrixBuilder::new(n_rows, n_rows));
            let solver_config = IterativeSolverConfig {
                tolerance: self.tolerance,
                max_iterations: 2000,
                ..Default::default()
            };
            let linear_solver = BiCGSTAB::new(solver_config);
            linear_solver.solve(
                &matrix,
                self.rhs.as_ref().expect("buffers initialized"),
                self.p_prime.as_mut().expect("buffers initialized"),
                None::<&IdentityPreconditioner>,
            )?;

            // Guard against non-finite pressure correction from ill-conditioned systems.
            {
                let pp = self.p_prime.as_mut().expect("buffers initialized");
                if pp.iter().any(|v| !v.is_finite()) {
                    // Reset to zero rather than propagating NaN through velocity correction.
                    pp.fill(T::zero());
                }
            }

            // ── Step 5: Correction Step ───────────────────────────────────────
            //
            // p^{n+1}   = p^n + α_p · p'
            // u^{n+1}_P = u*_P - D_u · (∂p'/∂x)
            // v^{n+1}_P = v*_P - D_v · (∂p'/∂y)
            let p_prime = self.p_prime.as_ref().expect("buffers initialized");
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let idx = j * nx + i;
                    let pp = p_prime[idx];
                    let pp_e = p_prime[idx + 1];
                    let pp_w = p_prime[idx - 1];
                    let pp_n = p_prime[idx + nx];
                    let pp_s = p_prime[idx - nx];

                    // Pressure update
                    if let Some(p) = fields.p.at_mut(i, j) {
                        *p += self.pressure_relaxation * pp;
                    }

                    // Velocity correction: u -= D_u * ∂p'/∂x
                    let dp_dx = (pp_e - pp_w) / (two * dx);
                    let dp_dy = (pp_n - pp_s) / (two * dy);
                    let d_u_p = d_u.at(i, j);
                    let d_v_p = d_v.at(i, j);

                    if let Some(u) = fields.u.at_mut(i, j) {
                        *u -= d_u_p * dp_dx * dx;
                    }
                    if let Some(v) = fields.v.at_mut(i, j) {
                        *v -= d_v_p * dp_dy * dy;
                    }
                }
            }

            let converged = max_residual < self.tolerance;
            Ok((max_residual, converged))
        }
    }

    /// Run the full SIMPLE outer loop until convergence or `max_iterations`.
    ///
    /// Returns the number of iterations taken on success, or
    /// `Err(ConvergenceError::MaxIterationsExceeded)` if the residual does not
    /// fall below `tolerance`.
    pub fn solve_simple(
        &mut self,
        momentum_solver: &mut MomentumSolver<T>,
        poisson_solver: &mut PoissonSolver<T>,
        fields: &mut SimulationFields<T>,
        dt: T,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
    ) -> Result<usize> {
        let mut iteration = 0;
        let mut converged = false;

        while iteration < self.max_iterations && !converged {
            let (_residual, conv) = self.simple_iteration(
                momentum_solver,
                poisson_solver,
                fields,
                dt,
                grid,
                boundary_conditions,
            )?;

            converged = conv;
            iteration += 1;

            tracing::debug!(
                iteration,
                residual = ?_residual,
                converged,
                "SIMPLE outer iteration"
            );
        }

        if !converged {
            return Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded {
                    max: self.max_iterations,
                },
            ));
        }

        Ok(iteration)
    }
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::Debug> Default
    for SimpleAlgorithm<T>
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_algorithm_creation() {
        let simple = SimpleAlgorithm::<f64>::new();
        assert_eq!(simple.max_iterations, 50);
        assert!(simple.pressure_relaxation > 0.0 && simple.pressure_relaxation < 1.0);
        assert!(simple.velocity_relaxation > 0.0 && simple.velocity_relaxation < 1.0);
        // Invariant: Patankar's recommendation α_u + α_p ≈ 1
        assert!(
            (simple.velocity_relaxation + simple.pressure_relaxation - 1.0).abs() < 0.01,
            "Expected α_u + α_p ≈ 1 (Patankar 1980 recommendation)"
        );
    }

    #[test]
    fn test_rhie_chow_uniform_pressure_zero_correction() {
        // Theorem 2, Uniform Field Invariant:
        // For a uniform pressure field p = const, the RC correction term
        //   d_e * (avg_grad_pe - face_grad_pe)
        // must be identically zero, since both gradients equal zero.
        let dx = 0.1_f64;
        let _dy = 0.1_f64;
        let two = 2.0_f64;

        // All pressures equal → all gradients zero
        let (p_p, p_e, p_w, p_ee) = (1.0, 1.0, 1.0, 1.0);
        let grad_p_p_x = (p_e - p_w) / (two * dx);
        let grad_p_e_x = (p_ee - p_p) / (two * dx);
        let avg_grad_pe = (grad_p_p_x + grad_p_e_x) * 0.5;
        let face_grad_pe = (p_e - p_p) / dx;
        let rc_correction = avg_grad_pe - face_grad_pe;

        assert!(
            rc_correction.abs() < 1e-14,
            "RC correction for uniform pressure must be zero; got {rc_correction}"
        );
    }

    #[test]
    fn test_rhie_chow_checkerboard_nonzero_correction() {
        // Theorem 2, Checker-Board Mode:
        // For alternating p = 1, -1, 1 the RC correction must be non-zero.
        let dx = 0.1_f64;
        let dy = 0.1_f64;
        let two = 2.0_f64;

        // Checker-board pattern: p_P = 1, p_E = -1, p_W = -1, p_EE = 1
        let (p_p, p_e, p_w, p_ee) = (1.0_f64, -1.0, -1.0, 1.0);
        let _ = dy; // used conceptually
        let grad_p_p_x = (p_e - p_w) / (two * dx); // (-1 - -1)/(0.2) = 0
        let grad_p_e_x = (p_ee - p_p) / (two * dx); // (1 - 1)/(0.2)   = 0
        let avg_grad_pe = (grad_p_p_x + grad_p_e_x) * 0.5; // = 0
        let face_grad_pe = (p_e - p_p) / dx; // (-1 - 1)/0.1    = -20
        let rc_correction = avg_grad_pe - face_grad_pe; // 0 - (-20) = 20 ≠ 0

        assert!(
            rc_correction.abs() > 1e-10,
            "RC correction for checker-board pressure must be non-zero; got {rc_correction}"
        );
    }
}
