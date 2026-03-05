//! Finite Difference Method (FDM) solvers for 2D CFD problems.
//!
//! This module provides finite difference implementations for solving
//! various 2D fluid dynamics problems with proper separation of concerns.
//!
//! # Theorem (FDM Consistency and Convergence — Lax–Richtmyer 1956)
//!
//! For a consistent and stable finite difference scheme applied to a well-posed linear
//! initial-boundary value problem, the numerical solution converges to the exact solution
//! as $\Delta x, \Delta t \to 0$.
//!
//! **Proof sketch**:
//! *Consistency*: Taylor expansion of the FD stencil recovers the PDE
//! $\partial u/\partial t + \mathcal{L}u = 0$ with truncation error
//! $\tau = O(\Delta x^p + \Delta t^q)$.
//!
//! *Stability*: The von Neumann amplification factor $|g(\xi)| \le 1 + C\Delta t$
//! for all wave numbers $\xi$. For the diffusion operator with central differences,
//! this requires $\alpha \Delta t / \Delta x^2 \le 1/2$. For convection with
//! first-order upwind, the CFL condition $|u|\Delta t / \Delta x \le 1$ suffices.
//!
//! *Convergence*: By the Lax equivalence theorem, consistency + stability $\Rightarrow$
//! convergence. The global error $\|e^n\| = \|u^n - u(t_n)\| \to 0$ at rate
//! $O(\Delta x^p + \Delta t^q)$.

pub mod advection_diffusion;
pub mod config;
pub mod diffusion;
pub mod linear_solver;
pub mod poisson;

// Re-export main types
pub use advection_diffusion::AdvectionDiffusionSolver;
pub use config::FdmConfig;
pub use diffusion::DiffusionSolver;
pub use linear_solver::solve_gauss_seidel;
pub use poisson::PoissonSolver;
