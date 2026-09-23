//! # Discontinuous Galerkin (DG) Methods
//!
//! This module provides a high-performance implementation of Discontinuous Galerkin (DG) methods
//! for solving partial differential equations (PDEs). DG methods combine the flexibility of
//! finite volume methods with the high-order accuracy of spectral methods, making them
//! particularly well-suited for problems with sharp gradients, discontinuities, and complex
//! geometries.
//!
//! ## Features
//!
//! - **High-order accuracy**: Support for arbitrary polynomial orders
//! - **Flexible time integration**: Multiple explicit, implicit, and IMEX time integration schemes
//! - **Shock-capturing**: Built-in limiters for handling discontinuities
//! - **Adaptive mesh refinement**: Support for h- and p-adaptivity
//! - **Parallel computing**: Designed for efficient parallel execution
//!
//! ## Quick Start
//!
//! ```no_run
//! use cfd_math::high_order::dg::*;
//! use cfd_math::error::Result;
//! use leto::{Array1, Array2};
//!
//! // Create a DG operator
//! let order = 3;
//! let num_components = 1;
//! let params = DGOperatorParams::new()
//!     .with_volume_flux(FluxType::Central)
//!     .with_surface_flux(FluxType::LaxFriedrichs)
//!     .with_limiter(LimiterType::Minmod);
//!
//! let dg_op = DGOperator::new(order, num_components, Some(params)).expect("expected value");
//!
//! // Create a time integrator
//! let integrator = TimeIntegratorFactory::create(TimeIntegration::SSPRK3);
//!
//! // Set up the solver
//! let t_final = 1.0;
//! let solver_params = TimeIntegrationParams::new(TimeIntegration::SSPRK3)
//!     .with_t_final(t_final)
//!     .with_dt(0.01)
//!     .with_verbose(true);
//!
//! let mut solver = DGSolver::new(dg_op, integrator, solver_params);
//!
//! // Set initial condition
//! let u0 = |x: f64| Array1::from_shape_vec([1], vec![x.sin()]).expect("expected value");
//! solver.initialize(u0).expect("expected value");
//!
//! // Define the right-hand side function
//! fn f(_t: f64, u: &Array2<f64>) -> Result<Array2<f64>> {
//!     // For the linear advection equation: du/dt = -du/dx
//!     // This is a simplified example; in practice, you would use the DG operator
//!     // to compute the spatial derivatives
//!     Ok(Array2::from_shape_fn(u.shape(), |idx| -u[idx]))
//! }
//!
//! // Run the solver
//! solver.solve(f, None::<fn(f64, &Array2<f64>) -> Result<Array2<f64>>>).expect("expected value");
//!
//! // Evaluate the solution
//! let x = 0.5;
//! let u = solver.evaluate(x);
//! ```
//!
//! ## Modules
//!
//! ## Modules
//!
//! - `basis`: Basis functions and quadrature rules
//! - `flux`: Numerical flux functions
//! - `linalg`: `leto` array shims shared by the operator, limiter and solver leaves
//! - `limiter`: Slope limiters for shock capturing
//! - `operators`: Core DG operators and discretization
//! - `solution`: The `DGSolution` carrier and the `DGMethod` trait
//! - `solver`: Time integration and solution algorithms
//!
//! This module root is a passthrough facade: it declares the concern leaves and
//! re-exports their public surface, so every path that existed before the split
//! (`dg::DGSolution`, `dg::matrix_cols`, ...) still resolves.


#![warn(missing_docs)]

mod basis;
mod flux;
mod linalg;
mod limiter;
mod operators;
mod solution;
mod solver;

pub use basis::*;
pub use flux::*;
pub use limiter::*;
pub use operators::*;
pub use solution::*;
pub use solver::*;

// The array shims were declared in this module root before the split, so they
// stay reachable at `dg::<shim>` for the sibling leaves and for in-crate callers.
pub(crate) use linalg::*;

#[cfg(test)]
mod tests;
