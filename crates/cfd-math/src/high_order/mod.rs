//! # High-Order Numerical Methods for CFD
//!
//! This module provides high-performance implementations of advanced numerical methods
//! for solving partial differential equations (PDEs) in computational fluid dynamics (CFD).
//! The implementations are designed for accuracy, efficiency, and ease of use.
//!
//! ## Features
//!
//! - **Discontinuous Galerkin (DG) Methods**: High-order accurate, flexible, and robust
//!   - Arbitrary polynomial orders
//!   - Support for explicit, implicit, and IMEX time integration
//!   - Built-in limiters for shock capturing
//!   - Adaptive mesh refinement (h- and p-adaptivity)
//!   - Parallel computing support
//!
//! - **Spectral Element Methods (SEM)**:
//!   - High-order accuracy with exponential convergence
//!   - Efficient quadrature rules
//!   - Tensor-product elements for multi-dimensional problems
//!
//! - **Weighted Essentially Non-Oscillatory (WENO) Schemes**:
//!   - High-order accurate shock-capturing schemes
//!   - Multiple reconstruction variants
//!   - Adaptive order selection
//!
//! ## Getting Started
//!
//! ### DG Method Example
//! ```no_run
//! use cfd_math::high_order::*;
//! use nalgebra::DVector;
//!
//! // Create a DG operator
//! let order = 3;
//! let num_components = 1;
//! let params = DGOperatorParams::new()
//!     .with_volume_flux(FluxType::Central)
//!     .with_surface_flux(FluxType::LaxFriedrichs)
//!     .with_limiter(LimiterType::Minmod);
//! 
//! let dg_op = DGOperator::new(order, num_components, Some(params)).unwrap();
//! 
//! // Set up and run the solver
//! let mut solver = DGSolver::new(
//!     dg_op,
//!     TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
//!     TimeIntegrationParams::new(TimeIntegration::SSPRK3)
//!         .with_t_final(1.0)
//!         .with_dt(0.01)
//! );
//! 
//! // Initialize and solve
//! solver.initialize(|x| DVector::from_vec(vec![x.sin()])).unwrap();
//! solver.solve(|_t, u| Ok(-u.clone()), None::<fn(_, _) -> _>).unwrap();
//! ```
//!
//! ## Performance Considerations
//!
//! - **Vectorization**: The implementations are designed to take advantage of
//!   SIMD instructions when available.
//! - **Memory Layout**: Data structures are optimized for cache efficiency.
//! - **Parallelism**: The code is designed to be easily parallelizable using
//!   Rust's concurrency primitives or external crates like `rayon`.
//!
//! ## References
//!
//! 1. Hesthaven, J. S., & Warburton, T. (2008). Nodal Discontinuous Galerkin Methods.
//! 2. Canuto, C., et al. (2006). Spectral Methods: Fundamentals in Single Domains.
//! 3. Shu, C.-W. (2009). High Order Weighted Essentially Nonoscillatory Schemes for
//!    Convection Dominated Problems.

#![warn(missing_docs)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::float_cmp)]

pub mod dg;
pub mod spectral;
pub mod weno;

/// Re-export commonly used items
pub use dg::*;
pub use spectral::*;
pub use weno::*;
