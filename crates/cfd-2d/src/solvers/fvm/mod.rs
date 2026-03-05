//! Finite Volume Method (FVM) solver for 2D CFD simulations.
//!
//! # Theorem (FVM Discrete Conservation — Godunov 1959, LeVeque 2002)
//!
//! A finite volume discretization satisfies exact discrete conservation:
//! the integral of the conserved quantity over the domain changes only
//! through boundary fluxes, not through internal discretization error.
//!
//! **Proof sketch**:
//! The FVM integrates the divergence form of the conservation law
//! $\partial \phi / \partial t + \nabla \cdot \mathbf{F} = S$ over each control volume
//! $\Omega_i$, yielding $\frac{d}{dt}\int_{\Omega_i}\phi\,dV = -\oint_{\partial\Omega_i}\mathbf{F}\cdot\mathbf{n}\,dA + \int_{\Omega_i}S\,dV$.
//!
//! Summing over all internal cells, each interior face flux appears exactly twice with
//! opposite signs (once for each adjacent cell), so all internal fluxes cancel
//! telescopically. The global change therefore equals the sum of boundary fluxes
//! plus source terms — i.e., exact discrete conservation regardless of mesh
//! resolution or flux scheme accuracy.

pub mod config;
pub mod flux;
pub mod geometry;
pub mod solver;

// Re-export main types
pub use config::FvmConfig;
pub use flux::{FluxScheme, FluxSchemeFactory};
pub use geometry::Face;
pub use solver::FvmSolver;
