//! Finite Element Method (FEM) solver for 3D incompressible flows.
//!
//! This module implements a mixed finite element formulation for the Stokes
//! and Navier-Stokes equations with stabilization.
//!
//! # Theorem — Discrete Inf-Sup Stability (Taylor–Hood)
//!
//! For the P2-P1 velocity-pressure pair used by this module, the discrete
//! inf-sup condition holds on shape-regular tetrahedral meshes:
//!
//! ```text
//! inf_{q_h} sup_{v_h} (div v_h, q_h) / (||v_h||_1 ||q_h||_0) >= β > 0
//! ```
//!
//! ensuring unique solvability and pressure stability of the mixed system.
//!
//! ## Literature References
//! - Hughes, T.J.R., Franca, L.P., Balestra, M. (1986). "A new finite element formulation
//!   for computational fluid dynamics: V. Circumventing the Babuška-Brezzi condition:
//!   A stable Petrov-Galerkin formulation of the Stokes problem accommodating equal-order interpolations"
//! - Donea, J., Huerta, A. (2003). "Finite Element Methods for Flow Problems"

pub mod boundary_classifier;
pub mod config;
pub mod constants;
pub mod element;
pub mod fluid;
pub mod mesh_utils;
pub mod mid_node_cache;
pub mod problem;
pub mod projection_solver;
pub mod quadrature;
pub mod shape_functions;
pub mod solution;
pub mod solver;
pub mod stabilization;
pub mod stress;

// Re-export main types for convenience
pub use boundary_classifier::{AxialBoundaryClassifier, BoundaryFaceSets};
pub use config::FemConfig;
pub use element::{ElementMatrices, FluidElement};
pub use fluid::FluidProperties;
pub use mid_node_cache::MidNodeCache;
pub use problem::StokesFlowProblem;
pub use projection_solver::ProjectionSolver;
pub use solution::StokesFlowSolution;
pub use solver::FemSolver;
