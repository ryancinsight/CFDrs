//! Finite Element Method (FEM) solver for 3D incompressible flows
//!
//! Reference implementations based on:
//! - Hughes et al. (1986) "A new finite element formulation for computational fluid dynamics"
//! - Brooks & Hughes (1982) "Streamline upwind/Petrov-Galerkin formulations"

pub mod config;
pub mod elements;
pub mod assembly;
pub mod stabilization;
pub mod boundary;
pub mod solver;

pub use config::FemConfig;
pub use solver::FemSolver;
pub use elements::{Element, ElementType};
pub use stabilization::{StabilizationType, StabilizationMethod};