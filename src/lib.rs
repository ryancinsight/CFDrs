//! Computational Fluid Dynamics Simulation Suite
//!
//! A comprehensive CFD simulation framework in Rust supporting 1D, 2D, and 3D simulations
//! with a plugin-based architecture.

#![warn(missing_docs)]

// Re-export all sub-crates
pub use cfd_core as core;
pub use cfd_math as math;
pub use cfd_io as io;
pub use cfd_mesh as mesh;
pub use cfd_1d as d1;
pub use cfd_2d as d2;
pub use cfd_3d as d3;
pub use cfd_validation as validation;

/// Prelude module for convenient imports with explicit re-exports to avoid ambiguity
pub mod prelude {
    // Core functionality
    pub use cfd_core::{Fluid, Error, Result};
    pub use cfd_core::BoundaryCondition as CoreBoundaryCondition;

    // Math functionality
    pub use cfd_math::{SparseMatrix, SparseMatrixBuilder};
    pub use cfd_math::integration::{GaussQuadrature, AdaptiveQuadrature};

    // I/O functionality
    pub use cfd_io::{VtkWriter, CsvWriter, JsonWriter};

    // Mesh functionality
    pub use cfd_mesh::Mesh;
    // TODO: Add MeshBuilder, Element, Node, Edge when implemented in cfd-mesh

    // 1D CFD functionality
    pub use cfd_1d::{Network, NetworkBuilder, NetworkSolver, SolverConfig};
    pub use cfd_1d::{Node as NetworkNode, Edge as NetworkEdge};
    pub use cfd_1d::BoundaryCondition as NetworkBoundaryCondition;
    pub use cfd_1d::{Component, RectangularChannel, CircularChannel, Micropump, Microvalve};

    // Validation functionality
    pub use cfd_validation::AnalyticalSolution;
    // TODO: Add ValidationSuite when implemented
}