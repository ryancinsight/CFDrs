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
    // Core functionality - fundamental types used across all modules
    pub use cfd_core::{Fluid, Error, Result, BoundaryCondition, Domain, Problem, Solver};
    pub use cfd_core::{Plugin, PluginRegistry, SimulationPlugin};

    // Math functionality - essential numerical methods
    pub use cfd_math::{SparseMatrix, SparseMatrixBuilder, LinearSolver, ConjugateGradient};
    pub use cfd_math::{GaussQuadrature, FiniteDifference};
    pub use cfd_math::integration::AdaptiveQuadrature;

    // I/O functionality - common file operations
    pub use cfd_io::{VtkWriter, CsvWriter, JsonWriter, CheckpointManager};
    #[cfg(feature = "hdf5-support")]
    pub use cfd_io::{Hdf5Writer, Hdf5Reader, DatasetMetadata, DataChunk};

    // Mesh functionality - geometric operations
    pub use cfd_mesh::{Mesh, Vertex, Face, Cell, MeshTopology};

    // 1D CFD functionality - microfluidic networks
    pub use cfd_1d::{Network, NetworkBuilder, NetworkSolver, SolverConfig};
    pub use cfd_1d::{Component, RectangularChannel, CircularChannel, Micropump, Microvalve};
    pub use cfd_1d::{NetworkAnalyzer, FlowAnalysis, PressureAnalysis};

    // 2D CFD functionality - grid-based methods
    pub use cfd_2d::{StructuredGrid2D, PoissonSolver, FvmSolver, LbmSolver};

    // 3D CFD functionality - advanced solvers
    pub use cfd_3d::{FemSolver, SpectralSolver, MeshAdapter, StlAdapter};

    // Validation functionality - testing and verification
    pub use cfd_validation::{AnalyticalSolution, ErrorMetric, ConvergenceAnalysis};
}