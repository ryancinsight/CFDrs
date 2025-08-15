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

/// Unified prelude module - Single Source of Truth for all common CFD functionality
///
/// This module provides a comprehensive set of imports for CFD simulations across all dimensions.
/// It follows the SSOT (Single Source of Truth) principle by being the primary interface for users.
///
/// # Usage
/// ```rust
/// use cfd_suite::prelude::*;
///
/// // Now you have access to all commonly used CFD functionality
/// let fluid = Fluid::<f64>::water();
/// let network = NetworkBuilder::<f64>::new();
/// let grid = StructuredGrid2D::<f64>::new(10, 10, 1.0, 1.0, 0.0, 1.0);
/// ```
pub mod prelude {
    // === Core Abstractions ===
    // Fundamental types used across all CFD simulations
    pub use cfd_core::{
        Fluid, Error, Result, BoundaryCondition, Domain, Problem, Solver,
        Plugin, PluginRegistry, SimulationPlugin, SimulationState, TimeIntegrator,
        SolverConfiguration, Configurable, Validatable, IterativeSolver, DirectSolver, SolverConfig, LinearSolverConfig, NetworkSolverConfig,
        ConcreteSolverFactory, FactoryCapability, SolverFactoryRegistry, Builder, ConfigurationBuilder, DynamicFactory, DynamicSolver, TypeErasedFactory, TypeErasedSolver
    };

    // === Mathematical Utilities ===
    // Essential numerical methods and linear algebra
    pub use cfd_math::{
                 SparseMatrix, SparseMatrixBuilder, LinearSolver, ConjugateGradient, BiCGSTAB,
GaussQuadrature, FiniteDifference, Interpolation, LinearInterpolation, CubicSplineInterpolation,
        MathIteratorExt, VectorOps, SliceOps, CfdIteratorChain, VectorizedOps, StencilOps
    };
    pub use cfd_math::integration::{AdaptiveQuadrature, VariableQuadrature};

    // === I/O Operations ===
    // File input/output for all supported formats
    pub use cfd_io::{VtkWriter, CsvWriter, JsonWriter, CheckpointManager, VtkMesh, VtkMeshBuilder};
    #[cfg(feature = "hdf5")]
    pub use cfd_io::{Hdf5Writer, Hdf5Reader, DatasetMetadata, DataChunk};

    // === Mesh Operations ===
    // Geometric operations and mesh handling
    pub use cfd_mesh::{Mesh, Vertex, Face, Cell, MeshTopology};

    // === 1D CFD Simulations ===
    // Microfluidic networks and pipe flow
    pub use cfd_1d::{
        Network, NetworkBuilder, NetworkSolver,
        Component, RectangularChannel, CircularChannel, Micropump, Microvalve,
        NetworkAnalyzer, FlowAnalysis, PressureAnalysis, ResistanceModel
    };

    // === 2D CFD Simulations ===
    // Grid-based methods for 2D flows
    pub use cfd_2d::{
        StructuredGrid2D, Grid2D, BoundaryType,
        PoissonSolver, FvmSolver, LbmSolver, PressureVelocityCouplerSolver,
        FdmConfig, FvmConfig, LbmConfig, PressureVelocityCouplingConfig
    };

    // === 3D CFD Simulations ===
    // Advanced 3D methods with CSGrs integration
    pub use cfd_3d::{
        FemSolver, FemConfig, SpectralSolver, SpectralConfig, SpectralBasis,
        IbmSolver, IbmConfig, LagrangianPoint,
        LevelSetSolver, LevelSetConfig,
        VofSolver, VofConfig
    };

    // === Validation Framework ===
    // Analytical solutions and error analysis
    pub use cfd_validation::{
        AnalyticalSolution, PoiseuilleFlow, CouetteFlow, StokesFlow, TaylorGreenVortex,
        ErrorMetric, L2Norm, L1Norm, LInfNorm, ErrorStatistics, ErrorAnalysis,
        ConvergenceAnalysis, ConvergenceStudy, RichardsonExtrapolation
    };
}