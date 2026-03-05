//! # cfd-1d — 1D Lumped-Network Microfluidic Solver
//!
//! `cfd-1d` solves **pressure-driven flow in microfluidic channel networks**
//! using the electrical circuit analogy: each channel is a lumped hydraulic
//! resistance and the network is solved as a linear system (Kirchhoff's laws).
//!
//! ## Physical Model
//!
//! The 1D model makes the following assumptions:
//! - Flow is **fully developed** in every cross-section (Hagen-Poiseuille).
//! - Channels are **straight** or have analytically correctable geometry.
//! - Inertial effects are negligible (low Reynolds number, Re ≪ 1).
//! - The fluid state at a node is described by a single scalar: **pressure**.
//!
//! Under these assumptions the governing equation per edge is:
//!
//! ```text
//!   Q = ΔP / R     where R = 128·μ·L / (π·D⁴)   (circular)
//!                        R = 12·μ·L / (w·h³·φ)   (rectangular, φ = shape factor)
//! ```
//!
//! This collapses the full 3D Navier-Stokes problem to a sparse linear system
//! solvable in microseconds, making it ideal for rapid design-space exploration.
//!
//! ## Relationship with `cfd-schematics`
//!
//! `cfd-schematics` is the **topology and geometry authority**. `cfd-1d` does
//! not define its own network types — it consumes:
//! - `NodeSpec` → `Node<T>` (via `From<&NodeSpec>`)
//! - `ChannelSpec` → `EdgeProperties<T>` (via `From<&ChannelSpec>`)
//! - `CrossSectionSpec` → `CrossSection<T>` + `ChannelGeometry<T>`
//!
//! The 2D schematic rendered by `cfd-schematics` is a **layout diagram** of
//! the network topology — it is not a simulation domain. The schematic's x-y
//! coordinates place channels on a chip diagram; the 1D solver only sees
//! channel length and cross-section, not position.
//!
//! ## Relationship with `cfd-2d`
//!
//! | Aspect              | `cfd-1d`                          | `cfd-2d`                              |
//! |---------------------|-----------------------------------|---------------------------------------|
//! | Governing equations | Kirchhoff / Hagen-Poiseuille      | Navier-Stokes (u, v, p fields)        |
//! | Spatial resolution  | Per-channel scalar (Q, ΔP)        | Per-cell vector/scalar field          |
//! | Grid                | Graph (nodes + edges)             | Structured / unstructured 2D grid     |
//! | Solvers             | Sparse linear system              | FDM, FVM, LBM, SIMPLE, PISO          |
//! | Turbulence          | Not applicable (Re ≪ 1)           | k-ε, k-ω SST, DES, Reynolds stress   |
//! | Cost                | Microseconds                      | Seconds to hours                      |
//! | Use case            | Network design, flow distribution | Detailed velocity/pressure fields     |
//!
//! `cfd-1d` results (pressure at nodes, flow rate per channel) can seed
//! boundary conditions for a `cfd-2d` simulation of a single channel segment.
//!
//! ## Modules
//! - **network**: Graph, `NetworkBuilder`, `EdgeProperties`, `Network` wrapper
//! - **solver**: `NetworkSolver`, `NetworkProblem`, `SolverConfig`
//! - **channel**: `CrossSection`, `ChannelGeometry`, `SurfaceProperties`
//! - **analysis**: `NetworkAnalyzerOrchestrator`, `FlowAnalysis`, `BloodShearLimits`
//! - **resistance**: Hydraulic resistance calculators
//! - **components**: Valves, pumps, filters
//! - **junctions**: T-junction, Y-junction, bifurcation models

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// 1D CFD network simulation allows
#![allow(clippy::similar_names)] // Network variables (p1,p2; q1,q2) often similar
#![allow(clippy::cast_precision_loss)] // Acceptable in flow calculations
#![allow(clippy::cast_possible_truncation)] // Network sizes typically small
#![allow(clippy::unused_self)] // Analyzer trait consistency
#![allow(clippy::must_use_candidate)] // Network analysis utilities
#![allow(clippy::missing_errors_doc)] // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)] // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)] // Signed to unsigned casts common in CFD indexing
#![allow(clippy::cast_possible_wrap)] // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)] // CFD functions often need many physical parameters
#![allow(clippy::float_cmp)] // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)] // Result types maintained for API consistency
#![allow(clippy::items_after_statements)] // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)] // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)] // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)] // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)] // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)] // Builder patterns used internally
#![allow(clippy::ptr_arg)] // &Vec used for API compatibility
#![allow(clippy::should_implement_trait)] // CFD-specific trait implementations

/// Domain modeling: topology, components, and channel geometries
pub mod domain;
/// Physical constraints, resistance laws, and biological models
pub mod physics;
/// Application orchestration: solvers and post-processors
pub mod solver;

// Export haemolysis models
pub use physics::hemolysis::{
    cavitation_amplified_hi, giersiepen_hi, HemolysisExposure,
    CAVITATION_HI_SLOPE, GIERSIEPEN_ALPHA, GIERSIEPEN_BETA, GIERSIEPEN_C,
};

// Export three-population inertial focusing model
pub use physics::cell_separation::{three_population_equilibria, ThreePopEquilibria};
pub use physics::cell_separation::{
    cascade_junction_separation, cascade_junction_separation_cross_junction,
    cascade_junction_separation_from_qfracs,
    cif_pretri_stage_center_fracs, cif_pretri_stage_q_fracs,
    incremental_filtration_separation, incremental_filtration_separation_cross_junction,
    incremental_filtration_separation_from_qfracs,
    incremental_filtration_separation_staged, tri_center_q_frac,
    tri_center_q_frac_cross_junction, CascadeJunctionResult,
    IncrementalFiltrationResult,
};

// Export network functionality
pub use domain::network::{
    BoundaryCondition, ChannelProperties, ComponentType, Edge, EdgeProperties, EdgeType,
    Network, NetworkBuilder, NetworkBuilderSink, NetworkGraph, NetworkMetadata, Node,
    NodeProperties, NodeType,
};

// Export solver functionality
pub use solver::core::{
    ConvergenceChecker, LinearSystemSolver, MatrixAssembler, NetworkDomain, NetworkProblem,
    NetworkSolver, NetworkState, SolverConfig, CompositionState, EdgeFlowEvent,
    InletCompositionEvent, MixtureComposition, PressureBoundaryEvent, SimulationTimeConfig,
    TransientCompositionSimulator, ChannelOccupancy, DropletBoundary, DropletInjection,
    DropletPosition, DropletSnapshot, DropletSplitPolicy, DropletState,
    DropletTrackingState, SplitMode, TransientDropletSimulator,
};

// Export component functionality
pub use domain::components::{
    CircularChannel, Component, ComponentFactory, FlowSensor, Micromixer, Micropump, Microvalve,
    MixerType, OrganCompartment, PorousMembrane, PumpType, RectangularChannel, SensorType,
    ValveType,
};

// Export analysis functionality
pub use solver::analysis::{
    BloodShearLimits, FlowAnalysis, NetworkAnalysisResult, NetworkAnalyzerOrchestrator,
    PerformanceMetrics, PressureAnalysis, ResistanceAnalysis, ShearLimitViolation,
};

// Export channel functionality
pub use domain::channel::{
    Channel, ChannelGeometry, ChannelType, CrossSection, FlowRegime, FlowState,
    KN_SLIP_MIN, NumericalParameters, SurfaceProperties, Wettability,
};

// Export resistance functionality
pub use physics::resistance::{
    BendType, ChannelGeometry as ResistanceChannelGeometry, CombinationMethod,
    DarcyWeisbachModel, ExpansionType, FlowConditions, HagenPoiseuilleModel,
    JunctionFlowDirection, JunctionLossModel, JunctionType, MembranePoreModel,
    RectangularChannelModel, ResistanceCalculator, ResistanceModel, ResistanceModelFactory,
    SerpentineAnalysis, SerpentineCrossSection, SerpentineModel, VenturiAnalysis,
    VenturiGeometry, VenturiModel,
};

// Export hierarchical junction functionality
pub use domain::junctions::branching::{
    BranchingNetworkConfig, BranchingNetworkSolver, BranchingValidationResult,
    BranchingValidator, ThreeWayBranchJunction, ThreeWayBranchSolution,
    TwoWayBranchJunction, TwoWayBranchSolution,
};

/// 1D CFD domain-specific prelude for microfluidic networks
///
/// This prelude exports 1D-specific functionality not available in the main prelude.
/// Use this when working extensively with microfluidic networks and need access to
/// specialized components, analysis tools, or scheme integration features.
///
/// For basic 1D functionality, prefer `cfd_suite::prelude::*`.
pub mod prelude {
    // === Extended Network Components ===
    // Specialized components not in main prelude
    pub use crate::{
        solver::analysis::{
            BloodShearLimits, NetworkAnalysisResult, NetworkAnalyzerOrchestrator,
            PerformanceMetrics, ShearLimitViolation,
        },
        domain::channel::{
            Channel, ChannelGeometry, ChannelType, CrossSection, FlowRegime, FlowState,
            SurfaceProperties, Wettability,
        },
        domain::components::{
            ComponentFactory, FlowSensor, Micromixer, MixerType, OrganCompartment,
            PorousMembrane, PumpType, SensorType, ValveType,
        },
        domain::network::{Edge, EdgeProperties, EdgeType, Node, NodeProperties, NodeType},
        domain::junctions::branching::{
            BranchingNetworkConfig, BranchingNetworkSolver, BranchingValidationResult,
            BranchingValidator, ThreeWayBranchJunction, ThreeWayBranchSolution,
            TwoWayBranchJunction, TwoWayBranchSolution,
        },
        physics::resistance::{
            DarcyWeisbachModel, FlowConditions, ResistanceCalculator, ResistanceModelFactory,
        },
        solver::core::{
            ChannelOccupancy, CompositionState, DropletBoundary, DropletInjection,
            DropletPosition, DropletSnapshot, DropletSplitPolicy, DropletState,
            DropletTrackingState, EdgeFlowEvent, InletCompositionEvent, MixtureComposition, PressureBoundaryEvent, SimulationTimeConfig, SplitMode,
            TransientCompositionSimulator, TransientDropletSimulator,
        },
    };
}
