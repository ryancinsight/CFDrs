//! 1D CFD simulations and microfluidic network solvers.
//!
//! This crate provides 1D computational fluid dynamics solvers specifically
//! designed for microfluidic and millifluidic networks. It includes support
//! for:
//! - Channel-based flow networks
//! - Electrical circuit analogy solvers
//! - Pressure-driven flow simulations
//! - Component-based microfluidic devices
//!
//! ## 2D Schematic Support
//!
//! This crate is designed to integrate with the `scheme` library for
//! 2D schematic representation of 1D microfluidic networks, similar to
//! how `cfd-3d` integrates with `csgrs` for 3D mesh handling.

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

pub mod analysis;
pub mod channel;
pub mod components;
pub mod junctions;
pub mod network;
pub mod resistance;

pub mod solver;
pub mod vascular;

// Export network functionality
pub use network::{
    BoundaryCondition, ChannelProperties, Edge, EdgeProperties, EdgeType, Network, NetworkBuilder,
    NetworkGraph, NetworkMetadata, Node, NodeProperties, NodeType,
};

// Export solver functionality
pub use solver::{
    ConvergenceChecker, LinearSystemSolver, MatrixAssembler, NetworkDomain, NetworkProblem,
    NetworkSolver, NetworkState, CompositionState, EdgeFlowEvent, InletCompositionEvent, MixtureComposition,
    PressureBoundaryEvent, SimulationTimeConfig, TransientCompositionSimulator, ChannelOccupancy, DropletBoundary, DropletInjection,
    DropletPosition, DropletSnapshot, DropletSplitPolicy, DropletState,
    DropletTrackingState, SplitMode, TransientDropletSimulator,
};

// Export component functionality
pub use components::{
    CircularChannel, Component, ComponentFactory, FlowSensor, Micromixer, Micropump, Microvalve,
    MixerType, OrganCompartment, PorousMembrane, PumpType, RectangularChannel, SensorType,
    ValveType,
};

// Export analysis functionality
pub use analysis::{
    BloodShearLimits, FlowAnalysis, NetworkAnalysisResult, NetworkAnalyzerOrchestrator,
    PerformanceMetrics, PressureAnalysis, ResistanceAnalysis, ShearLimitViolation,
};

// Export channel functionality
pub use channel::{
    Channel, ChannelGeometry, ChannelType, CrossSection, FlowRegime, FlowState,
    NumericalParameters, SurfaceProperties, Wettability,
};

// Export resistance functionality
pub use resistance::{
    BendType, ChannelGeometry as ResistanceChannelGeometry, CombinationMethod,
    DarcyWeisbachModel, ExpansionType, FlowConditions, HagenPoiseuilleModel,
    MembranePoreModel, RectangularChannelModel, ResistanceCalculator, ResistanceModel,
    ResistanceModelFactory, SerpentineAnalysis, SerpentineCrossSection, SerpentineModel,
    VenturiAnalysis, VenturiGeometry, VenturiModel,
};

// Export hierarchical junction functionality
pub use junctions::branching::{
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
        analysis::{
            BloodShearLimits, NetworkAnalysisResult, NetworkAnalyzerOrchestrator,
            PerformanceMetrics, ShearLimitViolation,
        },
        channel::{
            Channel, ChannelGeometry, ChannelType, CrossSection, FlowRegime, FlowState,
            SurfaceProperties, Wettability,
        },
        components::{
            ComponentFactory, FlowSensor, Micromixer, MixerType, OrganCompartment,
            PorousMembrane, PumpType, SensorType, ValveType,
        },
        network::{Edge, EdgeProperties, EdgeType, Node, NodeProperties, NodeType},
        junctions::branching::{
            BranchingNetworkConfig, BranchingNetworkSolver, BranchingValidationResult,
            BranchingValidator, ThreeWayBranchJunction, ThreeWayBranchSolution,
            TwoWayBranchJunction, TwoWayBranchSolution,
        },
        resistance::{
            DarcyWeisbachModel, FlowConditions, ResistanceCalculator, ResistanceModelFactory,
        },
        solver::{
            ChannelOccupancy, CompositionState, DropletBoundary, DropletInjection,
            DropletPosition, DropletSnapshot, DropletSplitPolicy, DropletState,
            DropletTrackingState, EdgeFlowEvent, InletCompositionEvent, MixtureComposition, PressureBoundaryEvent, SimulationTimeConfig, SplitMode,
            TransientCompositionSimulator, TransientDropletSimulator,
        },
    };
}
