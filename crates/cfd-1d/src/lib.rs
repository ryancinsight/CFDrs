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

pub mod analysis;
pub mod channel;
pub mod components;
pub mod network;
pub mod resistance;
pub mod solver;

// Export network functionality
pub use network::{
    BoundaryCondition, ChannelProperties, Edge, EdgeProperties, EdgeType, Network, NetworkBuilder,
    NetworkGraph, NetworkMetadata, Node, NodeProperties, NodeType,
};

// Export solver functionality
pub use solver::{
    ConvergenceChecker, LinearSystemSolver, MatrixAssembler, NetworkDomain, NetworkProblem,
    NetworkSolver, NetworkState,
};

// Export component functionality
pub use components::{
    CircularChannel, Component, ComponentFactory, FlowSensor, Micromixer, Micropump, Microvalve,
    MixerType, PumpType, RectangularChannel, SensorType, ValveType,
};

// Export analysis functionality
pub use analysis::{
    FlowAnalysis, NetworkAnalysisResult, NetworkAnalyzer, PerformanceMetrics, PressureAnalysis,
    ResistanceAnalysis,
};

// Export channel functionality
pub use channel::{
    Channel, ChannelGeometry, ChannelType, CrossSection, FlowRegime, FlowState,
    NumericalParameters, SurfaceProperties, Wettability,
};

// Export resistance functionality
pub use resistance::{
    ChannelGeometry as ResistanceChannelGeometry, CombinationMethod, DarcyWeisbachModel,
    FlowConditions, HagenPoiseuilleModel, RectangularChannelModel, ResistanceCalculator,
    ResistanceModel, ResistanceModelFactory,
};

/// 1D CFD domain-specific prelude for microfluidic networks
///
/// This prelude exports 1D-specific functionality not available in the main prelude.
/// Use this when working extensively with microfluidic networks and need access to
/// specialized components, analysis tools, or scheme integration features.
///
/// For basic 1D functionality, prefer `cfd_suite::prelude::*`.
pub mod prelude {
    // === Scheme Integration (2D Visualization) ===

    // === Extended Network Components ===
    // Specialized components not in main prelude
    pub use crate::{
        analysis::{NetworkAnalysisResult, PerformanceMetrics},
        channel::{
            Channel, ChannelGeometry, ChannelType, CrossSection, FlowRegime, FlowState,
            SurfaceProperties, Wettability,
        },
        components::{
            ComponentFactory, FlowSensor, Micromixer, MixerType, PumpType, SensorType, ValveType,
        },
        network::{Edge, EdgeProperties, EdgeType, Node, NodeProperties, NodeType},
        resistance::{
            DarcyWeisbachModel, FlowConditions, ResistanceCalculator, ResistanceModelFactory,
        },
    };
}
