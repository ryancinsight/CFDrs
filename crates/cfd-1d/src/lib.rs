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

pub mod constants;
pub mod network;
pub mod channel;
pub mod components;
pub mod solver;
pub mod resistance;
pub mod analysis;

#[cfg(feature = "scheme-integration")]
pub mod scheme_integration;

// When scheme integration is not enabled, provide a stub module
#[cfg(not(feature = "scheme-integration"))]
pub mod scheme_integration;

// Export network functionality
pub use network::{
    Network, Node, Edge, NetworkBuilder, NetworkGraph, NodeType, EdgeType,
    BoundaryCondition, NodeProperties, EdgeProperties, NetworkMetadata,
};

// Export solver functionality
pub use solver::{
    NetworkSolver, SolverConfig, SolutionResult,
};

// Export component functionality
pub use components::{
    Component, RectangularChannel, CircularChannel, Micropump, Microvalve,
    FlowSensor, Micromixer, ComponentFactory, PumpType, ValveType, SensorType, MixerType,
};

// Export analysis functionality
pub use analysis::{
    NetworkAnalyzer, FlowAnalysis, PressureAnalysis, ResistanceAnalysis,
    PerformanceMetrics, NetworkAnalysisResult,
};

// Export channel functionality
pub use channel::{
    Channel, ChannelGeometry, ChannelType, CrossSection, SurfaceProperties,
    FlowState, FlowRegime, NumericalParameters, Wettability,
};

// Export resistance functionality
pub use resistance::{
    ResistanceModel, FlowConditions, HagenPoiseuilleModel, RectangularChannelModel,
    DarcyWeisbachModel, ResistanceModelFactory, ResistanceCalculator,
    ChannelGeometry as ResistanceChannelGeometry, CombinationMethod,
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
    // Tools for converting 1D networks to 2D schematics
    pub use crate::scheme_integration::{SchemeConversion, ComponentType, SchematicLayout};

    #[cfg(feature = "scheme-integration")]
    pub use crate::scheme_integration::{JunctionType, ChannelPathType, helpers};

    // === Advanced Network Components ===
    // Specialized components not in main prelude
    pub use crate::{
        network::{Node, Edge, NodeType, EdgeType, NodeProperties, EdgeProperties},
        components::{ComponentFactory, FlowSensor, Micromixer, PumpType, ValveType, SensorType, MixerType},
        analysis::{PerformanceMetrics, NetworkAnalysisResult},
        channel::{Channel, ChannelGeometry, ChannelType, CrossSection, SurfaceProperties, FlowState, FlowRegime, Wettability},
        resistance::{ResistanceModelFactory, ResistanceCalculator, FlowConditions, DarcyWeisbachModel},
    };
}