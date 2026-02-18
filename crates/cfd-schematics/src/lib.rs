//! cfd-schematics - Design-time schematic generation for CFDrs
//!
//! This crate is the single authoritative source for all geometry generation,
//! network design, and schematic exports in the CFDrs workspace.
//!
//! # Modules
//! - **config** / **config_constants**: Channel and geometry configuration types
//! - **error**: Domain error types
//! - **geometry**: Core geometric types, generation, and optimisation logic
//! - **state_management**: Parameter management and bilateral-symmetry helpers
//! - **visualizations**: 2-D schematic rendering and SVG/PNG export
//! - **network**: CFDrs topology abstraction, graph builder, and preset factories

// ── Geometry & visualisation layer (inlined from former `scheme` crate) ──────
pub mod config;
pub mod config_constants;
pub mod error;
pub mod geometry;
pub mod state_management;
pub mod visualizations;

// ── CFDrs network / topology layer ───────────────────────────────────────────
pub mod application;
pub mod domain;
pub mod infrastructure;
pub mod interface;

// ── Flat convenience re-exports ───────────────────────────────────────────────
pub use application::use_cases::NetworkGenerationService;
pub use domain::model::{ChannelSpec, EdgeKind, NetworkBlueprint, NodeKind, NodeSpec};
pub use error::{
    ConfigurationError, GeometryError, SchemeError, SchemeResult, StrategyError,
    VisualizationError,
};
pub use infrastructure::adapters::{build_design_graph, DesignGraph, PetgraphGraphSink};
pub use interface::presets::{
    serpentine_chain, symmetric_bifurcation, symmetric_trifurcation, venturi_chain,
};
pub use state_management::{
    ConfigurableParameter, ConstraintError, ParameterConstraints, ParameterError,
    ParameterManager, ParameterRegistry, StateManagementError, StateManagementResult,
};
pub use visualizations::schematic::plot_geometry;

/// Unified network / design API — the single import point for downstream crates.
pub mod network {
    pub use crate::application::ports::GraphSink;
    pub use crate::application::use_cases::NetworkGenerationService;
    pub use crate::domain::model::{
        ChannelSpec, EdgeId, EdgeKind, NetworkBlueprint, NodeId, NodeKind, NodeSpec,
    };
    pub use crate::domain::rules::BlueprintValidator;
    pub use crate::infrastructure::adapters::{
        build_design_graph, DesignGraph, PetgraphGraphSink,
    };
    pub use crate::interface::presets::{
        serpentine_chain, symmetric_bifurcation, symmetric_trifurcation, venturi_chain,
    };
}
