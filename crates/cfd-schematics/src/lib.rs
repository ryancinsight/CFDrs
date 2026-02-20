//! # cfd-schematics — Design-Time Topology & Visualization Layer
//!
//! `cfd-schematics` is the **single authoritative source** for microfluidic
//! network topology, geometry generation, and 2D schematic visualization in
//! the CFDrs workspace. It is a *design-time* crate — it does not solve PDEs.
//!
//! ## Role in the CFDrs Architecture
//!
//! ```text
//!                  ┌─────────────────────────────────────┐
//!                  │          cfd-schematics              │
//!                  │  (topology · geometry · 2D render)  │
//!                  └──────────────┬──────────────────────┘
//!                                 │ NodeSpec / ChannelSpec
//!                    ┌────────────┴────────────┐
//!                    ▼                         ▼
//!             ┌─────────────┐          ┌─────────────┐
//!             │   cfd-1d    │          │   cfd-2d    │
//!             │ (1D network │          │ (2D PDE /   │
//!             │  solver)    │          │  FVM/FDM)   │
//!             └─────────────┘          └─────────────┘
//! ```
//!
//! ### What "2D schematic" means here
//!
//! The 2D coordinate system in `cfd-schematics` is a **layout plane** — it
//! describes where channels and nodes are positioned on a chip diagram (x, y
//! in millimetres). This is **not** the same as the 2D spatial domain solved
//! by `cfd-2d`. The schematic is always flat (z = 0); it is a *topology map*,
//! not a simulation domain.
//!
//! ### Relationship with `cfd-1d`
//!
//! `cfd-1d` consumes `NodeSpec` and `ChannelSpec` from this crate to build its
//! resistance-network graph. The 1D solver treats each channel as a lumped
//! element (Hagen-Poiseuille resistance), collapsing the 3D cross-section into
//! a single scalar per edge. `cfd-schematics` provides:
//! - `NodeSpec` / `ChannelSpec` — topology and cross-section geometry
//! - `CrossSectionSpec` — circular or rectangular geometry for resistance calc
//! - `AnalysisOverlay` — typed CFD field visualization on the schematic
//!
//! ### Relationship with `cfd-2d`
//!
//! `cfd-2d` solves the **full 2D Navier-Stokes equations** on a structured or
//! unstructured grid (FDM / FVM / LBM). It does not currently depend on
//! `cfd-schematics`, because its domain is a continuous PDE field (u, v, p at
//! every grid cell), not a lumped network graph.
//!
//! However, `cfd-schematics` *can* serve `cfd-2d` in two ways:
//! 1. **Geometry seeding**: `ChannelSystem` geometry can be used to define the
//!    bounding box and channel centrelines that seed a 2D structured grid.
//! 2. **Result overlay**: `AnalysisOverlay` can visualize 2D solver output
//!    (e.g. cross-section-averaged pressure or velocity) projected back onto
//!    the schematic for design-level inspection.
//!
//! # Modules
//! - **config** / **config_constants**: Channel and geometry configuration types
//! - **error**: Domain error types
//! - **geometry**: Core geometric types, generation, and optimisation logic
//! - **state_management**: Parameter management and bilateral-symmetry helpers
//! - **visualizations**: 2-D schematic rendering, SVG/PNG export, `AnalysisOverlay`
//! - **domain**: `NodeSpec`, `ChannelSpec`, `CrossSectionSpec`, `NetworkBlueprint`
//! - **application**: `NetworkGenerationService` use-case
//! - **infrastructure**: `PetgraphGraphSink`, `DesignGraph` adapters
//! - **interface**: Preset factories (bifurcation, trifurcation, serpentine, venturi)


// ── Geometry & visualisation layer (inlined from former `scheme` crate) ──────
pub mod config;
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
