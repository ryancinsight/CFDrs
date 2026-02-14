//! Fluidic network generation and design orchestration.
//!
//! This crate separates design/generation concerns from CFD solving concerns.
//! - `domain`: design models and invariants
//! - `application`: use-cases and ports
//! - `infrastructure`: adapters to concrete crates (`cfd-1d`)
//! - `interface`: user-facing presets/facades

pub mod application;
pub mod domain;
pub mod infrastructure;
pub mod interface;

pub use application::use_cases::NetworkGenerationService;
pub use domain::model::{ChannelSpec, EdgeKind, NetworkBlueprint, NodeKind, NodeSpec};
pub use infrastructure::adapters::Cfd1dGraphSink;
pub use interface::facade::FluidicDesigner;
pub use interface::presets::{
	serpentine_chain, symmetric_bifurcation, symmetric_trifurcation, venturi_chain,
};
