//! Core orchestrator for `BlueprintTopologyFactory`.
//!
//! This module provides the canonical entry point for converting declarative
//! [`BlueprintTopologySpec`]s into [`NetworkBlueprint`] graphs.  Instead of
//! maintaining its own geometry builders (SSOT violation), it delegates to the
//! canonical [`create_geometry`](crate::geometry::generator::create_geometry)
//! pipeline via [`GeometryGeneratorBuilder`].

pub use mutations::BlueprintTopologyMutation;
pub use orchestrator::BlueprintTopologyFactory;

mod build_impl;
mod mutations;
mod mutation_impl;
mod orchestrator;
mod spec_analysis_impl;

#[cfg(test)]
mod tests;
