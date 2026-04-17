//! Main geometry generation logic for 2D microfluidic channel systems.
//!
//! This module orchestrates the creation of nodes and channels using the
//! strategy pattern for channel type generation.
//!
//! # Architecture
//!
//! The `GeometryGenerator` follows the Builder pattern to incrementally
//! construct complex channel systems. It delegates channel type generation
//! to strategy objects, promoting loose coupling and extensibility.

mod entrypoints;
mod generator_impl;
mod linear;
mod selective;
pub mod shell;
mod splits;

#[cfg(test)]
mod tests;

pub use self::entrypoints::{
    create_blueprint_geometry, create_blueprint_geometry_with_metadata, create_geometry,
    create_geometry_with_metadata, GeometryGeneratorBuilder, MetadataConfig,
};
pub use self::linear::{create_parallel_geometry_from_spec, create_series_geometry_from_spec};
pub use self::selective::{
    create_primitive_selective_tree_geometry, create_primitive_selective_tree_geometry_from_spec,
    create_selective_tree_geometry, CenterSerpentinePathSpec, PrimitiveSelectiveSplitKind,
    PrimitiveSelectiveTreeRequest, SelectiveTreeRequest, SelectiveTreeTopology,
};
pub use self::shell::create_shell_cuboid;

// Private re-export: makes `GeometryGenerator` visible to submodules via `super::GeometryGenerator`.
use self::generator_impl::GeometryGenerator;

