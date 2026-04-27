//! Visualization abstraction traits and shared rendering style types.
//!
//! This module defines traits for visualization that abstract away the specific
//! plotting library implementation. This follows the Dependency Inversion
//! Principle by allowing the visualization logic to depend on abstractions
//! rather than concrete implementations.

mod interfaces;
mod styles;

pub use interfaces::{GeometricDrawer, SchematicRenderer, VisualizationEngine};
pub use styles::{
    ChannelTypeStyles, Color, LineStyle, OutputFormat, RenderConfig, TextStyle, VisualRoleStyles,
};
