//! visualizations/mod.rs - 2D Schematic Visualization
//!
//! This module provides visualization capabilities for 2D microfluidic schematics,
//! including plotting of channel layouts, bifurcation patterns, and trifurcation designs.
//!
//! # Architecture
//!
//! The visualization module follows the Dependency Inversion Principle by defining
//! abstract traits for rendering operations and providing concrete implementations
//! for specific backends (like plotters). This allows for easy extension with
//! new rendering backends without changing the core visualization logic.
//!
//! # Modules
//!
//! - `analysis_field`: Typed CFD field overlay (pressure, shear, velocity, etc.)
//! - `traits`: Abstract interfaces for visualization operations
//! - `plotters_backend`: Concrete implementation using the plotters library
//! - `schematic`: High-level schematic rendering functions
//! - `shared_utilities`: Common utilities for visualization operations

/// Typed CFD analysis field overlay for simulation result visualization
pub mod analysis_field;
pub mod annotations;
pub mod plotters_backend;
/// High-level schematic rendering functions
pub mod schematic;
/// Shared utilities for visualization operations
pub mod shared_utilities;
pub mod traits;

pub use analysis_field::{colorize, AnalysisField, AnalysisOverlay, ColormapKind};
pub use annotations::{
    center_biased_main_path, classify_node_roles, infer_terminal_nodes_by_x,
    project_markers_along_path, should_render_label, therapy_zone_presence,
    throat_count_from_blueprint_metadata, AnnotationMarker, AnnotationStyle, LabelDensity,
    MarkerRole, SchematicAnnotations,
};
pub use plotters_backend::{
    create_plotters_renderer, plot_geometry_with_plotters, plot_shell_cuboid, PlottersRenderer,
};
pub use schematic::{
    centerline_vertices, channel_system_from_blueprint, plot_blueprint,
    plot_blueprint_auto_annotated, plot_blueprint_with_annotations, plot_geometry,
    plot_geometry_with_annotations,
};
pub use traits::{
    ChannelTypeStyles, Color, LineStyle, OutputFormat, RenderConfig, SchematicRenderer, TextStyle,
};
