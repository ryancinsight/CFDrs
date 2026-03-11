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

mod linear;
mod selective;
pub mod shell;
mod splits;

pub use self::linear::{create_parallel_geometry_from_spec, create_series_geometry_from_spec};
pub use self::selective::{
    create_primitive_selective_tree_geometry, create_primitive_selective_tree_geometry_from_spec,
    create_selective_tree_geometry, CenterSerpentinePathSpec, PrimitiveSelectiveSplitKind,
    PrimitiveSelectiveTreeRequest, SelectiveTreeRequest, SelectiveTreeTopology,
};
pub use self::shell::create_shell_cuboid;

use super::builders::{ChannelBuilder, NodeBuilder};
use super::metadata::{
    BlueprintRenderHints, ChannelGeometryMetadata, GeometryAuthoringProvenance, MetadataContainer,
    OptimizationMetadata, PerformanceMetadata,
};
use super::strategies::ChannelTypeFactory;
use super::types::{polyline_length, ChannelType, Point2D, SplitType};
use crate::config::{ChannelTypeConfig, GeometryConfig};
use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::topology::{BlueprintTopologySpec, TopologyLineageMetadata};
use std::collections::HashMap;
use std::time::Instant;

/// Configuration for metadata generation
#[derive(Debug, Clone, Default)]
pub struct MetadataConfig {
    /// Whether to track performance metrics during generation
    pub track_performance: bool,
    /// Whether to track optimization metadata for serpentine channels
    pub track_optimization: bool,
    /// Optional design channel diameter (mm) used to enforce 3D-safe spacing.
    ///
    /// When set, split/merge spacing logic uses this value instead of
    /// `GeometryConfig::channel_width` for clearance calculations.
    pub channel_diameter_mm: Option<f64>,
}

impl MetadataConfig {
    /// Set a design channel diameter in millimeters for spacing calculations.
    #[must_use]
    pub const fn with_channel_diameter_mm(mut self, channel_diameter_mm: f64) -> Self {
        self.channel_diameter_mm = Some(channel_diameter_mm);
        self
    }
}

/// Internal geometry generator that builds channel systems incrementally
///
/// This struct follows the Builder pattern and uses the Strategy pattern
/// for channel type generation. It maintains state during the generation
/// process and produces a complete `NetworkBlueprint` when finalized.
///
/// Supports optional metadata tracking for performance analysis and optimization.
struct GeometryGenerator {
    box_dims: (f64, f64),
    nodes: Vec<NodeSpec>,
    channels: Vec<ChannelSpec>,
    node_counter: usize,
    channel_counter: usize,
    point_to_node_id: HashMap<(i64, i64), usize>,
    config: GeometryConfig,
    channel_type_config: ChannelTypeConfig,
    total_branches: usize,
    metadata_config: Option<MetadataConfig>,
    generation_start_time: Option<Instant>,
    render_hints: Option<BlueprintRenderHints>,
    topology: Option<BlueprintTopologySpec>,
    lineage: Option<TopologyLineageMetadata>,
    blueprint_metadata: Option<MetadataContainer>,
}

impl GeometryGenerator {
    fn channel_path(channel_type: &ChannelType) -> Vec<Point2D> {
        match channel_type {
            ChannelType::Straight => Vec::new(),
            ChannelType::SmoothStraight { path }
            | ChannelType::Serpentine { path }
            | ChannelType::Arc { path }
            | ChannelType::Frustum { path, .. } => path.clone(),
        }
    }

    fn channel_length_m(p1: Point2D, p2: Point2D, path: &[Point2D]) -> f64 {
        if path.len() >= 2 {
            polyline_length(path) * 1.0e-3
        } else {
            (p2.0 - p1.0).hypot(p2.1 - p1.1) * 1.0e-3
        }
    }

    fn estimate_local_bend_radius(path: &[Point2D], idx: usize) -> Option<f64> {
        if idx == 0 || idx + 1 >= path.len() {
            return None;
        }

        let a = path[idx - 1];
        let b = path[idx];
        let c = path[idx + 1];
        let ab = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
        let bc = ((c.0 - b.0).powi(2) + (c.1 - b.1).powi(2)).sqrt();
        let ac = ((c.0 - a.0).powi(2) + (c.1 - a.1).powi(2)).sqrt();
        if ab <= 1e-9 || bc <= 1e-9 || ac <= 1e-9 {
            return None;
        }

        let twice_area = ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)).abs();
        if twice_area <= 1e-12 {
            return Some(f64::INFINITY);
        }

        Some((ab * bc * ac) / (2.0 * twice_area))
    }

    fn infer_serpentine_shape(
        path: &[Point2D],
        p1: Point2D,
        p2: Point2D,
        channel_width: f64,
    ) -> ChannelShape {
        if path.len() < 3 {
            return ChannelShape::Serpentine {
                segments: 2,
                bend_radius_m: (channel_width * 0.5) * 1.0e-3,
            };
        }

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let length = dx.hypot(dy);
        if length <= 1e-9 {
            return ChannelShape::Serpentine {
                segments: 2,
                bend_radius_m: (channel_width * 0.5) * 1.0e-3,
            };
        }

        let nx = -dy / length;
        let ny = dx / length;
        let offsets: Vec<f64> = path
            .iter()
            .map(|point| ((point.0 - p1.0) * nx) + ((point.1 - p1.1) * ny))
            .collect();
        let extrema_threshold = channel_width * 0.15;
        let mut turns = 0usize;
        for idx in 1..offsets.len() - 1 {
            let prev_delta = offsets[idx] - offsets[idx - 1];
            let next_delta = offsets[idx + 1] - offsets[idx];
            if prev_delta.abs() <= 1e-9 || next_delta.abs() <= 1e-9 {
                continue;
            }
            if prev_delta.signum() != next_delta.signum() && offsets[idx].abs() > extrema_threshold
            {
                turns += 1;
            }
        }

        let min_radius_mm = path
            .iter()
            .enumerate()
            .filter_map(|(idx, _)| Self::estimate_local_bend_radius(path, idx))
            .filter(|radius| radius.is_finite())
            .fold(f64::INFINITY, f64::min);
        let bend_radius_mm = if min_radius_mm.is_finite() {
            min_radius_mm.max(channel_width * 0.5)
        } else {
            channel_width * 0.5
        };

        ChannelShape::Serpentine {
            segments: turns.saturating_add(1).max(2),
            bend_radius_m: bend_radius_mm * 1.0e-3,
        }
    }

    fn physical_shape_for_channel(
        channel_type: &ChannelType,
        path: &[Point2D],
        p1: Point2D,
        p2: Point2D,
        channel_width: f64,
    ) -> ChannelShape {
        match channel_type {
            ChannelType::Serpentine { .. } => {
                Self::infer_serpentine_shape(path, p1, p2, channel_width)
            }
            _ => ChannelShape::Straight,
        }
    }

    fn new(
        box_dims: (f64, f64),
        config: GeometryConfig,
        channel_type_config: ChannelTypeConfig,
        total_branches: usize,
    ) -> Self {
        Self {
            box_dims,
            nodes: Vec::new(),
            channels: Vec::new(),
            node_counter: 0,
            channel_counter: 0,
            point_to_node_id: HashMap::new(),
            config,
            channel_type_config,
            total_branches,
            metadata_config: None,
            generation_start_time: None,
            render_hints: None,
            topology: None,
            lineage: None,
            blueprint_metadata: None,
        }
    }

    fn new_with_metadata(
        box_dims: (f64, f64),
        config: GeometryConfig,
        channel_type_config: ChannelTypeConfig,
        total_branches: usize,
        metadata_config: MetadataConfig,
    ) -> Self {
        Self {
            box_dims,
            nodes: Vec::new(),
            channels: Vec::new(),
            node_counter: 0,
            channel_counter: 0,
            point_to_node_id: HashMap::new(),
            config,
            channel_type_config,
            total_branches,
            metadata_config: Some(metadata_config),
            generation_start_time: Some(Instant::now()),
            render_hints: None,
            topology: None,
            lineage: None,
            blueprint_metadata: None,
        }
    }

    /// Attach [`BlueprintRenderHints`] to the generated blueprint.
    pub(crate) fn with_render_hints(mut self, hints: BlueprintRenderHints) -> Self {
        self.render_hints = Some(hints);
        self
    }

    /// Attach a [`BlueprintTopologySpec`] to the generated blueprint.
    pub(crate) fn with_topology_spec(mut self, spec: BlueprintTopologySpec) -> Self {
        self.topology = Some(spec);
        self
    }

    /// Attach [`TopologyLineageMetadata`] to the generated blueprint.
    pub(crate) fn with_lineage(mut self, lineage: TopologyLineageMetadata) -> Self {
        self.lineage = Some(lineage);
        self
    }

    /// Attach blueprint-level [`MetadataContainer`] to the generated blueprint.
    pub(crate) fn with_blueprint_metadata(mut self, metadata: MetadataContainer) -> Self {
        self.blueprint_metadata = Some(metadata);
        self
    }

    fn point_to_key(p: Point2D) -> (i64, i64) {
        ((p.0 * 1e9) as i64, (p.1 * 1e9) as i64)
    }

    fn effective_channel_diameter(&self) -> f64 {
        self.metadata_config
            .as_ref()
            .and_then(|config| config.channel_diameter_mm)
            .filter(|diameter| diameter.is_finite() && *diameter > 0.0)
            .unwrap_or(self.config.channel_width)
    }

    fn get_or_create_node(&mut self, p: Point2D) -> usize {
        let key = Self::point_to_key(p);
        if let Some(id) = self.point_to_node_id.get(&key) {
            return *id;
        }

        let id = self.node_counter;

        // Create node with optional metadata
        let node = if let Some(ref metadata_config) = self.metadata_config {
            let mut node_builder = NodeBuilder::new(id, p);

            // Add performance metadata if enabled
            if metadata_config.track_performance {
                if let Some(start_time) = self.generation_start_time {
                    let perf_metadata = PerformanceMetadata {
                        generation_time_us: start_time.elapsed().as_micros() as u64,
                        memory_usage_bytes: std::mem::size_of::<NodeSpec>(),
                        path_points_count: 1, // Single point for node
                    };
                    node_builder = node_builder.with_metadata(perf_metadata);
                }
            }

            node_builder.build()
        } else {
            // Fast path for no metadata
            NodeSpec::new_at(format!("node_{}", id), NodeKind::Junction, p)
        };

        self.nodes.push(node);
        self.point_to_node_id.insert(key, id);
        self.node_counter += 1;
        id
    }

    fn determine_channel_type(
        &self,
        p1: Point2D,
        p2: Point2D,
        neighbor_info: Option<&[f64]>,
    ) -> ChannelType {
        let strategy =
            ChannelTypeFactory::create_strategy(&self.channel_type_config, p1, p2, self.box_dims);

        strategy.create_channel(
            p1,
            p2,
            &self.config,
            self.box_dims,
            self.total_branches,
            neighbor_info,
        )
    }

    fn add_channel_with_neighbors(
        &mut self,
        p1: Point2D,
        p2: Point2D,
        neighbor_y_coords: &[f64],
        width: Option<f64>,
    ) {
        let channel_type = self.determine_channel_type(p1, p2, Some(neighbor_y_coords));
        self.add_channel_with_type(p1, p2, Some(channel_type), width);
    }

    fn add_channel_with_type(
        &mut self,
        p1: Point2D,
        p2: Point2D,
        channel_type: Option<ChannelType>,
        width: Option<f64>,
    ) {
        let from_id = self.get_or_create_node(p1);
        let to_id = self.get_or_create_node(p2);
        let id = self.channel_counter;

        let final_channel_type =
            channel_type.unwrap_or_else(|| self.determine_channel_type(p1, p2, None));

        let channel_width = width.unwrap_or(self.config.channel_width);
        let mut channel_path = Self::channel_path(&final_channel_type);
        if channel_path.len() < 2 {
            channel_path = vec![p1, p2];
        }
        let physical_length_m = Self::channel_length_m(p1, p2, &channel_path);
        let physical_shape = Self::physical_shape_for_channel(
            &final_channel_type,
            &channel_path,
            p1,
            p2,
            channel_width,
        );

        // Create channel with optional metadata
        let channel = if let Some(ref metadata_config) = self.metadata_config {
            let mut channel_builder = ChannelBuilder::new(
                id,
                from_id,
                to_id,
                channel_width,
                self.config.channel_height,
                final_channel_type.clone(),
            );
            channel_builder = channel_builder
                .with_physical_length_m(physical_length_m)
                .with_physical_dims_m(channel_width * 1e-3, self.config.channel_height * 1e-3)
                .with_physical_shape(physical_shape);

            // Add performance metadata if enabled
            if metadata_config.track_performance {
                if let Some(start_time) = self.generation_start_time {
                    let path_points = match &final_channel_type {
                        ChannelType::Straight => 2,
                        ChannelType::SmoothStraight { path }
                        | ChannelType::Serpentine { path }
                        | ChannelType::Arc { path }
                        | ChannelType::Frustum { path, .. } => path.len(),
                    };

                    let perf_metadata = PerformanceMetadata {
                        generation_time_us: start_time.elapsed().as_micros() as u64,
                        memory_usage_bytes: std::mem::size_of::<ChannelSpec>()
                            + path_points * std::mem::size_of::<Point2D>(),
                        path_points_count: path_points,
                    };
                    channel_builder = channel_builder.with_metadata(perf_metadata);
                }
            }

            // Add optimization metadata if enabled and channel is serpentine
            if metadata_config.track_optimization {
                if let ChannelType::Serpentine { .. } = &final_channel_type {
                    let opt_metadata = OptimizationMetadata {
                        original_length: 0.0, // Will be updated by optimization system
                        optimized_length: 0.0,
                        improvement_percentage: 0.0,
                        iterations: 0,
                        optimization_time_ms: 0,
                        optimization_profile: "None".to_string(),
                    };
                    channel_builder = channel_builder.with_metadata(opt_metadata);
                }
            }

            if let Some(channel_diameter_mm) = metadata_config
                .channel_diameter_mm
                .filter(|diameter| diameter.is_finite() && *diameter > 0.0)
            {
                let geometry_metadata = ChannelGeometryMetadata {
                    channel_diameter_mm,
                };
                channel_builder = channel_builder.with_metadata(geometry_metadata);
            }

            channel_builder.build()
        } else {
            // Fast path for no metadata
            let mut channel = ChannelSpec::new_pipe_rect(
                format!("ch_{}", id),
                format!("node_{}", from_id),
                format!("node_{}", to_id),
                physical_length_m,
                channel_width * 1e-3,
                self.config.channel_height * 1e-3,
                0.0,
                0.0,
            );
            channel.path = channel_path;
            channel.channel_shape = physical_shape;
            channel
        };

        self.channels.push(channel);
        self.channel_counter += 1;
    }

    fn finalize(mut self) -> NetworkBlueprint {
        let (length, width) = self.box_dims;
        let box_outline = vec![
            ((0.0, 0.0), (length, 0.0)),
            ((length, 0.0), (length, width)),
            ((length, width), (0.0, width)),
            ((0.0, width), (0.0, 0.0)),
        ];

        let mut in_deg: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
        let mut out_deg: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
        for ch in &self.channels {
            *out_deg.entry(ch.from.as_str()).or_default() += 1;
            *in_deg.entry(ch.to.as_str()).or_default() += 1;
        }

        let mut terminals = Vec::new();
        for (idx, node) in self.nodes.iter().enumerate() {
            let deg_in = in_deg.get(node.id.as_str()).copied().unwrap_or(0);
            let deg_out = out_deg.get(node.id.as_str()).copied().unwrap_or(0);
            if deg_in + deg_out == 1 {
                terminals.push(idx);
            }
        }

        let inlet_idx = terminals
            .iter()
            .copied()
            .min_by(|&a, &b| self.nodes[a].point.0.total_cmp(&self.nodes[b].point.0));
        let outlet_idx = terminals
            .iter()
            .copied()
            .max_by(|&a, &b| self.nodes[a].point.0.total_cmp(&self.nodes[b].point.0));

        for (idx, node) in self.nodes.iter_mut().enumerate() {
            if node.kind == NodeKind::Junction {
                node.kind = if Some(idx) == inlet_idx {
                    NodeKind::Inlet
                } else if Some(idx) == outlet_idx {
                    NodeKind::Outlet
                } else {
                    NodeKind::Junction
                };
            }
        }

        let mut blueprint = NetworkBlueprint {
            name: "Generated Schematic".to_string(),
            box_dims: self.box_dims,
            nodes: self.nodes,
            channels: self.channels,
            box_outline,
            render_hints: self.render_hints,
            topology: self.topology,
            lineage: self.lineage,
            metadata: self.blueprint_metadata,
            geometry_authored: true,
        };
        blueprint.insert_metadata(GeometryAuthoringProvenance::create_geometry());
        blueprint
    }
}

/// Creates a complete 2D microfluidic channel system
///
/// This is the main entry point for generating microfluidic geometries.
/// It creates a blueprint with the specified split pattern and
/// channel types within the given bounding box.
///
/// # Arguments
///
/// * `box_dims` - Dimensions of the containing box (width, height)
/// * `splits` - Array of split types defining the branching pattern
/// * `config` - Geometry configuration (channel dimensions, clearances)
/// * `channel_type_config` - Configuration for channel type generation
///
/// # Returns
///
/// A complete `NetworkBlueprint` containing all nodes, channels, and boundary information
///
/// # Examples
///
/// ```rust
/// use cfd_schematics::{
///     geometry::{generator::create_geometry, SplitType},
///     config::{GeometryConfig, ChannelTypeConfig},
/// };
///
/// let blueprint = create_geometry(
///     (200.0, 100.0),
///     &[SplitType::Bifurcation],
///     &GeometryConfig::default(),
///     &ChannelTypeConfig::AllStraight,
/// );
/// ```
#[must_use]
pub fn create_geometry(
    box_dims: (f64, f64),
    splits: &[SplitType],
    config: &GeometryConfig,
    channel_type_config: &ChannelTypeConfig,
) -> NetworkBlueprint {
    let total_branches = splits
        .iter()
        .map(SplitType::branch_count)
        .product::<usize>()
        .max(1);
    GeometryGenerator::new(box_dims, *config, *channel_type_config, total_branches).generate(splits)
}

/// Creates a canonical split-tree [`NetworkBlueprint`] for branching geometry generation.
#[must_use]
pub fn create_blueprint_geometry(
    box_dims: (f64, f64),
    splits: &[SplitType],
    config: &GeometryConfig,
    channel_type_config: &ChannelTypeConfig,
) -> NetworkBlueprint {
    create_geometry(box_dims, splits, config, channel_type_config)
}

/// Creates a complete 2D microfluidic channel system with metadata support
///
/// This function provides the same functionality as `create_geometry` but with
/// optional metadata tracking for performance analysis and optimization.
///
/// # Arguments
///
/// * `box_dims` - Dimensions of the containing box (width, height)
/// * `splits` - Array of split types defining the branching pattern
/// * `config` - Geometry configuration (channel dimensions, clearances)
/// * `channel_type_config` - Configuration for channel type generation
/// * `metadata_config` - Configuration for metadata tracking
///
/// # Returns
///
/// A complete `NetworkBlueprint` containing all nodes, channels, boundary information,
/// and optional metadata based on the configuration.
///
/// # Examples
///
/// ```rust
/// use cfd_schematics::{
///     geometry::{generator::{create_geometry_with_metadata, MetadataConfig}, SplitType},
///     config::{GeometryConfig, ChannelTypeConfig},
/// };
///
/// let metadata_config = MetadataConfig {
///     track_performance: true,
///     track_optimization: true,
///     channel_diameter_mm: None,
/// };
///
/// let blueprint = create_geometry_with_metadata(
///     (200.0, 100.0),
///     &[SplitType::Bifurcation],
///     &GeometryConfig::default(),
///     &ChannelTypeConfig::AllStraight,
///     &metadata_config,
/// );
/// ```
#[must_use]
pub fn create_geometry_with_metadata(
    box_dims: (f64, f64),
    splits: &[SplitType],
    config: &GeometryConfig,
    channel_type_config: &ChannelTypeConfig,
    metadata_config: &MetadataConfig,
) -> NetworkBlueprint {
    let total_branches = splits
        .iter()
        .map(SplitType::branch_count)
        .product::<usize>()
        .max(1);
    GeometryGenerator::new_with_metadata(
        box_dims,
        *config,
        *channel_type_config,
        total_branches,
        metadata_config.clone(),
    )
    .generate(splits)
}

/// Creates a canonical split-tree [`NetworkBlueprint`] with metadata support.
#[must_use]
pub fn create_blueprint_geometry_with_metadata(
    box_dims: (f64, f64),
    splits: &[SplitType],
    config: &GeometryConfig,
    channel_type_config: &ChannelTypeConfig,
    metadata_config: &MetadataConfig,
) -> NetworkBlueprint {
    create_geometry_with_metadata(
        box_dims,
        splits,
        config,
        channel_type_config,
        metadata_config,
    )
}

/// Fluent builder that exposes all [`NetworkBlueprint`] capabilities when
/// generating geometry via [`create_geometry`].
///
/// Use this when you need to attach rendering hints, topology metadata,
/// lineage records, or a blueprint-level [`MetadataContainer`] to a
/// programmatically generated schematic.
///
/// # Example
///
/// ```rust
/// use cfd_schematics::geometry::{
///     generator::{GeometryGeneratorBuilder, MetadataConfig},
///     SplitType,
/// };
/// use cfd_schematics::config::{GeometryConfig, ChannelTypeConfig};
///
/// let blueprint = GeometryGeneratorBuilder::new(
///     (200.0, 100.0),
///     &[SplitType::Bifurcation],
///     &GeometryConfig::default(),
///     &ChannelTypeConfig::AllStraight,
/// )
/// .build();
/// ```
pub struct GeometryGeneratorBuilder<'s> {
    box_dims: (f64, f64),
    splits: &'s [SplitType],
    config: GeometryConfig,
    channel_type_config: ChannelTypeConfig,
    metadata_config: Option<MetadataConfig>,
    render_hints: Option<BlueprintRenderHints>,
    topology: Option<BlueprintTopologySpec>,
    lineage: Option<TopologyLineageMetadata>,
    blueprint_metadata: Option<MetadataContainer>,
}

impl<'s> GeometryGeneratorBuilder<'s> {
    /// Create a new builder with the minimum required geometry parameters.
    #[must_use]
    pub fn new(
        box_dims: (f64, f64),
        splits: &'s [SplitType],
        config: &GeometryConfig,
        channel_type_config: &ChannelTypeConfig,
    ) -> Self {
        Self {
            box_dims,
            splits,
            config: *config,
            channel_type_config: *channel_type_config,
            metadata_config: None,
            render_hints: None,
            topology: None,
            lineage: None,
            blueprint_metadata: None,
        }
    }

    /// Enable per-node / per-channel metadata tracking during generation.
    #[must_use]
    pub fn with_metadata_config(mut self, cfg: MetadataConfig) -> Self {
        self.metadata_config = Some(cfg);
        self
    }

    /// Attach [`BlueprintRenderHints`] to the produced blueprint.
    #[must_use]
    pub fn with_render_hints(mut self, hints: BlueprintRenderHints) -> Self {
        self.render_hints = Some(hints);
        self
    }

    /// Attach a [`BlueprintTopologySpec`] to the produced blueprint.
    #[must_use]
    pub fn with_topology_spec(mut self, spec: BlueprintTopologySpec) -> Self {
        self.topology = Some(spec);
        self
    }

    /// Attach [`TopologyLineageMetadata`] to the produced blueprint.
    #[must_use]
    pub fn with_lineage(mut self, lineage: TopologyLineageMetadata) -> Self {
        self.lineage = Some(lineage);
        self
    }

    /// Attach a blueprint-level [`MetadataContainer`] to the produced blueprint.
    #[must_use]
    pub fn with_blueprint_metadata(mut self, metadata: MetadataContainer) -> Self {
        self.blueprint_metadata = Some(metadata);
        self
    }

    /// Generate the [`NetworkBlueprint`], applying all configured options.
    #[must_use]
    pub fn build(self) -> NetworkBlueprint {
        let total_branches = self
            .splits
            .iter()
            .map(SplitType::branch_count)
            .product::<usize>()
            .max(1);

        let mut gen = if let Some(mc) = self.metadata_config {
            GeometryGenerator::new_with_metadata(
                self.box_dims,
                self.config,
                self.channel_type_config,
                total_branches,
                mc,
            )
        } else {
            GeometryGenerator::new(
                self.box_dims,
                self.config,
                self.channel_type_config,
                total_branches,
            )
        };

        if let Some(hints) = self.render_hints {
            gen = gen.with_render_hints(hints);
        }
        if let Some(spec) = self.topology {
            gen = gen.with_topology_spec(spec);
        }
        if let Some(lineage) = self.lineage {
            gen = gen.with_lineage(lineage);
        }
        if let Some(meta) = self.blueprint_metadata {
            gen = gen.with_blueprint_metadata(meta);
        }

        gen.generate(self.splits)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::SerpentineConfig;
    use crate::domain::model::ChannelShape;
    use crate::geometry::builders::{ChannelExt, NodeExt};
    use crate::geometry::metadata::{ChannelGeometryMetadata, PerformanceMetadata};

    #[test]
    fn test_generator_with_performance_metadata() {
        let metadata_config = MetadataConfig {
            track_performance: true,
            track_optimization: false,
            channel_diameter_mm: None,
        };

        let system = create_geometry_with_metadata(
            (100.0, 50.0),
            &[],
            &GeometryConfig::default(),
            &ChannelTypeConfig::AllStraight,
            &metadata_config,
        );

        // Check that channels have performance metadata
        for channel in &system.channels {
            assert!(channel.has_metadata::<PerformanceMetadata>());
            let perf_data = channel
                .get_metadata::<PerformanceMetadata>()
                .expect("Performance metadata should be available after creation");
            assert!(perf_data.generation_time_us > 0);
            assert!(perf_data.memory_usage_bytes > 0);
        }

        // Check that nodes have performance metadata
        for node in &system.nodes {
            assert!(node.has_metadata::<PerformanceMetadata>());
        }
    }

    #[test]
    fn test_channel_diameter_metadata_adjusts_split_spacing() {
        let box_dims = (200.0, 100.0);
        let splits = [SplitType::Bifurcation];
        let config = GeometryConfig::default();
        let channel_type_config = ChannelTypeConfig::AllStraight;

        let baseline_metadata = MetadataConfig::default();
        let baseline_system = create_geometry_with_metadata(
            box_dims,
            &splits,
            &config,
            &channel_type_config,
            &baseline_metadata,
        );

        let diameter_metadata = MetadataConfig::default().with_channel_diameter_mm(40.0);
        let diameter_system = create_geometry_with_metadata(
            box_dims,
            &splits,
            &config,
            &channel_type_config,
            &diameter_metadata,
        );

        let baseline_min_y = baseline_system
            .nodes
            .iter()
            .map(|node| node.point.1)
            .fold(f64::INFINITY, f64::min);
        let diameter_min_y = diameter_system
            .nodes
            .iter()
            .map(|node| node.point.1)
            .fold(f64::INFINITY, f64::min);

        assert!(
            diameter_min_y > baseline_min_y + 1.0,
            "Larger diameter metadata should push branches farther from walls (baseline_min_y={baseline_min_y:.3}, diameter_min_y={diameter_min_y:.3})"
        );

        for channel in &diameter_system.channels {
            let metadata = channel.get_metadata::<ChannelGeometryMetadata>().expect(
                "channel should include ChannelGeometryMetadata when diameter is configured",
            );
            assert!((metadata.channel_diameter_mm - 40.0).abs() < 1e-10);
        }
    }

    #[test]
    fn generated_serpentine_channels_persist_physical_length_and_shape() {
        let system = create_geometry(
            (200.0, 100.0),
            &[SplitType::Bifurcation],
            &GeometryConfig::default(),
            &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
        );

        let serpentine_channels: Vec<_> = system
            .channels
            .iter()
            .filter(|channel| channel.path.len() > 2)
            .collect();
        assert!(
            !serpentine_channels.is_empty(),
            "expected generated geometry to include routed serpentine channels"
        );

        for channel in serpentine_channels {
            let expected_length_m = polyline_length(&channel.path) * 1.0e-3;
            assert!(channel.length_m > 0.0);
            assert!(
                (channel.length_m - expected_length_m).abs() < 1e-9,
                "channel {:?} length should follow stored polyline length",
                channel.id
            );

            match channel.channel_shape {
                ChannelShape::Serpentine {
                    segments,
                    bend_radius_m,
                } => {
                    assert!(
                        segments >= 2,
                        "channel {:?} should expose at least one serpentine bend",
                        channel.id
                    );
                    assert!(
                        bend_radius_m > 0.0,
                        "channel {:?} should expose a positive bend radius",
                        channel.id
                    );
                }
                ref shape => panic!(
                    "channel {:?} should be marked serpentine for 1D modeling, got {:?}",
                    channel.id, shape
                ),
            }
        }
    }
}
