//! geometry/generator.rs - Main Geometry Generation Logic
//!
//! This module contains the core geometry generation logic for creating
//! 2D microfluidic channel systems. It orchestrates the creation of nodes
//! and channels using the strategy pattern for channel type generation.
//!
//! # Architecture
//!
//! The `GeometryGenerator` follows the Builder pattern to incrementally
//! construct complex channel systems. It delegates channel type generation
//! to strategy objects, promoting loose coupling and extensibility.

use super::builders::{ChannelBuilder, NodeBuilder};
use super::metadata::{ChannelGeometryMetadata, OptimizationMetadata, PerformanceMetadata};
use super::strategies::ChannelTypeFactory;
use super::types::{Channel, ChannelSystem, ChannelType, Node, Point2D, SplitType};
use crate::config::{ChannelTypeConfig, GeometryConfig};
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
/// process and produces a complete `ChannelSystem` when finalized.
///
/// Supports optional metadata tracking for performance analysis and optimization.
struct GeometryGenerator {
    box_dims: (f64, f64),
    nodes: Vec<Node>,
    channels: Vec<Channel>,
    node_counter: usize,
    channel_counter: usize,
    point_to_node_id: HashMap<(i64, i64), usize>,
    config: GeometryConfig,
    channel_type_config: ChannelTypeConfig,
    total_branches: usize,
    metadata_config: Option<MetadataConfig>,
    generation_start_time: Option<Instant>,
}

impl GeometryGenerator {
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
        }
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
                        memory_usage_bytes: std::mem::size_of::<Node>(),
                        path_points_count: 1, // Single point for node
                    };
                    node_builder = node_builder.with_metadata(perf_metadata);
                }
            }

            node_builder.build()
        } else {
            // Fast path for no metadata
            Node {
                id,
                point: p,
                metadata: None,
            }
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

    fn add_channel_with_neighbors(&mut self, p1: Point2D, p2: Point2D, neighbor_y_coords: &[f64]) {
        let channel_type = self.determine_channel_type(p1, p2, Some(neighbor_y_coords));
        self.add_channel_with_type(p1, p2, Some(channel_type));
    }

    fn add_channel_with_type(
        &mut self,
        p1: Point2D,
        p2: Point2D,
        channel_type: Option<ChannelType>,
    ) {
        let from_id = self.get_or_create_node(p1);
        let to_id = self.get_or_create_node(p2);
        let id = self.channel_counter;

        let final_channel_type =
            channel_type.unwrap_or_else(|| self.determine_channel_type(p1, p2, None));

        // Create channel with optional metadata
        let channel = if let Some(ref metadata_config) = self.metadata_config {
            let mut channel_builder = ChannelBuilder::new(
                id,
                from_id,
                to_id,
                self.config.channel_width,
                self.config.channel_height,
                final_channel_type.clone(),
            );

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
                        memory_usage_bytes: std::mem::size_of::<Channel>()
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
            Channel {
                id,
                from_node: from_id,
                to_node: to_id,
                width: self.config.channel_width,
                height: self.config.channel_height,
                channel_type: final_channel_type,
                metadata: None,
            }
        };

        self.channels.push(channel);
        self.channel_counter += 1;
    }

    fn generate(mut self, splits: &[SplitType]) -> ChannelSystem {
        let (length, width) = self.box_dims;

        if splits.is_empty() {
            let p1 = (0.0, width / 2.0);
            let p2 = (length, width / 2.0);
            // For single channel, pass empty neighbor list so it uses box boundaries
            self.add_channel_with_neighbors(p1, p2, &[]);
            return self.finalize();
        }

        let first_half_lines = self.generate_first_half(splits);

        // Collect y-coordinates for dynamic amplitude calculation using iterator combinators
        let y_coords_for_amplitude: Vec<f64> = first_half_lines
            .iter()
            .map(|(p1, p2)| f64::midpoint(p1.1, p2.1))
            .collect();

        for (p1, p2) in &first_half_lines {
            self.add_channel_with_neighbors(*p1, *p2, &y_coords_for_amplitude);
        }

        // Generate the second half with proper merge pattern (inverse of splits)
        self.generate_second_half(splits);

        self.finalize()
    }

    fn generate_first_half(&self, splits: &[SplitType]) -> Vec<(Point2D, Point2D)> {
        let (length, width) = self.box_dims;
        let effective_width = (-2.0f64).mul_add(self.config.wall_clearance, width);
        let half_l = length / 2.0;
        let num_splits = splits.len() as u32;
        let num_segments_per_half = f64::from(num_splits).mul_add(2.0, 1.0);
        let dx = half_l / num_segments_per_half;

        let mut y_coords: Vec<f64> = vec![width / 2.0];
        let mut y_ranges: Vec<f64> = vec![effective_width];
        let mut current_x = 0.0;
        let mut lines = Vec::new();

        for split_type in splits {
            for y in &y_coords {
                lines.push(((current_x, *y), (current_x + dx, *y)));
            }
            current_x += dx;

            let (next_y_coords, next_y_ranges, new_lines) =
                self.apply_split(*split_type, &y_coords, &y_ranges, current_x, dx);

            y_coords = next_y_coords;
            y_ranges = next_y_ranges;
            lines.extend(new_lines);

            current_x += dx;
        }

        for y in &y_coords {
            lines.push(((current_x, *y), (half_l, *y)));
        }

        lines
    }

    fn apply_split(
        &self,
        split_type: SplitType,
        y_coords: &[f64],
        y_ranges: &[f64],
        current_x: f64,
        dx: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<(Point2D, Point2D)>) {
        let mut next_y_coords = Vec::new();
        let mut next_y_ranges = Vec::new();
        let mut new_lines = Vec::new();

        for (j, y_center) in y_coords.iter().enumerate() {
            let y_range = y_ranges[j];
            let n_branches = split_type.branch_count();
            let effective_channel_diameter = self.effective_channel_diameter();
            let edge_padding = self.calculate_edge_padding(y_range);
            let usable_range = (y_range - 2.0 * edge_padding).max(0.0);

            let lower_bound = y_center - y_range / 2.0 + edge_padding;
            let slot_size = if n_branches > 0 {
                usable_range / n_branches as f64
            } else {
                0.0
            };
            let child_range = (y_range / n_branches as f64).max(effective_channel_diameter);

            for i in 0..n_branches {
                let y_new = if n_branches <= 1 {
                    *y_center
                } else {
                    lower_bound + slot_size * (i as f64 + 0.5)
                };

                new_lines.push(((current_x, *y_center), (current_x + dx, y_new)));
                next_y_coords.push(y_new);
                next_y_ranges.push(child_range);
            }
        }
        (next_y_coords, next_y_ranges, new_lines)
    }

    fn apply_merge(
        &self,
        split_type: SplitType,
        y_coords: &[f64],
        y_ranges: &[f64],
        current_x: f64,
        dx: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<(Point2D, Point2D)>) {
        let mut next_y_coords = Vec::new();
        let mut next_y_ranges = Vec::new();
        let mut new_lines = Vec::new();

        let n_branches = split_type.branch_count();

        // Group the y_coords by n_branches to create merges
        for chunk in y_coords.chunks(n_branches) {
            // Calculate the center y-coordinate for this merge group
            let y_center = chunk.iter().sum::<f64>() / chunk.len() as f64;

            // Calculate the range for the merged group
            let y_range = if chunk.len() > 1 {
                let min_y = chunk.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                let max_y = chunk.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                max_y - min_y
            } else {
                y_ranges[0]
            };

            // Create merge lines from each branch to the center
            for &y_branch in chunk {
                new_lines.push(((current_x, y_branch), (current_x + dx, y_center)));
            }

            next_y_coords.push(y_center);
            next_y_ranges.push(y_range.max(self.effective_channel_diameter()));
        }

        (next_y_coords, next_y_ranges, new_lines)
    }

    /// Determine if this geometry configuration should reserve extra split spacing
    /// for curved channels (serpentine/arcs).
    fn requires_curved_spacing(&self) -> bool {
        matches!(
            self.channel_type_config,
            ChannelTypeConfig::AllSerpentine(_)
                | ChannelTypeConfig::AllArcs(_)
                | ChannelTypeConfig::MixedByPosition { .. }
                | ChannelTypeConfig::Adaptive { .. }
                | ChannelTypeConfig::SmoothSerpentineWithTransitions { .. }
        )
    }

    /// Compute edge padding inside a split range to keep branches away from
    /// local boundaries and preserve room for curved channels.
    fn calculate_edge_padding(&self, y_range: f64) -> f64 {
        let effective_channel_diameter = self.effective_channel_diameter();

        let curvature_padding = if self.requires_curved_spacing() {
            effective_channel_diameter.mul_add(1.5, self.config.wall_clearance * 0.5)
        } else {
            effective_channel_diameter * 0.5
        };

        curvature_padding.min(y_range * 0.2)
    }

    fn generate_second_half(&mut self, splits: &[SplitType]) {
        let (length, width) = self.box_dims;
        let effective_width = (-2.0f64).mul_add(self.config.wall_clearance, width);
        let half_l = length / 2.0;
        let num_splits = splits.len() as u32;
        let num_segments_per_half = f64::from(num_splits).mul_add(2.0, 1.0);
        let dx = half_l / num_segments_per_half;

        // Calculate the final y-coordinates at the center (end of first half)
        // Start with the initial state
        let mut y_coords = vec![width / 2.0];
        let mut y_ranges = vec![effective_width];

        // Apply all splits to get the final state
        for split_type in splits {
            let (next_y_coords, next_y_ranges, _) =
                self.apply_split(*split_type, &y_coords, &y_ranges, 0.0, dx);
            y_coords = next_y_coords;
            y_ranges = next_y_ranges;
        }

        // Now generate the second half by reversing the splits (creating merges)
        let mut current_x = half_l;
        let mut lines = Vec::new();

        // Process splits in reverse order to create merges
        for split_type in splits.iter().rev() {
            // Add horizontal segments from current position
            for y in &y_coords {
                lines.push(((current_x, *y), (current_x + dx, *y)));
            }
            current_x += dx;

            // Apply merge (reverse of split)
            let (next_y_coords, next_y_ranges, new_lines) =
                self.apply_merge(*split_type, &y_coords, &y_ranges, current_x, dx);

            y_coords = next_y_coords;
            y_ranges = next_y_ranges;
            lines.extend(new_lines);

            current_x += dx;
        }

        // Final horizontal segments to the right edge
        for y in &y_coords {
            lines.push(((current_x, *y), (length, *y)));
        }

        // Collect y-coordinates for amplitude calculation
        let mut y_coords_for_amplitude: Vec<f64> = Vec::with_capacity(lines.len());
        for (p1, p2) in &lines {
            y_coords_for_amplitude.push(f64::midpoint(p1.1, p2.1));
        }

        // Add all the second half channels
        for (p1, p2) in &lines {
            self.add_channel_with_neighbors(*p1, *p2, &y_coords_for_amplitude);
        }
    }

    fn finalize(self) -> ChannelSystem {
        let (length, width) = self.box_dims;
        let box_outline = vec![
            ((0.0, 0.0), (length, 0.0)),
            ((length, 0.0), (length, width)),
            ((length, width), (0.0, width)),
            ((0.0, width), (0.0, 0.0)),
        ];
        ChannelSystem {
            box_dims: self.box_dims,
            nodes: self.nodes,
            channels: self.channels,
            box_outline,
        }
    }
}

/// Creates a complete 2D microfluidic channel system
///
/// This is the main entry point for generating microfluidic geometries.
/// It creates a channel system with the specified split pattern and
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
/// A complete `ChannelSystem` containing all nodes, channels, and boundary information
///
/// # Examples
///
/// ```rust
/// use scheme::{
///     geometry::{generator::create_geometry, SplitType},
///     config::{GeometryConfig, ChannelTypeConfig},
/// };
///
/// let system = create_geometry(
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
) -> ChannelSystem {
    let total_branches = splits
        .iter()
        .map(SplitType::branch_count)
        .product::<usize>()
        .max(1);
    GeometryGenerator::new(box_dims, *config, *channel_type_config, total_branches).generate(splits)
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
/// A complete `ChannelSystem` containing all nodes, channels, boundary information,
/// and optional metadata based on the configuration.
///
/// # Examples
///
/// ```rust
/// use scheme::{
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
/// let system = create_geometry_with_metadata(
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
) -> ChannelSystem {
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

#[cfg(test)]
mod tests {
    use super::*;
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
}
