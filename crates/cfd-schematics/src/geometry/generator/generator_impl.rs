//! Internal `GeometryGenerator` builder — constructs channel systems incrementally.
//!
//! This module is private to the `generator` module family. Public entry-points
//! live in `entrypoints.rs`; split-sequence generation lives in `splits.rs`.

use super::super::builders::{ChannelBuilder, NodeBuilder};
use super::super::metadata::{
    BlueprintRenderHints, ChannelGeometryMetadata, GeometryAuthoringProvenance, MetadataContainer,
    OptimizationMetadata, PerformanceMetadata,
};
use super::super::strategies::ChannelTypeFactory;
use super::super::types::{polyline_length, ChannelType, Point2D};
use crate::config::{ChannelTypeConfig, GeometryConfig};
use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::topology::{BlueprintTopologySpec, TopologyLineageMetadata};
use std::collections::HashMap;
use std::time::Instant;

// `MetadataConfig` is re-exported from the parent module (entrypoints.rs → mod.rs → here).
use super::MetadataConfig;

/// Internal geometry generator that builds channel systems incrementally.
///
/// Follows the Builder pattern and delegates channel-type creation to strategy
/// objects (see `ChannelTypeFactory`).  Produces a `NetworkBlueprint` via
/// [`finalize`](Self::finalize).
///
/// Supports optional metadata tracking for performance analysis and optimization.
pub(super) struct GeometryGenerator {
    pub(super) box_dims: (f64, f64),
    pub(super) nodes: Vec<NodeSpec>,
    pub(super) channels: Vec<ChannelSpec>,
    pub(super) node_counter: usize,
    pub(super) channel_counter: usize,
    pub(super) point_to_node_id: HashMap<(i64, i64), usize>,
    pub(super) config: GeometryConfig,
    pub(super) channel_type_config: ChannelTypeConfig,
    pub(super) total_branches: usize,
    pub(super) metadata_config: Option<MetadataConfig>,
    pub(super) generation_start_time: Option<Instant>,
    pub(super) render_hints: Option<BlueprintRenderHints>,
    pub(super) topology: Option<BlueprintTopologySpec>,
    pub(super) lineage: Option<TopologyLineageMetadata>,
    pub(super) blueprint_metadata: Option<MetadataContainer>,
}

impl GeometryGenerator {
    pub(super) fn channel_path(channel_type: &ChannelType) -> Vec<Point2D> {
        match channel_type {
            ChannelType::Straight => Vec::new(),
            ChannelType::SmoothStraight { path }
            | ChannelType::Serpentine { path }
            | ChannelType::Arc { path }
            | ChannelType::Frustum { path, .. } => path.clone(),
        }
    }

    pub(super) fn channel_length_m(p1: Point2D, p2: Point2D, path: &[Point2D]) -> f64 {
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
                wave_type: crate::topology::SerpentineWaveType::Sine,
            };
        }

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let length = dx.hypot(dy);
        if length <= 1e-9 {
            return ChannelShape::Serpentine {
                segments: 2,
                bend_radius_m: (channel_width * 0.5) * 1.0e-3,
                wave_type: crate::topology::SerpentineWaveType::Sine,
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
            wave_type: crate::topology::SerpentineWaveType::Sine,
        }
    }

    pub(super) fn physical_shape_for_channel(
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

    pub(super) fn new(
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

    pub(super) fn new_with_metadata(
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

    pub(super) fn point_to_key(p: Point2D) -> (i64, i64) {
        ((p.0 * 1e9) as i64, (p.1 * 1e9) as i64)
    }

    pub(super) fn effective_channel_diameter(&self) -> f64 {
        self.metadata_config
            .as_ref()
            .and_then(|config| config.channel_diameter_mm)
            .filter(|diameter| diameter.is_finite() && *diameter > 0.0)
            .unwrap_or(self.config.channel_width)
    }

    pub(super) fn get_or_create_node(&mut self, p: Point2D) -> usize {
        let key = Self::point_to_key(p);
        if let Some(id) = self.point_to_node_id.get(&key) {
            return *id;
        }

        let id = self.node_counter;

        let node = if let Some(ref metadata_config) = self.metadata_config {
            let mut node_builder = NodeBuilder::new(id, p);

            if metadata_config.track_performance {
                if let Some(start_time) = self.generation_start_time {
                    let perf_metadata = PerformanceMetadata {
                        generation_time_us: start_time.elapsed().as_micros() as u64,
                        memory_usage_bytes: std::mem::size_of::<NodeSpec>(),
                        path_points_count: 1,
                    };
                    node_builder = node_builder.with_metadata(perf_metadata);
                }
            }

            node_builder.build()
        } else {
            NodeSpec::new_at(format!("node_{id}"), NodeKind::Junction, p)
        };

        self.nodes.push(node);
        self.point_to_node_id.insert(key, id);
        self.node_counter += 1;
        id
    }

    pub(super) fn determine_channel_type(
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

    pub(super) fn add_channel_with_neighbors(
        &mut self,
        p1: Point2D,
        p2: Point2D,
        neighbor_y_coords: &[f64],
        width: Option<f64>,
    ) {
        let channel_type = self.determine_channel_type(p1, p2, Some(neighbor_y_coords));
        self.add_channel_with_type(p1, p2, Some(channel_type), width);
    }

    pub(super) fn add_channel_with_type(
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

            if metadata_config.track_optimization {
                if let ChannelType::Serpentine { .. } = &final_channel_type {
                    let opt_metadata = OptimizationMetadata {
                        original_length: 0.0,
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
            let mut channel = ChannelSpec::new_pipe_rect(
                format!("ch_{id}"),
                format!("node_{from_id}"),
                format!("node_{to_id}"),
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

    pub(super) fn finalize(mut self) -> NetworkBlueprint {
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
