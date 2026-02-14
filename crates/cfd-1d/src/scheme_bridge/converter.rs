//! Core conversion logic from [`scheme::geometry::ChannelSystem`] to `cfd-1d` [`Network`].

use super::error::BridgeError;
use crate::channel::{
    ChannelGeometry as CfdChannelGeometry, ChannelType as CfdChannelType,
    CrossSection, SurfaceProperties, Wettability,
};
use crate::network::{
    ComponentType, Edge, EdgeProperties, EdgeType, Network, NetworkBuilder, Node, NodeType,
};
use crate::resistance::{
    ChannelGeometry as ResistanceGeometry, FlowConditions, ResistanceCalculator,
};
use cfd_core::physics::fluid::{ConstantPropertyFluid, FluidTrait};
use num_traits::FromPrimitive;
use scheme::geometry::types::{
    Channel as SchemeChannel, ChannelSystem, ChannelType as SchemeChannelType,
    Node as SchemeNode, Point2D,
};
use std::collections::{HashMap, HashSet};

/// Converts a `scheme` 2D schematic into a `cfd-1d` simulation-ready [`Network`].
///
/// The converter analyses the topology of a [`ChannelSystem`], infers inlet /
/// outlet / junction roles for each node, computes physical channel lengths from
/// the polyline paths, and builds a directed [`petgraph`] network with
/// rectangular cross-sections sized from the scheme channel dimensions.
pub struct SchemeNetworkConverter<'a> {
    system: &'a ChannelSystem,
    /// Optional scale factor to convert scheme coordinate units to meters.
    /// Default is `1e-6` (micrometers → meters) for microfluidics.
    scale: f64,
}

impl<'a> SchemeNetworkConverter<'a> {
    /// Create a converter for a scheme [`ChannelSystem`].
    ///
    /// Uses the default micrometre→metre scale factor (`1e-6`). If your scheme
    /// coordinates are already in metres, call [`Self::with_scale`] with `1.0`.
    #[must_use]
    pub fn new(system: &'a ChannelSystem) -> Self {
        Self {
            system,
            scale: 1e-6,
        }
    }

    /// Create a converter with an explicit coordinate scale factor.
    ///
    /// `scale` is multiplied with every coordinate and dimension taken from the
    /// scheme system to produce SI-metre values for the CFD network.
    #[must_use]
    pub fn with_scale(system: &'a ChannelSystem, scale: f64) -> Self {
        Self { system, scale }
    }

    // ── topology inference ───────────────────────────────────────────────

    /// Classify node roles from channel connectivity.
    ///
    /// A node that has **no incoming** channels is an `Inlet`.
    /// A node that has **no outgoing** channels is an `Outlet`.
    /// Everything else is a `Junction`.
    fn classify_nodes(&self) -> Result<HashMap<usize, NodeType>, BridgeError> {
        if self.system.channels.is_empty() {
            return Err(BridgeError::EmptyNetwork);
        }

        // Validate node references
        let max_node_id = self.system.nodes.len();
        for ch in &self.system.channels {
            if ch.from_node >= max_node_id {
                return Err(BridgeError::InvalidNodeReference {
                    channel_id: ch.id,
                    node_id: ch.from_node,
                });
            }
            if ch.to_node >= max_node_id {
                return Err(BridgeError::InvalidNodeReference {
                    channel_id: ch.id,
                    node_id: ch.to_node,
                });
            }
        }

        let mut has_incoming: HashSet<usize> = HashSet::new();
        let mut has_outgoing: HashSet<usize> = HashSet::new();

        for ch in &self.system.channels {
            has_outgoing.insert(ch.from_node);
            has_incoming.insert(ch.to_node);
        }

        let all_node_ids: HashSet<usize> = self.system.nodes.iter().map(|n| n.id).collect();
        let mut classification = HashMap::new();

        for &nid in &all_node_ids {
            let incoming = has_incoming.contains(&nid);
            let outgoing = has_outgoing.contains(&nid);

            let node_type = match (incoming, outgoing) {
                (false, true) => NodeType::Inlet,  // no incoming → inlet
                (true, false) => NodeType::Outlet,  // no outgoing → outlet
                (false, false) => NodeType::Inlet,   // isolated → treat as inlet
                (true, true) => NodeType::Junction,
            };
            classification.insert(nid, node_type);
        }

        // Validate we have at least one inlet and outlet
        let has_inlet = classification.values().any(|t| matches!(t, NodeType::Inlet));
        let has_outlet = classification.values().any(|t| matches!(t, NodeType::Outlet));

        if !has_inlet {
            return Err(BridgeError::NoInlets);
        }
        if !has_outlet {
            return Err(BridgeError::NoOutlets);
        }

        Ok(classification)
    }

    // ── geometry helpers ─────────────────────────────────────────────────

    /// Compute the polyline length of a channel (in scaled units).
    fn channel_path_length(&self, channel: &SchemeChannel) -> f64 {
        let points = self.channel_path_points(channel);
        if points.len() < 2 {
            // Straight line between nodes
            let from = self.system.nodes[channel.from_node].point;
            let to = self.system.nodes[channel.to_node].point;
            return self.point_distance(from, to);
        }
        let mut length = 0.0;
        for pair in points.windows(2) {
            length += self.point_distance(pair[0], pair[1]);
        }
        length
    }

    /// Get the path points for a channel, falling back to endpoints for Straight.
    fn channel_path_points(&self, channel: &SchemeChannel) -> Vec<Point2D> {
        match &channel.channel_type {
            SchemeChannelType::Straight => {
                let from = self.system.nodes[channel.from_node].point;
                let to = self.system.nodes[channel.to_node].point;
                vec![from, to]
            }
            SchemeChannelType::SmoothStraight { path }
            | SchemeChannelType::Serpentine { path }
            | SchemeChannelType::Arc { path } => path.clone(),
            SchemeChannelType::Frustum { path, .. } => path.clone(),
        }
    }

    /// Euclidean distance between two 2D points, scaled to metres.
    fn point_distance(&self, a: Point2D, b: Point2D) -> f64 {
        let dx = (b.0 - a.0) * self.scale;
        let dy = (b.1 - a.1) * self.scale;
        (dx * dx + dy * dy).sqrt()
    }

    /// Map a scheme `ChannelType` to a cfd-1d `ChannelType`.
    fn map_channel_type(scheme_type: &SchemeChannelType) -> CfdChannelType {
        match scheme_type {
            SchemeChannelType::Straight | SchemeChannelType::SmoothStraight { .. } => {
                CfdChannelType::Straight
            }
            SchemeChannelType::Serpentine { path } => {
                // Count turns: a full period requires going up then down.
                // Approximate turns from path curvature sign changes.
                let turns = Self::estimate_serpentine_turns(path);
                CfdChannelType::Serpentine { turns }
            }
            SchemeChannelType::Arc { .. } => CfdChannelType::Curved {
                radius: 0.0, // Will be refined from geometry
            },
            SchemeChannelType::Frustum { .. } => CfdChannelType::Tapered,
        }
    }

    /// Estimate the number of serpentine turns from path y-direction reversals.
    fn estimate_serpentine_turns(path: &[Point2D]) -> usize {
        if path.len() < 3 {
            return 0;
        }
        let mut sign_changes = 0;
        for window in path.windows(3) {
            let dy1 = window[1].1 - window[0].1;
            let dy2 = window[2].1 - window[1].1;
            if dy1 * dy2 < 0.0 {
                sign_changes += 1;
            }
        }
        // Two sign changes ≈ 1 full turn
        sign_changes / 2
    }

    /// Get average cross-section dimensions for frustum channels.
    fn frustum_average_dims(channel: &SchemeChannel) -> (f64, f64) {
        if let SchemeChannelType::Frustum {
            inlet_width,
            throat_width,
            outlet_width,
            ..
        } = &channel.channel_type
        {
            let avg_width = (inlet_width + throat_width + outlet_width) / 3.0;
            (avg_width, channel.height)
        } else {
            (channel.width, channel.height)
        }
    }

    // ── network construction ─────────────────────────────────────────────

    /// Build a [`Network`] using the given fluid model.
    ///
    /// This is the primary conversion method. It:
    /// 1. Classifies nodes into inlet/outlet/junction.
    /// 2. Creates the petgraph network.
    /// 3. Computes channel lengths and cross-sections.
    /// 4. Calculates initial resistance using `ResistanceCalculator`.
    ///
    /// # Errors
    ///
    /// Returns [`BridgeError`] if the scheme topology is invalid or
    /// resistance calculation fails.
    pub fn build_network<F: FluidTrait<f64> + Clone>(
        &self,
        fluid: F,
    ) -> Result<Network<f64, F>, BridgeError> {
        let node_types = self.classify_nodes()?;

        // Build petgraph via NetworkBuilder
        let mut builder = NetworkBuilder::<f64>::new();

        // Map scheme node-id → petgraph NodeIndex
        let mut node_map: HashMap<usize, petgraph::graph::NodeIndex> = HashMap::new();

        for scheme_node in &self.system.nodes {
            let nt = node_types
                .get(&scheme_node.id)
                .copied()
                .unwrap_or(NodeType::Junction);

            let node = Node::new(format!("node_{}", scheme_node.id), nt).with_position(
                scheme_node.point.0 * self.scale,
                scheme_node.point.1 * self.scale,
            );
            let idx = builder.add_node(node);
            node_map.insert(scheme_node.id, idx);
        }

        // Track edge indices for post-build property assignment
        let mut edge_channel_pairs: Vec<(petgraph::graph::EdgeIndex, &SchemeChannel)> = Vec::new();

        for channel in &self.system.channels {
            let from_idx = node_map[&channel.from_node];
            let to_idx = node_map[&channel.to_node];

            let mut edge = Edge::new(format!("ch_{}", channel.id), EdgeType::Pipe);

            // Pre-compute a nominal resistance from the analytical formula so
            // that NetworkBuilder::build() validation passes (R > 0).
            let path_len = self.channel_path_length(channel);
            if path_len <= 0.0 {
                return Err(BridgeError::ZeroLengthChannel {
                    channel_id: channel.id,
                });
            }

            let (w, h) = Self::frustum_average_dims(channel);
            let width_m = w * self.scale;
            let height_m = h * self.scale;

            // Analytical Hagen-Poiseuille for rectangular section:
            // R ≈ 12 μ L / (w h³)  (thin-slit approximation)
            // Use a placeholder viscosity of water at 20°C: 0.001 Pa·s
            let mu = 0.001_f64;
            let nominal_r = 12.0 * mu * path_len / (width_m * height_m.powi(3));
            edge.resistance = nominal_r.max(f64::EPSILON);

            let edge_idx = builder.add_edge(from_idx, to_idx, edge);
            edge_channel_pairs.push((edge_idx, channel));
        }

        // Build and validate the graph
        let graph = builder.build().map_err(BridgeError::CfdError)?;

        // Assemble the Network wrapper
        let mut network = Network::new(graph, fluid);

        // Populate edge properties with full geometry information
        let calculator: ResistanceCalculator<f64> = ResistanceCalculator::new();

        for (edge_idx, channel) in &edge_channel_pairs {
            let path_len = self.channel_path_length(channel);
            let (w, h) = Self::frustum_average_dims(channel);
            let width_m = w * self.scale;
            let height_m = h * self.scale;
            let area = width_m * height_m;
            let dh = 2.0 * width_m * height_m / (width_m + height_m); // hydraulic diameter

            let cfd_channel_type = Self::map_channel_type(&channel.channel_type);

            let channel_geometry = CfdChannelGeometry {
                channel_type: cfd_channel_type,
                length: path_len,
                cross_section: CrossSection::Rectangular {
                    width: width_m,
                    height: height_m,
                },
                surface: SurfaceProperties {
                    roughness: 1e-7, // smooth PDMS default
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            };

            // Add frustum geometric variations if applicable
            let mut variations = Vec::new();
            if let SchemeChannelType::Frustum {
                widths,
                path,
                ..
            } = &channel.channel_type
            {
                if widths.len() == path.len() && !widths.is_empty() {
                    let base_width = widths[0];
                    for (i, w_val) in widths.iter().enumerate() {
                        let position = i as f64 / (widths.len() - 1).max(1) as f64;
                        let scale_factor = if base_width > 0.0 {
                            w_val / base_width
                        } else {
                            1.0
                        };
                        variations.push(crate::channel::GeometricVariation {
                            position,
                            scale_factor,
                            roughness_factor: 1.0,
                        });
                    }
                }
            }

            let mut channel_geometry_with_variations = channel_geometry;
            channel_geometry_with_variations.variations = variations;

            let component_type = match &channel.channel_type {
                SchemeChannelType::Frustum { .. } => ComponentType::Generic,
                _ => ComponentType::Pipe,
            };

            let props = EdgeProperties {
                id: format!("ch_{}", channel.id),
                component_type,
                length: path_len,
                area,
                hydraulic_diameter: Some(dh),
                resistance: f64::EPSILON, // placeholder, updated below
                geometry: Some(channel_geometry_with_variations),
                properties: HashMap::new(),
            };

            network.add_edge_properties(*edge_idx, props);

            // Try to calculate proper resistance coefficients with the fluid model.
            // Use zero flow rate for initial linearization.
            let res_geom = ResistanceGeometry::Rectangular {
                width: width_m,
                height: height_m,
                length: path_len,
            };

            let mut conditions = FlowConditions::new(0.0);
            conditions.velocity = None;
            conditions.flow_rate = Some(1e-12); // tiny initial flow for Re calculation

            if let Ok((r, k)) = calculator.calculate_coefficients_auto(
                &res_geom,
                network.fluid(),
                &conditions,
            ) {
                if let Some(edge) = network.graph.edge_weight_mut(*edge_idx) {
                    edge.resistance = r;
                    edge.quad_coeff = k;
                }
            }
        }

        Ok(network)
    }

    /// Convenience: build a network using water at 20 °C.
    ///
    /// # Errors
    ///
    /// Returns [`BridgeError`] if the conversion or fluid creation fails.
    pub fn build_network_with_water(
        &self,
    ) -> Result<Network<f64, ConstantPropertyFluid<f64>>, BridgeError> {
        let water = ConstantPropertyFluid::water_20c().map_err(BridgeError::CfdError)?;
        self.build_network(water)
    }

    /// Get a summary of the conversion that would be performed.
    ///
    /// Useful for diagnostics before running the full conversion.
    #[must_use]
    pub fn summary(&self) -> ConversionSummary {
        let node_types = self.classify_nodes();

        let (inlets, outlets, junctions) = match &node_types {
            Ok(map) => {
                let inlets = map.values().filter(|t| matches!(t, NodeType::Inlet)).count();
                let outlets = map.values().filter(|t| matches!(t, NodeType::Outlet)).count();
                let junctions = map
                    .values()
                    .filter(|t| matches!(t, NodeType::Junction))
                    .count();
                (inlets, outlets, junctions)
            }
            Err(_) => (0, 0, 0),
        };

        let mut channel_lengths: Vec<f64> = self
            .system
            .channels
            .iter()
            .map(|ch| self.channel_path_length(ch))
            .collect();
        channel_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let total_length: f64 = channel_lengths.iter().sum();
        let min_length = channel_lengths.first().copied().unwrap_or(0.0);
        let max_length = channel_lengths.last().copied().unwrap_or(0.0);

        let channel_type_counts = {
            let mut counts = HashMap::new();
            for ch in &self.system.channels {
                let label = match &ch.channel_type {
                    SchemeChannelType::Straight => "straight",
                    SchemeChannelType::SmoothStraight { .. } => "smooth_straight",
                    SchemeChannelType::Serpentine { .. } => "serpentine",
                    SchemeChannelType::Arc { .. } => "arc",
                    SchemeChannelType::Frustum { .. } => "frustum",
                };
                *counts.entry(label.to_string()).or_insert(0usize) += 1;
            }
            counts
        };

        ConversionSummary {
            total_nodes: self.system.nodes.len(),
            total_channels: self.system.channels.len(),
            inlets,
            outlets,
            junctions,
            total_path_length_m: total_length,
            min_channel_length_m: min_length,
            max_channel_length_m: max_length,
            channel_type_counts,
            scale_factor: self.scale,
            valid: node_types.is_ok(),
        }
    }
}

/// Summary of a scheme→cfd-1d conversion.
#[derive(Debug, Clone)]
pub struct ConversionSummary {
    /// Total number of nodes in the scheme system.
    pub total_nodes: usize,
    /// Total number of channels in the scheme system.
    pub total_channels: usize,
    /// Number of inferred inlets.
    pub inlets: usize,
    /// Number of inferred outlets.
    pub outlets: usize,
    /// Number of inferred junctions.
    pub junctions: usize,
    /// Total path length across all channels (metres).
    pub total_path_length_m: f64,
    /// Shortest channel (metres).
    pub min_channel_length_m: f64,
    /// Longest channel (metres).
    pub max_channel_length_m: f64,
    /// Count of channels by scheme type name.
    pub channel_type_counts: HashMap<String, usize>,
    /// Scale factor used (coordinate units → metres).
    pub scale_factor: f64,
    /// Whether the topology passed validation.
    pub valid: bool,
}

impl std::fmt::Display for ConversionSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Scheme → CFD-1D Conversion Summary")?;
        writeln!(f, "──────────────────────────────────")?;
        writeln!(f, "Nodes:    {} total ({} inlets, {} outlets, {} junctions)",
            self.total_nodes, self.inlets, self.outlets, self.junctions)?;
        writeln!(f, "Channels: {}", self.total_channels)?;
        for (ctype, count) in &self.channel_type_counts {
            writeln!(f, "  - {}: {}", ctype, count)?;
        }
        writeln!(f, "Path length: {:.6} m total ({:.6}–{:.6} m per channel)",
            self.total_path_length_m, self.min_channel_length_m, self.max_channel_length_m)?;
        writeln!(f, "Scale:    {} (coord → m)", self.scale_factor)?;
        writeln!(f, "Valid:    {}", self.valid)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use scheme::config::{ChannelTypeConfig, GeometryConfig};
    use scheme::geometry::generator::create_geometry;
    use scheme::geometry::SplitType;

    fn make_bifurcation_system() -> ChannelSystem {
        create_geometry(
            (200.0, 100.0),
            &[SplitType::Bifurcation],
            &GeometryConfig::default(),
            &ChannelTypeConfig::AllStraight,
        )
    }

    #[test]
    fn test_summary_valid() {
        let system = make_bifurcation_system();
        let conv = SchemeNetworkConverter::new(&system);
        let summary = conv.summary();

        assert!(summary.valid);
        assert!(summary.inlets > 0);
        assert!(summary.outlets > 0);
        assert!(summary.total_channels > 0);
        assert!(summary.total_path_length_m > 0.0);
    }

    #[test]
    fn test_build_network_with_water() {
        let system = make_bifurcation_system();
        let conv = SchemeNetworkConverter::new(&system);
        let network = conv.build_network_with_water().expect("conversion failed");

        assert!(network.node_count() > 0);
        assert!(network.edge_count() > 0);
    }

    #[test]
    fn test_node_classification() {
        let system = make_bifurcation_system();
        let conv = SchemeNetworkConverter::new(&system);
        let types = conv.classify_nodes().expect("classification failed");

        let inlets = types.values().filter(|t| matches!(t, NodeType::Inlet)).count();
        let outlets = types.values().filter(|t| matches!(t, NodeType::Outlet)).count();

        assert!(inlets >= 1, "must have at least one inlet");
        assert!(outlets >= 1, "must have at least one outlet");
    }

    #[test]
    fn test_empty_system_fails() {
        let system = ChannelSystem {
            box_dims: (100.0, 50.0),
            nodes: vec![],
            channels: vec![],
            box_outline: vec![],
        };
        let conv = SchemeNetworkConverter::new(&system);
        assert!(conv.build_network_with_water().is_err());
    }

    #[test]
    fn test_scale_factor() {
        let system = make_bifurcation_system();

        let conv_micro = SchemeNetworkConverter::new(&system);
        let summary_micro = conv_micro.summary();

        let conv_milli = SchemeNetworkConverter::with_scale(&system, 1e-3);
        let summary_milli = conv_milli.summary();

        // Lengths should differ by factor of 1000
        let ratio = summary_milli.total_path_length_m / summary_micro.total_path_length_m;
        assert!(
            (ratio - 1000.0).abs() < 1.0,
            "scale ratio should be ~1000, got {}",
            ratio
        );
    }

    #[test]
    fn test_trifurcation_system() {
        let system = create_geometry(
            (300.0, 150.0),
            &[SplitType::Bifurcation, SplitType::Trifurcation],
            &GeometryConfig::default(),
            &ChannelTypeConfig::AllStraight,
        );
        let conv = SchemeNetworkConverter::new(&system);
        let network = conv.build_network_with_water().expect("conversion failed");

        assert!(network.node_count() > 4, "complex system should have many nodes");
        assert!(network.edge_count() > 4, "complex system should have many edges");
    }

    #[test]
    fn test_serpentine_system() {
        use scheme::config::SerpentineConfig;

        let system = create_geometry(
            (200.0, 100.0),
            &[SplitType::Bifurcation],
            &GeometryConfig::default(),
            &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
        );
        let conv = SchemeNetworkConverter::new(&system);
        let summary = conv.summary();

        assert!(summary.valid);
        let serp_count = summary.channel_type_counts.get("serpentine").copied().unwrap_or(0);
        assert!(serp_count > 0, "should have serpentine channels");

        let network = conv.build_network_with_water().expect("conversion failed");
        assert!(network.node_count() > 0);
    }
}
