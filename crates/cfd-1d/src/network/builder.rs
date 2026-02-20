//! Network builder for constructing fluid networks

use super::{Edge, EdgeProperties, EdgeType, NetworkGraph, Node, NodeType};
use cfd_core::error::Result;
use cfd_schematics::domain::model::NetworkBlueprint;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

use cfd_core::physics::fluid::FluidTrait;

/// Build a solver-ready [`Network`](super::super::network::wrapper::Network) from a
/// [`NetworkBlueprint`].
///
/// This is the **canonical entry-point** for constructing a `cfd-1d` network.
/// Callers obtain a blueprint from [`cfd_schematics::geometry::types::ChannelSystem::to_blueprint`]
/// or from the `cfd_schematics::interface::presets` factory functions.
///
/// # Type parameters
/// - `T` : numeric field type (e.g. `f64`)
/// - `F` : fluid model implementing [`FluidTrait<T>`]
///
/// # Errors
/// Returns [`cfd_core::error::Error`] when:
/// - the blueprint has no nodes or no channels
/// - node ID references are inconsistent
/// - resistance coefficients are invalid (negatives, all-zero)
pub fn network_from_blueprint<T, F>(
    blueprint: &NetworkBlueprint,
    fluid: F,
) -> Result<super::super::network::wrapper::Network<T, F>>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T> + Clone,
{
    use crate::{
        channel::{ChannelGeometry, ChannelType, CrossSection, SurfaceProperties, Wettability},
        resistance::{ChannelGeometry as ResGeometry, FlowConditions, ResistanceCalculator},
    };
    use cfd_schematics::domain::model::CrossSectionSpec;

    if blueprint.nodes.is_empty() {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "NetworkBlueprint has no nodes".to_string(),
        ));
    }
    if blueprint.channels.is_empty() {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "NetworkBlueprint has no channels".to_string(),
        ));
    }

    let mut builder = NetworkBuilder::<T>::new();
    // Map NodeId string → petgraph NodeIndex
    let mut node_map: HashMap<String, petgraph::graph::NodeIndex> = HashMap::new();

    for node_spec in &blueprint.nodes {
        let node = Node::new(node_spec.id.as_str().to_string(), node_spec.kind);
        let idx = builder.add_node(node);
        node_map.insert(node_spec.id.as_str().to_string(), idx);
    }

    // Collect (EdgeIndex, ChannelSpec ref) for post-build property population.
    let mut edge_specs: Vec<(petgraph::graph::EdgeIndex, &cfd_schematics::domain::model::ChannelSpec)> =
        Vec::new();

    for ch_spec in &blueprint.channels {
        let from_id = ch_spec.from.as_str();
        let to_id   = ch_spec.to.as_str();

        let from_idx = node_map.get(from_id).copied().ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration(format!(
                "Channel '{}' references missing node '{}'",
                ch_spec.id.as_str(),
                from_id
            ))
        })?;
        let to_idx = node_map.get(to_id).copied().ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration(format!(
                "Channel '{}' references missing node '{}'",
                ch_spec.id.as_str(),
                to_id
            ))
        })?;

        let seed_r = T::from_f64(ch_spec.resistance.max(f64::EPSILON))
            .unwrap_or_else(T::default_epsilon);
        let mut edge = Edge::new(ch_spec.id.as_str().to_string(), ch_spec.kind);
        edge.resistance = seed_r;

        let edge_idx = builder.add_edge(from_idx, to_idx, edge);
        edge_specs.push((edge_idx, ch_spec));
    }

    let graph = builder.build()?;
    let mut network = super::super::network::wrapper::Network::new(graph, fluid);

    let calculator: ResistanceCalculator<T> = ResistanceCalculator::new();

    for (edge_idx, ch_spec) in &edge_specs {
        let length     = T::from_f64(ch_spec.length_m).unwrap_or(T::zero());
        let (area, dh, cross_section, res_geom) = match ch_spec.cross_section {
            CrossSectionSpec::Circular { diameter_m } => {
                let d  = T::from_f64(diameter_m).unwrap_or(T::zero());
                let a  = T::from_f64(std::f64::consts::PI * (diameter_m / 2.0).powi(2))
                    .unwrap_or(T::zero());
                (
                    a,
                    Some(d),
                    CrossSection::Circular { diameter: d },
                    ResGeometry::Circular { diameter: T::from_f64(diameter_m).unwrap_or(T::zero()), length },
                )
            }
            CrossSectionSpec::Rectangular { width_m, height_m } => {
                let w  = T::from_f64(width_m).unwrap_or(T::zero());
                let h  = T::from_f64(height_m).unwrap_or(T::zero());
                let a  = T::from_f64(width_m * height_m).unwrap_or(T::zero());
                let dh = T::from_f64(2.0 * width_m * height_m / (width_m + height_m))
                    .unwrap_or(T::zero());
                (
                    a,
                    Some(dh),
                    CrossSection::Rectangular { width: w, height: h },
                    ResGeometry::Rectangular {
                        width: T::from_f64(width_m).unwrap_or(T::zero()),
                        height: T::from_f64(height_m).unwrap_or(T::zero()),
                        length,
                    },
                )
            }
        };

        let channel_geometry = ChannelGeometry {
            channel_type: ChannelType::Straight,
            length,
            cross_section,
            surface: SurfaceProperties {
                roughness: T::from_f64(1e-7).unwrap_or(T::zero()),
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        };

        let props = EdgeProperties {
            id: ch_spec.id.as_str().to_string(),
            component_type: super::ComponentType::Pipe,
            length,
            area,
            hydraulic_diameter: dh,
            resistance: T::from_f64(ch_spec.resistance).unwrap_or(T::default_epsilon()),
            geometry: Some(channel_geometry),
            properties: HashMap::new(),
        };
        network.add_edge_properties(*edge_idx, props);

        // Compute refined resistance with the fluid model at near-zero flow.
        let mut conds = FlowConditions::new(T::zero());
        conds.flow_rate = Some(T::from_f64(1e-12).unwrap_or(T::zero()));
        if let Ok((r, k)) = calculator.calculate_coefficients_auto(
            &res_geom,
            network.fluid(),
            &conds,
        ) {
            if let Some(edge) = network.graph.edge_weight_mut(*edge_idx) {
                edge.resistance = r;
                edge.quad_coeff = k;
            }
        }
    }

    Ok(network)
}

/// Builder for constructing network graphs directly (internal / advanced use).
///
/// Prefer [`network_from_blueprint`] for any network derived from a
/// [`cfd_schematics`] topology. This builder is retained for low-level
/// graph construction in blood-vessel models and other domain-specific
/// applications that build the graph programmatically.
pub struct NetworkBuilder<T: RealField + Copy> {
    graph: NetworkGraph<T>,
}

impl<T: RealField + Copy> NetworkBuilder<T> {
    /// Create a new network builder
    #[must_use]
    pub fn new() -> Self {
        Self {
            graph: NetworkGraph::new(),
        }
    }

    /// Add a node to the network
    pub fn add_node(&mut self, node: Node<T>) -> petgraph::graph::NodeIndex {
        self.graph.add_node(node)
    }

    /// Add an edge between two nodes
    pub fn add_edge(
        &mut self,
        from: petgraph::graph::NodeIndex,
        to: petgraph::graph::NodeIndex,
        edge: Edge<T>,
    ) -> petgraph::graph::EdgeIndex {
        self.graph.add_edge(from, to, edge)
    }

    /// Create an inlet node
    pub fn add_inlet(&mut self, id: String) -> petgraph::graph::NodeIndex {
        self.add_node(Node::new(id, NodeType::Inlet))
    }

    /// Create an outlet node
    pub fn add_outlet(&mut self, id: String) -> petgraph::graph::NodeIndex {
        self.add_node(Node::new(id, NodeType::Outlet))
    }

    /// Create a junction node
    pub fn add_junction(&mut self, id: String) -> petgraph::graph::NodeIndex {
        self.add_node(Node::new(id, NodeType::Junction))
    }

    /// Connect two nodes with a pipe
    pub fn connect_with_pipe(
        &mut self,
        from: petgraph::graph::NodeIndex,
        to: petgraph::graph::NodeIndex,
        id: String,
    ) -> petgraph::graph::EdgeIndex {
        self.add_edge(from, to, Edge::new(id, EdgeType::Pipe))
    }

    /// Build the final network with validation
    pub fn build(self) -> Result<NetworkGraph<T>> {
        if self.graph.node_count() == 0 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network has no nodes".to_string(),
            ));
        }

        let has_inlet = self
            .graph
            .node_weights()
            .any(|n| matches!(n.node_type, NodeType::Inlet));
        let has_outlet = self
            .graph
            .node_weights()
            .any(|n| matches!(n.node_type, NodeType::Outlet));

        if !has_inlet {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network requires at least one inlet".to_string(),
            ));
        }
        if !has_outlet {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network requires at least one outlet".to_string(),
            ));
        }

        use petgraph::algo::connected_components;
        if connected_components(&self.graph) > 1 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network has disconnected components".to_string(),
            ));
        }

        let eps = T::default_epsilon();
        for edge_ref in self.graph.edge_references() {
            let idx = edge_ref.id();
            let w = edge_ref.weight();
            let r = w.resistance;
            let k = w.quad_coeff;
            if r < T::zero() {
                return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                    "Edge {} has negative resistance: {}",
                    idx.index(),
                    r
                )));
            }
            if k < T::zero() {
                return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                    "Edge {} has negative quadratic coefficient: {}",
                    idx.index(),
                    k
                )));
            }
            if r.abs() < eps && k.abs() < eps {
                return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                    "Edge {} has zero resistance and zero quadratic coefficient",
                    idx.index()
                )));
            }
        }

        Ok(self.graph)
    }
}

impl<T: RealField + Copy> Default for NetworkBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}
