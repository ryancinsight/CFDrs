//! Network builder for constructing fluid networks

use super::{
    junction_losses::apply_blueprint_junction_losses, Edge, EdgeProperties, EdgeType, NetworkGraph,
    Node, NodeType,
};
use cfd_core::error::Result;
use cfd_schematics::domain::model::{ChannelShape, NetworkBlueprint};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

use cfd_core::physics::fluid::FluidTrait;

/// Build a solver-ready [`Network`](crate::domain::network::wrapper::Network) from a
/// [`NetworkBlueprint`].
///
/// This is the **canonical entry-point** for constructing a `cfd-1d` network.
/// Callers obtain a blueprint from [`cfd_schematics::geometry::types::ChannelSystem::to_blueprint`]
/// or from the `cfd_schematics::interface::presets` factory functions.
///
/// Fidelity boundary: the 1D solve consumes blueprint lengths, cross-sections,
/// serpentine bend metadata, venturi throat metadata, and junction-angle
/// metadata. It does not resolve 2D/3D separation zones, recirculation pockets,
/// or secondary vortical structure beyond the reduced resistance models.
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
#[allow(clippy::too_many_lines)]
pub fn network_from_blueprint<T, F>(
    blueprint: &NetworkBlueprint,
    fluid: F,
) -> Result<crate::domain::network::wrapper::Network<T, F>>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T> + Clone,
{
    use crate::domain::channel::{
        ChannelGeometry, ChannelType, CrossSection, SurfaceProperties, Wettability,
    };
    use crate::physics::resistance::{
        ChannelGeometry as ResGeometry, FlowConditions, ResistanceCalculator,
    };
    use cfd_schematics::domain::model::CrossSectionSpec;
    use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;

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
    let mut edge_specs: Vec<(
        petgraph::graph::EdgeIndex,
        &cfd_schematics::domain::model::ChannelSpec,
    )> = Vec::new();

    for ch_spec in &blueprint.channels {
        let from_id = ch_spec.from.as_str();
        let to_id = ch_spec.to.as_str();

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

        if ch_spec.resistance < 0.0 || !ch_spec.resistance.is_finite() {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Channel '{}' has invalid resistance: {}",
                ch_spec.id.as_str(),
                ch_spec.resistance
            )));
        }
        let seed_r = T::from_f64(ch_spec.resistance.max(f64::EPSILON)).ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration(
                "Numeric conversion failed for resistance".to_string(),
            )
        })?;
        let mut edge = Edge::new(ch_spec.id.as_str().to_string(), ch_spec.kind);
        edge.resistance = seed_r;
        edge.area = match ch_spec.cross_section {
            CrossSectionSpec::Circular { diameter_m } => {
                if diameter_m <= 0.0 || !diameter_m.is_finite() {
                    return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                        "Channel '{}' has invalid circular diameter: {}",
                        ch_spec.id.as_str(),
                        diameter_m
                    )));
                }
                T::from_f64(std::f64::consts::PI * (diameter_m / 2.0).powi(2)).ok_or_else(|| {
                    cfd_core::error::Error::InvalidConfiguration(
                        "Numeric conversion failed for area".to_string(),
                    )
                })?
            }
            CrossSectionSpec::Rectangular { width_m, height_m } => {
                if width_m <= 0.0
                    || !width_m.is_finite()
                    || height_m <= 0.0
                    || !height_m.is_finite()
                {
                    return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                        "Channel '{}' has invalid rectangular dimensions: {}x{}",
                        ch_spec.id.as_str(),
                        width_m,
                        height_m
                    )));
                }
                T::from_f64(width_m * height_m).ok_or_else(|| {
                    cfd_core::error::Error::InvalidConfiguration(
                        "Numeric conversion failed for area".to_string(),
                    )
                })?
            }
        };

        let edge_idx = builder.add_edge(from_idx, to_idx, edge);
        edge_specs.push((edge_idx, ch_spec));
    }

    let graph = builder.build()?;
    let mut network = crate::domain::network::wrapper::Network::new(graph, fluid);

    let calculator: ResistanceCalculator<T> = ResistanceCalculator::new();

    for (edge_idx, ch_spec) in &edge_specs {
        if ch_spec.length_m <= 0.0 || !ch_spec.length_m.is_finite() {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Channel '{}' has invalid length: {}",
                ch_spec.id.as_str(),
                ch_spec.length_m
            )));
        }
        let length = T::from_f64(ch_spec.length_m).ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration(
                "Numeric conversion failed for length".to_string(),
            )
        })?;

        let (area, dh, cross_section, res_geom) =
            match ch_spec.cross_section {
                CrossSectionSpec::Circular { diameter_m } => {
                    let d = T::from_f64(diameter_m).ok_or_else(|| {
                        cfd_core::error::Error::InvalidConfiguration(
                            "Numeric conversion failed for diameter".to_string(),
                        )
                    })?;
                    let a = T::from_f64(std::f64::consts::PI * (diameter_m / 2.0).powi(2))
                        .ok_or_else(|| {
                            cfd_core::error::Error::InvalidConfiguration(
                                "Numeric conversion failed for area".to_string(),
                            )
                        })?;
                    (
                        a,
                        Some(d),
                        CrossSection::Circular { diameter: d },
                        ResGeometry::Circular {
                            diameter: d,
                            length,
                        },
                    )
                }
                CrossSectionSpec::Rectangular { width_m, height_m } => {
                    let w = T::from_f64(width_m).ok_or_else(|| {
                        cfd_core::error::Error::InvalidConfiguration(
                            "Numeric conversion failed for width".to_string(),
                        )
                    })?;
                    let h = T::from_f64(height_m).ok_or_else(|| {
                        cfd_core::error::Error::InvalidConfiguration(
                            "Numeric conversion failed for height".to_string(),
                        )
                    })?;
                    let a = T::from_f64(width_m * height_m).ok_or_else(|| {
                        cfd_core::error::Error::InvalidConfiguration(
                            "Numeric conversion failed for area".to_string(),
                        )
                    })?;
                    let dh = T::from_f64(2.0 * width_m * height_m / (width_m + height_m))
                        .ok_or_else(|| {
                            cfd_core::error::Error::InvalidConfiguration(
                                "Numeric conversion failed for hydraulic diameter".to_string(),
                            )
                        })?;
                    (
                        a,
                        Some(dh),
                        CrossSection::Rectangular {
                            width: w,
                            height: h,
                        },
                        ResGeometry::Rectangular {
                            width: w,
                            height: h,
                            length,
                        },
                    )
                }
            };

        let channel_type = match ch_spec.channel_shape {
            ChannelShape::Serpentine { segments, .. } => ChannelType::Serpentine {
                turns: segments.saturating_sub(1),
            },
            ChannelShape::Straight => ChannelType::Straight,
        };

        let channel_geometry = ChannelGeometry {
            channel_type,
            length,
            cross_section,
            surface: SurfaceProperties {
                roughness: T::from_f64(1e-7).expect("Mathematical constant conversion compromised"),
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
            resistance: T::from_f64(ch_spec.resistance)
                .expect("Mathematical constant conversion compromised"),
            geometry: Some(channel_geometry),
            properties: HashMap::new(),
        };
        network.add_edge_properties(*edge_idx, props);

        // Compute refined resistance with the fluid model at near-zero flow.
        let mut refined = None;
        if let Some(venturi) = ch_spec.venturi_geometry.as_ref().or_else(|| {
            ch_spec
                .metadata
                .as_ref()
                .and_then(|meta| meta.get::<VenturiGeometryMetadata>())
        }) {
            refined = venturi_coefficients::<T, F>(venturi, network.fluid()).ok();
        }

        let mut conds = FlowConditions::new(T::zero());
        conds.flow_rate =
            Some(T::from_f64(1e-12).expect("Mathematical constant conversion compromised"));
        if refined.is_none() {
            refined = calculator
                .calculate_coefficients_auto(&res_geom, network.fluid(), &conds)
                .ok();
        }

        if let Some((r, k)) = refined {
            if let Some(edge) = network.graph.edge_weight_mut(*edge_idx) {
                edge.resistance = r;
                edge.quad_coeff = k;
            }
        }
    }

    let _junction_stats = apply_blueprint_junction_losses(&mut network, blueprint);

    // ── Dean-flow correction post-pass ───────────────────────────────────────
    //
    // For channels flagged as `ChannelShape::Serpentine`, replace the base
    // rectangular/circular resistance with the full `SerpentineModel` which
    // includes:
    //   1. Shah-London straight-section friction (same as before)
    //   2. Dean-number curvature enhancement to friction
    //   3. Bend K-factor minor losses at each 180° U-turn
    //
    // This gives physically correct pressure-drop predictions for serpentine
    // channels instead of treating them as straight ducts.
    {
        use crate::physics::resistance::models::FlowConditions as FC;
        use crate::physics::resistance::models::ResistanceModel as RM;
        use crate::physics::resistance::models::{SerpentineCrossSection, SerpentineModel};

        for (edge_idx, ch_spec) in &edge_specs {
            let ChannelShape::Serpentine {
                segments: total_segments,
                bend_radius_m,
            } = ch_spec.channel_shape
            else {
                continue;
            };

            let length_t = T::from_f64(ch_spec.length_m).ok_or_else(|| {
                cfd_core::error::Error::InvalidConfiguration(format!(
                    "Numeric conversion failed for length {}",
                    ch_spec.id.as_str()
                ))
            })?;
            let bend_r = if bend_radius_m <= 0.0 || !bend_radius_m.is_finite() {
                return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                    "Channel '{}' has invalid serpentine bend radius: {}",
                    ch_spec.id.as_str(),
                    bend_radius_m
                )));
            } else {
                T::from_f64(bend_radius_m).ok_or_else(|| {
                    cfd_core::error::Error::InvalidConfiguration(
                        "Numeric conversion failed for bend_radius".to_string(),
                    )
                })?
            };

            let cross_section = match ch_spec.cross_section {
                CrossSectionSpec::Circular { diameter_m } => SerpentineCrossSection::Circular {
                    diameter: diameter_m,
                },
                CrossSectionSpec::Rectangular { width_m, height_m } => {
                    SerpentineCrossSection::Rectangular {
                        width: width_m,
                        height: height_m,
                    }
                }
            };

            // Each channel spec is ONE segment; bend count = (total_segments - 1)
            // distributed evenly across segments. Each segment "owns" the bend
            // at its downstream end, so segment i sees
            //   bends_owned = if i < total_segments-1 { 1 } else { 0 }
            // In aggregate the SerpentineModel with num_segments=1 produces
            // num_bends=0 for the last segment. To average the bend losses
            // evenly, we construct a model for ONE segment with a fractional
            // bend mass: total model = 1 segment, bend losses from
            // (total_segments - 1) / total_segments bends, scaled by 1 segment.
            //
            // Simpler: build a SerpentineModel for the whole serpentine, then
            // divide per-segment to get the per-channel coefficients.

            let total_length = length_t * T::from_usize(total_segments).unwrap_or(T::one());

            let serp_model =
                SerpentineModel::new(total_length, total_segments, cross_section, bend_r);

            let mut conds: FC<T> = FC::new(T::zero());
            conds.flow_rate =
                Some(T::from_f64(1e-12).expect("Mathematical constant conversion compromised"));

            if let Ok((r_total, k_total)) =
                serp_model.calculate_coefficients(network.fluid(), &conds)
            {
                // Per-segment share
                let n_segs = T::from_usize(total_segments.max(1)).unwrap_or(T::one());
                let r_seg = r_total / n_segs;
                let k_seg = k_total / n_segs;

                if let Some(edge) = network.graph.edge_weight_mut(*edge_idx) {
                    edge.resistance = r_seg;
                    edge.quad_coeff = k_seg;
                }
            }
        }
    }

    Ok(network)
}

fn venturi_coefficients<T, F>(
    metadata: &cfd_schematics::geometry::metadata::VenturiGeometryMetadata,
    fluid: &F,
) -> Result<(T, T)>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T> + Clone,
{
    use crate::physics::resistance::models::{
        ExpansionType, FlowConditions as ModelFlowConditions, ResistanceModel, VenturiGeometry,
        VenturiModel,
    };
    use crate::{discharge_coefficient_from_convergent_half_angle_deg, venturi_taper_length_m};

    let hydraulic_diameter =
        |width_m: f64, height_m: f64| 2.0 * width_m * height_m / (width_m + height_m).max(1e-18);
    let inlet_d = hydraulic_diameter(metadata.inlet_width_m, metadata.throat_height_m);
    let throat_d = hydraulic_diameter(metadata.throat_width_m, metadata.throat_height_m);
    let outlet_d = hydraulic_diameter(metadata.outlet_width_m, metadata.throat_height_m);
    let conv_len = venturi_taper_length_m(
        metadata.inlet_width_m,
        metadata.throat_width_m,
        metadata.convergent_half_angle_deg,
    );
    let diff_len = venturi_taper_length_m(
        metadata.outlet_width_m,
        metadata.throat_width_m,
        metadata.divergent_half_angle_deg,
    );
    let total_length = metadata.throat_length_m + conv_len + diff_len;

    let inlet_d_t = T::from_f64(inlet_d).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration("inlet diameter conversion failed".to_string())
    })?;
    let throat_d_t = T::from_f64(throat_d).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration(
            "throat diameter conversion failed".to_string(),
        )
    })?;
    let outlet_d_t = T::from_f64(outlet_d).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration(
            "outlet diameter conversion failed".to_string(),
        )
    })?;
    let throat_len_t = T::from_f64(total_length).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration("throat length conversion failed".to_string())
    })?;
    let total_len_t = T::from_f64(total_length).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration(
            "venturi total length conversion failed".to_string(),
        )
    })?;

    let mut model = VenturiModel::new(inlet_d_t, throat_d_t, outlet_d_t, throat_len_t, total_len_t)
        .with_geometry(VenturiGeometry::Custom {
            discharge_coefficient: discharge_coefficient_from_convergent_half_angle_deg(
                metadata.convergent_half_angle_deg,
            ),
        })
        .with_expansion(ExpansionType::Gradual {
            half_angle_deg: metadata.divergent_half_angle_deg.clamp(1.0, 45.0),
        });
    model.throat_roughness =
        T::from_f64(1e-7).expect("Mathematical constant conversion compromised");

    let mut conds = ModelFlowConditions::new(T::zero());
    conds.flow_rate =
        Some(T::from_f64(1e-12).expect("Mathematical constant conversion compromised"));
    model.calculate_coefficients(fluid, &conds)
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
