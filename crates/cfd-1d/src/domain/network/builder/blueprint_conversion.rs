//! Blueprint-to-network conversion pipeline.
//!
//! Transforms a [`NetworkBlueprint`] from `cfd-schematics` into a solver-ready
//! [`Network`] with physically refined resistance coefficients.

use super::super::{
    blueprint_validation::validate_blueprint_for_1d_solve,
    junction_losses::apply_blueprint_junction_losses, Edge, EdgeProperties, Node,
    ResistanceUpdatePolicy,
};
use super::network_builder::NetworkBuilder;
use super::venturi_coefficients::venturi_coefficients;
use crate::physics::resistance::models::BendType;
use cfd_core::error::{Error, Result};
use cfd_schematics::domain::model::{ChannelShape, EdgeKind, NetworkBlueprint};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;
use std::hash::BuildHasher;

use crate::domain::network::wrapper::{
    blood_microchannel_apparent_viscosity, EDGE_PROPERTY_HEMATOCRIT,
    EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S,
};
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
    validate_blueprint_for_1d_solve(blueprint)?;

    use crate::domain::channel::{
        ChannelGeometry, ChannelType, CrossSection, SurfaceProperties, Wettability,
    };
    use crate::physics::resistance::{
        ChannelGeometry as ResGeometry, FlowConditions, ResistanceCalculator,
    };
    use cfd_schematics::domain::model::CrossSectionSpec;
    use cfd_schematics::geometry::metadata::{ChannelVenturiSpec, VenturiGeometryMetadata};

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

    let mut node_degree: HashMap<&str, usize> = HashMap::with_capacity(blueprint.nodes.len());
    for channel in &blueprint.channels {
        *node_degree.entry(channel.from.as_str()).or_default() += 1;
        *node_degree.entry(channel.to.as_str()).or_default() += 1;
    }
    let expects_explicit_junction_geometry = blueprint
        .nodes
        .iter()
        .any(|node| node.junction_geometry.is_some());
    if expects_explicit_junction_geometry {
        for node in &blueprint.nodes {
            let degree = node_degree.get(node.id.as_str()).copied().unwrap_or(0);
            if degree >= 3 && node.junction_geometry.is_none() {
                return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                    "Branching node '{}' is missing explicit junction metadata",
                    node.id.as_str()
                )));
            }
        }
    }

    let mut builder = NetworkBuilder::<T>::new();
    let mut node_map: HashMap<&str, petgraph::graph::NodeIndex> =
        HashMap::with_capacity(blueprint.nodes.len());

    for node_spec in &blueprint.nodes {
        let node = Node::new(node_spec.id.as_str().to_string(), node_spec.kind);
        let idx = builder.add_node(node);
        node_map.insert(node_spec.id.as_str(), idx);
    }

    let mut edge_specs: Vec<(
        petgraph::graph::EdgeIndex,
        &cfd_schematics::domain::model::ChannelSpec,
    )> = Vec::with_capacity(blueprint.channels.len());

    for ch_spec in &blueprint.channels {
        if matches!(ch_spec.kind, EdgeKind::Pump) {
            return Err(Error::InvalidConfiguration(format!(
                "Pump edge '{}' is not supported by the 1D passive-network solver yet",
                ch_spec.id.as_str()
            )));
        }

        let serial_venturi_throats = ch_spec
            .metadata
            .as_ref()
            .and_then(|meta| meta.get::<ChannelVenturiSpec>())
            .map_or(0, |spec| spec.n_throats);
        if serial_venturi_throats > 0
            && ch_spec.venturi_geometry.is_none()
            && ch_spec
                .metadata
                .as_ref()
                .and_then(|meta| meta.get::<VenturiGeometryMetadata>())
                .is_none()
        {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Channel '{}' declares venturi throats but has no explicit venturi geometry",
                ch_spec.id.as_str()
            )));
        }

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
    let is_blood_like = network
        .fluid()
        .name()
        .to_ascii_lowercase()
        .contains("blood");
    let blood_state = network.fluid().properties_at(
        T::from_f64(cfd_core::physics::constants::physics::thermo::T_STANDARD)
            .unwrap_or_else(T::one),
        T::from_f64(cfd_core::physics::constants::physics::thermo::P_ATM).unwrap_or_else(T::zero),
    )?;
    let default_hematocrit = T::from_f64(0.45).unwrap_or_else(T::zero);
    let default_plasma_viscosity =
        blood_state.dynamic_viscosity / T::from_f64(3.2).unwrap_or_else(T::one);

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

        let mut edge_property_overrides = HashMap::new();
        if is_blood_like {
            edge_property_overrides
                .insert(EDGE_PROPERTY_HEMATOCRIT.to_string(), default_hematocrit);
            edge_property_overrides.insert(
                EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S.to_string(),
                default_plasma_viscosity,
            );
        }

        let props = EdgeProperties {
            id: ch_spec.id.as_str().to_string(),
            component_type: super::super::ComponentType::Pipe,
            length,
            area,
            hydraulic_diameter: dh,
            resistance: T::from_f64(ch_spec.resistance)
                .expect("Mathematical constant conversion compromised"),
            geometry: Some(channel_geometry),
            resistance_update_policy: ResistanceUpdatePolicy::FlowDependent,
            properties: edge_property_overrides,
        };
        network.add_edge_properties(*edge_idx, props);

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
            let mut resistance_scale = T::one();
            if let Some(dh_val) = dh {
                let dh_f64 = nalgebra::try_convert::<T, f64>(dh_val).unwrap_or(1e-3);
                let l_f64 = nalgebra::try_convert::<T, f64>(length).unwrap_or(0.01);
                let l_over_dh = l_f64 / dh_f64.max(1e-12);

                if l_over_dh < 50.0 {
                    let re_estimate = 10.0_f64;
                    let multiplier =
                        crate::physics::resistance::models::durst_resistance_multiplier(
                            re_estimate,
                            l_over_dh,
                        );
                    resistance_scale *= T::from_f64(multiplier).unwrap_or(T::one());
                }

                if is_blood_like && dh_f64 < 300.0e-6 {
                    let q_seed_t = conds.flow_rate.unwrap_or_else(T::zero);
                    let q_seed = nalgebra::try_convert::<T, f64>(q_seed_t).unwrap_or(0.0);
                    let area_f64 = nalgebra::try_convert::<T, f64>(area).unwrap_or(0.0);
                    if let Some(target_mu) = blood_microchannel_apparent_viscosity(
                        dh_val,
                        q_seed_t,
                        area,
                        default_hematocrit,
                        default_plasma_viscosity,
                    ) {
                        let local_shear_rate = if dh_f64 > 0.0 && area_f64 > 0.0 {
                            8.0 * (q_seed / area_f64).abs() / dh_f64
                        } else {
                            0.0
                        };
                        let local_mu = network
                            .fluid()
                            .viscosity_at_shear(
                                T::from_f64(local_shear_rate).unwrap_or_else(T::zero),
                                conds.temperature,
                                conds.pressure,
                            )
                            .unwrap_or(blood_state.dynamic_viscosity);
                        if local_mu > T::default_epsilon() {
                            resistance_scale *= target_mu / local_mu;
                        }
                    }
                }
            }

            if let Some(edge) = network.graph.edge_weight_mut(*edge_idx) {
                edge.resistance = r * resistance_scale;
                edge.quad_coeff = k;
            }
        }
    }

    let _junction_stats = apply_blueprint_junction_losses(&mut network, blueprint);
    let eps = T::default_epsilon();
    for (edge_idx, ch_spec) in &edge_specs {
        let has_venturi = ch_spec.venturi_geometry.is_some()
            || ch_spec
                .metadata
                .as_ref()
                .and_then(|meta| meta.get::<VenturiGeometryMetadata>())
                .is_some();
        let update_policy = match ch_spec.channel_shape {
            ChannelShape::Serpentine { .. } => ResistanceUpdatePolicy::FlowDependent,
            ChannelShape::Straight if has_venturi => ResistanceUpdatePolicy::FlowDependent,
            ChannelShape::Straight => network.graph.edge_weight(*edge_idx).map_or(
                ResistanceUpdatePolicy::FlowDependent,
                |edge| {
                    if edge.quad_coeff.abs() <= eps {
                        ResistanceUpdatePolicy::FlowInvariant
                    } else {
                        ResistanceUpdatePolicy::FlowDependent
                    }
                },
            ),
        };
        if let Some(props) = network.properties.get_mut(edge_idx) {
            props.resistance_update_policy = update_policy;
        }
    }

    // ── Dean-flow correction post-pass ───────────────────────────────────────
    {
        use crate::physics::resistance::models::FlowConditions as FC;
        use crate::physics::resistance::models::ResistanceModel;
        use crate::physics::resistance::models::{SerpentineCrossSection, SerpentineModel};

        for (edge_idx, ch_spec) in &edge_specs {
            let ChannelShape::Serpentine {
                segments: total_segments,
                bend_radius_m,
                wave_type,
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

            let total_length = length_t * T::from_usize(total_segments).unwrap_or(T::one());

            let serp_model = {
                use cfd_schematics::SerpentineWaveType;
                let base =
                    SerpentineModel::new(total_length, total_segments, cross_section, bend_r);
                match wave_type {
                    SerpentineWaveType::Square => base.with_bend_type(BendType::Sharp),
                    SerpentineWaveType::Sine | SerpentineWaveType::Triangular => base,
                }
            };

            let mut conds: FC<T> = FC::new(T::zero());
            conds.flow_rate =
                Some(T::from_f64(1e-12).expect("Mathematical constant conversion compromised"));

            if let Ok((r_total, k_total)) =
                serp_model.calculate_coefficients(network.fluid(), &conds)
            {
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

    for edge_ref in network.graph.edge_references() {
        let edge = edge_ref.weight();
        if !edge.resistance.is_finite() || edge.resistance <= T::zero() {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Edge '{}' has invalid effective linear resistance after refinement",
                edge.id
            )));
        }
        if !edge.quad_coeff.is_finite() || edge.quad_coeff < T::zero() {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Edge '{}' has invalid effective quadratic coefficient after refinement",
                edge.id
            )));
        }
    }

    Ok(network)
}

/// Apply schematics-authored branch boundary metadata to a solver network.
///
/// Explicit [`BranchBoundaryMetadata`] overrides the legacy inlet/outlet
/// defaults derived from [`cfd_schematics::domain::model::NodeKind`].
pub fn apply_blueprint_boundary_conditions<T, F, S>(
    network: &mut crate::domain::network::wrapper::Network<T, F>,
    blueprint: &NetworkBlueprint,
    node_indices: &HashMap<String, NodeIndex, S>,
    inlet_pressure: T,
    outlet_pressure: T,
) -> Result<()>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
    S: BuildHasher,
{
    use cfd_schematics::geometry::metadata::{BranchBoundaryMetadata, BranchBoundarySpecification};

    for node in &blueprint.nodes {
        let Some(node_idx) = node_indices.get(node.id.as_str()).copied() else {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Blueprint node '{}' is missing from the solver network",
                node.id.as_str()
            )));
        };

        if let Some(boundary) = node
            .metadata
            .as_ref()
            .and_then(|metadata| metadata.get::<BranchBoundaryMetadata>())
        {
            match boundary.boundary {
                BranchBoundarySpecification::Pressure { pressure_pa } => {
                    let pressure = T::from_f64(pressure_pa).ok_or_else(|| {
                        cfd_core::error::Error::InvalidConfiguration(format!(
                            "Blueprint node '{}' pressure boundary could not be converted",
                            node.id.as_str()
                        ))
                    })?;
                    network.set_pressure(node_idx, pressure);
                }
                BranchBoundarySpecification::FlowRate { flow_rate_m3_s } => {
                    let flow_rate = T::from_f64(flow_rate_m3_s).ok_or_else(|| {
                        cfd_core::error::Error::InvalidConfiguration(format!(
                            "Blueprint node '{}' flow boundary could not be converted",
                            node.id.as_str()
                        ))
                    })?;
                    network.set_neumann_flow(node_idx, flow_rate);
                }
            }
            continue;
        }

        match node.kind {
            cfd_schematics::domain::model::NodeKind::Inlet => {
                network.set_pressure(node_idx, inlet_pressure);
            }
            cfd_schematics::domain::model::NodeKind::Outlet => {
                network.set_pressure(node_idx, outlet_pressure);
            }
            cfd_schematics::domain::model::NodeKind::Reservoir
            | cfd_schematics::domain::model::NodeKind::Junction => {}
        }
    }

    Ok(())
}
