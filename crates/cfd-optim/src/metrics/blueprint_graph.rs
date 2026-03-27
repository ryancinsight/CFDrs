use std::collections::{HashMap, HashSet};

use cfd_1d::domain::network::network_from_blueprint;
use cfd_1d::{NetworkProblem, NetworkSolver, SolvePathStatus};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_schematics::domain::model::{ChannelShape, CrossSectionSpec, NetworkBlueprint};
use cfd_schematics::domain::therapy_metadata::TherapyZone;
use cfd_schematics::geometry::metadata::{
    ChannelPathMetadata, ChannelVenturiSpec, ChannelVisualRole, JunctionFamily,
    VenturiGeometryMetadata,
};
use petgraph::visit::EdgeRef;

use crate::constraints::BLOOD_DENSITY_KG_M3;
use crate::domain::BlueprintCandidate;
use crate::error::OptimError;

#[derive(Debug, Clone)]
pub struct BlueprintSolveSample<'a> {
    pub id: &'a str,
    pub from_node: &'a str,
    pub to_node: &'a str,
    pub length_m: f64,
    pub cross_section: CrossSectionSpec,
    pub channel_shape: ChannelShape,
    pub flow_m3_s: f64,
    pub from_pressure_pa: f64,
    pub is_treatment_channel: bool,
    pub is_venturi_channel: bool,
    pub venturi_throat_count: u8,
}

#[derive(Debug, Clone)]
pub struct BlueprintSolveSummary<'a> {
    pub inlet_pressure_pa: f64,
    pub inlet_flow_m3_s: f64,
    pub flow_uniformity: f64,
    pub venturi_flow_fraction: f64,
    pub mean_residence_time_s: f64,
    pub solve_path_status: SolvePathStatus,
    pub remerge_loss_pa: f64,
    pub channel_samples: Vec<BlueprintSolveSample<'a>>,
}

pub fn solve_blueprint_candidate(
    candidate: &BlueprintCandidate,
) -> Result<BlueprintSolveSummary<'_>, OptimError> {
    let blood = CassonBlood::<f64>::normal_blood();
    let mut network = network_from_blueprint(&candidate.blueprint, blood).map_err(|error| {
        OptimError::PhysicsError {
            id: candidate.id.clone(),
            reason: format!("network_from_blueprint failed: {error}"),
        }
    })?;

    let mut inlet_nodes = Vec::new();
    let mut outlet_nodes = Vec::new();
    for index in network.graph.node_indices() {
        if let Some(node) = network.graph.node_weight(index) {
            if node.node_type == cfd_1d::NodeType::Inlet {
                inlet_nodes.push(index);
            } else if node.node_type == cfd_1d::NodeType::Outlet {
                outlet_nodes.push(index);
            }
        }
    }
    if inlet_nodes.is_empty() || outlet_nodes.is_empty() {
        return Err(OptimError::PhysicsError {
            id: candidate.id.clone(),
            reason: "blueprint is missing inlet/outlet nodes".to_string(),
        });
    }

    let q_per_inlet = candidate.operating_point.flow_rate_m3_s / inlet_nodes.len() as f64;
    for inlet in &inlet_nodes {
        network.set_neumann_flow(*inlet, q_per_inlet);
    }
    for outlet in &outlet_nodes {
        network.set_pressure(*outlet, 0.0);
    }

    let config = cfd_1d::solver::core::SolverConfig {
        tolerance: 1.0e-5,
        max_iterations: if cfg!(debug_assertions) { 50 } else { 500 },
    };
    let solver = NetworkSolver::<f64, CassonBlood<f64>>::with_config(config);
    let primary_problem = NetworkProblem::new(network);
    let (solved, solve_path_status, _diagnostics) = match solver
        .solve_network_with_diagnostics(&primary_problem)
    {
        Ok((solved, diagnostics)) => (solved, SolvePathStatus::PrimaryConverged, diagnostics),
        Err(primary_error) => {
            let mut linear_network = primary_problem.network.clone();
            let edge_indices: Vec<_> = linear_network.graph.edge_indices().collect();
            for edge_index in edge_indices {
                linear_network.set_flow_rate(edge_index, q_per_inlet);
            }
            let _ = linear_network.update_resistances();
            for properties in linear_network.properties.values_mut() {
                properties.geometry = None;
            }
            for edge_index in linear_network.graph.edge_indices() {
                if let Some(edge) = linear_network.graph.edge_weight_mut(edge_index) {
                    edge.quad_coeff = 0.0;
                    edge.resistance = edge.resistance.max(1.0e-12);
                }
            }
            for inlet in &inlet_nodes {
                linear_network.set_pressure(*inlet, 1.0);
            }
            for outlet in &outlet_nodes {
                linear_network.set_pressure(*outlet, 0.0);
            }

            let fallback_problem = NetworkProblem::new(linear_network);
            let mut fallback_solved =
                solver
                    .solve_network(&fallback_problem)
                    .map_err(|fallback_error| OptimError::PhysicsError {
                        id: candidate.id.clone(),
                        reason: format!(
                        "NetworkSolver failed: {primary_error}; fallback failed: {fallback_error}"
                    ),
                    })?;

            let mut inlet_flow_unit = 0.0_f64;
            for edge_ref in fallback_solved.graph.edge_references() {
                let edge_index = edge_ref.id();
                let Some((from_index, _)) = fallback_solved.graph.edge_endpoints(edge_index) else {
                    continue;
                };
                if inlet_nodes.contains(&from_index) {
                    inlet_flow_unit += fallback_solved
                        .flow_rates
                        .get(edge_index.index())
                        .copied()
                        .unwrap_or(edge_ref.weight().flow_rate)
                        .abs();
                }
            }
            if inlet_flow_unit <= 1.0e-18 {
                return Err(OptimError::PhysicsError {
                    id: candidate.id.clone(),
                    reason: "fallback solve produced zero inlet flow".to_string(),
                });
            }
            let scale = candidate.operating_point.flow_rate_m3_s / inlet_flow_unit;
            for pressure in &mut fallback_solved.pressures {
                *pressure *= scale;
            }
            let edge_indices: Vec<_> = fallback_solved.graph.edge_indices().collect();
            for edge_index in edge_indices {
                let idx = edge_index.index();
                if let Some(flow) = fallback_solved.flow_rates.get_mut(idx) {
                    *flow *= scale;
                    if let Some(edge) = fallback_solved.graph.edge_weight_mut(edge_index) {
                        edge.flow_rate = *flow;
                    }
                } else if let Some(edge) = fallback_solved.graph.edge_weight_mut(edge_index) {
                    edge.flow_rate *= scale;
                }
            }
            (
                fallback_solved,
                SolvePathStatus::RecoveredLinearized {
                    reason: primary_error.reason,
                },
                primary_error.diagnostics,
            )
        }
    };

    let inlet_pressure_pa = inlet_nodes
        .iter()
        .filter_map(|index| solved.pressures.get(index.index()).copied())
        .sum::<f64>()
        / inlet_nodes.len() as f64;

    let mut edge_by_id: HashMap<&str, (f64, f64)> =
        HashMap::with_capacity(solved.graph.edge_count());
    for edge_ref in solved.graph.edge_references() {
        let edge_index = edge_ref.id();
        let Some((from_index, _)) = solved.graph.edge_endpoints(edge_index) else {
            continue;
        };
        let from_pressure_pa = solved
            .pressures
            .get(from_index.index())
            .copied()
            .unwrap_or_default();
        let flow_m3_s = solved
            .flow_rates
            .get(edge_index.index())
            .copied()
            .unwrap_or(edge_ref.weight().flow_rate);
        edge_by_id.insert(edge_ref.weight().id.as_str(), (flow_m3_s, from_pressure_pa));
    }

    let treatment_ids: HashSet<String> = candidate.treatment_channel_ids().into_iter().collect();
    let mut channel_samples = Vec::with_capacity(candidate.blueprint.channels.len());
    let mut total_volume_m3 = 0.0_f64;
    let mut venturi_flows = Vec::new();
    for channel in &candidate.blueprint.channels {
        let (flow_m3_s, from_pressure_pa) = edge_by_id
            .get(channel.id.as_str())
            .copied()
            .unwrap_or((0.0, 0.0));
        let is_venturi_channel = channel.venturi_geometry.is_some()
            || channel.visual_role == Some(ChannelVisualRole::VenturiThroat)
            || channel.metadata.as_ref().is_some_and(|metadata| {
                metadata.contains::<VenturiGeometryMetadata>()
                    || metadata
                        .get::<ChannelPathMetadata>()
                        .is_some_and(|path| path.visual_role == ChannelVisualRole::VenturiThroat)
            });
        let venturi_throat_count = channel
            .metadata
            .as_ref()
            .and_then(|metadata| metadata.get::<ChannelVenturiSpec>())
            .map_or(u8::from(is_venturi_channel), |spec| spec.n_throats.max(1));
        let is_treatment_channel = channel.therapy_zone == Some(TherapyZone::CancerTarget)
            || treatment_ids.contains(channel.id.as_str())
            || treatment_ids
                .iter()
                .any(|target| channel.id.as_str().starts_with(target));
        if is_venturi_channel {
            venturi_flows.push(flow_m3_s.abs());
        }
        total_volume_m3 += channel.length_m * channel.cross_section.area();
        channel_samples.push(BlueprintSolveSample {
            id: channel.id.as_str(),
            from_node: channel.from.as_str(),
            to_node: channel.to.as_str(),
            length_m: channel.length_m,
            cross_section: channel.cross_section,
            channel_shape: channel.channel_shape,
            flow_m3_s,
            from_pressure_pa,
            is_treatment_channel,
            is_venturi_channel,
            venturi_throat_count,
        });
    }

    let venturi_flow_fraction = (venturi_flows.iter().sum::<f64>()
        / candidate.operating_point.flow_rate_m3_s.max(1.0e-12))
    .clamp(0.0, 1.0);
    let flow_uniformity = compute_blueprint_flow_uniformity(&candidate.blueprint, &channel_samples);
    let mean_residence_time_s =
        total_volume_m3 / candidate.operating_point.flow_rate_m3_s.max(1.0e-12);
    let remerge_loss_pa = compute_blueprint_remerge_loss(&candidate.blueprint, &channel_samples);

    Ok(BlueprintSolveSummary {
        inlet_pressure_pa,
        inlet_flow_m3_s: candidate.operating_point.flow_rate_m3_s,
        flow_uniformity,
        venturi_flow_fraction,
        mean_residence_time_s,
        solve_path_status,
        remerge_loss_pa,
        channel_samples,
    })
}

fn compute_blueprint_flow_uniformity(
    blueprint: &NetworkBlueprint,
    samples: &[BlueprintSolveSample<'_>],
) -> f64 {
    let selected_ids = blueprint.topology_spec().map_or_else(Vec::new, |topology| {
        if !topology.venturi_placements.is_empty() {
            topology
                .venturi_placements
                .iter()
                .map(|placement| placement.target_channel_id.clone())
                .collect()
        } else if let Some(last_stage) = topology.split_stages.last() {
            last_stage
                .branches
                .iter()
                .map(|branch| {
                    cfd_schematics::BlueprintTopologySpec::branch_channel_id(
                        &last_stage.stage_id,
                        &branch.label,
                    )
                })
                .collect()
        } else {
            Vec::new()
        }
    });

    let selected_flows = if selected_ids.is_empty() {
        let from_nodes: HashSet<&str> = samples.iter().map(|sample| sample.from_node).collect();
        samples
            .iter()
            .filter(|sample| !from_nodes.contains(sample.to_node))
            .map(|sample| sample.flow_m3_s.abs())
            .collect::<Vec<_>>()
    } else {
        samples
            .iter()
            .filter(|sample| {
                selected_ids
                    .iter()
                    .any(|target| sample.id == target || sample.id.starts_with(target))
            })
            .map(|sample| sample.flow_m3_s.abs())
            .collect::<Vec<_>>()
    };

    if selected_flows.len() <= 1 {
        return 1.0;
    }
    let mean = selected_flows.iter().sum::<f64>() / selected_flows.len() as f64;
    if mean <= 1.0e-15 {
        return 0.0;
    }
    let variance = selected_flows
        .iter()
        .map(|flow| (flow - mean).powi(2))
        .sum::<f64>()
        / selected_flows.len() as f64;
    (1.0 / (1.0 + variance.sqrt() / mean)).clamp(0.0, 1.0)
}

fn compute_blueprint_remerge_loss(
    blueprint: &NetworkBlueprint,
    samples: &[BlueprintSolveSample<'_>],
) -> f64 {
    let mut incoming: HashMap<&str, (usize, f64)> = HashMap::new();
    let mut outgoing_area: HashMap<&str, f64> = HashMap::new();
    for sample in samples {
        let entry = incoming.entry(sample.to_node).or_insert((0, 0.0));
        entry.0 += 1;
        entry.1 += sample.flow_m3_s.abs();
        let out_entry = outgoing_area.entry(sample.from_node).or_insert(0.0);
        *out_entry = out_entry.max(sample.cross_section.area());
    }

    let mut total_loss_pa = 0.0_f64;
    for node in &blueprint.nodes {
        let Some((incoming_count, q_combined)) = incoming.get(node.id.as_str()).copied() else {
            continue;
        };
        if incoming_count < 2 || q_combined <= 1.0e-18 {
            continue;
        }
        let k_merge = match node
            .junction_geometry
            .as_ref()
            .map(|geometry| geometry.junction_family)
        {
            Some(JunctionFamily::Merge) => 1.0,
            Some(JunctionFamily::Tee) => 1.5,
            Some(JunctionFamily::Bifurcation) => 1.0,
            Some(JunctionFamily::Trifurcation) => 0.9,
            Some(JunctionFamily::Cross) => 1.2,
            None => 1.0,
        };
        let a_out = outgoing_area
            .get(node.id.as_str())
            .copied()
            .filter(|area| *area > 1.0e-18)
            .unwrap_or_else(|| {
                samples
                    .iter()
                    .filter(|sample| sample.to_node == node.id.as_str())
                    .map(|sample| sample.cross_section.area())
                    .sum::<f64>()
                    / incoming_count as f64
            })
            .max(1.0e-18);
        let velocity_m_s = q_combined / a_out;
        total_loss_pa += k_merge * 0.5 * BLOOD_DENSITY_KG_M3 * velocity_m_s * velocity_m_s;
    }
    total_loss_pa
}
