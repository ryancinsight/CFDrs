use std::collections::HashMap;

use cfd_1d::domain::network::{
    apply_blueprint_boundary_conditions, network_from_blueprint,
};
use cfd_1d::{
    NetworkProblem, NetworkSolver, PrimarySolveDiagnostics, SolvePathStatus, SolverConfig,
};
use cfd_1d::BoundaryCondition;
use cfd_core::error::{Error, Result as CfdResult};
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_schematics::domain::model::{NetworkBlueprint, NodeKind};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use petgraph::graph::NodeIndex;

/// Per-node diagnostic data from the authoritative cfd-1d reference solve.
#[derive(Debug, Clone)]
pub struct NodeReferenceTrace<T> {
    /// Blueprint node identifier.
    pub node_id: String,
    /// Blueprint node classification.
    pub node_kind: NodeKind,
    /// Scaled nodal pressure from the reference solve [Pa].
    pub pressure_pa: T,
    /// Sum of incoming channel flow rates [m^3/s].
    pub incoming_flow_m3_s: T,
    /// Sum of outgoing channel flow rates [m^3/s].
    pub outgoing_flow_m3_s: T,
    /// Net boundary source or sink applied at the node [m^3/s].
    pub prescribed_boundary_flow_m3_s: T,
    /// Continuity residual after accounting for internal and boundary flows [m^3/s].
    pub continuity_residual_m3_s: T,
}

/// Per-channel diagnostic data from the authoritative cfd-1d reference solve.
#[derive(Debug, Clone)]
pub struct ChannelReferenceTrace<T> {
    /// Blueprint channel identifier.
    pub channel_id: String,
    /// Upstream blueprint node identifier.
    pub from_node_id: String,
    /// Downstream blueprint node identifier.
    pub to_node_id: String,
    /// Scaled volumetric flow rate [m^3/s].
    pub flow_rate_m3_s: T,
    /// Scaled channel pressure drop [Pa].
    pub pressure_drop_pa: T,
    /// Mean velocity based on the true blueprint area [m/s].
    pub mean_velocity_m_s: T,
    /// Hydraulic resistance carried by the blueprint channel contract [Pa·s/m^3].
    pub resistance_pa_s_per_m3: T,
    /// True blueprint cross-sectional area [m^2].
    pub cross_section_area_m2: T,
    /// Blueprint channel length [m].
    pub length_m: T,
}

/// Complete cfd-1d reference trace used to configure a cfd-2d network build.
#[derive(Debug, Clone)]
pub struct NetworkReferenceTrace<T> {
    /// Solve-path classification reported by the reference solver.
    pub solve_path_status: SolvePathStatus,
    /// Primary-path diagnostics from the cfd-1d solve.
    pub diagnostics: PrimarySolveDiagnostics,
    /// Normalized inlet pressure used before scaling [Pa].
    pub normalized_inlet_pressure_pa: T,
    /// Inlet pressure after scaling to the requested total flow [Pa].
    pub scaled_inlet_pressure_pa: T,
    /// Pressure multiplier applied to the normalized solve.
    pub pressure_scale_factor: T,
    /// Total inlet flow in the scaled reference solve [m^3/s].
    pub total_inlet_flow_m3_s: T,
    /// Total outlet flow in the scaled reference solve [m^3/s].
    pub total_outlet_flow_m3_s: T,
    /// Per-node continuity and pressure diagnostics.
    pub node_traces: Vec<NodeReferenceTrace<T>>,
    /// Per-channel flow and pressure diagnostics.
    pub channel_traces: Vec<ChannelReferenceTrace<T>>,
}

/// Solve a normalized cfd-1d reference network and scale it to the requested total flow.
pub fn solve_reference_trace<T>(
    blueprint: &NetworkBlueprint,
    density_kg_m3: f64,
    viscosity_pa_s: f64,
    target_total_flow_m3_s: f64,
) -> CfdResult<NetworkReferenceTrace<T>>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
{
    let fluid = reference_fluid::<T>(density_kg_m3, viscosity_pa_s)?;
    let mut network = network_from_blueprint::<T, _>(blueprint, fluid).map_err(|e| {
        Error::InvalidInput(format!(
            "Network2DSolver reference trace failed to build cfd-1d network: {e}"
        ))
    })?;

    let node_indices: HashMap<String, NodeIndex> = network
        .graph
        .node_indices()
        .filter_map(|idx| {
            network
                .graph
                .node_weight(idx)
                .map(|node| (node.id.clone(), idx))
        })
        .collect();

    let normalized_inlet_pressure_pa = T::one();
    apply_blueprint_boundary_conditions(
        &mut network,
        blueprint,
        &node_indices,
        normalized_inlet_pressure_pa,
        T::zero(),
    )?;

    let solver = NetworkSolver::<T, ConstantPropertyFluid<T>>::with_config(SolverConfig::<T> {
        tolerance: T::from_f64(1e-12).unwrap_or_else(T::default_epsilon),
        max_iterations: 500,
    });
    let (solved, diagnostics) = solver
        .solve_network_with_diagnostics(&NetworkProblem::<T, ConstantPropertyFluid<T>>::new(
            network,
        ))
        .map_err(|e| {
            Error::InvalidInput(format!(
                "Network2DSolver reference trace failed during cfd-1d solve: {e}"
            ))
        })?;

    build_reference_trace_from_solved_network(
        blueprint,
        &solved,
        diagnostics,
        normalized_inlet_pressure_pa,
        target_total_flow_m3_s,
    )
}

pub(crate) fn build_reference_trace_from_solved_network<T>(
    blueprint: &NetworkBlueprint,
    solved: &cfd_1d::domain::network::Network<T, ConstantPropertyFluid<T>>,
    diagnostics: PrimarySolveDiagnostics,
    normalized_inlet_pressure_pa: T,
    target_total_flow_m3_s: f64,
) -> CfdResult<NetworkReferenceTrace<T>>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
{
    let edge_flow_by_id: HashMap<String, T> = solved
        .edges_with_properties()
        .into_iter()
        .map(|edge| (edge.id, edge.flow_rate))
        .collect();

    let solved_pressure_by_id: HashMap<String, T> = solved
        .graph
        .node_indices()
        .filter_map(|idx| {
            solved.graph.node_weight(idx).map(|node| {
                (
                    node.id.clone(),
                    *solved.pressures().get(idx.index()).unwrap_or(&T::zero()),
                )
            })
        })
        .collect();
    let node_index_by_id: HashMap<&str, NodeIndex> = solved
        .graph
        .node_indices()
        .filter_map(|idx| solved.graph.node_weight(idx).map(|node| (node.id.as_str(), idx)))
        .collect();
    let boundary_conditions = solved.boundary_conditions();
    let scale_factor = T::one();

    let mut min_edge_resistance = f64::INFINITY;
    let mut max_edge_resistance = 0.0_f64;
    for edge in solved.edges_with_properties() {
        let resistance = edge.properties.resistance.to_f64().unwrap_or(0.0).abs();
        min_edge_resistance = min_edge_resistance.min(resistance);
        max_edge_resistance = max_edge_resistance.max(resistance);
    }

    let mut channel_traces = Vec::with_capacity(blueprint.channels.len());
    let mut incoming_by_node: HashMap<String, T> = HashMap::with_capacity(blueprint.nodes.len());
    let mut outgoing_by_node: HashMap<String, T> = HashMap::with_capacity(blueprint.nodes.len());
    for node in &blueprint.nodes {
        incoming_by_node.insert(node.id.as_str().to_owned(), T::zero());
        outgoing_by_node.insert(node.id.as_str().to_owned(), T::zero());
    }

    for channel in &blueprint.channels {
        let flow = edge_flow_by_id
            .get(channel.id.as_str())
            .copied()
            .expect("analytical constant conversion");
        let p_from = solved_pressure_by_id
            .get(channel.from.as_str())
            .copied()
            .expect("analytical constant conversion")
            * scale_factor;
        let p_to = solved_pressure_by_id
            .get(channel.to.as_str())
            .copied()
            .expect("analytical constant conversion")
            * scale_factor;
        let area = T::from_f64(channel.cross_section.area()).ok_or_else(|| {
            Error::InvalidInput(format!(
                "Network2DSolver reference trace area conversion failed for '{}'",
                channel.id.as_str()
            ))
        })?;
        let length = T::from_f64(channel.length_m).ok_or_else(|| {
            Error::InvalidInput(format!(
                "Network2DSolver reference trace length conversion failed for '{}'",
                channel.id.as_str()
            ))
        })?;
        let resistance = T::from_f64(channel.resistance).expect("analytical constant conversion");
        let mean_velocity = if area > T::zero() {
            flow / area
        } else {
            T::zero()
        };

        *outgoing_by_node
            .entry(channel.from.as_str().to_owned())
            .or_insert_with(T::zero) += flow;
        *incoming_by_node
            .entry(channel.to.as_str().to_owned())
            .or_insert_with(T::zero) += flow;

        channel_traces.push(ChannelReferenceTrace {
            channel_id: channel.id.as_str().to_owned(),
            from_node_id: channel.from.as_str().to_owned(),
            to_node_id: channel.to.as_str().to_owned(),
            flow_rate_m3_s: flow,
            pressure_drop_pa: p_from - p_to,
            mean_velocity_m_s: mean_velocity,
            resistance_pa_s_per_m3: resistance,
            cross_section_area_m2: area,
            length_m: length,
        });
    }

    let mut normalized_inlet_flow = T::zero();
    let mut normalized_outlet_flow = T::zero();
    for node in &blueprint.nodes {
        let incoming = incoming_by_node
            .get(node.id.as_str())
            .copied()
            .expect("analytical constant conversion")
            * scale_factor;
        let outgoing = outgoing_by_node
            .get(node.id.as_str())
            .copied()
            .expect("analytical constant conversion")
            * scale_factor;

        match node.kind {
            NodeKind::Inlet => normalized_inlet_flow += outgoing - incoming,
            NodeKind::Outlet => normalized_outlet_flow += incoming - outgoing,
            NodeKind::Reservoir | NodeKind::Junction => {}
        }
    }

    let normalized_inlet_flow_f64 = normalized_inlet_flow.to_f64().unwrap_or(0.0).abs();
    if normalized_inlet_flow_f64 <= 1e-30 {
        let boundary_summary = blueprint
            .nodes
            .iter()
            .filter(|node| matches!(node.kind, NodeKind::Inlet | NodeKind::Outlet))
            .map(|node| {
                let incoming = incoming_by_node
                    .get(node.id.as_str())
                    .copied()
                    .unwrap_or_else(T::zero)
                    .to_f64()
                    .unwrap_or(0.0);
                let outgoing = outgoing_by_node
                    .get(node.id.as_str())
                    .copied()
                    .unwrap_or_else(T::zero)
                    .to_f64()
                    .unwrap_or(0.0);
                format!(
                    "{}({:?}):in={incoming:.3e},out={outgoing:.3e}",
                    node.id.as_str(),
                    node.kind
                )
            })
            .collect::<Vec<_>>()
            .join("; ");
        return Err(Error::InvalidInput(
            format!(
                "Network2DSolver reference trace produced zero inlet flow from the cfd-1d solve (normalized_inlet_flow={normalized_inlet_flow_f64:.3e}, normalized_outlet_flow={:.3e}, edge_resistance_range=[{min_edge_resistance:.3e}, {max_edge_resistance:.3e}], boundary_nodes=[{}])",
                normalized_outlet_flow.to_f64().unwrap_or(0.0).abs(),
                boundary_summary
            ),
        ));
    }

    let scale_factor =
        T::from_f64(target_total_flow_m3_s / normalized_inlet_flow_f64).ok_or_else(|| {
            Error::InvalidInput(
                "Network2DSolver reference trace scale conversion failed".to_string(),
            )
        })?;

    for channel_trace in &mut channel_traces {
        channel_trace.flow_rate_m3_s *= scale_factor;
        channel_trace.pressure_drop_pa *= scale_factor;
        channel_trace.mean_velocity_m_s *= scale_factor;
    }

    let mut node_traces = Vec::with_capacity(blueprint.nodes.len());
    for node in &blueprint.nodes {
        let node_idx = node_index_by_id
            .get(node.id.as_str())
            .copied()
            .expect("analytical constant conversion");
        let incoming = incoming_by_node
            .get(node.id.as_str())
            .copied()
            .expect("analytical constant conversion")
            * scale_factor;
        let outgoing = outgoing_by_node
            .get(node.id.as_str())
            .copied()
            .expect("analytical constant conversion")
            * scale_factor;
        let pressure = solved_pressure_by_id
            .get(node.id.as_str())
            .copied()
            .expect("analytical constant conversion")
            * scale_factor;
        let boundary_condition = boundary_conditions.get(&node_idx);
        let prescribed_boundary_flow = match boundary_condition {
            Some(BoundaryCondition::Neumann { gradient }) => *gradient * scale_factor,
            _ => match node.kind {
                NodeKind::Inlet => outgoing - incoming,
                NodeKind::Outlet => incoming - outgoing,
                NodeKind::Reservoir | NodeKind::Junction => T::zero(),
            },
        };
        let continuity_residual = match boundary_condition {
            Some(BoundaryCondition::Neumann { gradient }) => {
                incoming + *gradient * scale_factor - outgoing
            }
            _ => match node.kind {
                NodeKind::Inlet => incoming + prescribed_boundary_flow - outgoing,
                NodeKind::Outlet => incoming - outgoing - prescribed_boundary_flow,
                NodeKind::Reservoir | NodeKind::Junction => incoming - outgoing,
            },
        };

        node_traces.push(NodeReferenceTrace {
            node_id: node.id.as_str().to_owned(),
            node_kind: node.kind,
            pressure_pa: pressure,
            incoming_flow_m3_s: incoming,
            outgoing_flow_m3_s: outgoing,
            prescribed_boundary_flow_m3_s: prescribed_boundary_flow,
            continuity_residual_m3_s: continuity_residual,
        });
    }

    Ok(NetworkReferenceTrace {
        solve_path_status: SolvePathStatus::PrimaryConverged,
        diagnostics,
        normalized_inlet_pressure_pa,
        scaled_inlet_pressure_pa: normalized_inlet_pressure_pa * scale_factor,
        pressure_scale_factor: scale_factor,
        total_inlet_flow_m3_s: normalized_inlet_flow * scale_factor,
        total_outlet_flow_m3_s: normalized_outlet_flow * scale_factor,
        node_traces,
        channel_traces,
    })
}

pub(crate) fn reference_fluid<T>(
    density_kg_m3: f64,
    viscosity_pa_s: f64,
) -> CfdResult<ConstantPropertyFluid<T>>
where
    T: RealField + Copy + Float + FromPrimitive,
{
    let fluid = ConstantPropertyFluid::new(
        "Blood reference fluid".to_string(),
        T::from_f64(density_kg_m3).ok_or_else(|| {
            Error::InvalidInput("Network2DSolver reference density conversion failed".to_string())
        })?,
        T::from_f64(viscosity_pa_s).ok_or_else(|| {
            Error::InvalidInput("Network2DSolver reference viscosity conversion failed".to_string())
        })?,
        T::from_f64(3617.0).expect("analytical constant conversion"),
        T::from_f64(0.52).expect("analytical constant conversion"),
        T::from_f64(1570.0).expect("analytical constant conversion"),
    );
    fluid.validate()?;
    Ok(fluid)
}
