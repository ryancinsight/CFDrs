//! Solved-network extraction utilities for `compute_metrics`.

use std::collections::HashMap;

use cfd_1d::{network::network_from_blueprint, NetworkProblem, NetworkSolver};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;
use cfd_schematics::domain::model::CrossSectionSpec;
use petgraph::visit::EdgeRef;

use crate::design::DesignTopology;
use crate::error::OptimError;

/// Per-channel solved-flow sample with local node pressures.
#[derive(Debug, Clone)]
pub(super) struct ChannelSolveSample {
    pub id: String,
    pub from_node: String,
    pub to_node: String,
    pub length_m: f64,
    pub cross_section: CrossSectionSpec,
    pub is_venturi_throat: bool,
    pub flow_m3_s: f64,
    pub from_pressure_pa: f64,
}

/// Solved-network summary used by metric computation.
#[derive(Debug, Clone)]
pub(super) struct NetworkSolveSummary {
    pub inlet_pressure_pa: f64,
    pub inlet_flow_m3_s: f64,
    pub flow_uniformity: f64,
    pub venturi_flow_fraction: f64,
    pub mean_residence_time_s: f64,
    /// Estimated additional pressure loss at merge junctions where branched
    /// flows recombine pre-outlet [Pa]. Uses momentum-conserving T/Y-junction
    /// model: ΔP_merge = K_merge · ½ρv²_combined, K_merge ∈ [0.5, 1.5].
    pub remerge_loss_pa: f64,
    pub channel_samples: Vec<ChannelSolveSample>,
}

pub(super) fn solve_blueprint_network(
    candidate_id: &str,
    topology: DesignTopology,
    blueprint: &cfd_schematics::NetworkBlueprint,
    blood: &CassonBlood<f64>,
    inlet_flow_m3_s: f64,
) -> Result<NetworkSolveSummary, OptimError> {
    let mut network =
        network_from_blueprint(blueprint, *blood).map_err(|e| OptimError::PhysicsError {
            id: candidate_id.to_string(),
            reason: format!("network_from_blueprint failed: {e}"),
        })?;

    let mut inlet_nodes = Vec::new();
    let mut outlet_nodes = Vec::new();
    for idx in network.graph.node_indices() {
        if let Some(node) = network.graph.node_weight(idx) {
            if node.node_type == cfd_1d::NodeType::Inlet {
                inlet_nodes.push(idx);
            } else if node.node_type == cfd_1d::NodeType::Outlet {
                outlet_nodes.push(idx);
            }
        }
    }

    if inlet_nodes.is_empty() || outlet_nodes.is_empty() {
        return Err(OptimError::PhysicsError {
            id: candidate_id.to_string(),
            reason: "network is missing inlet/outlet nodes".to_string(),
        });
    }

    let q_per_inlet = inlet_flow_m3_s / inlet_nodes.len() as f64;
    for inlet in &inlet_nodes {
        network.set_neumann_flow(*inlet, q_per_inlet);
    }
    for outlet in &outlet_nodes {
        network.set_pressure(*outlet, 0.0);
    }

    let solver = NetworkSolver::<f64, CassonBlood<f64>>::new();
    let primary_problem = NetworkProblem::new(network.clone());
    let solved = match solver.solve_network(&primary_problem) {
        Ok(solved) => solved,
        Err(primary_err) => {
            // Fallback for short-channel blueprints where dynamic resistance
            // updates can violate entrance-length assumptions.
            //
            // Strategy:
            // 1) freeze geometry-dependent updates and drop quadratic terms
            // 2) solve a linear Dirichlet-normalized problem (Pin=1, Pout=0)
            // 3) scale the solved flow/pressure field to the target inlet flow
            let mut linear_network = network.clone();
            for props in linear_network.properties.values_mut() {
                props.geometry = None;
            }
            for edge_idx in linear_network.graph.edge_indices() {
                if let Some(edge) = linear_network.graph.edge_weight_mut(edge_idx) {
                    edge.quad_coeff = 0.0;
                    edge.resistance = edge.resistance.max(1e-12);
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
                    .map_err(|fallback_err| OptimError::PhysicsError {
                        id: candidate_id.to_string(),
                        reason: format!(
                            "NetworkSolver failed: {primary_err}; fallback failed: {fallback_err}"
                        ),
                    })?;

            let mut inlet_flow_unit = 0.0_f64;
            for edge_ref in fallback_solved.graph.edge_references() {
                let edge_idx = edge_ref.id();
                let Some((from_idx, _)) = fallback_solved.graph.edge_endpoints(edge_idx) else {
                    continue;
                };
                if inlet_nodes.contains(&from_idx) {
                    inlet_flow_unit += fallback_solved
                        .flow_rates
                        .get(&edge_idx)
                        .copied()
                        .unwrap_or(edge_ref.weight().flow_rate)
                        .abs();
                }
            }
            if inlet_flow_unit <= 1e-18 {
                return Err(OptimError::PhysicsError {
                    id: candidate_id.to_string(),
                    reason: format!(
                        "NetworkSolver failed: {primary_err}; fallback produced zero inlet flow"
                    ),
                });
            }

            let scale = inlet_flow_m3_s / inlet_flow_unit;
            for pressure in fallback_solved.pressures.values_mut() {
                *pressure *= scale;
            }
            let edge_indices: Vec<_> = fallback_solved.graph.edge_indices().collect();
            for edge_idx in edge_indices {
                if let Some(q) = fallback_solved.flow_rates.get_mut(&edge_idx) {
                    *q *= scale;
                    if let Some(edge) = fallback_solved.graph.edge_weight_mut(edge_idx) {
                        edge.flow_rate = *q;
                    }
                } else if let Some(edge) = fallback_solved.graph.edge_weight_mut(edge_idx) {
                    edge.flow_rate *= scale;
                }
            }

            fallback_solved
        }
    };

    let inlet_pressure_pa = inlet_nodes
        .iter()
        .filter_map(|idx| solved.pressures.get(idx).copied())
        .sum::<f64>()
        / inlet_nodes.len() as f64;

    let mut edge_by_id: HashMap<String, (f64, f64)> = HashMap::new();
    for edge_ref in solved.graph.edge_references() {
        let edge_idx = edge_ref.id();
        let edge = edge_ref.weight();
        let (from_idx, _) =
            solved
                .graph
                .edge_endpoints(edge_idx)
                .ok_or_else(|| OptimError::PhysicsError {
                    id: candidate_id.to_string(),
                    reason: "missing solved edge endpoints".to_string(),
                })?;
        let from_p = solved.pressures.get(&from_idx).copied().unwrap_or(0.0);
        let q = solved
            .flow_rates
            .get(&edge_idx)
            .copied()
            .unwrap_or(edge.flow_rate);
        edge_by_id.insert(edge.id.clone(), (q, from_p));
    }

    let mut channel_samples = Vec::with_capacity(blueprint.channels.len());
    let mut total_volume_m3 = 0.0_f64;
    let mut throat_flows = Vec::new();

    for ch in &blueprint.channels {
        let (q, from_p) = edge_by_id
            .get(ch.id.as_str())
            .copied()
            .unwrap_or((0.0, 0.0));
        let is_venturi_throat = ch.metadata.as_ref().is_some_and(|meta| {
            meta.contains::<VenturiGeometryMetadata>()
        }) || is_venturi_throat_id(ch.id.as_str());

        if is_venturi_throat {
            throat_flows.push(q.abs());
        }
        total_volume_m3 += ch.length_m * cross_section_area(ch.cross_section);
        channel_samples.push(ChannelSolveSample {
            id: ch.id.as_str().to_string(),
            from_node: ch.from.as_str().to_string(),
            to_node: ch.to.as_str().to_string(),
            length_m: ch.length_m,
            cross_section: ch.cross_section,
            is_venturi_throat,
            flow_m3_s: q,
            from_pressure_pa: from_p,
        });
    }

    let venturi_flow_fraction =
        (throat_flows.iter().sum::<f64>() / inlet_flow_m3_s.max(1e-12)).clamp(0.0, 1.0);
    let flow_uniformity = compute_flow_uniformity(topology, &channel_samples);
    let mean_residence_time_s = total_volume_m3 / inlet_flow_m3_s.max(1e-12);

    // Estimate remerge junction pressure loss for branching topologies.
    // When N branches merge into a single outlet trunk, the momentum-conserving
    // loss is K_merge · ½ρv²_combined per merge level.
    // K_merge = 1.0 for symmetric T-junction (Idelchik 1994).
    let remerge_loss_pa = compute_remerge_loss(topology, &channel_samples, inlet_flow_m3_s);

    Ok(NetworkSolveSummary {
        inlet_pressure_pa,
        inlet_flow_m3_s,
        flow_uniformity,
        venturi_flow_fraction,
        mean_residence_time_s,
        remerge_loss_pa,
        channel_samples,
    })
}

fn is_venturi_throat_id(id: &str) -> bool {
    id.as_bytes()
        .windows("throat".len())
        .any(|window| window.eq_ignore_ascii_case(b"throat"))
}

fn cross_section_area(cs: CrossSectionSpec) -> f64 {
    match cs {
        CrossSectionSpec::Rectangular { width_m, height_m } => width_m * height_m,
        CrossSectionSpec::Circular { diameter_m } => {
            std::f64::consts::PI * 0.25 * diameter_m * diameter_m
        }
    }
}

fn compute_flow_uniformity(topology: DesignTopology, samples: &[ChannelSolveSample]) -> f64 {
    let mut selected = Vec::new();

    if topology.has_venturi() {
        for sample in samples {
            if sample.is_venturi_throat {
                selected.push(sample.flow_m3_s.abs());
            }
        }
    } else {
        for sample in samples {
            let id = sample.id.as_str();
            if id == "segment_1"
                || id == "arm2_seg_1"
                || id == "arm3_seg_1"
                || id == "wide_seg_1"
                || id == "narrow_seg_1"
                || id.starts_with("ch_")
            {
                selected.push(sample.flow_m3_s.abs());
            }
        }
    }

    if selected.len() <= 1 {
        return 1.0;
    }
    let mean = selected.iter().sum::<f64>() / selected.len() as f64;
    if mean <= 1e-15 {
        return 0.0;
    }
    let variance = selected
        .iter()
        .map(|v| (v - mean) * (v - mean))
        .sum::<f64>()
        / selected.len() as f64;
    let cv = variance.sqrt() / mean;
    (1.0 / (1.0 + cv)).clamp(0.0, 1.0)
}

/// Estimate additional pressure loss at merge junctions where branching
/// topologies recombine flow before the outlet.
///
/// Uses the Idelchik (1994) momentum-conserving merge model:
///   ΔP = K · ½ρv²_combined per merge level
///
/// K values are topology-dependent:
/// - **T-junction merge** (90° perpendicular branch): K = 1.5 (Idelchik §7.23)
/// - **Y-junction merge** (30–60° angled branch):     K = 0.9 (Idelchik §7.27)
/// - **Symmetric T-merge** (equal branches):           K = 1.0 (Idelchik §7.20)
///
/// For CCT/CIF cascaded topologies each merge level is a T-junction where the
/// center arm re-joins the main trunk at ~90°.
fn compute_remerge_loss(
    topology: DesignTopology,
    samples: &[ChannelSolveSample],
    inlet_flow_m3_s: f64,
) -> f64 {
    use crate::constraints::BLOOD_DENSITY_KG_M3;

    // (n_merges, k_per_merge) from topology geometry.
    let (n_merges, k_merge): (usize, f64) = match topology {
        // CCT cascades: each level has a T-junction merge on the return path.
        DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => (n_levels as usize, 1.5),
        // CIF staged: pre-tri stages + terminal tri + terminal bi → all T-merges.
        DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => {
            (n_pretri as usize + 2, 1.5)
        }
        // Symmetric bifurcation merges: equal-flow T-merge (K ≈ 1.0).
        DesignTopology::BifurcationVenturi
        | DesignTopology::BifurcationSerpentine
        | DesignTopology::AsymmetricBifurcationSerpentine => (1, 1.0),
        // Trifurcation families: two merges per level, Y-junction geometry.
        DesignTopology::TrifurcationSerpentine
        | DesignTopology::TrifurcationBifurcationVenturi
        | DesignTopology::TripleTrifurcationVenturi
        | DesignTopology::TrifurcationBifurcationBifurcationVenturi
        | DesignTopology::QuadTrifurcationVenturi => (2, 0.9),
        _ => (0, 0.0),
    };

    if n_merges == 0 {
        return 0.0;
    }

    // Find the outlet trunk area (last channel before outlet).
    // Fall back to the largest non-throat channel area.
    let outlet_area = samples
        .iter()
        .filter(|s| !s.is_venturi_throat)
        .map(|s| cross_section_area(s.cross_section))
        .fold(0.0_f64, f64::max)
        .max(1e-18);

    let v_outlet = inlet_flow_m3_s / outlet_area;
    let dp_per_merge = k_merge * 0.5 * BLOOD_DENSITY_KG_M3 * v_outlet * v_outlet;

    n_merges as f64 * dp_per_merge
}
