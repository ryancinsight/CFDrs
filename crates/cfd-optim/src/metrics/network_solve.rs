//! Solved-network extraction utilities for `compute_metrics`.

use std::collections::HashMap;

use cfd_1d::{network::network_from_blueprint, NetworkProblem, NetworkSolver};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_schematics::domain::model::CrossSectionSpec;
use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;
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
        let is_venturi_throat = ch
            .metadata
            .as_ref()
            .is_some_and(|meta| meta.contains::<VenturiGeometryMetadata>())
            || is_venturi_throat_id(ch.id.as_str());

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
/// # Per-Junction Physics
///
/// For each identified merge node (a node where ≥2 channels converge), the
/// Idelchik (1994) momentum-conserving merge loss is computed from the **actual
/// solved flow velocities** at that node:
///
/// ```text
/// ΔP_j = K_j · ½ρv²_combined_at_j
/// ```
///
/// where:
/// - `v_combined = ΣQ_incoming / A_out` is the velocity in the run channel
///   downstream of the merge (Idelchik §7.20, combining flow reference branch)
/// - `K_j` is looked up by junction type:
///   - T-junction combining (90°): K = 1.5 (Idelchik §7.23)
///   - Y-junction combining (30–60°): K = 0.9 (Idelchik §7.27)
///   - Symmetric T-merge: K = 1.0 (Idelchik §7.20)
///
/// ## Fallback
///
/// When merge junctions cannot be resolved from the solved samples (e.g. simple
/// in-out topologies), the aggregate fixed-K estimate is used as a lower bound.
fn compute_remerge_loss(
    topology: DesignTopology,
    samples: &[ChannelSolveSample],
    inlet_flow_m3_s: f64,
) -> f64 {
    use crate::constraints::BLOOD_DENSITY_KG_M3;
    use std::collections::HashMap;

    // ── Step 1: Build a node → (list of incoming Q, list of outgoing Q) map ──
    // Each channel in `samples` goes from `from_node` to `to_node`.
    // A merge node is one with ≥2 incoming channels.
    let mut incoming: HashMap<&str, Vec<f64>> = HashMap::new();
    let mut outgoing_area: HashMap<&str, f64> = HashMap::new();

    for s in samples {
        incoming
            .entry(s.to_node.as_str())
            .or_default()
            .push(s.flow_m3_s.abs());
        // The outgoing area from a node is the area of channels leaving that node.
        // When multiple channels share a from_node, use the largest (trunk) area.
        let a = cross_section_area(s.cross_section);
        let entry = outgoing_area.entry(s.from_node.as_str()).or_insert(0.0);
        if a > *entry {
            *entry = a;
        }
    }

    // ── Step 2: Classify topology K-factor for merge junctions ──────────────
    // K = 1.5 → T-junction combining (CCT/CIF cascades, 90° angle)
    // K = 1.0 → symmetric equal-flow T-merge (bifurcation families)
    // K = 0.9 → Y-junction combining (30–60°, trifurcation families)
    let k_merge: f64 = match topology {
        DesignTopology::CascadeCenterTrifurcationSeparator { .. }
        | DesignTopology::IncrementalFiltrationTriBiSeparator { .. } => 1.5,
        DesignTopology::BifurcationVenturi
        | DesignTopology::BifurcationSerpentine
        | DesignTopology::AsymmetricBifurcationSerpentine => 1.0,
        DesignTopology::TrifurcationSerpentine
        | DesignTopology::TrifurcationBifurcationVenturi
        | DesignTopology::TripleTrifurcationVenturi
        | DesignTopology::TrifurcationBifurcationBifurcationVenturi
        | DesignTopology::QuadTrifurcationVenturi => 0.9,
        _ => 1.0,
    };

    // ── Step 3: Accumulate per-junction ΔP at merge nodes ───────────────────
    let mut total_merge_loss_pa = 0.0_f64;
    let mut merge_count = 0usize;

    for (node_id, q_in_list) in &incoming {
        if q_in_list.len() < 2 {
            // Not a merge node — skip.
            continue;
        }
        // Combined flow entering this merge node [m³/s].
        let q_combined: f64 = q_in_list.iter().sum();
        if q_combined < 1e-18 {
            continue;
        }

        // Area of the outgoing (run) channel from this node.
        // If node is "outlet" or has no outgoing channel, fall back to
        // the mean incoming-channel area.
        let a_out: f64 = outgoing_area
            .get(node_id)
            .copied()
            .filter(|&a| a > 1e-30)
            .unwrap_or_else(|| {
                // Mean area of channels arriving at this node.
                samples
                    .iter()
                    .filter(|s| s.to_node.as_str() == *node_id)
                    .map(|s| cross_section_area(s.cross_section))
                    .sum::<f64>()
                    / q_in_list.len() as f64
            })
            .max(1e-18);

        let v_combined = q_combined / a_out;
        let dp = k_merge * 0.5 * BLOOD_DENSITY_KG_M3 * v_combined * v_combined;
        total_merge_loss_pa += dp;
        merge_count += 1;
    }

    // ── Step 4: Fallback aggregate estimate when no merge nodes found ────────
    if merge_count == 0 {
        // Use the original topology-table estimate as a lower bound.
        let n_merges_fallback: usize = match topology {
            DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => n_levels as usize,
            DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => {
                n_pretri as usize + 2
            }
            DesignTopology::BifurcationVenturi
            | DesignTopology::BifurcationSerpentine
            | DesignTopology::AsymmetricBifurcationSerpentine => 1,
            DesignTopology::TrifurcationSerpentine
            | DesignTopology::TrifurcationBifurcationVenturi
            | DesignTopology::TripleTrifurcationVenturi
            | DesignTopology::TrifurcationBifurcationBifurcationVenturi
            | DesignTopology::QuadTrifurcationVenturi => 2,
            _ => 0,
        };
        if n_merges_fallback == 0 {
            return 0.0;
        }
        let outlet_area = samples
            .iter()
            .filter(|s| !s.is_venturi_throat)
            .map(|s| cross_section_area(s.cross_section))
            .fold(0.0_f64, f64::max)
            .max(1e-18);
        let v_outlet = inlet_flow_m3_s / outlet_area;
        return n_merges_fallback as f64
            * k_merge
            * 0.5
            * BLOOD_DENSITY_KG_M3
            * v_outlet
            * v_outlet;
    }

    total_merge_loss_pa
}
