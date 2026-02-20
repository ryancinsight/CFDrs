//! Physics-based metric computation for SDT millifluidic design candidates.
//!
//! Each [`DesignCandidate`] is evaluated against the 1-D physics models from
//! `cfd-1d` (venturi + serpentine) and the Casson blood rheology from
//! `cfd-core`.  The resulting [`SdtMetrics`] struct gathers all values needed
//! for multi-objective scoring.
//!
//! # Cavitation physics
//! The cavitation number is defined as:
//! ```text
//! σ = (P_inlet − P_vapor) / (½ ρ V_throat²)
//! ```
//! Cavitation onset occurs when **σ < 1**.  In practical millifluidic devices
//! driven at 1–3 bar gauge by syringe pumps, σ < 1 is achievable only with
//! very small throat diameters (50–150 μm) at the higher flow rates (5–10 mL/min).
//!
//! # Haemolysis model
//! The combined haemolysis index (HI) per device pass is computed via the
//! Giersiepen (1990) model summed over each flow segment:
//! ```text
//! ΔHb/Hb = C · t^α · τ_wall^β
//! ```
//! where *t* is the transit time through a segment and *τ_wall* is the local
//! wall shear stress.  A per-pass HI < 0.001 (0.1 %) is the acceptability
//! threshold used here.


use cfd_core::physics::fluid::blood::CassonBlood;
use serde::{Deserialize, Serialize};

use crate::constraints::*;
use crate::design::{DesignCandidate, DesignTopology};
use crate::error::OptimError;

// ── SDT metrics struct ───────────────────────────────────────────────────────

/// All physics-derived metrics for one [`DesignCandidate`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SdtMetrics {
    // ── Cavitation ──
    /// Cavitation number σ at the primary venturi throat.
    /// `< 1` implies active cavitation inception; lower values = more intense.
    /// `f64::INFINITY` when the topology has no venturi.
    pub cavitation_number: f64,

    /// Dimensionless cavitation potential: `max(0, 1 − σ)`.
    /// `0` = no cavitation; approaching `1` = very strong cavitation.
    pub cavitation_potential: f64,

    /// Wall shear rate at the venturi throat [1/s].
    pub throat_shear_rate_inv_s: f64,

    /// Estimated wall shear stress at the venturi throat [Pa].
    /// May greatly exceed the FDA 150 Pa limit; transit time is < 1 μs–1 ms.
    pub throat_shear_pa: f64,

    /// Whether the throat shear stress exceeds FDA 150 Pa.
    /// Expected `true` for effective SDT; cumulative haemolysis is tracked
    /// separately via [`hemolysis_index_per_pass`].
    pub throat_exceeds_fda: bool,

    // ── Main-channel safety ──
    /// Maximum wall shear stress in sustained-contact (non-throat) channels [Pa].
    /// Must be < 150 Pa for FDA compliance.
    pub max_main_channel_shear_pa: f64,

    /// FDA compliance of main distribution / serpentine channels.
    /// `true` iff `max_main_channel_shear_pa < 150 Pa`.
    pub fda_main_compliant: bool,

    // ── Cumulative haemolysis ──
    /// Giersiepen (1990) haemolysis index per device pass (dimensionless fraction).
    /// Sum of throat + main-channel contributions.
    pub hemolysis_index_per_pass: f64,

    // ── Distribution / exposure ──
    /// Flow uniformity across parallel outlet branches: `1 − CV(Q_outlets)`.
    /// `1.0` = perfectly uniform.  Symmetric designs reach 1.0 by construction.
    pub flow_uniformity: f64,

    /// Fraction of the 36 treatment wells covered by channel passes.
    /// `1.0` = all 36 wells are in the channel path.
    pub well_coverage_fraction: f64,

    /// Mean fluid residence time in the treatment zone [s].
    pub mean_residence_time_s: f64,

    // ── System-level ──
    /// Total pressure drop across the device [Pa].
    pub total_pressure_drop_pa: f64,

    /// Total channel path length [mm].
    pub total_path_length_mm: f64,

    /// `true` if the total ΔP < available gauge pressure (syringe pump can
    /// drive the flow without stalling).
    pub pressure_feasible: bool,

    // ── Cell separation ──
    /// Inertial focusing separation efficiency for MCF-7 cancer cells vs RBCs.
    ///
    /// Defined as `|x̃_cancer − x̃_rbc|` ∈ [0, 1] where `x̃` is the normalised
    /// lateral equilibrium position (0 = center, 1 = wall).
    ///
    /// `0.0` for topologies without a cell separation stage.
    pub cell_separation_efficiency: f64,

    /// Fraction of MCF-7 cancer cells collected in the center channel.
    ///
    /// `0.0` for topologies without a cell separation stage.
    pub cancer_center_fraction: f64,

    /// Fraction of RBCs collected in the peripheral (wall) channels.
    ///
    /// `0.0` for topologies without a cell separation stage.
    pub rbc_peripheral_fraction: f64,
}

// ── Public entry point ───────────────────────────────────────────────────────

/// Evaluate all physics-based metrics for a single design candidate.
///
/// Uses [`CassonBlood::normal_blood`] (Casson non-Newtonian blood at 37 °C)
/// and the `cfd-1d` network solver mapped to `VenturiModel` / `SerpentineModel`.
///
/// # Errors
/// Returns [`OptimError::PhysicsError`] if network solving or a 1-D model fails.
pub fn compute_metrics(candidate: &DesignCandidate) -> Result<SdtMetrics, OptimError> {
    use cfd_1d::network::network_from_blueprint;
    use cfd_1d::solver::{NetworkProblem, NetworkSolver};
    use cfd_1d::{FlowConditions, SerpentineModel, VenturiModel};
    use cfd_1d::cell_separation::{CellProperties, CellSeparationModel};

    let bp = generate_network(candidate)?;
    let blood = CassonBlood::<f64>::normal_blood();

    let mut network = network_from_blueprint(&bp, blood.clone()).map_err(|e| {
        OptimError::PhysicsError {
            id: candidate.id.clone(),
            reason: format!("Network build failed: {}", e),
        }
    })?;

    // Apply Boundary Conditions
    let inlet_node = network
        .graph
        .node_indices()
        .find(|i| network.graph.node_weight(*i).unwrap().id == "inlet")
        .ok_or_else(|| OptimError::PhysicsError {
            id: candidate.id.clone(),
            reason: "Missing inlet node".into(),
        })?;

    network.set_neumann_flow(inlet_node, candidate.flow_rate_m3_s);

    for idx in network.graph.node_indices() {
        if network.graph.node_weight(idx).unwrap().node_type
            == cfd_1d::network::NodeType::Outlet
        {
            network.set_pressure(idx, 0.0);
        }
    }

    let solver = NetworkSolver::<f64, cfd_core::physics::fluid::blood::CassonBlood<f64>>::new();
    let problem = NetworkProblem::new(network);

    let solved_network = solver.solve_network(&problem).map_err(|e| {
        OptimError::PhysicsError {
            id: candidate.id.clone(),
            reason: format!("Network solve failed: {}", e),
        }
    })?;

    let mut max_main_shear_pa = 0.0_f64;
    let mut total_hi = 0.0_f64;
    let mut venturi_sigma = f64::INFINITY;
    let mut venturi_shear_rate: f64 = 0.0;
    let mut venturi_shear_pa: f64 = 0.0;
    let mut pressure_deficit: f64 = 0.0;
    
    let mut distribution_path_len_m = 0.0;
    let mut total_path_len_m = 0.0;
    let mut residence_time_s = 0.0;

    for edge in solved_network.edges_with_properties() {
        let q = edge.flow_rate.abs();
        let p_from = solved_network.pressures[&edge.nodes.0];
        let p_to = solved_network.pressures[&edge.nodes.1];
        let dp_graph = (p_from - p_to).abs();
        
        // Accumulate paths (use specific edges that form a single chain from inlet to outlet)
        let is_primary_path = match edge.id.as_str() {
            "trunk" | "b1" | "b_peri" | "venturi_1" | "serpentine" | "trunk_v" | "trunk_s" => true,
            _ => false,
        };

        if edge.id.starts_with("venturi") {
            if candidate.topology.has_venturi() {
                let d_inlet = candidate.inlet_diameter_m;
                let d_throat = candidate.throat_diameter_m;
                let l_throat = candidate.throat_length_m.max(d_throat * 2.0);

                let model = VenturiModel::<f64>::millifluidic(d_inlet, d_throat, l_throat);
                let mut cond = FlowConditions::<f64>::new(0.0);
                cond.flow_rate = Some(q);
                cond.pressure = p_from;

                if let Ok(analysis) = model.analyze(&blood, &cond) {
                    let v_throat = analysis.throat_velocity;
                    let shear_rate = analysis.throat_shear_rate;
                    let shear_pa = blood.apparent_viscosity(shear_rate) * shear_rate;

                    let dynamic_pressure = 0.5 * BLOOD_DENSITY_KG_M3 * v_throat * v_throat;
                    let sigma = if v_throat > 0.0 {
                        (p_from + P_ATM_PA - BLOOD_VAPOR_PRESSURE_PA) / dynamic_pressure
                    } else {
                        f64::INFINITY
                    };

                    venturi_sigma = venturi_sigma.min(sigma);
                    venturi_shear_rate = venturi_shear_rate.max(shear_rate);
                    venturi_shear_pa = venturi_shear_pa.max(shear_pa);

                    if is_primary_path {
                        pressure_deficit += (analysis.dp_total - dp_graph).max(0.0);
                        let t_throat = l_throat / v_throat.max(1e-9);
                        total_hi += giersiepen_hi(shear_pa, t_throat);
                        // Approximate venturi path len
                        let v_len = d_inlet * 10.0 + l_throat;
                        total_path_len_m += v_len;
                    }
                }
            }
        } else if edge.id.starts_with("serpentine") {
            let w = candidate.channel_width_m;
            let h = candidate.channel_height_m;
            let v_ch = q / (w * h);
            let model = cfd_1d::SerpentineModel::<f64>::millifluidic_rectangular(
                w, h, candidate.segment_length_m, candidate.serpentine_segments, candidate.bend_radius_m,
            );
            let mut cond = FlowConditions::<f64>::new(v_ch);
            cond.flow_rate = Some(q);
            if let Ok(analysis) = model.analyze(&blood, &cond) {
                let shear_rate = analysis.wall_shear_rate;
                let shear_pa = blood.apparent_viscosity(shear_rate) * shear_rate;

                max_main_shear_pa = max_main_shear_pa.max(shear_pa);

                if is_primary_path {
                    pressure_deficit += (analysis.dp_total - dp_graph).max(0.0);
                    let n = candidate.serpentine_segments as f64;
                    let bends = (n - 1.0).max(0.0);
                    let path_m = n * candidate.segment_length_m + bends * std::f64::consts::PI * candidate.bend_radius_m;
                    let res = path_m * w * h / q.max(1e-12);
                    total_hi += giersiepen_hi(shear_pa, res);
                    distribution_path_len_m += path_m;
                    total_path_len_m += path_m;
                    residence_time_s += res;
                }
            }
        } else {
            // Normal trunk/branch rectangular distribution channel
            let w = candidate.channel_width_m;
            let h = candidate.channel_height_m;
            let v_ch = q / (w * h);
            let shear_rate = 6.0 * v_ch / h;
            let shear_pa = blood.apparent_viscosity(shear_rate) * shear_rate;

            max_main_shear_pa = max_main_shear_pa.max(shear_pa);

            if is_primary_path {
                let path_len = edge.properties.length;
                let res = path_len * w * h / q.max(1e-12);
                total_hi += giersiepen_hi(shear_pa, res);
                distribution_path_len_m += path_len;
                total_path_len_m += path_len;
                residence_time_s += res;
            }
        }
    }

    let mut outlet_flows = Vec::new();
    for idx in solved_network.graph.node_indices() {
        if solved_network.graph.node_weight(idx).unwrap().node_type == cfd_1d::network::NodeType::Outlet {
            // Find edge connecting to this node
            for edge in solved_network.edges_with_properties() {
                if edge.nodes.1 == idx {
                    outlet_flows.push(edge.flow_rate.abs());
                }
            }
        }
    }
    
    let flow_uniformity = if outlet_flows.len() > 1 {
        let mean = outlet_flows.iter().sum::<f64>() / outlet_flows.len() as f64;
        let var = outlet_flows.iter().map(|f| (f - mean).powi(2)).sum::<f64>() / outlet_flows.len() as f64;
        let std_dev = var.sqrt();
        (1.0 - (std_dev / mean.max(1e-12))).max(0.0_f64).min(1.0_f64)
    } else {
        1.0_f64
    };

    let total_dp = solved_network.pressures[&inlet_node] + pressure_deficit;

    let cav_potential = if venturi_sigma < SIGMA_CRIT && venturi_sigma >= 0.0 {
        (1.0 - venturi_sigma / SIGMA_CRIT).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let well_coverage_fraction = candidate.topology.nominal_well_coverage();
    let pressure_feasible = total_dp <= candidate.inlet_gauge_pa;

    // Cell separation
    let sep_metrics = if candidate.topology == DesignTopology::CellSeparationVenturi {
        let cancer = cfd_1d::cell_separation::CellProperties::mcf7_breast_cancer();
        let rbc = cfd_1d::cell_separation::CellProperties::red_blood_cell();
        let model = cfd_1d::cell_separation::CellSeparationModel::new(
            candidate.channel_width_m,
            candidate.channel_height_m,
            Some(candidate.bend_radius_m),
        );
        let mean_v = candidate.flow_rate_m3_s / (candidate.channel_width_m * candidate.channel_height_m);
        let shear_est = 6.0 * mean_v / candidate.channel_height_m;
        match model.analyze(&cancer, &rbc, blood.density, blood.apparent_viscosity(shear_est), mean_v) {
            Some(a) => (a.separation_efficiency, a.target_center_fraction, a.background_peripheral_fraction),
            None => (0.0, 0.0, 0.0),
        }
    } else {
        (0.0, 0.0, 0.0)
    };

    Ok(SdtMetrics {
        cavitation_number: venturi_sigma,
        cavitation_potential: cav_potential,
        throat_shear_rate_inv_s: venturi_shear_rate,
        throat_shear_pa: venturi_shear_pa,
        throat_exceeds_fda: venturi_shear_pa > FDA_MAX_WALL_SHEAR_PA,
        max_main_channel_shear_pa: max_main_shear_pa,
        fda_main_compliant: max_main_shear_pa <= FDA_MAX_WALL_SHEAR_PA,
        hemolysis_index_per_pass: total_hi,
        flow_uniformity,
        well_coverage_fraction,
        mean_residence_time_s: residence_time_s,
        total_pressure_drop_pa: total_dp,
        total_path_length_mm: total_path_len_m * 1000.0,
        pressure_feasible,
        cell_separation_efficiency: sep_metrics.0,
        cancer_center_fraction: sep_metrics.1,
        rbc_peripheral_fraction: sep_metrics.2,
    })
}

fn generate_network(
    candidate: &DesignCandidate,
) -> Result<cfd_schematics::domain::model::NetworkBlueprint, OptimError> {
    use cfd_schematics::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
    let mut bp = NetworkBlueprint::new(candidate.id.clone());

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));

    match candidate.topology {
        DesignTopology::SingleVenturi => {
            bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
            bp.add_channel(ChannelSpec::new_pipe(
                "venturi_1", "inlet", "outlet",
                candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0,
            ));
        }
        DesignTopology::BifurcationVenturi => {
            bp.add_node(NodeSpec::new("split1", NodeKind::Junction));
            bp.add_node(NodeSpec::new("split2_1", NodeKind::Junction));
            bp.add_node(NodeSpec::new("split2_2", NodeKind::Junction));

            for i in 1..=4 {
                bp.add_node(NodeSpec::new(format!("outlet_{i}"), NodeKind::Outlet));
            }

            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;

            bp.add_channel(ChannelSpec::new_pipe_rect("trunk", "inlet", "split1", trunk_len, candidate.channel_width_m, candidate.channel_height_m, 1.0, 0.0));
            bp.add_channel(ChannelSpec::new_pipe_rect("b1", "split1", "split2_1", branch_len, candidate.channel_width_m, candidate.channel_height_m, 1.0, 0.0));
            bp.add_channel(ChannelSpec::new_pipe_rect("b2", "split1", "split2_2", branch_len, candidate.channel_width_m, candidate.channel_height_m, 1.0, 0.0));

            bp.add_channel(ChannelSpec::new_pipe("venturi_1", "split2_1", "outlet_1", candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0));
            bp.add_channel(ChannelSpec::new_pipe("venturi_2", "split2_1", "outlet_2", candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0));
            bp.add_channel(ChannelSpec::new_pipe("venturi_3", "split2_2", "outlet_3", candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0));
            bp.add_channel(ChannelSpec::new_pipe("venturi_4", "split2_2", "outlet_4", candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0));
        }
        DesignTopology::TrifurcationVenturi => {
            bp.add_node(NodeSpec::new("split1", NodeKind::Junction));
            for i in 1..=3 {
                bp.add_node(NodeSpec::new(format!("outlet_{i}"), NodeKind::Outlet));
            }

            let trunk_len = TREATMENT_HEIGHT_MM * 0.5e-3;
            bp.add_channel(ChannelSpec::new_pipe_rect("trunk", "inlet", "split1", trunk_len, candidate.channel_width_m, candidate.channel_height_m, 1.0, 0.0));

            bp.add_channel(ChannelSpec::new_pipe("venturi_1", "split1", "outlet_1", candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0));
            bp.add_channel(ChannelSpec::new_pipe("venturi_2", "split1", "outlet_2", candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0));
            bp.add_channel(ChannelSpec::new_pipe("venturi_3", "split1", "outlet_3", candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0));
        }
        DesignTopology::VenturiSerpentine => {
            bp.add_node(NodeSpec::new("venturi_out", NodeKind::Junction));
            bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

            bp.add_channel(ChannelSpec::new_pipe("venturi_1", "inlet", "venturi_out", candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0));
            bp.add_channel(ChannelSpec::new_pipe_rect("serpentine", "venturi_out", "outlet", candidate.segment_length_m, candidate.channel_width_m, candidate.channel_height_m, 1.0, 0.0));
        }
        DesignTopology::SerpentineGrid => {
            bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
            bp.add_channel(ChannelSpec::new_pipe_rect("serpentine", "inlet", "outlet", candidate.segment_length_m, candidate.channel_width_m, candidate.channel_height_m, 1.0, 0.0));
        }
        DesignTopology::CellSeparationVenturi => {
            bp.add_node(NodeSpec::new("split1", NodeKind::Junction));
            bp.add_node(NodeSpec::new("outlet_center", NodeKind::Outlet));
            bp.add_node(NodeSpec::new("outlet_peri", NodeKind::Outlet));

            let trunk_len = TREATMENT_HEIGHT_MM * 0.5e-3;
            bp.add_channel(ChannelSpec::new_pipe_rect("trunk", "inlet", "split1", trunk_len, candidate.channel_width_m, candidate.channel_height_m, 1.0, 0.0));

            bp.add_channel(ChannelSpec::new_pipe("venturi_1", "split1", "outlet_center", candidate.throat_length_m, candidate.inlet_diameter_m, 1.0, 0.0));
            bp.add_channel(ChannelSpec::new_pipe_rect("b_peri", "split1", "outlet_peri", trunk_len, candidate.channel_width_m, candidate.channel_height_m, 1.0, 0.0));
        }
    }

    Ok(bp)
}

// ── Physics helpers ──────────────────────────────────────────────────────────

/// Giersiepen (1990) haemolysis index contribution:
/// `HI = C · t^α · τ^β`
///
/// Returns 0 for non-positive inputs (no damage).
#[inline]
pub fn giersiepen_hi(tau_pa: f64, t_s: f64) -> f64 {
    if tau_pa <= 0.0 || t_s <= 0.0 {
        return 0.0;
    }
    GIERSIEPEN_C * t_s.powf(GIERSIEPEN_ALPHA) * tau_pa.powf(GIERSIEPEN_BETA)
}
