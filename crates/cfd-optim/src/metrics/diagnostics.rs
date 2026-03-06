//! Compact 1D audit reporting for a single candidate.

use std::collections::HashMap;
use std::time::Instant;

use cfd_1d::domain::network::network_from_blueprint;
use cfd_1d::{evaluate_venturi_screening, VenturiScreeningInput};
use cfd_core::physics::fluid::blood::CassonBlood;
use serde::Serialize;

use crate::constraints::{
    BLOOD_DENSITY_KG_M3, BLOOD_VAPOR_PRESSURE_PA, DIFFUSER_DISCHARGE_COEFF, P_ATM_PA, VENTURI_CC,
};
use crate::design::DesignCandidate;
use crate::error::OptimError;

use super::{compute_metrics, network_solve::solve_blueprint_network, SdtMetrics};

/// Per-edge 1D hydraulic audit details.
#[derive(Debug, Clone, Serialize)]
pub struct OneDimensionalEdgeAudit {
    pub id: String,
    pub from_node: String,
    pub to_node: String,
    pub length_m: f64,
    pub area_m2: f64,
    pub hydraulic_diameter_m: f64,
    pub linear_resistance: f64,
    pub quadratic_coeff: f64,
    pub solved_flow_m3_s: f64,
    pub upstream_pressure_pa: f64,
    pub wall_shear_rate_per_s: f64,
    pub wall_shear_stress_pa: f64,
    pub is_venturi_throat: bool,
    pub is_bypass_channel: bool,
    pub per_channel_throat_count: u8,
}

/// Audit report for one solved 1D candidate.
#[derive(Debug, Clone, Serialize)]
pub struct OneDimensionalAuditReport {
    pub candidate_id: String,
    pub topology: String,
    pub build_duration_ms: f64,
    pub solve_duration_ms: f64,
    pub metrics_duration_ms: f64,
    pub fallback_used: bool,
    pub branching_junctions_without_geometry: usize,
    pub venturi_channels_without_geometry: usize,
    pub inlet_pressure_pa: f64,
    pub flow_uniformity: f64,
    pub venturi_flow_fraction: f64,
    pub mean_residence_time_s: f64,
    pub total_pressure_drop_pa: f64,
    pub min_cavitation_number: Option<f64>,
    pub max_wall_shear_pa: f64,
    pub cancer_center_fraction: f64,
    pub wbc_center_fraction: f64,
    pub rbc_peripheral_fraction: f64,
    pub metrics: SdtMetrics,
    pub edges: Vec<OneDimensionalEdgeAudit>,
}

/// Build, solve, and summarize the current 1D interpretation of one candidate.
pub fn audit_candidate_1d(
    candidate: &DesignCandidate,
) -> Result<OneDimensionalAuditReport, OptimError> {
    let blood = CassonBlood::<f64>::normal_blood();
    let blueprint = candidate.to_blueprint();

    let branching_junctions_without_geometry = blueprint
        .nodes
        .iter()
        .filter(|node| matches!(node.kind, cfd_schematics::domain::model::NodeKind::Junction))
        .filter(|node| {
            let degree = blueprint
                .channels
                .iter()
                .filter(|channel| channel.from == node.id || channel.to == node.id)
                .count();
            degree >= 3 && node.junction_geometry.is_none()
        })
        .count();
    let venturi_channels_without_geometry = blueprint
        .channels
        .iter()
        .filter(|channel| channel.id.as_str().contains("throat"))
        .filter(|channel| channel.venturi_geometry.is_none())
        .count();

    let build_start = Instant::now();
    let network =
        network_from_blueprint(&blueprint, blood).map_err(|e| OptimError::PhysicsError {
            id: candidate.id.clone(),
            reason: format!("network_from_blueprint failed during audit: {e}"),
        })?;
    let build_duration_ms = build_start.elapsed().as_secs_f64() * 1e3;

    let solve_start = Instant::now();
    let solved = solve_blueprint_network(
        &candidate.id,
        candidate.topology,
        &blueprint,
        &blood,
        candidate.flow_rate_m3_s,
    )?;
    let solve_duration_ms = solve_start.elapsed().as_secs_f64() * 1e3;

    let metrics_start = Instant::now();
    let metrics = compute_metrics(candidate)?;
    let metrics_duration_ms = metrics_start.elapsed().as_secs_f64() * 1e3;

    let mut edge_state: HashMap<&str, (f64, f64)> = HashMap::new();
    for edge_ref in network.graph.edge_references() {
        edge_state.insert(
            edge_ref.weight().id.as_str(),
            (edge_ref.weight().resistance, edge_ref.weight().quad_coeff),
        );
    }

    let mut min_sigma = None::<f64>;
    let mut max_wall_shear_pa = 0.0_f64;
    let mut edges = Vec::with_capacity(blueprint.channels.len());

    for sample in &solved.channel_samples {
        let (linear_resistance, quadratic_coeff) =
            edge_state.get(sample.id).copied().unwrap_or((0.0, 0.0));
        let solved_flow = sample.flow_m3_s;
        let upstream_pressure_pa = sample.from_pressure_pa;
        let area = sample.cross_section.area().max(1e-18);
        let velocity = solved_flow.abs() / area;
        let gamma = sample.cross_section.wall_shear_rate(velocity);
        let mu = blood.apparent_viscosity(gamma);
        let tau = mu * gamma;
        max_wall_shear_pa = max_wall_shear_pa.max(tau);

        if sample.is_venturi_throat {
            let upstream_velocity = solved
                .channel_samples
                .iter()
                .find(|other| other.to_node == sample.from_node && !other.is_venturi_throat)
                .map_or_else(
                    || {
                        candidate.flow_rate_m3_s
                            / (candidate.channel_width_m * candidate.channel_height_m).max(1e-18)
                    },
                    |other| other.flow_m3_s.abs() / other.cross_section.area().max(1e-18),
                );
            let venturi = evaluate_venturi_screening(VenturiScreeningInput {
                upstream_pressure_pa: upstream_pressure_pa.max(0.0) + P_ATM_PA,
                upstream_velocity_m_s: upstream_velocity,
                throat_velocity_m_s: velocity,
                throat_hydraulic_diameter_m: sample.cross_section.hydraulic_diameter(),
                throat_length_m: sample.length_m,
                density_kg_m3: BLOOD_DENSITY_KG_M3,
                viscosity_pa_s: mu.max(1e-18),
                vapor_pressure_pa: BLOOD_VAPOR_PRESSURE_PA,
                vena_contracta_coeff: VENTURI_CC,
                diffuser_recovery_coeff: DIFFUSER_DISCHARGE_COEFF,
            });
            min_sigma = Some(min_sigma.map_or(venturi.cavitation_number, |sigma| {
                sigma.min(venturi.cavitation_number)
            }));
        }

        edges.push(OneDimensionalEdgeAudit {
            id: sample.id.to_string(),
            from_node: sample.from_node.to_string(),
            to_node: sample.to_node.to_string(),
            length_m: sample.length_m,
            area_m2: area,
            hydraulic_diameter_m: sample.cross_section.hydraulic_diameter(),
            linear_resistance,
            quadratic_coeff,
            solved_flow_m3_s: solved_flow,
            upstream_pressure_pa,
            wall_shear_rate_per_s: gamma,
            wall_shear_stress_pa: tau,
            is_venturi_throat: sample.is_venturi_throat,
            is_bypass_channel: sample.is_bypass_channel,
            per_channel_throat_count: sample.per_channel_throat_count,
        });
    }

    Ok(OneDimensionalAuditReport {
        candidate_id: candidate.id.clone(),
        topology: candidate.topology.name().to_string(),
        build_duration_ms,
        solve_duration_ms,
        metrics_duration_ms,
        fallback_used: solved.fallback_used,
        branching_junctions_without_geometry,
        venturi_channels_without_geometry,
        inlet_pressure_pa: solved.inlet_pressure_pa,
        flow_uniformity: solved.flow_uniformity,
        venturi_flow_fraction: solved.venturi_flow_fraction,
        mean_residence_time_s: solved.mean_residence_time_s,
        total_pressure_drop_pa: metrics.total_pressure_drop_pa,
        min_cavitation_number: min_sigma,
        max_wall_shear_pa,
        cancer_center_fraction: metrics.cancer_center_fraction,
        wbc_center_fraction: metrics.wbc_center_fraction,
        rbc_peripheral_fraction: metrics.rbc_peripheral_fraction_three_pop,
        metrics,
        edges,
    })
}
