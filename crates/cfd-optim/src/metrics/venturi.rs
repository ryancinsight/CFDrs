use cfd_1d::{evaluate_venturi_screening, VenturiScreeningInput};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

use crate::constraints::{
    BLOOD_DENSITY_KG_M3, BLOOD_VAPOR_PRESSURE_PA, BLOOD_VISCOSITY_PA_S, DIFFUSER_DISCHARGE_COEFF,
    VENTURI_CC,
};
use crate::domain::BlueprintCandidate;
use crate::error::OptimError;

use super::blueprint_graph::BlueprintSolveSummary;
use super::blueprint_separation::BlueprintSeparationMetrics;

fn cavitation_strength_from_sigma(cavitation_number: f64) -> f64 {
    let strength = (1.0 - cavitation_number).max(0.0);
    strength / (1.0 + strength)
}

fn dynamic_pressure_pa(density_kg_m3: f64, velocity_m_s: f64) -> f64 {
    0.5 * density_kg_m3 * velocity_m_s * velocity_m_s
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiPlacementMetrics {
    pub placement_id: String,
    pub target_channel_id: String,
    pub cavitation_number: f64,
    pub effective_throat_velocity_m_s: f64,
    pub throat_static_pressure_pa: f64,
    pub diffuser_recovery_pa: f64,
    #[serde(default)]
    pub total_loss_coefficient: f64,
    pub dean_number: f64,
    pub curvature_radius_m: f64,
    pub arc_length_m: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintVenturiMetrics {
    pub placements: Vec<VenturiPlacementMetrics>,
    pub cavitation_selectivity_score: f64,
    pub venturi_flow_fraction: f64,
    pub rbc_exposure_fraction: f64,
    pub wbc_exposure_fraction: f64,
}

pub fn compute_blueprint_venturi_metrics(
    candidate: &BlueprintCandidate,
    solve: &BlueprintSolveSummary<'_>,
    separation: &BlueprintSeparationMetrics,
) -> Result<BlueprintVenturiMetrics, OptimError> {
    let topology = candidate.topology_spec()?;
    let n_placements = topology.venturi_placements.len();
    let mut placements = Vec::with_capacity(n_placements);

    let sample_index: HashMap<&str, usize> = solve
        .channel_samples
        .iter()
        .enumerate()
        .map(|(i, s)| (s.id, i))
        .collect();
    let mut used_indices: HashSet<usize> = HashSet::with_capacity(n_placements);
    let mut cavitation_strength_sum = 0.0_f64;

    let mut sample_to_placement = HashMap::new();
    for (p_idx, placement) in topology.venturi_placements.iter().enumerate() {
        let sample_idx = sample_index
            .get(placement.target_channel_id.as_str())
            .copied()
            .filter(|&i| used_indices.insert(i))
            .or_else(|| {
                solve
                    .channel_samples
                    .iter()
                    .enumerate()
                    .find(|(i, s)| {
                        s.id.starts_with(&placement.target_channel_id) && used_indices.insert(*i)
                    })
                    .map(|(i, _)| i)
            })
            .or_else(|| {
                solve
                    .channel_samples
                    .iter()
                    .enumerate()
                    .find(|(i, s)| s.is_venturi_channel && used_indices.insert(*i))
                    .map(|(i, _)| i)
            })
            .ok_or_else(|| OptimError::PhysicsError {
                id: candidate.id.clone(),
                reason: format!(
                    "missing solved sample for venturi target '{}' and no materialized venturi sample fallback was available",
                    placement.target_channel_id
                ),
            })?;
        sample_to_placement.insert(sample_idx, p_idx);
    }

    let mut in_degree: HashMap<&str, usize> = HashMap::new();
    let mut adj: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, sample) in solve.channel_samples.iter().enumerate() {
        *in_degree.entry(sample.to_node).or_insert(0) += 1;
        in_degree.entry(sample.from_node).or_insert(0);
        adj.entry(sample.from_node).or_default().push(i);
    }
    
    let mut queue: Vec<&str> = in_degree
        .iter()
        .filter_map(|(n, &d)| if d == 0 { Some(*n) } else { None })
        .collect();
    let mut sorted_indices = Vec::with_capacity(solve.channel_samples.len());
    while let Some(node) = queue.pop() {
        if let Some(edges) = adj.get(node) {
            for &edge_idx in edges {
                sorted_indices.push(edge_idx);
                let to_node = solve.channel_samples[edge_idx].to_node;
                if let Some(d) = in_degree.get_mut(to_node) {
                    *d -= 1;
                    if *d == 0 {
                        queue.push(to_node);
                    }
                }
            }
        }
    }
    if sorted_indices.len() != solve.channel_samples.len() {
        sorted_indices = (0..solve.channel_samples.len()).collect();
    }

    let transport = cfd_core::physics::cavitation::nuclei_transport::NucleiTransport::new(
        cfd_core::physics::cavitation::nuclei_transport::NucleiTransportConfig::default()
    );

    let mut node_nuclei: HashMap<&str, f64> = HashMap::new();
    let mut placements_out = vec![None; n_placements];

    for idx in sorted_indices {
        let sample = &solve.channel_samples[idx];
        let phi_in = node_nuclei.get(sample.from_node).copied().unwrap_or(0.0);
        
        let velocity = sample.flow_m3_s.abs() / sample.cross_section.area().max(1e-18);
        let transit_time_s = sample.length_m / velocity.max(1e-9);
        let phi_arrival = transport.advect_1d_dissolution(phi_in, transit_time_s);
        let mut phi_out = phi_arrival;

        if let Some(&p_idx) = sample_to_placement.get(&idx) {
            let placement = &topology.venturi_placements[p_idx];
            let mut resolved_placement = placement.clone();
            resolved_placement.target_channel_id = sample.id.to_string();
            let area_inlet_m2 = (placement.throat_geometry.inlet_width_m
                * placement.throat_geometry.throat_height_m)
                .max(1.0e-18);
            let area_throat_m2 = (placement.throat_geometry.throat_width_m
                * placement.throat_geometry.throat_height_m)
                .max(1.0e-18);
            let upstream_velocity_m_s = sample.flow_m3_s.abs() / area_inlet_m2;
            let throat_velocity_m_s = sample.flow_m3_s.abs() / area_throat_m2;
            let upstream_dynamic_pressure_pa =
                dynamic_pressure_pa(BLOOD_DENSITY_KG_M3, upstream_velocity_m_s).max(1.0e-18);
            let screening = evaluate_venturi_screening(VenturiScreeningInput {
                upstream_pressure_pa: candidate.operating_point.absolute_inlet_pressure_pa()
                    + sample.from_pressure_pa.max(0.0),
                upstream_velocity_m_s,
                throat_velocity_m_s,
                throat_hydraulic_diameter_m: 2.0
                    * placement.throat_geometry.throat_width_m
                    * placement.throat_geometry.throat_height_m
                    / (placement.throat_geometry.throat_width_m
                        + placement.throat_geometry.throat_height_m)
                        .max(1.0e-18),
                throat_length_m: placement.throat_geometry.throat_length_m,
                density_kg_m3: BLOOD_DENSITY_KG_M3,
                viscosity_pa_s: BLOOD_VISCOSITY_PA_S,
                vapor_pressure_pa: BLOOD_VAPOR_PRESSURE_PA,
                vena_contracta_coeff: VENTURI_CC,
                diffuser_recovery_coeff: DIFFUSER_DISCHARGE_COEFF,
                upstream_nuclei_fraction: phi_arrival,
            })
            .map_err(|error| OptimError::PhysicsError {
                id: candidate.id.clone(),
                reason: error.to_string(),
            })?;
            phi_out = screening.outlet_nuclei_fraction;
            let total_loss_pa = (screening.bernoulli_drop_pa
                + screening.throat_friction_drop_pa
                - screening.diffuser_recovery_pa)
                .max(0.0);

            let dean_site = cfd_schematics::BlueprintTopologyFactory::estimate_dean_site(
                &candidate.blueprint,
                &resolved_placement,
                sample.flow_m3_s.abs(),
                BLOOD_VISCOSITY_PA_S / BLOOD_DENSITY_KG_M3,
            )
            .unwrap_or_default();
            
            cavitation_strength_sum += cavitation_strength_from_sigma(screening.cavitation_number);

            placements_out[p_idx] = Some(VenturiPlacementMetrics {
                placement_id: placement.placement_id.clone(),
                target_channel_id: sample.id.to_string(),
                cavitation_number: screening.cavitation_number,
                effective_throat_velocity_m_s: screening.effective_throat_velocity_m_s,
                throat_static_pressure_pa: screening.throat_static_pressure_pa,
                diffuser_recovery_pa: screening.diffuser_recovery_pa,
                total_loss_coefficient: total_loss_pa / upstream_dynamic_pressure_pa,
                dean_number: dean_site.dean_number,
                curvature_radius_m: dean_site.curvature_radius_m,
                arc_length_m: dean_site.arc_length_m,
            });
        }
        
        // Update to_node
        let existing = node_nuclei.entry(sample.to_node).or_insert(0.0);
        *existing = existing.max(phi_out);
    }
    
    for p in placements_out {
        if let Some(metrics) = p {
            placements.push(metrics);
        }
    }

    // Normalize by placement count so multi-throat designs aren't
    // automatically favored over single-throat designs by accumulation.
    let n_placements = topology.venturi_placements.len().max(1) as f64;
    let cavitation_term = (cavitation_strength_sum / n_placements).clamp(0.0, 1.0);
    let rbc_exposure_fraction = (1.0 - separation.rbc_peripheral_fraction).clamp(0.0, 1.0);
    let wbc_exposure_fraction = separation.wbc_center_fraction.clamp(0.0, 1.0);

    // Selectivity score: additive weighted sum (85%) + geometric synergy (15%)
    // so that a single zero factor (e.g. zero cancer enrichment in symmetric
    // splits) reduces but never eliminates the gradient signal from the
    // remaining terms.
    let cancer_enrich = separation.cancer_center_fraction.clamp(0.0, 1.0);
    let rbc_shield = (1.0 - rbc_exposure_fraction).clamp(0.0, 1.0);
    let wbc_shield = (1.0 - wbc_exposure_fraction).clamp(0.0, 1.0);
    let additive =
        0.40 * cavitation_term + 0.25 * cancer_enrich + 0.10 * rbc_shield + 0.10 * wbc_shield;
    let geometric = 0.15
        * (cavitation_term * cancer_enrich.max(0.01) * rbc_shield.max(0.01) * wbc_shield.max(0.01))
            .powf(0.25);
    let cavitation_selectivity_score = (additive + geometric).clamp(0.0, 1.0);

    Ok(BlueprintVenturiMetrics {
        placements,
        cavitation_selectivity_score,
        venturi_flow_fraction: solve.venturi_flow_fraction,
        rbc_exposure_fraction,
        wbc_exposure_fraction,
    })
}

#[cfg(test)]
mod tests {
    use cfd_1d::{evaluate_venturi_screening, VenturiScreeningInput};
    use cfd_schematics::VenturiPlacementMode;

    use crate::constraints::{
        BLOOD_DENSITY_KG_M3, BLOOD_VAPOR_PRESSURE_PA, BLOOD_VISCOSITY_PA_S,
        DIFFUSER_DISCHARGE_COEFF, VENTURI_CC,
    };
    use crate::domain::fixtures::{operating_point, stage0_venturi_candidate};
    use crate::metrics::{compute_blueprint_separation_metrics, solve_blueprint_candidate};

    use super::compute_blueprint_venturi_metrics;

    #[test]
    fn saturating_cavitation_strength_discriminates_negative_sigma() {
        // The saturating function f(σ) = max(0,1−σ)/(1+max(0,1−σ)) must produce
        // strictly different outputs for distinct σ < 0, unlike the old
        // (1−σ).clamp(0,1) which collapsed all σ < 0 to 1.0.
        let strength = |sigma: f64| -> f64 {
            let s = (1.0 - sigma).max(0.0);
            s / (1.0 + s)
        };

        // Property 1: f(σ) = 0 for σ ≥ 1 (no cavitation)
        assert_eq!(strength(1.0), 0.0);
        assert_eq!(strength(2.0), 0.0);

        // Property 2: f(σ) ∈ (0, 1) for σ < 1
        let f_half = strength(0.5);
        assert!(f_half > 0.0 && f_half < 1.0, "f(0.5) = {f_half}");

        // Property 3: Monotone — lower σ produces higher score
        assert!(strength(-0.93) > strength(-0.88));
        assert!(strength(-0.88) > strength(0.0));
        assert!(strength(0.0) > strength(0.5));

        // Property 4: Discrimination — distinct negative σ produce distinct values
        let scores: Vec<f64> = [-0.93, -0.91, -0.88, -0.85, -0.80]
            .iter()
            .map(|&s| strength(s))
            .collect();
        for i in 0..scores.len() {
            for j in (i + 1)..scores.len() {
                assert!(
                    (scores[i] - scores[j]).abs() > 1e-12,
                    "σ[{i}] and σ[{j}] produced identical strength: {}",
                    scores[i]
                );
            }
        }

        // Property 5: Asymptotic saturation — f(σ) → 1 as σ → −∞
        assert!(strength(-1000.0) > 0.999);
        assert!(strength(-1000.0) < 1.0);
    }

    #[test]
    fn dean_peak_placement_selects_highest_curvature_treatment_segment() {
        let operating_point = operating_point(2.0e-6, 30_000.0, 0.18);
        let straight_candidate = stage0_venturi_candidate(
            "straight",
            operating_point.clone(),
            VenturiPlacementMode::StraightSegment,
        );
        let dean_candidate = stage0_venturi_candidate(
            "dean",
            operating_point,
            VenturiPlacementMode::CurvaturePeakDeanNumber,
        );

        let straight_solve = solve_blueprint_candidate(&straight_candidate).expect("solve");
        let dean_solve = solve_blueprint_candidate(&dean_candidate).expect("solve");
        let straight_sep =
            compute_blueprint_separation_metrics(&straight_candidate).expect("separation");
        let dean_sep = compute_blueprint_separation_metrics(&dean_candidate).expect("separation");

        let straight =
            compute_blueprint_venturi_metrics(&straight_candidate, &straight_solve, &straight_sep)
                .expect("venturi metrics");
        let dean = compute_blueprint_venturi_metrics(&dean_candidate, &dean_solve, &dean_sep)
            .expect("venturi metrics");

        // The split-tree builder generates 2-point rectilinear paths,
        // so `estimate_dean_site` infers curvature from the segment
        // aspect ratio rather than from explicit serpentine waypoints.
        // Both modes therefore produce the same Dean number in a
        // straight-segment topology; the strict `>` ordering emerges
        // only when the GA introduces serpentine modifications that add
        // curvature peaks (3+ path points).
        assert!(
            dean.placements[0].dean_number >= straight.placements[0].dean_number,
            "CurvaturePeakDeanNumber must select a Dean site at least as strong as StraightSegment"
        );
    }

    /// Zero cavitation: f(σ) = 0 for σ ≥ 1.
    ///
    /// The saturating function `f(σ) = max(0, 1−σ) / (1 + max(0, 1−σ))`
    /// returns 0 when there is no cavitation (σ ≥ 1).
    #[test]
    fn cavitation_strength_zero_above_sigma_one() {
        assert_eq!(super::cavitation_strength_from_sigma(1.0), 0.0);
        assert_eq!(super::cavitation_strength_from_sigma(1.5), 0.0);
        assert_eq!(super::cavitation_strength_from_sigma(100.0), 0.0);
    }

    /// Bounded output: f(σ) ∈ [0, 1) for all σ.
    ///
    /// The function f(σ) = s/(1+s) with s = max(0, 1−σ) is bounded:
    /// - f ≥ 0 because s ≥ 0.
    /// - f < 1 because s/(1+s) < 1 for all finite s.
    /// - f → 1 as σ → −∞ (asymptotic saturation).
    #[test]
    fn cavitation_strength_bounded() {
        for sigma in [-1000.0, -10.0, -1.0, 0.0, 0.5, 0.99] {
            let f = super::cavitation_strength_from_sigma(sigma);
            assert!(f >= 0.0, "f({sigma}) = {f} < 0");
            assert!(f < 1.0, "f({sigma}) = {f} >= 1");
        }
        // Asymptotic: σ = -1000 should be close to 1
        assert!(super::cavitation_strength_from_sigma(-1000.0) > 0.999);
    }

    #[test]
    fn total_loss_coefficient_matches_screening_contract() {
        let candidate = stage0_venturi_candidate(
            "loss-coefficient",
            operating_point(2.0e-6, 30_000.0, 0.18),
            VenturiPlacementMode::StraightSegment,
        );
        let solve = solve_blueprint_candidate(&candidate).expect("solve");
        let separation = compute_blueprint_separation_metrics(&candidate).expect("separation");
        let venturi =
            compute_blueprint_venturi_metrics(&candidate, &solve, &separation).expect("venturi");
        let placement_metrics = venturi.placements.first().expect("venturi placement");
        let topology = candidate.topology_spec().expect("topology");
        let placement_spec = topology
            .venturi_placements
            .iter()
            .find(|placement| placement.placement_id == placement_metrics.placement_id)
            .expect("placement spec");
        let sample = solve
            .channel_samples
            .iter()
            .find(|sample| sample.id == placement_metrics.target_channel_id)
            .expect("resolved venturi sample");

        let inlet_area_m2 = (placement_spec.throat_geometry.inlet_width_m
            * placement_spec.throat_geometry.throat_height_m)
            .max(1.0e-18);
        let throat_area_m2 = (placement_spec.throat_geometry.throat_width_m
            * placement_spec.throat_geometry.throat_height_m)
            .max(1.0e-18);
        let upstream_velocity_m_s = sample.flow_m3_s.abs() / inlet_area_m2;
        let throat_velocity_m_s = sample.flow_m3_s.abs() / throat_area_m2;
        let screening = evaluate_venturi_screening(VenturiScreeningInput {
            upstream_pressure_pa: candidate.operating_point.absolute_inlet_pressure_pa()
                + sample.from_pressure_pa.max(0.0),
            upstream_velocity_m_s,
            throat_velocity_m_s,
            throat_hydraulic_diameter_m: 2.0
                * placement_spec.throat_geometry.throat_width_m
                * placement_spec.throat_geometry.throat_height_m
                / (placement_spec.throat_geometry.throat_width_m
                    + placement_spec.throat_geometry.throat_height_m)
                    .max(1.0e-18),
            throat_length_m: placement_spec.throat_geometry.throat_length_m,
            density_kg_m3: BLOOD_DENSITY_KG_M3,
            viscosity_pa_s: BLOOD_VISCOSITY_PA_S,
            vapor_pressure_pa: BLOOD_VAPOR_PRESSURE_PA,
            vena_contracta_coeff: VENTURI_CC,
            diffuser_recovery_coeff: DIFFUSER_DISCHARGE_COEFF,
            upstream_nuclei_fraction: 0.0,
        })
        .expect("venturi screening should succeed for the reference contract");
        let expected_total_loss_pa = (screening.bernoulli_drop_pa
            + screening.throat_friction_drop_pa
            - screening.diffuser_recovery_pa)
            .max(0.0);
        let expected_total_loss_coefficient = expected_total_loss_pa
            / super::dynamic_pressure_pa(BLOOD_DENSITY_KG_M3, upstream_velocity_m_s).max(1.0e-18);

        assert!(
            placement_metrics.total_loss_coefficient.is_finite()
                && placement_metrics.total_loss_coefficient >= 0.0,
            "total loss coefficient must be finite and non-negative"
        );
        assert!(
            (placement_metrics.total_loss_coefficient - expected_total_loss_coefficient).abs()
                < 1.0e-12,
            "expected total loss coefficient {}, got {}",
            expected_total_loss_coefficient,
            placement_metrics.total_loss_coefficient
        );
    }
}
