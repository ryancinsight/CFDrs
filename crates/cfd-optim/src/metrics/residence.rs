use serde::{Deserialize, Serialize};

use crate::domain::BlueprintCandidate;

use super::blueprint_graph::BlueprintSolveSummary;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResidenceMetrics {
    pub treatment_volume_m3: f64,
    pub treatment_residence_time_s: f64,
    pub treatment_flow_fraction: f64,
    pub mean_treatment_velocity_m_s: f64,
}

pub fn compute_residence_metrics(
    _candidate: &BlueprintCandidate,
    solve: &BlueprintSolveSummary<'_>,
) -> ResidenceMetrics {
    let mut treatment_volume_m3 = 0.0_f64;
    let mut treatment_flow_fraction = 0.0_f64;
    let mut velocities = Vec::new();
    let mut per_channel_residence = Vec::new();

    for sample in &solve.channel_samples {
        if !sample.is_treatment_channel {
            continue;
        }
        let area_m2 = sample.cross_section.area().max(1.0e-18);
        let volume_m3 = sample.length_m * area_m2;
        treatment_volume_m3 += volume_m3;
        // Per-channel residence time: V / Q for this channel segment.
        let channel_res_s = volume_m3 / sample.flow_m3_s.abs().max(1.0e-18);
        per_channel_residence.push(channel_res_s);
        treatment_flow_fraction +=
            (sample.flow_m3_s.abs() / solve.inlet_flow_m3_s.max(1.0e-18)).clamp(0.0, 1.0);
        velocities.push(sample.flow_m3_s.abs() / area_m2);
    }
    // Treatment residence time is the MAXIMUM per-channel time, not the sum.
    // In a parallel split tree, cells traverse one path (the longest).
    // For serial channels on the same path, their times should be summed,
    // but for parallel branches only the slowest path matters.
    //
    // As an approximation (exact would require path tracing through the
    // network graph), use the flow-weighted average: each channel's
    // contribution is weighted by the fraction of total treatment flow
    // it carries.  This gives the expected residence time for a randomly
    // selected cell entering the treatment zone.
    let total_treatment_flow = solve
        .channel_samples
        .iter()
        .filter(|s| s.is_treatment_channel)
        .map(|s| s.flow_m3_s.abs())
        .sum::<f64>()
        .max(1e-18);
    let treatment_residence_time_s = solve
        .channel_samples
        .iter()
        .filter(|s| s.is_treatment_channel)
        .map(|s| {
            let area = s.cross_section.area().max(1e-18);
            let vol = s.length_m * area;
            let res = vol / s.flow_m3_s.abs().max(1e-18);
            let weight = s.flow_m3_s.abs() / total_treatment_flow;
            res * weight
        })
        .sum::<f64>();

    let mean_treatment_velocity_m_s = if velocities.is_empty() {
        0.0
    } else {
        velocities.iter().sum::<f64>() / velocities.len() as f64
    };

    ResidenceMetrics {
        treatment_volume_m3,
        treatment_residence_time_s,
        treatment_flow_fraction: treatment_flow_fraction.clamp(0.0, 1.0),
        mean_treatment_velocity_m_s,
    }
}
