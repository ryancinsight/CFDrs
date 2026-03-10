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
    let mut treatment_residence_time_s = 0.0_f64;
    let mut treatment_flow_fraction = 0.0_f64;
    let mut velocities = Vec::new();

    for sample in &solve.channel_samples {
        if !sample.is_treatment_channel {
            continue;
        }
        let area_m2 = sample.cross_section.area().max(1.0e-18);
        let volume_m3 = sample.length_m * area_m2;
        treatment_volume_m3 += volume_m3;
        treatment_residence_time_s += volume_m3 / sample.flow_m3_s.abs().max(1.0e-18);
        // Accumulate: for multi-branch split trees, treatment channels are
        // parallel and each carries a fraction of inlet flow.  Sum gives the
        // total fraction routed through treatment zones.
        treatment_flow_fraction +=
            (sample.flow_m3_s.abs() / solve.inlet_flow_m3_s.max(1.0e-18)).clamp(0.0, 1.0);
        velocities.push(sample.flow_m3_s.abs() / area_m2);
    }

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
