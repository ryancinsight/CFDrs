use serde::{Deserialize, Serialize};

use crate::domain::BlueprintCandidate;
use aequitas::systems::si::quantities::{Area, Length, Time, Velocity, Volume, VolumetricFlowRate};

use super::blueprint_graph::BlueprintSolveSummary;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResidenceMetrics {
    pub treatment_volume_m3: f64,
    pub treatment_residence_time_s: f64,
    pub treatment_flow_fraction: f64,
    pub mean_treatment_velocity_m_s: f64,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct TypedResidenceMetrics {
    pub(crate) treatment_volume: Volume,
    pub(crate) treatment_residence_time: Time,
    pub(crate) treatment_flow_fraction: f64,
    pub(crate) mean_treatment_velocity: Velocity,
}

impl TypedResidenceMetrics {
    fn into_serialized(self) -> ResidenceMetrics {
        ResidenceMetrics {
            treatment_volume_m3: self.treatment_volume.into_base(),
            treatment_residence_time_s: self.treatment_residence_time.into_base(),
            treatment_flow_fraction: self.treatment_flow_fraction,
            mean_treatment_velocity_m_s: self.mean_treatment_velocity.into_base(),
        }
    }
}

pub(crate) fn compute_typed_residence_metrics(
    _candidate: &BlueprintCandidate,
    solve: &BlueprintSolveSummary<'_>,
) -> TypedResidenceMetrics {
    let mut treatment_volume = Volume::from_base(0.0);
    let mut treatment_flow_fraction = 0.0_f64;
    let mut velocities = Vec::new();

    for sample in &solve.channel_samples {
        if !sample.is_treatment_channel {
            continue;
        }
        let area_m2 = sample.cross_section.area().max(1.0e-18);
        let volume = Length::from_base(sample.length_m) * Area::from_base(area_m2);
        let flow = VolumetricFlowRate::from_base(sample.flow_m3_s.abs().max(1.0e-18));
        treatment_volume += volume;
        treatment_flow_fraction +=
            (sample.flow_m3_s.abs() / solve.inlet_flow_m3_s.max(1.0e-18)).clamp(0.0, 1.0);
        velocities.push(flow / Area::from_base(area_m2));
    }

    let total_treatment_flow = solve
        .channel_samples
        .iter()
        .filter(|s| s.is_treatment_channel)
        .map(|s| s.flow_m3_s.abs())
        .sum::<f64>()
        .max(1.0e-18);
    let treatment_residence_time = solve
        .channel_samples
        .iter()
        .filter(|s| s.is_treatment_channel)
        .map(|s| {
            let area = s.cross_section.area().max(1.0e-18);
            let volume = Length::from_base(s.length_m) * Area::from_base(area);
            let flow = VolumetricFlowRate::from_base(s.flow_m3_s.abs().max(1.0e-18));
            let residence = volume / flow;
            residence * (s.flow_m3_s.abs() / total_treatment_flow)
        })
        .fold(Time::from_base(0.0), |total, residence| total + residence);

    let mean_treatment_velocity = if velocities.is_empty() {
        Velocity::from_base(0.0)
    } else {
        Velocity::from_base(
            velocities
                .iter()
                .map(|velocity| velocity.into_base())
                .sum::<f64>()
                / velocities.len() as f64,
        )
    };

    TypedResidenceMetrics {
        treatment_volume,
        treatment_residence_time,
        treatment_flow_fraction: treatment_flow_fraction.clamp(0.0, 1.0),
        mean_treatment_velocity,
    }
}

pub fn compute_residence_metrics(
    candidate: &BlueprintCandidate,
    solve: &BlueprintSolveSummary<'_>,
) -> ResidenceMetrics {
    compute_typed_residence_metrics(candidate, solve).into_serialized()
}
