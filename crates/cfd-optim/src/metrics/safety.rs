use serde::{Deserialize, Serialize};

use crate::constraints::{BLOOD_VISCOSITY_PA_S, FDA_MAX_WALL_SHEAR_PA, FDA_TRANSIENT_TIME_S};
use crate::domain::BlueprintCandidate;

use super::blueprint_graph::BlueprintSolveSummary;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintSafetyMetrics {
    pub max_main_channel_shear_pa: f64,
    pub max_venturi_shear_pa: f64,
    pub pressure_drop_pa: f64,
    pub pressure_feasible: bool,
    pub main_channel_margin: f64,
    pub cavitation_safety_margin: f64,
    pub mean_device_residence_time_s: f64,
}

pub fn compute_blueprint_safety_metrics(
    candidate: &BlueprintCandidate,
    solve: &BlueprintSolveSummary<'_>,
) -> BlueprintSafetyMetrics {
    let mut max_main_channel_shear_pa = 0.0_f64;
    let mut max_venturi_shear_pa = 0.0_f64;
    let mut max_venturi_transit_time_s = 0.0_f64;

    for sample in &solve.channel_samples {
        let area_m2 = sample.cross_section.area().max(1.0e-18);
        let mean_velocity_m_s = sample.flow_m3_s.abs() / area_m2;
        let shear_rate_inv_s = sample.cross_section.wall_shear_rate(mean_velocity_m_s);
        let shear_pa = BLOOD_VISCOSITY_PA_S * shear_rate_inv_s;
        let transit_time_s = sample.length_m / mean_velocity_m_s.max(1.0e-18);

        if sample.is_venturi_channel {
            max_venturi_shear_pa = max_venturi_shear_pa.max(shear_pa);
            max_venturi_transit_time_s = max_venturi_transit_time_s.max(transit_time_s);
        } else {
            max_main_channel_shear_pa = max_main_channel_shear_pa.max(shear_pa);
        }
    }

    let pressure_drop_pa = solve.inlet_pressure_pa.max(0.0) + solve.remerge_loss_pa.max(0.0);
    let pressure_feasible =
        pressure_drop_pa <= candidate.operating_point.inlet_gauge_pa.max(0.0) + 1.0e-9;
    let main_channel_margin =
        (1.0 - max_main_channel_shear_pa / FDA_MAX_WALL_SHEAR_PA).clamp(0.0, 1.0);
    let cavitation_safety_margin = if max_venturi_transit_time_s <= 1.0e-18 {
        main_channel_margin
    } else {
        (FDA_TRANSIENT_TIME_S / max_venturi_transit_time_s).clamp(0.0, 1.0)
    };

    BlueprintSafetyMetrics {
        max_main_channel_shear_pa,
        max_venturi_shear_pa,
        pressure_drop_pa,
        pressure_feasible,
        main_channel_margin,
        cavitation_safety_margin,
        mean_device_residence_time_s: solve.mean_residence_time_s,
    }
}
