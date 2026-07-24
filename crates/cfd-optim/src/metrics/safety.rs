use serde::{Deserialize, Serialize};

use crate::constraints::{FDA_MAX_WALL_SHEAR_PA, FDA_TRANSIENT_TIME_S};
use crate::domain::BlueprintCandidate;
use aequitas::systems::si::quantities::{Length, Pressure, Time, Velocity};
use cfd_2d::physics::non_newtonian::CarreauYasudaModel;

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

#[derive(Debug, Clone, Copy)]
pub(crate) struct TypedBlueprintSafetyMetrics {
    pub(crate) max_main_channel_shear: Pressure,
    pub(crate) max_venturi_shear: Pressure,
    pub(crate) pressure_drop: Pressure,
    pub(crate) pressure_feasible: bool,
    pub(crate) main_channel_margin: f64,
    pub(crate) cavitation_safety_margin: f64,
    pub(crate) mean_device_residence_time: Time,
}

impl TypedBlueprintSafetyMetrics {
    fn into_serialized(self) -> BlueprintSafetyMetrics {
        BlueprintSafetyMetrics {
            max_main_channel_shear_pa: self.max_main_channel_shear.into_base(),
            max_venturi_shear_pa: self.max_venturi_shear.into_base(),
            pressure_drop_pa: self.pressure_drop.into_base(),
            pressure_feasible: self.pressure_feasible,
            main_channel_margin: self.main_channel_margin,
            cavitation_safety_margin: self.cavitation_safety_margin,
            mean_device_residence_time_s: self.mean_device_residence_time.into_base(),
        }
    }
}

pub(crate) fn compute_typed_blueprint_safety_metrics(
    candidate: &BlueprintCandidate,
    solve: &BlueprintSolveSummary<'_>,
) -> TypedBlueprintSafetyMetrics {
    let mut max_main_channel_shear = Pressure::from_base(0.0);
    let mut max_venturi_shear = Pressure::from_base(0.0);
    let mut max_venturi_transit_time = Time::from_base(0.0);

    let cy_model = CarreauYasudaModel::<f64>::typical_blood();

    for sample in &solve.channel_samples {
        let area_m2 = sample.cross_section.area().max(1.0e-18);
        let mean_velocity = Velocity::from_base(sample.flow_m3_s.abs() / area_m2);
        let shear_rate_inv_s = sample
            .cross_section
            .wall_shear_rate(mean_velocity.into_base());
        let apparent_viscosity_pa_s = cy_model.apparent_viscosity(shear_rate_inv_s);
        let shear = Pressure::from_base(apparent_viscosity_pa_s * shear_rate_inv_s);
        let transit_time = Length::from_base(sample.length_m)
            / Velocity::from_base(mean_velocity.into_base().max(1.0e-18));

        if sample.is_venturi_channel {
            if shear > max_venturi_shear {
                max_venturi_shear = shear;
            }
            if transit_time > max_venturi_transit_time {
                max_venturi_transit_time = transit_time;
            }
        } else if shear > max_main_channel_shear {
            max_main_channel_shear = shear;
        }
    }

    let pressure_drop =
        Pressure::from_base(solve.inlet_pressure_pa.max(0.0) + solve.remerge_loss_pa.max(0.0));
    let pressure_feasible =
        pressure_drop.into_base() <= candidate.operating_point.inlet_gauge_pa.max(0.0) + 1.0e-9;
    let main_channel_margin =
        (1.0 - max_main_channel_shear.into_base() / FDA_MAX_WALL_SHEAR_PA).clamp(0.0, 1.0);
    let cavitation_safety_margin = if max_venturi_transit_time.into_base() <= 1.0e-18 {
        main_channel_margin
    } else {
        (FDA_TRANSIENT_TIME_S / max_venturi_transit_time.into_base()).clamp(0.0, 1.0)
    };

    TypedBlueprintSafetyMetrics {
        max_main_channel_shear,
        max_venturi_shear,
        pressure_drop,
        pressure_feasible,
        main_channel_margin,
        cavitation_safety_margin,
        mean_device_residence_time: Time::from_base(solve.mean_residence_time_s),
    }
}

pub fn compute_blueprint_safety_metrics(
    candidate: &BlueprintCandidate,
    solve: &BlueprintSolveSummary<'_>,
) -> BlueprintSafetyMetrics {
    compute_typed_blueprint_safety_metrics(candidate, solve).into_serialized()
}
