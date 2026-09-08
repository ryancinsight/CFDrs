use crate::domain::OperatingPoint;
use aequitas::systems::si::quantities::{Pressure, VolumetricFlowRate};

pub(crate) fn operating_point(
    flow_rate_m3_s: f64,
    inlet_gauge_pa: f64,
    feed_hematocrit: f64,
) -> OperatingPoint {
    OperatingPoint {
        flow_rate_m3_s: VolumetricFlowRate::from_base(flow_rate_m3_s),
        inlet_gauge_pa: Pressure::from_base(inlet_gauge_pa),
        feed_hematocrit,
        patient_context: None,
    }
}
