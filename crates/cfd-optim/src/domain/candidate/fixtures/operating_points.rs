use crate::domain::OperatingPoint;

pub(crate) fn operating_point(
    flow_rate_m3_s: f64,
    inlet_gauge_pa: f64,
    feed_hematocrit: f64,
) -> OperatingPoint {
    OperatingPoint {
        flow_rate_m3_s,
        inlet_gauge_pa,
        feed_hematocrit,
        patient_context: None,
    }
}
