use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PatientContext {
    pub label: String,
    pub weight_kg: f64,
    pub blood_volume_ml: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct OperatingPoint {
    pub flow_rate_m3_s: f64,
    pub inlet_gauge_pa: f64,
    pub feed_hematocrit: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub patient_context: Option<PatientContext>,
}

impl OperatingPoint {
    #[must_use]
    pub fn absolute_inlet_pressure_pa(&self) -> f64 {
        crate::constraints::P_ATM_PA + self.inlet_gauge_pa.max(0.0)
    }
}
