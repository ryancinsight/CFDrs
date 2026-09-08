use aequitas::systems::si::quantities::{Pressure, VolumetricFlowRate};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PatientContext {
    pub label: String,
    pub weight_kg: f64,
    pub blood_volume_ml: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OperatingPoint {
    pub flow_rate_m3_s: VolumetricFlowRate,
    pub inlet_gauge_pa: Pressure,
    pub feed_hematocrit: f64,
    pub patient_context: Option<PatientContext>,
}

#[derive(Debug, Deserialize, Serialize)]
struct OperatingPointRepr {
    flow_rate_m3_s: f64,
    inlet_gauge_pa: f64,
    feed_hematocrit: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    patient_context: Option<PatientContext>,
}

impl Default for OperatingPoint {
    fn default() -> Self {
        Self {
            flow_rate_m3_s: VolumetricFlowRate::from_base(0.0),
            inlet_gauge_pa: Pressure::from_base(0.0),
            feed_hematocrit: 0.0,
            patient_context: None,
        }
    }
}

impl Serialize for OperatingPoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        OperatingPointRepr {
            flow_rate_m3_s: self.flow_rate_m3_s.into_base(),
            inlet_gauge_pa: self.inlet_gauge_pa.into_base(),
            feed_hematocrit: self.feed_hematocrit,
            patient_context: self.patient_context.clone(),
        }
        .serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for OperatingPoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let repr = OperatingPointRepr::deserialize(deserializer)?;
        Ok(Self {
            flow_rate_m3_s: VolumetricFlowRate::from_base(repr.flow_rate_m3_s),
            inlet_gauge_pa: Pressure::from_base(repr.inlet_gauge_pa),
            feed_hematocrit: repr.feed_hematocrit,
            patient_context: repr.patient_context,
        })
    }
}

impl OperatingPoint {
    #[must_use]
    pub fn absolute_inlet_pressure_pa(&self) -> Pressure {
        Pressure::from_base(crate::constraints::P_ATM_PA + self.inlet_gauge_pa.into_base().max(0.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serde_roundtrip_preserves_typed_operating_values() {
        let point = OperatingPoint {
            flow_rate_m3_s: VolumetricFlowRate::from_base(2.0e-6),
            inlet_gauge_pa: Pressure::from_base(30_000.0),
            feed_hematocrit: 0.45,
            patient_context: None,
        };

        let encoded = serde_json::to_string(&point).expect("operating point should serialize");
        let decoded: OperatingPoint =
            serde_json::from_str(&encoded).expect("operating point should deserialize");

        assert_eq!(
            decoded.flow_rate_m3_s.into_base().to_bits(),
            point.flow_rate_m3_s.into_base().to_bits()
        );
        assert_eq!(
            decoded.inlet_gauge_pa.into_base().to_bits(),
            point.inlet_gauge_pa.into_base().to_bits()
        );
        assert_eq!(
            decoded.feed_hematocrit.to_bits(),
            point.feed_hematocrit.to_bits()
        );
        assert_eq!(
            decoded.absolute_inlet_pressure_pa().into_base().to_bits(),
            (crate::constraints::P_ATM_PA + 30_000.0).to_bits()
        );
    }
}
