//! Network metadata and properties

use aequitas::systems::si::quantities::{Pressure, ThermodynamicTemperature, Volume};
use cfd_core::CfdScalar;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Metadata for the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkMetadata<T: CfdScalar + Copy> {
    /// Name of the network
    pub name: String,
    /// Description
    pub description: Option<String>,
    /// Total volume of the network
    pub total_volume: Option<Volume<T>>,
    /// Operating pressure range
    pub pressure_range: Option<(Pressure<T>, Pressure<T>)>,
    /// Temperature range
    pub temperature_range: Option<(ThermodynamicTemperature<T>, ThermodynamicTemperature<T>)>,
    /// Additional properties
    pub properties: HashMap<String, String>,
}

impl<T: CfdScalar + Copy> Default for NetworkMetadata<T> {
    fn default() -> Self {
        Self {
            name: "Unnamed Network".to_string(),
            description: None,
            total_volume: None,
            pressure_range: None,
            temperature_range: None,
            properties: HashMap::new(),
        }
    }
}

impl<T: CfdScalar + Copy> NetworkMetadata<T> {
    /// Create new metadata with a name
    #[must_use]
    pub fn new(name: String) -> Self {
        Self {
            name,
            ..Default::default()
        }
    }

    /// Set the description
    pub fn with_description(mut self, description: String) -> Self {
        self.description = Some(description);
        self
    }

    /// Set the pressure range
    pub fn with_pressure_range(mut self, min: Pressure<T>, max: Pressure<T>) -> Self {
        self.pressure_range = Some((min, max));
        self
    }

    /// Set the total network volume.
    #[must_use]
    pub fn with_total_volume(mut self, total_volume: Volume<T>) -> Self {
        self.total_volume = Some(total_volume);
        self
    }

    /// Set the operating temperature range.
    #[must_use]
    pub fn with_temperature_range(
        mut self,
        min: ThermodynamicTemperature<T>,
        max: ThermodynamicTemperature<T>,
    ) -> Self {
        self.temperature_range = Some((min, max));
        self
    }
}

#[cfg(test)]
mod tests {
    use super::NetworkMetadata;
    use aequitas::systems::si::quantities::{Pressure, ThermodynamicTemperature, Volume};

    #[test]
    fn builder_preserves_typed_network_metrics() {
        let metadata = NetworkMetadata::<f64>::new("vascular".to_string())
            .with_total_volume(Volume::from_base(2.5e-6))
            .with_pressure_range(Pressure::from_base(8.0e3), Pressure::from_base(1.2e4))
            .with_temperature_range(
                ThermodynamicTemperature::from_base(293.15),
                ThermodynamicTemperature::from_base(310.15),
            );

        assert_eq!(metadata.total_volume.unwrap().into_base(), 2.5e-6);
        let (pressure_min, pressure_max) = metadata.pressure_range.unwrap();
        assert_eq!(pressure_min.into_base(), 8.0e3);
        assert_eq!(pressure_max.into_base(), 1.2e4);
        let (temperature_min, temperature_max) = metadata.temperature_range.unwrap();
        assert_eq!(temperature_min.into_base(), 293.15);
        assert_eq!(temperature_max.into_base(), 310.15);
    }

    #[test]
    fn default_metadata_has_no_assumed_metric_values() {
        let metadata = NetworkMetadata::<f64>::default();

        assert!(metadata.total_volume.is_none());
        assert!(metadata.pressure_range.is_none());
        assert!(metadata.temperature_range.is_none());
    }
}
