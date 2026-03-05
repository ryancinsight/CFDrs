//! Network metadata and properties

use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Metadata for the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkMetadata<T: RealField + Copy> {
    /// Name of the network
    pub name: String,
    /// Description
    pub description: Option<String>,
    /// Total volume of the network
    pub total_volume: Option<T>,
    /// Operating pressure range
    pub pressure_range: Option<(T, T)>,
    /// Temperature range
    pub temperature_range: Option<(T, T)>,
    /// Additional properties
    pub properties: HashMap<String, String>,
}

impl<T: RealField + Copy> Default for NetworkMetadata<T> {
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

impl<T: RealField + Copy> NetworkMetadata<T> {
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
    pub fn with_pressure_range(mut self, min: T, max: T) -> Self {
        self.pressure_range = Some((min, max));
        self
    }
}
