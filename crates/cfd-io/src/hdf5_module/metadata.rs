//! Metadata structures for HDF5 datasets.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Metadata for CFD datasets
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetMetadata {
    /// Dataset name
    pub name: String,
    /// Data type description
    pub data_type: String,
    /// Dimensions
    pub dimensions: Vec<usize>,
    /// Physical units
    pub units: Option<String>,
    /// Time step (if time-dependent)
    pub time_step: Option<f64>,
    /// Additional attributes
    pub attributes: HashMap<String, String>,
}

impl DatasetMetadata {
    /// Create metadata for a scalar field
    pub fn scalar_field(name: impl Into<String>, dimensions: Vec<usize>) -> Self {
        Self {
            name: name.into(),
            data_type: "scalar".to_string(),
            dimensions,
            units: None,
            time_step: None,
            attributes: HashMap::new(),
        }
    }

    /// Create metadata for a vector field
    pub fn vector_field(name: impl Into<String>, dimensions: Vec<usize>) -> Self {
        Self {
            name: name.into(),
            data_type: "vector".to_string(),
            dimensions,
            units: None,
            time_step: None,
            attributes: HashMap::new(),
        }
    }

    /// Set units for the dataset
    pub fn with_units(mut self, units: impl Into<String>) -> Self {
        self.units = Some(units.into());
        self
    }

    /// Set time step for time-dependent data
    pub fn with_time_step(mut self, time_step: f64) -> Self {
        self.time_step = Some(time_step);
        self
    }

    /// Add an attribute
    pub fn with_attribute(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.attributes.insert(key.into(), value.into());
        self
    }
}
