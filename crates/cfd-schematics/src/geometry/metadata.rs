//! Extensible metadata system for channels and nodes
//!
//! This module provides a flexible metadata system that allows for easy addition
//! of new tracking variables without requiring changes to core data structures.
//! It uses trait-based extensibility with type-safe metadata storage.

use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::fmt::Debug;

/// Base trait for all metadata types
///
/// This trait provides the foundation for type-safe metadata storage.
/// All metadata types must implement this trait to be stored in the system.
pub trait Metadata: Any + Debug + Send + Sync {
    /// Returns a unique name for this metadata type
    fn metadata_type_name(&self) -> &'static str;

    /// Clone the metadata as a boxed trait object
    fn clone_metadata(&self) -> Box<dyn Metadata>;

    /// Convert to Any for downcasting
    fn as_any(&self) -> &dyn Any;

    /// Convert to mutable Any for downcasting
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

/// Metadata storage container
///
/// This container provides type-safe storage and retrieval of metadata
/// using `TypeId` as keys for efficient lookup.
#[derive(Debug)]
pub struct MetadataContainer {
    data: HashMap<TypeId, Box<dyn Metadata>>,
}

impl MetadataContainer {
    /// Create a new empty metadata container
    #[must_use]
    pub fn new() -> Self {
        Self {
            data: HashMap::new(),
        }
    }

    /// Insert metadata of a specific type
    pub fn insert<T: Metadata + Clone + 'static>(&mut self, metadata: T) {
        self.data.insert(TypeId::of::<T>(), Box::new(metadata));
    }

    /// Get metadata of a specific type
    #[must_use]
    pub fn get<T: Metadata + 'static>(&self) -> Option<&T> {
        self.data
            .get(&TypeId::of::<T>())
            .and_then(|boxed| boxed.as_any().downcast_ref::<T>())
    }

    /// Get mutable metadata of a specific type
    pub fn get_mut<T: Metadata + 'static>(&mut self) -> Option<&mut T> {
        self.data
            .get_mut(&TypeId::of::<T>())
            .and_then(|boxed| boxed.as_any_mut().downcast_mut::<T>())
    }

    /// Remove metadata of a specific type
    pub fn remove<T: Metadata + 'static>(&mut self) -> Option<Box<dyn Metadata>> {
        self.data.remove(&TypeId::of::<T>())
    }

    /// Check if metadata of a specific type exists
    #[must_use]
    pub fn contains<T: Metadata + 'static>(&self) -> bool {
        self.data.contains_key(&TypeId::of::<T>())
    }

    /// Get all metadata type names (for debugging)
    #[must_use]
    pub fn metadata_types(&self) -> Vec<&'static str> {
        self.data
            .values()
            .map(|metadata| metadata.metadata_type_name())
            .collect()
    }

    /// Check if container is empty
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get number of metadata entries
    #[must_use]
    pub fn len(&self) -> usize {
        self.data.len()
    }
}

impl Clone for MetadataContainer {
    fn clone(&self) -> Self {
        let mut new_container = Self::new();
        for (type_id, metadata) in &self.data {
            new_container
                .data
                .insert(*type_id, metadata.clone_metadata());
        }
        new_container
    }
}

impl Default for MetadataContainer {
    fn default() -> Self {
        Self::new()
    }
}

/// Flow-related metadata for channels
#[derive(Debug, Clone, PartialEq)]
pub struct FlowMetadata {
    /// Flow rate in μL/min
    pub flow_rate: f64,
    /// Pressure drop in Pa
    pub pressure_drop: f64,
    /// Reynolds number
    pub reynolds_number: f64,
    /// Velocity in m/s
    pub velocity: f64,
}

impl Metadata for FlowMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "FlowMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Thermal metadata for channels
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalMetadata {
    /// Temperature in Celsius
    pub temperature: f64,
    /// Heat transfer coefficient in W/(m²·K)
    pub heat_transfer_coefficient: f64,
    /// Thermal conductivity in W/(m·K)
    pub thermal_conductivity: f64,
}

impl Metadata for ThermalMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "ThermalMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Manufacturing tolerance metadata
#[derive(Debug, Clone, PartialEq)]
pub struct ManufacturingMetadata {
    /// Width tolerance in micrometers
    pub width_tolerance: f64,
    /// Height tolerance in micrometers
    pub height_tolerance: f64,
    /// Surface roughness in micrometers
    pub surface_roughness: f64,
    /// Manufacturing method
    pub manufacturing_method: String,
}

impl Metadata for ManufacturingMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "ManufacturingMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Channel geometry metadata for downstream 3D-aware workflows.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelGeometryMetadata {
    /// Design channel diameter in millimeters used for spacing and reconstruction.
    pub channel_diameter_mm: f64,
}

impl Metadata for ChannelGeometryMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "ChannelGeometryMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Optimization history metadata
#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationMetadata {
    /// Original channel length before optimization
    pub original_length: f64,
    /// Optimized channel length
    pub optimized_length: f64,
    /// Length improvement percentage
    pub improvement_percentage: f64,
    /// Optimization iterations used
    pub iterations: usize,
    /// Optimization time in milliseconds
    pub optimization_time_ms: u64,
    /// Optimization profile used
    pub optimization_profile: String,
}

impl Metadata for OptimizationMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "OptimizationMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Runtime performance metadata
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PerformanceMetadata {
    /// Generation time in microseconds
    pub generation_time_us: u64,
    /// Memory usage in bytes
    pub memory_usage_bytes: usize,
    /// Number of path points generated
    pub path_points_count: usize,
}

impl Metadata for PerformanceMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "PerformanceMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Convenience macro for implementing Metadata trait
#[macro_export]
macro_rules! impl_metadata {
    ($type:ty, $name:expr) => {
        impl Metadata for $type {
            fn metadata_type_name(&self) -> &'static str {
                $name
            }

            fn clone_metadata(&self) -> Box<dyn Metadata> {
                Box::new(self.clone())
            }

            fn as_any(&self) -> &dyn Any {
                self
            }

            fn as_any_mut(&mut self) -> &mut dyn Any {
                self
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_metadata_container_basic_operations() {
        let mut container = MetadataContainer::new();

        // Test insertion and retrieval
        let flow_data = FlowMetadata {
            flow_rate: 10.0,
            pressure_drop: 1000.0,
            reynolds_number: 0.1,
            velocity: 0.001,
        };

        container.insert(flow_data.clone());

        let retrieved = container.get::<FlowMetadata>().unwrap();
        assert_eq!(retrieved, &flow_data);

        // Test contains
        assert!(container.contains::<FlowMetadata>());
        assert!(!container.contains::<ThermalMetadata>());

        // Test removal
        let removed = container.remove::<FlowMetadata>();
        assert!(removed.is_some());
        assert!(!container.contains::<FlowMetadata>());
    }

    #[test]
    fn test_multiple_metadata_types() {
        let mut container = MetadataContainer::new();

        let flow_data = FlowMetadata {
            flow_rate: 10.0,
            pressure_drop: 1000.0,
            reynolds_number: 0.1,
            velocity: 0.001,
        };

        let thermal_data = ThermalMetadata {
            temperature: 25.0,
            heat_transfer_coefficient: 100.0,
            thermal_conductivity: 0.6,
        };

        container.insert(flow_data.clone());
        container.insert(thermal_data.clone());

        assert_eq!(container.len(), 2);
        assert!(container.contains::<FlowMetadata>());
        assert!(container.contains::<ThermalMetadata>());

        let retrieved_flow = container.get::<FlowMetadata>().unwrap();
        let retrieved_thermal = container.get::<ThermalMetadata>().unwrap();

        assert_eq!(retrieved_flow, &flow_data);
        assert_eq!(retrieved_thermal, &thermal_data);
    }
}
