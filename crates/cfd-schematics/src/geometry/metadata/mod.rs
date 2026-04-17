//! Extensible metadata system for channels and nodes.

use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::fmt::Debug;

mod blueprint;
mod common;
mod schematic;
mod therapy;

pub use self::{blueprint::*, common::*, schematic::*, therapy::*};

/// Base trait for all metadata types.
pub trait Metadata: Any + Debug + Send + Sync {
    fn metadata_type_name(&self) -> &'static str;
    fn clone_metadata(&self) -> Box<dyn Metadata>;
    fn as_any(&self) -> &dyn Any;
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

/// Metadata storage container.
#[derive(Debug)]
pub struct MetadataContainer {
    data: HashMap<TypeId, Box<dyn Metadata>>,
}

impl MetadataContainer {
    #[must_use]
    pub fn new() -> Self {
        Self {
            data: HashMap::new(),
        }
    }

    pub fn insert<T: Metadata + Clone + 'static>(&mut self, metadata: T) {
        self.data.insert(TypeId::of::<T>(), Box::new(metadata));
    }

    #[must_use]
    pub fn get<T: Metadata + 'static>(&self) -> Option<&T> {
        self.data
            .get(&TypeId::of::<T>())
            .and_then(|boxed| boxed.as_any().downcast_ref::<T>())
    }

    pub fn get_mut<T: Metadata + 'static>(&mut self) -> Option<&mut T> {
        self.data
            .get_mut(&TypeId::of::<T>())
            .and_then(|boxed| boxed.as_any_mut().downcast_mut::<T>())
    }

    pub fn remove<T: Metadata + 'static>(&mut self) -> Option<Box<dyn Metadata>> {
        self.data.remove(&TypeId::of::<T>())
    }

    #[must_use]
    pub fn contains<T: Metadata + 'static>(&self) -> bool {
        self.data.contains_key(&TypeId::of::<T>())
    }

    #[must_use]
    pub fn metadata_types(&self) -> Vec<&'static str> {
        self.data
            .values()
            .map(|metadata| metadata.metadata_type_name())
            .collect()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.data.len()
    }
}

impl Clone for MetadataContainer {
    fn clone(&self) -> Self {
        let mut new_container = Self::new();
        for (type_id, metadata) in &self.data {
            new_container.data.insert(*type_id, metadata.clone_metadata());
        }
        new_container
    }
}

impl Default for MetadataContainer {
    fn default() -> Self {
        Self::new()
    }
}

/// Convenience macro for implementing the Metadata trait.
#[macro_export]
macro_rules! impl_metadata {
    ($type:ty, $name:expr) => {
        impl $crate::geometry::metadata::Metadata for $type {
            fn metadata_type_name(&self) -> &'static str {
                $name
            }

            fn clone_metadata(&self) -> Box<dyn $crate::geometry::metadata::Metadata> {
                Box::new(self.clone())
            }

            fn as_any(&self) -> &dyn std::any::Any {
                self
            }

            fn as_any_mut(&mut self) -> &mut dyn std::any::Any {
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

        let flow_data = FlowMetadata {
            flow_rate: 10.0,
            pressure_drop: 1000.0,
            reynolds_number: 0.1,
            velocity: 0.001,
        };

        container.insert(flow_data.clone());

        let retrieved = container
            .get::<FlowMetadata>()
            .expect("structural invariant");
        assert_eq!(retrieved, &flow_data);
        assert!(container.contains::<FlowMetadata>());
        assert!(!container.contains::<ThermalMetadata>());

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

        assert_eq!(container.get::<FlowMetadata>(), Some(&flow_data));
        assert_eq!(container.get::<ThermalMetadata>(), Some(&thermal_data));
        assert_eq!(container.len(), 2);
    }
}