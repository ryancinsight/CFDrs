//! Builder patterns for creating channels and nodes with metadata
//!
//! This module provides convenient builder patterns for creating channels and nodes
//! with optional metadata, making it easy to add tracking variables without
//! breaking existing code.

use super::metadata::{Metadata, MetadataContainer};
use super::types::{Channel, ChannelType, Node, Point2D};
use crate::domain::model::{ChannelShape, NodeKind};
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::{ChannelVisualRole, JunctionGeometryMetadata, VenturiGeometryMetadata};

/// Builder for creating nodes with optional metadata
#[derive(Debug)]
pub struct NodeBuilder {
    id: usize,
    name: Option<String>,
    point: Point2D,
    kind: Option<NodeKind>,
    junction_geometry: Option<JunctionGeometryMetadata>,
    metadata: Option<MetadataContainer>,
}

impl NodeBuilder {
    /// Create a new node builder
    #[must_use]
    pub const fn new(id: usize, point: Point2D) -> Self {
        Self {
            id,
            name: None,
            point,
            kind: None,
            junction_geometry: None,
            metadata: None,
        }
    }

    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    #[must_use]
    pub const fn with_kind(mut self, kind: NodeKind) -> Self {
        self.kind = Some(kind);
        self
    }

    #[must_use]
    pub fn with_junction_geometry(mut self, geometry: JunctionGeometryMetadata) -> Self {
        self.junction_geometry = Some(geometry);
        self
    }

    /// Add metadata to the node
    ///
    /// # Panics
    ///
    /// Panics if the metadata container is in an invalid state (should not happen in normal usage).
    #[must_use]
    pub fn with_metadata<T: Metadata + Clone + 'static>(mut self, metadata: T) -> Self {
        if self.metadata.is_none() {
            self.metadata = Some(MetadataContainer::new());
        }
        self.metadata.as_mut().unwrap().insert(metadata);
        self
    }

    /// Add multiple metadata entries
    #[must_use]
    pub fn with_metadata_container(mut self, container: MetadataContainer) -> Self {
        self.metadata = Some(container);
        self
    }

    /// Build the node
    #[must_use]
    pub fn build(self) -> Node {
        Node {
            id: self.id,
            name: self.name,
            point: self.point,
            kind: self.kind,
            junction_geometry: self.junction_geometry,
            metadata: self.metadata,
        }
    }
}

/// Builder for creating channels with optional metadata
#[derive(Debug)]
pub struct ChannelBuilder {
    id: usize,
    name: Option<String>,
    from_node: usize,
    to_node: usize,
    width: f64,
    height: f64,
    channel_type: ChannelType,
    visual_role: Option<ChannelVisualRole>,
    physical_length_m: Option<f64>,
    physical_width_m: Option<f64>,
    physical_height_m: Option<f64>,
    physical_shape: Option<ChannelShape>,
    therapy_zone: Option<TherapyZone>,
    venturi_geometry: Option<VenturiGeometryMetadata>,
    metadata: Option<MetadataContainer>,
}

impl ChannelBuilder {
    /// Create a new channel builder
    #[must_use]
    pub const fn new(
        id: usize,
        from_node: usize,
        to_node: usize,
        width: f64,
        height: f64,
        channel_type: ChannelType,
    ) -> Self {
        Self {
            id,
            name: None,
            from_node,
            to_node,
            width,
            height,
            channel_type,
            visual_role: None,
            physical_length_m: None,
            physical_width_m: None,
            physical_height_m: None,
            physical_shape: None,
            therapy_zone: None,
            venturi_geometry: None,
            metadata: None,
        }
    }

    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    #[must_use]
    pub const fn with_visual_role(mut self, role: ChannelVisualRole) -> Self {
        self.visual_role = Some(role);
        self
    }

    #[must_use]
    pub const fn with_physical_length_m(mut self, length_m: f64) -> Self {
        self.physical_length_m = Some(length_m);
        self
    }

    #[must_use]
    pub const fn with_physical_dims_m(mut self, width_m: f64, height_m: f64) -> Self {
        self.physical_width_m = Some(width_m);
        self.physical_height_m = Some(height_m);
        self
    }

    #[must_use]
    pub const fn with_physical_shape(mut self, shape: ChannelShape) -> Self {
        self.physical_shape = Some(shape);
        self
    }

    #[must_use]
    pub const fn with_therapy_zone(mut self, zone: TherapyZone) -> Self {
        self.therapy_zone = Some(zone);
        self
    }

    #[must_use]
    pub fn with_venturi_geometry(mut self, geometry: VenturiGeometryMetadata) -> Self {
        self.venturi_geometry = Some(geometry);
        self
    }

    /// Add metadata to the channel
    ///
    /// # Panics
    ///
    /// This method will panic if the metadata container is in an invalid state.
    /// This should never happen under normal usage.
    #[must_use]
    pub fn with_metadata<T: Metadata + Clone + 'static>(mut self, metadata: T) -> Self {
        if self.metadata.is_none() {
            self.metadata = Some(MetadataContainer::new());
        }
        self.metadata.as_mut().unwrap().insert(metadata);
        self
    }

    /// Add multiple metadata entries
    #[must_use]
    pub fn with_metadata_container(mut self, container: MetadataContainer) -> Self {
        self.metadata = Some(container);
        self
    }

    /// Build the channel
    #[must_use]
    pub fn build(self) -> Channel {
        Channel {
            id: self.id,
            name: self.name,
            from_node: self.from_node,
            to_node: self.to_node,
            width: self.width,
            height: self.height,
            channel_type: self.channel_type,
            visual_role: self.visual_role,
            physical_length_m: self.physical_length_m,
            physical_width_m: self.physical_width_m,
            physical_height_m: self.physical_height_m,
            physical_shape: self.physical_shape,
            therapy_zone: self.therapy_zone,
            venturi_geometry: self.venturi_geometry,
            metadata: self.metadata,
        }
    }
}

/// Extension trait for Node to provide convenient metadata access
pub trait NodeExt {
    /// Get metadata of a specific type
    fn get_metadata<T: Metadata + 'static>(&self) -> Option<&T>;

    /// Get mutable metadata of a specific type
    fn get_metadata_mut<T: Metadata + 'static>(&mut self) -> Option<&mut T>;

    /// Add metadata to the node
    fn add_metadata<T: Metadata + Clone + 'static>(&mut self, metadata: T);

    /// Check if node has metadata of a specific type
    fn has_metadata<T: Metadata + 'static>(&self) -> bool;

    /// Remove metadata of a specific type
    fn remove_metadata<T: Metadata + 'static>(&mut self) -> bool;

    /// Get all metadata type names
    fn metadata_types(&self) -> Vec<&'static str>;
}

impl NodeExt for Node {
    fn get_metadata<T: Metadata + 'static>(&self) -> Option<&T> {
        self.metadata.as_ref()?.get::<T>()
    }

    fn get_metadata_mut<T: Metadata + 'static>(&mut self) -> Option<&mut T> {
        self.metadata.as_mut()?.get_mut::<T>()
    }

    fn add_metadata<T: Metadata + Clone + 'static>(&mut self, metadata: T) {
        if self.metadata.is_none() {
            self.metadata = Some(MetadataContainer::new());
        }
        if let Some(container) = self.metadata.as_mut() {
            container.insert(metadata);
        }
    }

    fn has_metadata<T: Metadata + 'static>(&self) -> bool {
        self.metadata
            .as_ref()
            .is_some_and(super::metadata::MetadataContainer::contains::<T>)
    }

    fn remove_metadata<T: Metadata + 'static>(&mut self) -> bool {
        self.metadata
            .as_mut()
            .is_some_and(|m| m.remove::<T>().is_some())
    }

    fn metadata_types(&self) -> Vec<&'static str> {
        self.metadata.as_ref().map_or(
            Vec::new(),
            super::metadata::MetadataContainer::metadata_types,
        )
    }
}

/// Extension trait for Channel to provide convenient metadata access
pub trait ChannelExt {
    /// Get metadata of a specific type
    fn get_metadata<T: Metadata + 'static>(&self) -> Option<&T>;

    /// Get mutable metadata of a specific type
    fn get_metadata_mut<T: Metadata + 'static>(&mut self) -> Option<&mut T>;

    /// Add metadata to the channel
    fn add_metadata<T: Metadata + Clone + 'static>(&mut self, metadata: T);

    /// Check if channel has metadata of a specific type
    fn has_metadata<T: Metadata + 'static>(&self) -> bool;

    /// Remove metadata of a specific type
    fn remove_metadata<T: Metadata + 'static>(&mut self) -> bool;

    /// Get all metadata type names
    fn metadata_types(&self) -> Vec<&'static str>;
}

impl ChannelExt for Channel {
    fn get_metadata<T: Metadata + 'static>(&self) -> Option<&T> {
        self.metadata.as_ref()?.get::<T>()
    }

    fn get_metadata_mut<T: Metadata + 'static>(&mut self) -> Option<&mut T> {
        self.metadata.as_mut()?.get_mut::<T>()
    }

    fn add_metadata<T: Metadata + Clone + 'static>(&mut self, metadata: T) {
        if self.metadata.is_none() {
            self.metadata = Some(MetadataContainer::new());
        }
        if let Some(container) = self.metadata.as_mut() {
            container.insert(metadata);
        }
    }

    fn has_metadata<T: Metadata + 'static>(&self) -> bool {
        self.metadata
            .as_ref()
            .is_some_and(super::metadata::MetadataContainer::contains::<T>)
    }

    fn remove_metadata<T: Metadata + 'static>(&mut self) -> bool {
        self.metadata
            .as_mut()
            .is_some_and(|m| m.remove::<T>().is_some())
    }

    fn metadata_types(&self) -> Vec<&'static str> {
        self.metadata.as_ref().map_or(
            Vec::new(),
            super::metadata::MetadataContainer::metadata_types,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::metadata::{FlowMetadata, ThermalMetadata};

    #[test]
    fn test_node_builder() {
        let flow_data = FlowMetadata {
            flow_rate: 10.0,
            pressure_drop: 1000.0,
            reynolds_number: 0.1,
            velocity: 0.001,
        };

        let node = NodeBuilder::new(0, (5.0, 10.0))
            .with_metadata(flow_data.clone())
            .build();

        assert_eq!(node.id, 0);
        assert_eq!(node.point, (5.0, 10.0));
        assert!(node.has_metadata::<FlowMetadata>());

        let retrieved = node.get_metadata::<FlowMetadata>().unwrap();
        assert_eq!(retrieved, &flow_data);
    }

    #[test]
    fn test_channel_builder() {
        let thermal_data = ThermalMetadata {
            temperature: 25.0,
            heat_transfer_coefficient: 100.0,
            thermal_conductivity: 0.6,
        };

        let channel = ChannelBuilder::new(0, 0, 1, 1.0, 0.5, ChannelType::Straight)
            .with_metadata(thermal_data.clone())
            .build();

        assert_eq!(channel.id, 0);
        assert_eq!(channel.from_node, 0);
        assert_eq!(channel.to_node, 1);
        assert!(channel.has_metadata::<ThermalMetadata>());

        let retrieved = channel.get_metadata::<ThermalMetadata>().unwrap();
        assert_eq!(retrieved, &thermal_data);
    }

    #[test]
    fn test_extension_traits() {
        let mut node = Node {
            id: 0,
            point: (0.0, 0.0),
            metadata: None,
        };

        let flow_data = FlowMetadata {
            flow_rate: 5.0,
            pressure_drop: 500.0,
            reynolds_number: 0.05,
            velocity: 0.0005,
        };

        // Test adding metadata
        node.add_metadata(flow_data.clone());
        assert!(node.has_metadata::<FlowMetadata>());

        // Test getting metadata
        let retrieved = node.get_metadata::<FlowMetadata>().unwrap();
        assert_eq!(retrieved, &flow_data);

        // Test metadata types
        let types = node.metadata_types();
        assert_eq!(types, vec!["FlowMetadata"]);

        // Test removing metadata
        assert!(node.remove_metadata::<FlowMetadata>());
        assert!(!node.has_metadata::<FlowMetadata>());
    }
}
