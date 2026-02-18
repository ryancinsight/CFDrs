//! Network node definitions

use nalgebra::RealField;
use scheme::domain::model::NodeKind;
use serde::{Deserialize, Serialize};

/// Node in the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node<T: RealField + Copy> {
    /// Unique identifier
    pub id: String,
    /// Type of node
    pub node_type: NodeKind,
    /// Spatial position (x, y)
    pub position: (T, T),
}

impl<T: RealField + Copy> Node<T> {
    /// Create a new node
    #[must_use]
    pub fn new(id: String, node_type: NodeKind) -> Self {
        Self {
            id,
            node_type,
            position: (T::zero(), T::zero()),
        }
    }

    /// Set the position of the node
    pub fn with_position(mut self, x: T, y: T) -> Self {
        self.position = (x, y);
        self
    }
}

/// Properties associated with a node
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeProperties<T: RealField + Copy> {
    /// Pressure at the node
    pub pressure: T,
    /// Temperature at the node
    pub temperature: T,
    /// Additional properties
    pub metadata: std::collections::HashMap<String, T>,
}

impl<T: RealField + Copy> Default for NodeProperties<T> {
    fn default() -> Self {
        Self {
            pressure: T::zero(),
            temperature: T::from_f64(293.15).unwrap_or(T::zero()), // 20°C default
            metadata: std::collections::HashMap::new(),
        }
    }
}
