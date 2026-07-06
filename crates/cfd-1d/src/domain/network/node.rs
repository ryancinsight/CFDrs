//! Network node definitions

use crate::scalar::Cfd1dScalar;
use cfd_core::conversion::SafeFromF64;
use cfd_schematics::domain::model::NodeKind;
use serde::{Deserialize, Serialize};

/// Node in the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node<T: Cfd1dScalar + Copy> {
    /// Unique identifier
    pub id: String,
    /// Type of node
    pub node_type: NodeKind,
    /// Spatial position (x, y)
    pub position: (T, T),
}

use cfd_schematics::domain::model::NodeSpec;

impl<T: Cfd1dScalar + Copy> From<&NodeSpec> for Node<T> {
    fn from(spec: &NodeSpec) -> Self {
        Self {
            id: spec.id.as_str().to_string(),
            node_type: spec.kind,
            position: (T::zero(), T::zero()), // Positions are assigned later from layout
        }
    }
}

impl<T: Cfd1dScalar + Copy> Node<T> {
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
pub struct NodeProperties<T: Cfd1dScalar + Copy> {
    /// Pressure at the node
    pub pressure: T,
    /// Temperature at the node
    pub temperature: T,
    /// Additional properties
    pub metadata: std::collections::HashMap<String, T>,
}

impl<T: Cfd1dScalar + Copy + SafeFromF64> Default for NodeProperties<T> {
    fn default() -> Self {
        Self {
            pressure: T::zero(),
            temperature: T::from_f64_or_zero(293.15), // 20°C default
            metadata: std::collections::HashMap::new(),
        }
    }
}
