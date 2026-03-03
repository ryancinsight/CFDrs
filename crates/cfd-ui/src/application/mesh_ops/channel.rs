//! Channel mesh creation command.

use crate::domain::document::history::UndoableCommand;
use crate::domain::document::project::ProjectDocument;
use crate::domain::scene::graph::{MeshHandle, SceneEntity};

/// Specification for creating a millifluidic channel mesh.
#[derive(Clone, Debug)]
pub enum ChannelSpec {
    /// Venturi constriction channel.
    Venturi {
        inlet_radius: f64,
        throat_radius: f64,
        length: f64,
        segments: usize,
    },
    /// Serpentine (S-bend) channel.
    Serpentine {
        radius: f64,
        length: f64,
        n_bends: usize,
        segments: usize,
    },
}

impl ChannelSpec {
    /// Human-readable name for this channel type.
    #[must_use]
    pub fn type_name(&self) -> &str {
        match self {
            Self::Venturi { .. } => "Venturi",
            Self::Serpentine { .. } => "Serpentine",
        }
    }
}

/// Command to create a channel mesh.
pub struct CreateChannelCommand {
    spec: ChannelSpec,
    name: String,
    created_handle: Option<MeshHandle>,
    created_node: Option<usize>,
}

impl CreateChannelCommand {
    /// Create a new channel creation command.
    #[must_use]
    pub fn new(spec: ChannelSpec, name: String) -> Self {
        Self {
            spec,
            name,
            created_handle: None,
            created_node: None,
        }
    }
}

impl UndoableCommand for CreateChannelCommand {
    fn execute(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        // Channel builders depend on specific cfd-mesh channel API.
        // For now, create a placeholder cylinder to represent the channel.
        let mesh = match &self.spec {
            ChannelSpec::Venturi { inlet_radius, length, segments, .. } => {
                use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
                use nalgebra::Point3;
                (Cylinder { base_center: Point3::<f64>::origin(), radius: *inlet_radius, height: *length, segments: *segments })
                    .build()
                    .map_err(|e| anyhow::anyhow!("{e}"))?
            }
            ChannelSpec::Serpentine { radius, length, segments, .. } => {
                use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
                use nalgebra::Point3;
                (Cylinder { base_center: Point3::<f64>::origin(), radius: *radius, height: *length, segments: *segments })
                    .build()
                    .map_err(|e| anyhow::anyhow!("{e}"))?
            }
        };

        let handle = doc.add_mesh(mesh);
        let node_idx = doc.scene.add_node(
            self.name.clone(),
            SceneEntity::Mesh(handle),
        );
        self.created_handle = Some(handle);
        self.created_node = Some(node_idx);
        Ok(())
    }

    fn undo(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        if let Some(node_idx) = self.created_node {
            doc.scene.remove_node(node_idx);
        }
        if let Some(handle) = self.created_handle {
            doc.remove_mesh(handle);
        }
        self.created_handle = None;
        self.created_node = None;
        Ok(())
    }

    fn description(&self) -> &str {
        "Create Channel"
    }
}
