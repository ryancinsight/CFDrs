//! Channel mesh creation command.
//!
//! # Theorem
//! A UI-created channel mesh is admissible for downstream inspection only if it
//! is generated from the canonical schematic authority and the resulting fluid
//! surface is watertight, positively oriented, and a single connected
//! component.
//!
//! **Proof sketch**: the command first materializes a geometry-authored
//! `NetworkBlueprint` using canonical `cfd-schematics` presets, then delegates
//! meshing to `gaia`'s `BlueprintMeshPipeline`. Acceptance checks apply the same
//! surface-topology invariants used elsewhere in the workspace: watertightness
//! rejects holes and boundary leaks, positive signed volume rejects inverted
//! orientation, and a single connected component rejects phantom islands.

use anyhow::{anyhow, ensure, Context};
use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig};
use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::topology::connectivity::connected_components;
use cfd_mesh::domain::topology::AdjacencyGraph;
use cfd_mesh::IndexedMesh;
use cfd_schematics::{serpentine_chain, venturi_chain, NetworkBlueprint};

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

fn circular_segments(requested: usize) -> usize {
    requested.max(12)
}

fn axial_rings(requested: usize) -> usize {
    requested.max(8)
}

fn build_channel_blueprint(spec: &ChannelSpec, name: &str) -> anyhow::Result<NetworkBlueprint> {
    match spec {
        ChannelSpec::Venturi {
            inlet_radius,
            throat_radius,
            length,
            ..
        } => {
            ensure!(*inlet_radius > 0.0, "venturi inlet radius must be positive");
            ensure!(
                *throat_radius > 0.0,
                "venturi throat radius must be positive"
            );
            ensure!(
                throat_radius < inlet_radius,
                "venturi throat radius must be smaller than inlet radius"
            );
            ensure!(*length > 0.0, "venturi length must be positive");
            Ok(venturi_chain(
                name,
                *length,
                2.0 * inlet_radius,
                2.0 * throat_radius,
            ))
        }
        ChannelSpec::Serpentine {
            radius,
            length,
            n_bends,
            ..
        } => {
            ensure!(*radius > 0.0, "serpentine radius must be positive");
            ensure!(*length > 0.0, "serpentine length must be positive");
            let segment_count = n_bends.saturating_add(1).max(1);
            let segment_length_m = *length / segment_count as f64;
            Ok(serpentine_chain(
                name,
                segment_count,
                segment_length_m,
                2.0 * radius,
            ))
        }
    }
}

fn build_channel_mesh(spec: &ChannelSpec, name: &str) -> anyhow::Result<IndexedMesh<f64>> {
    let blueprint = build_channel_blueprint(spec, name)?;
    blueprint
        .validate()
        .map_err(|error| anyhow!("blueprint validation failed: {error}"))?;

    let segment_hint = match spec {
        ChannelSpec::Venturi { segments, .. } | ChannelSpec::Serpentine { segments, .. } => {
            *segments
        }
    };

    let mut output = BlueprintMeshPipeline::run(
        &blueprint,
        &PipelineConfig {
            circular_segments: circular_segments(segment_hint),
            axial_rings: axial_rings(segment_hint),
            include_chip_body: false,
            skip_diameter_constraint: true,
            ..PipelineConfig::default()
        },
    )
    .map_err(|error| anyhow!("blueprint meshing failed: {error}"))?;

    ensure_channel_mesh_acceptance(&mut output.fluid_mesh)?;
    Ok(output.fluid_mesh)
}

fn ensure_channel_mesh_acceptance(mesh: &mut IndexedMesh<f64>) -> anyhow::Result<()> {
    mesh.rebuild_edges();
    let edges = mesh
        .edges_ref()
        .context("channel mesh must rebuild edges before validation")?;
    let watertight = check_watertight(&mesh.vertices, &mesh.faces, edges);
    ensure!(
        watertight.is_watertight,
        "channel mesh must be watertight; boundary edges={}, non-manifold edges={}",
        watertight.boundary_edge_count,
        watertight.non_manifold_edge_count
    );

    let adjacency = AdjacencyGraph::build(&mesh.faces, edges);
    let component_count = connected_components(&mesh.faces, &adjacency).len();
    ensure!(
        component_count == 1,
        "channel mesh must be a single connected component, found {component_count}"
    );

    let signed_volume = mesh.signed_volume();
    ensure!(
        signed_volume.is_finite() && signed_volume > 0.0,
        "channel mesh must have positive finite signed volume, found {signed_volume}"
    );

    Ok(())
}

#[cfg(test)]
pub(crate) fn build_test_venturi_mesh() -> IndexedMesh<f64> {
    build_channel_mesh(
        &ChannelSpec::Venturi {
            inlet_radius: 1.0e-3,
            throat_radius: 4.0e-4,
            length: 12.0e-3,
            segments: 10,
        },
        "reference-test-venturi",
    )
    .expect("test venturi mesh build must succeed")
}

impl UndoableCommand for CreateChannelCommand {
    fn execute(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        let mesh = build_channel_mesh(&self.spec, &self.name)?;

        let handle = doc.add_mesh(mesh);
        let node_idx = doc
            .scene
            .add_node(self.name.clone(), SceneEntity::Mesh(handle));
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

    fn description(&self) -> &'static str {
        "Create Channel"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::scene::graph::MeshHandle;

    #[test]
    fn create_venturi_channel_command_builds_valid_mesh() {
        let mut command = CreateChannelCommand::new(
            ChannelSpec::Venturi {
                inlet_radius: 1.0e-3,
                throat_radius: 4.0e-4,
                length: 12.0e-3,
                segments: 10,
            },
            "test venturi".to_string(),
        );
        let mut document = ProjectDocument::new();

        command
            .execute(&mut document)
            .expect("venturi build must succeed");

        assert_eq!(document.mesh_count(), 1);
        let mut mesh = document
            .mesh(MeshHandle(0))
            .expect("created mesh should be stored in the document")
            .clone();
        ensure_channel_mesh_acceptance(&mut mesh).expect("stored venturi mesh must remain valid");
    }

    #[test]
    fn create_serpentine_channel_command_builds_valid_mesh() {
        let mut command = CreateChannelCommand::new(
            ChannelSpec::Serpentine {
                radius: 7.5e-4,
                length: 24.0e-3,
                n_bends: 3,
                segments: 12,
            },
            "test serpentine".to_string(),
        );
        let mut document = ProjectDocument::new();

        command
            .execute(&mut document)
            .expect("serpentine build must succeed");

        assert_eq!(document.mesh_count(), 1);
        let mut mesh = document
            .mesh(MeshHandle(0))
            .expect("created mesh should be stored in the document")
            .clone();
        ensure_channel_mesh_acceptance(&mut mesh)
            .expect("stored serpentine mesh must remain valid");
    }

    #[test]
    fn create_channel_command_rejects_nonphysical_venturi_geometry() {
        let error = build_channel_mesh(
            &ChannelSpec::Venturi {
                inlet_radius: 5.0e-4,
                throat_radius: 6.0e-4,
                length: 8.0e-3,
                segments: 8,
            },
            "invalid venturi",
        )
        .err()
        .expect("nonphysical venturi geometry must be rejected");

        assert!(error.to_string().contains("throat radius"));
    }
}
