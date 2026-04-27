//! CSG boolean operation commands.

use crate::domain::document::history::UndoableCommand;
use crate::domain::document::project::ProjectDocument;
use crate::domain::scene::graph::{MeshHandle, SceneEntity};

/// CSG boolean operation type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CsgOp {
    Union,
    Intersection,
    Difference,
}

/// Command to perform a CSG boolean operation on two meshes.
pub struct CsgBooleanCommand {
    op: CsgOp,
    mesh_a: MeshHandle,
    mesh_b: MeshHandle,
    name: String,
    result_handle: Option<MeshHandle>,
    result_node: Option<usize>,
}

impl CsgBooleanCommand {
    /// Create a new CSG boolean command.
    #[must_use]
    pub fn new(op: CsgOp, mesh_a: MeshHandle, mesh_b: MeshHandle, name: String) -> Self {
        Self {
            op,
            mesh_a,
            mesh_b,
            name,
            result_handle: None,
            result_node: None,
        }
    }
}

impl UndoableCommand for CsgBooleanCommand {
    fn execute(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        let a = doc
            .mesh(self.mesh_a)
            .ok_or_else(|| anyhow::anyhow!("mesh A not found"))?
            .clone();
        let b = doc
            .mesh(self.mesh_b)
            .ok_or_else(|| anyhow::anyhow!("mesh B not found"))?
            .clone();

        let csg_op = match self.op {
            CsgOp::Union => cfd_mesh::application::csg::boolean::operations::BooleanOp::Union,
            CsgOp::Intersection => {
                cfd_mesh::application::csg::boolean::operations::BooleanOp::Intersection
            }
            CsgOp::Difference => {
                cfd_mesh::application::csg::boolean::operations::BooleanOp::Difference
            }
        };

        let mut result = cfd_mesh::application::csg::boolean::indexed::csg_boolean(csg_op, &a, &b)
            .map_err(|e| anyhow::anyhow!("{e}"))?;

        result.orient_outward();
        result.retain_largest_component();

        let handle = doc.add_mesh(result);
        let node_idx = doc
            .scene
            .add_node(self.name.clone(), SceneEntity::Mesh(handle));
        self.result_handle = Some(handle);
        self.result_node = Some(node_idx);
        Ok(())
    }

    fn undo(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        if let Some(node_idx) = self.result_node {
            doc.scene.remove_node(node_idx);
        }
        if let Some(handle) = self.result_handle {
            doc.remove_mesh(handle);
        }
        self.result_handle = None;
        self.result_node = None;
        Ok(())
    }

    fn description(&self) -> &'static str {
        "CSG Boolean"
    }
}
