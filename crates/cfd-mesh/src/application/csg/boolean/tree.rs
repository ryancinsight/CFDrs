//! CsgNode expression tree

use nalgebra::Isometry3;
use crate::domain::core::error::MeshResult;
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::Real;
use crate::domain::mesh::IndexedMesh;
use super::indexed::csg_boolean_indexed;
use super::operations::BooleanOp;

/// A composable CSG expression tree over [`IndexedMesh`] operands.
pub enum CsgNode {
    /// A terminal mesh operand.
    Leaf(IndexedMesh),
    /// A ∪ B — mesh union.
    Union {
        /// Left operand.
        left: Box<CsgNode>,
        /// Right operand.
        right: Box<CsgNode>,
    },
    /// A ∩ B — mesh intersection.
    Intersection {
        /// Left operand.
        left: Box<CsgNode>,
        /// Right operand.
        right: Box<CsgNode>,
    },
    /// A \ B — subtract B from A.
    Difference {
        /// Minuend (mesh to subtract from).
        left: Box<CsgNode>,
        /// Subtrahend (mesh to subtract).
        right: Box<CsgNode>,
    },
    /// Apply a rigid-body transform to the sub-tree result.
    Transform {
        /// Sub-tree.
        node: Box<CsgNode>,
        /// Rigid-body isometry.
        iso: Isometry3<Real>,
    },
}

impl CsgNode {
    /// Evaluate the expression tree, consuming `self`.
    pub fn evaluate(self) -> MeshResult<IndexedMesh> {
        match self {
            CsgNode::Leaf(mesh) => Ok(mesh),
            CsgNode::Union { left, right } => {
                csg_boolean_indexed(BooleanOp::Union, &left.evaluate()?, &right.evaluate()?)
            }
            CsgNode::Intersection { left, right } => csg_boolean_indexed(
                BooleanOp::Intersection,
                &left.evaluate()?,
                &right.evaluate()?,
            ),
            CsgNode::Difference { left, right } => {
                csg_boolean_indexed(BooleanOp::Difference, &left.evaluate()?, &right.evaluate()?)
            }
            CsgNode::Transform { node, iso } => Ok(transform_mesh(node.evaluate()?, &iso)),
        }
    }
}

/// Apply a rigid-body `Isometry3` transform to all vertices of a mesh.
fn transform_mesh(mesh: IndexedMesh, iso: &Isometry3<Real>) -> IndexedMesh {
    use crate::domain::core::index::VertexId;
    let mut new_mesh = IndexedMesh::new();
    let mut remap: Vec<Option<VertexId>> = vec![None; mesh.vertices.len()];

    for face in mesh.faces.iter() {
        let mut new_verts = [VertexId::default(); 3];
        for (k, &vid) in face.vertices.iter().enumerate() {
            let idx = vid.as_usize();
            let new_id = match remap[idx] {
                Some(id) => id,
                None => {
                    let pos = iso * *mesh.vertices.position(vid);
                    let nrm = iso.rotation * *mesh.vertices.normal(vid);
                    let id = new_mesh.add_vertex(pos, nrm);
                    remap[idx] = Some(id);
                    id
                }
            };
            new_verts[k] = new_id;
        }
        if face.region == RegionId::INVALID {
            new_mesh.add_face(new_verts[0], new_verts[1], new_verts[2]);
        } else {
            new_mesh.add_face_with_region(new_verts[0], new_verts[1], new_verts[2], face.region);
        }
    }
    new_mesh
}
