//! `CsgNode` expression tree
//!
//! ## N-ary flattening optimisation
//!
//! When a flat chain of Union or Intersection nodes is detected (e.g.
//! `Union(Union(A, B), C)`), `evaluate()` collects all leaf operands and
//! dispatches them through the single-pass n-ary arrangement engine
//! (`csg_boolean_nary`).  This avoids O(N²) pairwise intermediate meshes
//! and the associated remap / reconstruct / repair overhead.
//!
//! The flattening is **associativity-preserving**: Union and Intersection
//! are associative, so `(A ∪ B) ∪ C = A ∪ B ∪ C`.  Difference is *not*
//! associative, so it is never flattened and always evaluated as a binary
//! operation.

use super::indexed::{csg_boolean, csg_boolean_nary};
use crate::application::csg::arrangement::boolean_csg::BooleanOp;
use crate::domain::core::error::MeshResult;
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::Real;
use crate::domain::mesh::IndexedMesh;
use nalgebra::Isometry3;

/// A composable CSG expression tree over [`IndexedMesh`] operands.
pub enum CsgNode {
    /// A terminal mesh operand.
    Leaf(Box<IndexedMesh>),
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
    ///
    /// ## N-ary flattening
    ///
    /// Flat chains of the same associative operation (Union or Intersection)
    /// are collected and dispatched through [`csg_boolean_nary`], which runs
    /// a single arrangement pass over all operands simultaneously.
    ///
    /// **Theorem** (associativity of Union/Intersection):
    ///
    /// For closed, orientable surfaces $S_1, S_2, \ldots, S_n$ in $\mathbb{R}^3$,
    /// $$\bigcup_{i=1}^{n} S_i = (\cdots((S_1 \cup S_2) \cup S_3) \cdots \cup S_n)$$
    /// and likewise for $\bigcap$.  Difference is *not* associative
    /// ($(A \setminus B) \setminus C \neq A \setminus (B \setminus C)$),
    /// so it is never flattened.
    ///
    /// **Proof sketch**: Both Union and Intersection are defined point-wise
    /// on the solid interior.  Point-wise OR and AND are associative and
    /// commutative.  The boundary operator commutes with finite Boolean
    /// combinations of closed regular sets (Requicha 1980, "Representations
    /// for Rigid Solids"), so the resulting boundary surface is independent
    /// of evaluation order.  ∎
    pub fn evaluate(self) -> MeshResult<IndexedMesh> {
        match self {
            CsgNode::Leaf(mesh) => Ok(*mesh),
            CsgNode::Union { left, right } => {
                let mut operands = Vec::new();
                Self::collect_flat_operands(*left, BooleanOp::Union, &mut operands)?;
                Self::collect_flat_operands(*right, BooleanOp::Union, &mut operands)?;
                if operands.len() > 2 {
                    csg_boolean_nary(BooleanOp::Union, &operands)
                } else {
                    csg_boolean(BooleanOp::Union, &operands[0], &operands[1])
                }
            }
            CsgNode::Intersection { left, right } => {
                let mut operands = Vec::new();
                Self::collect_flat_operands(*left, BooleanOp::Intersection, &mut operands)?;
                Self::collect_flat_operands(*right, BooleanOp::Intersection, &mut operands)?;
                if operands.len() > 2 {
                    csg_boolean_nary(BooleanOp::Intersection, &operands)
                } else {
                    csg_boolean(BooleanOp::Intersection, &operands[0], &operands[1])
                }
            }
            CsgNode::Difference { left, right } => {
                csg_boolean(BooleanOp::Difference, &left.evaluate()?, &right.evaluate()?)
            }
            CsgNode::Transform { node, iso } => Ok(transform_mesh(node.evaluate()?, &iso)),
        }
    }

    /// Recursively collect leaf operands from a flat chain of the same
    /// associative operation.  Non-matching or non-associative nodes are
    /// evaluated eagerly and pushed as a single operand.
    fn collect_flat_operands(
        node: CsgNode,
        target_op: BooleanOp,
        out: &mut Vec<IndexedMesh>,
    ) -> MeshResult<()> {
        match (target_op, node) {
            (BooleanOp::Union, CsgNode::Union { left, right }) => {
                Self::collect_flat_operands(*left, BooleanOp::Union, out)?;
                Self::collect_flat_operands(*right, BooleanOp::Union, out)?;
            }
            (BooleanOp::Intersection, CsgNode::Intersection { left, right }) => {
                Self::collect_flat_operands(*left, BooleanOp::Intersection, out)?;
                Self::collect_flat_operands(*right, BooleanOp::Intersection, out)?;
            }
            (_, other) => {
                out.push(other.evaluate()?);
            }
        }
        Ok(())
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
            let new_id = if let Some(id) = remap[idx] {
                id
            } else {
                let pos = iso * *mesh.vertices.position(vid);
                let nrm = iso.rotation * *mesh.vertices.normal(vid);
                let id = new_mesh.add_vertex(pos, nrm);
                remap[idx] = Some(id);
                id
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
