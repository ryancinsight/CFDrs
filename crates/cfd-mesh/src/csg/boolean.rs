//! Boolean operations: union, difference, intersection.
//!
//! This module exposes two API levels:
//!
//! ## Low-level face-soup API (backward compatible)
//!
//! [`csg_boolean`] operates on `Vec<FaceData>` + a shared `VertexPool`.
//! Internally uses the BSP-tree approach from `csgrs`.
//!
//! ## High-level `IndexedMesh` API (new in Phase 8)
//!
//! [`csg_boolean_indexed`] takes two `IndexedMesh` objects, merges their
//! vertex pools, runs the BSP pipeline, and reconstructs a fresh mesh.
//!
//! [`CsgNode`] provides a composable CSG tree that evaluates lazily.

use crate::core::error::{MeshError, MeshResult};
use crate::core::index::RegionId;
use crate::core::scalar::Real;
use crate::mesh::IndexedMesh;
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;
use crate::csg::bsp::BspNode;
use crate::csg::reconstruct;
use nalgebra::Isometry3;

/// CSG boolean operation type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BooleanOp {
    /// A ∪ B
    Union,
    /// A ∩ B
    Intersection,
    /// A \ B
    Difference,
}

/// Perform a CSG boolean operation on two sets of faces.
///
/// Both sets share the same `VertexPool`; new intersection vertices are
/// inserted into the pool during splitting.
pub fn csg_boolean(
    op: BooleanOp,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> MeshResult<Vec<FaceData>> {
    let result = match op {
        BooleanOp::Union => boolean_union(faces_a, faces_b, pool),
        BooleanOp::Intersection => boolean_intersection(faces_a, faces_b, pool),
        BooleanOp::Difference => boolean_difference(faces_a, faces_b, pool),
    };

    if result.is_empty() {
        return Err(MeshError::EmptyBooleanResult {
            op: format!("{:?}", op),
        });
    }

    Ok(result)
}

/// A ∪ B
///
/// Matches csgrs exactly:
/// ```text
/// a.clip_to(b); b.clip_to(a); b.invert(); b.clip_to(a); b.invert();
/// a.build(b.all()); return a.all()
/// ```
/// The invert/clip/invert around b de-duplicates coplanar faces that would
/// otherwise appear in both a and b, producing non-manifold geometry.
fn boolean_union(
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    let mut bsp_a = BspNode::build(faces_a, pool);
    let mut bsp_b = BspNode::build(faces_b, pool);

    bsp_a.clip_to(&bsp_b, pool);
    bsp_b.clip_to(&bsp_a, pool);
    bsp_b.invert();
    bsp_b.clip_to(&bsp_a, pool);
    bsp_b.invert();

    let b_faces = bsp_b.all_faces();
    bsp_a.add_faces(&b_faces, pool);
    bsp_a.all_faces()
}

/// A ∩ B
///
/// Matches csgrs exactly:
/// ```text
/// a.invert(); b.clip_to(a); b.invert(); a.clip_to(b); b.clip_to(a);
/// a.build(b.all()); a.invert(); return a.all()
/// ```
fn boolean_intersection(
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    let mut bsp_a = BspNode::build(faces_a, pool);
    let mut bsp_b = BspNode::build(faces_b, pool);

    bsp_a.invert();
    bsp_b.clip_to(&bsp_a, pool);
    bsp_b.invert();
    bsp_a.clip_to(&bsp_b, pool);
    bsp_b.clip_to(&bsp_a, pool);

    let b_faces = bsp_b.all_faces();
    bsp_a.add_faces(&b_faces, pool);
    bsp_a.invert();
    bsp_a.all_faces()
}

/// A \ B  (subtract B from A)
///
/// Matches csgrs exactly:
/// ```text
/// a.invert(); a.clip_to(b); b.clip_to(a); b.invert(); b.clip_to(a); b.invert();
/// a.build(b.all()); a.invert(); return a.all()
/// ```
fn boolean_difference(
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    let mut bsp_a = BspNode::build(faces_a, pool);
    let mut bsp_b = BspNode::build(faces_b, pool);

    // Difference: A \ B
    bsp_a.invert();
    bsp_a.clip_to(&bsp_b, pool);
    bsp_b.clip_to(&bsp_a, pool);
    bsp_b.invert();
    bsp_b.clip_to(&bsp_a, pool);
    bsp_b.invert();
    
    let b_faces = bsp_b.all_faces();
    bsp_a.add_faces(&b_faces, pool);
    bsp_a.invert();
    bsp_a.all_faces()
}

// ── IndexedMesh API (Phase 8) ─────────────────────────────────────────────────

/// High-level Boolean operation on two [`IndexedMesh`] objects.
///
/// Merges the vertex pools of both meshes, runs the BSP Boolean pipeline,
/// and reconstructs a fresh deduplicated `IndexedMesh`.
///
/// This is the recommended entry point for Phase 8 and beyond.  The
/// lower-level [`csg_boolean`] face-soup API remains available for
/// backward compatibility.
pub fn csg_boolean_indexed(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    use std::collections::HashMap;
    use crate::core::index::VertexId;

    // ── Build a combined vertex pool from both meshes ─────────────────────
    let mut combined = VertexPool::default_millifluidic();

    // Map each mesh's VertexId → combined-pool VertexId.
    let mut remap_a: HashMap<VertexId, VertexId> = HashMap::new();
    for (old_id, _) in mesh_a.vertices.iter() {
        let pos = *mesh_a.vertices.position(old_id);
        let nrm = *mesh_a.vertices.normal(old_id);
        remap_a.insert(old_id, combined.insert_or_weld(pos, nrm));
    }

    let mut remap_b: HashMap<VertexId, VertexId> = HashMap::new();
    for (old_id, _) in mesh_b.vertices.iter() {
        let pos = *mesh_b.vertices.position(old_id);
        let nrm = *mesh_b.vertices.normal(old_id);
        remap_b.insert(old_id, combined.insert_or_weld(pos, nrm));
    }

    // ── Remap face vertex IDs to the combined pool ────────────────────────
    let faces_a: Vec<FaceData> = mesh_a.faces.iter().map(|f| {
        let v = f.vertices.map(|vid| remap_a[&vid]);
        FaceData { vertices: v, region: f.region }
    }).collect();

    let faces_b: Vec<FaceData> = mesh_b.faces.iter().map(|f| {
        let v = f.vertices.map(|vid| remap_b[&vid]);
        FaceData { vertices: v, region: f.region }
    }).collect();

    // ── Run BSP Boolean on the remapped face soups ────────────────────────
    let result_faces = csg_boolean(op, &faces_a, &faces_b, &mut combined)?;

    // ── Reconstruct into a fresh, deduplicated IndexedMesh ────────────────
    Ok(reconstruct::reconstruct_mesh(&result_faces, &combined))
}

// ── CsgNode expression tree ───────────────────────────────────────────────────

/// A composable CSG expression tree over [`IndexedMesh`] operands.
///
/// Nodes are evaluated lazily (depth-first) when [`CsgNode::evaluate`] is
/// called.  Each non-leaf node allocates a new `IndexedMesh` for its result.
///
/// # Example
///
/// ```rust,ignore
/// use cfd_mesh::csg::{CsgNode, BooleanOp};
///
/// let cube   = CsgNode::Leaf(build_cube());
/// let sphere = CsgNode::Leaf(build_sphere());
///
/// let swiss_cheese = CsgNode::Difference {
///     left:  Box::new(cube),
///     right: Box::new(sphere),
/// };
/// let mesh = swiss_cheese.evaluate().expect("CSG evaluation");
/// ```
pub enum CsgNode {
    /// A terminal mesh operand.
    Leaf(IndexedMesh),
    /// A ∪ B — mesh union.
    Union {
        /// Left operand of the union.
        left: Box<CsgNode>,
        /// Right operand of the union.
        right: Box<CsgNode>,
    },
    /// A ∩ B — mesh intersection.
    Intersection {
        /// Left operand of the intersection.
        left: Box<CsgNode>,
        /// Right operand of the intersection.
        right: Box<CsgNode>,
    },
    /// A \ B — subtract B from A.
    Difference {
        /// Mesh to subtract from (minuend).
        left: Box<CsgNode>,
        /// Mesh to subtract (subtrahend).
        right: Box<CsgNode>,
    },
    /// Apply a rigid-body transform to the sub-tree before evaluation.
    Transform {
        /// Sub-tree whose result will be transformed.
        node: Box<CsgNode>,
        /// Rigid-body isometry to apply to all vertices.
        iso: Isometry3<Real>,
    },
}

impl CsgNode {
    /// Evaluate the CSG expression tree, consuming `self`.
    ///
    /// Returns an `IndexedMesh` representing the result, or a [`MeshError`]
    /// if any Boolean operation fails (e.g. two disjoint meshes intersected).
    pub fn evaluate(self) -> MeshResult<IndexedMesh> {
        match self {
            CsgNode::Leaf(mesh) => Ok(mesh),
            CsgNode::Union { left, right } => {
                let a = left.evaluate()?;
                let b = right.evaluate()?;
                csg_boolean_indexed(BooleanOp::Union, &a, &b)
            }
            CsgNode::Intersection { left, right } => {
                let a = left.evaluate()?;
                let b = right.evaluate()?;
                csg_boolean_indexed(BooleanOp::Intersection, &a, &b)
            }
            CsgNode::Difference { left, right } => {
                let a = left.evaluate()?;
                let b = right.evaluate()?;
                csg_boolean_indexed(BooleanOp::Difference, &a, &b)
            }
            CsgNode::Transform { node, iso } => {
                let mesh = node.evaluate()?;
                Ok(transform_mesh(mesh, &iso))
            }
        }
    }
}

// ── Transform helper ──────────────────────────────────────────────────────────

/// Apply a rigid-body `Isometry3` transform to all vertices of a mesh.
///
/// Positions are translated + rotated; normals are only rotated (no
/// translation for pseudo-vectors).  The output is a fresh `IndexedMesh`.
fn transform_mesh(mesh: IndexedMesh, iso: &Isometry3<Real>) -> IndexedMesh {
    use crate::core::index::VertexId;
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
                    let new_id = new_mesh.add_vertex(pos, nrm);
                    remap[idx] = Some(new_id);
                    new_id
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
