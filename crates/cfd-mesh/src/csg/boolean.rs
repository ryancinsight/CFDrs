//! Boolean operations: union, difference, intersection.

use crate::core::error::{MeshError, MeshResult};
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;
use crate::csg::bsp::BspNode;

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

/// A ∪ B = clip_a_to_b_complement ∪ clip_b_to_a_complement
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

    let mut result = bsp_a.all_faces();
    result.extend(bsp_b.all_faces());
    result
}

/// A ∩ B
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

    let mut result = bsp_a.all_faces();
    result.extend(bsp_b.all_faces());

    // Invert back
    for face in &mut result {
        face.flip();
    }
    result
}

/// A \ B (difference)
fn boolean_difference(
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

    let mut result = bsp_a.all_faces();
    result.extend(bsp_b.all_faces());

    for face in &mut result {
        face.flip();
    }
    result
}
