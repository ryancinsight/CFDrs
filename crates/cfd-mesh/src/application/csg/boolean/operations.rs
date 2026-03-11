//! Low-level face-soup boolean operations.

use crate::domain::core::error::MeshResult;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

pub use crate::application::csg::arrangement::boolean_csg::BooleanOp;

/// Perform a Boolean operation on two face soups sharing one vertex pool.
pub fn csg_boolean(
    op: BooleanOp,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> MeshResult<Vec<FaceData>> {
    crate::application::csg::arrangement::boolean_csg::csg_boolean(
        op,
        &[faces_a.to_vec(), faces_b.to_vec()],
        pool,
    )
}
