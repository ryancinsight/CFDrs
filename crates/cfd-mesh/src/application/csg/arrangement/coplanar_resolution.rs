//! Shared coplanar-group resolution for arrangement CSG.
//!
//! Resolves one exact coplanar face group across one or more meshes while
//! respecting face orientation relative to the group plane normal.

use super::super::boolean::BooleanOp;
use super::super::coplanar::basis::PlaneBasis;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

fn accumulate_oriented_faces(
    op: BooleanOp,
    faces_by_mesh: &[Vec<FaceData>],
    basis: &PlaneBasis,
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    let mut accum_opt: Option<Vec<FaceData>> = None;

    for (mesh_idx, mesh_faces) in faces_by_mesh.iter().enumerate() {
        if let Some(accum) = accum_opt.take() {
            accum_opt = Some(crate::application::csg::coplanar::boolean_coplanar(
                op, &accum, mesh_faces, pool, basis,
            ));
        } else if mesh_idx == 0 {
            accum_opt = Some(mesh_faces.clone());
        } else if op == BooleanOp::Union {
            accum_opt = Some(mesh_faces.clone());
        } else {
            accum_opt = Some(Vec::new());
            break;
        }
    }

    accum_opt.unwrap_or_default()
}

/// Resolve one coplanar group using the canonical orientation-aware policy.
///
/// Faces whose normals oppose the plane basis are flipped into the basis
/// orientation before 2-D clipping, then re-separated from co-directed faces so
/// same-side and opposite-side sheets do not cancel incorrectly.
pub(crate) fn resolve_oriented_coplanar_group(
    op: BooleanOp,
    faces_by_mesh: &[Vec<FaceData>],
    basis: &PlaneBasis,
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    let mut same_by_mesh = vec![Vec::new(); faces_by_mesh.len()];
    let mut opposite_by_mesh = vec![Vec::new(); faces_by_mesh.len()];

    for (mesh_idx, mesh_faces) in faces_by_mesh.iter().enumerate() {
        for face in mesh_faces {
            let p0 = *pool.position(face.vertices[0]);
            let p1 = *pool.position(face.vertices[1]);
            let p2 = *pool.position(face.vertices[2]);
            let normal = (p1 - p0).cross(&(p2 - p0));

            if normal.dot(&basis.normal) > 0.0 {
                same_by_mesh[mesh_idx].push(*face);
            } else {
                let mut flipped = *face;
                flipped.vertices.swap(1, 2);
                opposite_by_mesh[mesh_idx].push(flipped);
            }
        }
    }

    let same_faces = accumulate_oriented_faces(op, &same_by_mesh, basis, pool);
    let opposite_faces = accumulate_oriented_faces(op, &opposite_by_mesh, basis, pool);

    let mut result = if opposite_faces.is_empty() {
        same_faces.clone()
    } else if same_faces.is_empty() {
        Vec::new()
    } else {
        crate::application::csg::coplanar::boolean_coplanar(
            BooleanOp::Difference,
            &same_faces,
            &opposite_faces,
            pool,
            basis,
        )
    };

    let mut opposite_residual = if same_faces.is_empty() {
        opposite_faces.clone()
    } else if opposite_faces.is_empty() {
        Vec::new()
    } else {
        crate::application::csg::coplanar::boolean_coplanar(
            BooleanOp::Difference,
            &opposite_faces,
            &same_faces,
            pool,
            basis,
        )
    };

    for face in &mut opposite_residual {
        face.vertices.swap(1, 2);
    }

    result.extend(opposite_residual);
    result
}
