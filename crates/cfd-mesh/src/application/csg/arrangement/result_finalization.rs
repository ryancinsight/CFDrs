//! Boolean-result finalization for arrangement CSG.
//!
//! Applies the shared seam-stitching and hole-patching pass used by both the
//! binary arrangement pipeline and the canonical generalized Boolean engine.

use super::patch::patch_small_boundary_holes;
use super::seam::stitch_boundary_seams;
use super::stitch::fill_boundary_loops;
use crate::application::csg::reconstruct::reconstruct_mesh;
use crate::application::watertight::check::check_watertight;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Finalize a Boolean face soup into the canonical watertight output pass.
pub(crate) fn finalize_boolean_faces(result_faces: &mut Vec<FaceData>, pool: &mut VertexPool) {
    let mut preview = reconstruct_mesh(result_faces, pool);
    preview.rebuild_edges();
    let preview_report =
        check_watertight(&preview.vertices, &preview.faces, preview.edges_ref().unwrap());
    if preview_report.is_watertight {
        return;
    }

    stitch_boundary_seams(result_faces, pool);
    fill_boundary_loops(result_faces, pool);
    patch_small_boundary_holes(result_faces, pool);
}
