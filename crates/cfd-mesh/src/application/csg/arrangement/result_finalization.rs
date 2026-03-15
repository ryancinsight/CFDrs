//! Boolean-result finalization for arrangement CSG.
//!
//! Applies the shared seam-stitching and hole-patching pass used by both the
//! binary arrangement pipeline and the canonical generalized Boolean engine.

use super::curvature_refine::refine_high_curvature_faces;
use super::patch::patch_small_boundary_holes;
use super::seam::stitch_boundary_seams;
use super::stitch::fill_boundary_loops;
use crate::application::csg::reconstruct::reconstruct_mesh;
use crate::application::watertight::check::check_watertight;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Finalize a Boolean face soup into the canonical watertight output pass.
///
/// Pipeline order:
/// 1. Quick watertight check — skip stitching if already clean
/// 2. Seam stitching — resolve T-junctions at shared edges
/// 3. Boundary loop filling — CDT-fill remaining open loops
/// 4. Boundary hole patching — fan-fill small open holes
/// 5. Curvature-adaptive refinement — split high-curvature faces via centroid
///    insertion.  Applied only when the mesh is already watertight (either
///    before stitching or verified after), because centroid splitting of a
///    non-manifold face soup can amplify topology defects.
pub(crate) fn finalize_boolean_faces(result_faces: &mut Vec<FaceData>, pool: &mut VertexPool) {
    let mut preview = reconstruct_mesh(result_faces, pool);
    preview.rebuild_edges();
    let preview_report =
        check_watertight(&preview.vertices, &preview.faces, preview.edges_ref().unwrap());
    if preview_report.is_watertight {
        refine_high_curvature_faces(result_faces, pool);
        return;
    }

    stitch_boundary_seams(result_faces, pool);
    fill_boundary_loops(result_faces, pool);
    patch_small_boundary_holes(result_faces, pool);

    // Post-stitch curvature refinement: only apply if stitch/fill/patch achieved
    // watertightness.  Centroid splitting of a non-manifold soup amplifies edge
    // valence defects, so we guard with a second watertight check.
    let mut post_check = reconstruct_mesh(result_faces, pool);
    post_check.rebuild_edges();
    let post_report =
        check_watertight(&post_check.vertices, &post_check.faces, post_check.edges_ref().unwrap());
    if post_report.is_watertight {
        refine_high_curvature_faces(result_faces, pool);
    }
}
