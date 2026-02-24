//! IndexedMesh wrapper API for CSG Boolean operations

use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::mesh::IndexedMesh;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;
use super::operations::{BooleanOp, csg_boolean};
use crate::application::csg::reconstruct;

/// High-level Boolean operation on two [`IndexedMesh`] objects.
///
/// Merges the vertex pools, runs the Boolean pipeline with containment
/// detection, and reconstructs a fresh deduplicated `IndexedMesh`.
pub fn csg_boolean_indexed(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    use crate::domain::core::index::VertexId;
    use std::collections::HashMap;

    let mut combined = VertexPool::default_millifluidic();

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

    let faces_a: Vec<FaceData> = mesh_a
        .faces
        .iter()
        .map(|f| FaceData {
            vertices: f.vertices.map(|vid| remap_a[&vid]),
            region: f.region,
        })
        .collect();

    let faces_b: Vec<FaceData> = mesh_b
        .faces
        .iter()
        .map(|f| FaceData {
            vertices: f.vertices.map(|vid| remap_b[&vid]),
            region: f.region,
        })
        .collect();

    // Detect coplanar flat-surface operands before consuming the face soups.
    // Coplanar (2-D) operations produce open surfaces with zero signed volume;
    // the watertight invariant only applies to 3-D solid outputs.
    let is_coplanar = crate::application::csg::coplanar::detect_flat_plane(&faces_a, &combined).is_some()
        && crate::application::csg::coplanar::detect_flat_plane(&faces_b, &combined).is_some();

    let result_faces = csg_boolean(op, &faces_a, &faces_b, &mut combined)?;
    let mut mesh = reconstruct::reconstruct_mesh(&result_faces, &combined);

    // Recompute vertex normals from the final face geometry.
    // After boolean operations, B-sourced faces may have had their winding
    // flipped (e.g. for Difference) but retain the original vertex normals from
    // mesh B.  Recomputing from face geometry ensures alignment consistency.
    mesh.recompute_normals();

    if !is_coplanar {
        mesh.rebuild_edges();
        let report = crate::application::watertight::check::check_watertight(
            &mesh.vertices,
            &mesh.faces,
            mesh.edges_ref().unwrap(),
        );
        if !report.is_watertight {
            return Err(MeshError::NotWatertight {
                count: report.boundary_edge_count + report.non_manifold_edge_count,
            });
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::domain::geometry::primitives::{Cube, Cylinder, Disk, PrimitiveMesh, UvSphere};
    use crate::domain::core::scalar::Point3r;

    fn sphere() -> IndexedMesh {
        UvSphere {
            radius: 1.0,
            center: Point3r::origin(),
            segments: 16,
            stacks: 8,
        }
        .build()
        .expect("sphere build")
    }

    fn cylinder() -> IndexedMesh {
        Cylinder {
            base_center: Point3r::new(0.0, -1.5, 0.0),
            radius: 0.4,
            height: 3.0,
            segments: 16,
        }
        .build()
        .expect("cylinder build")
    }

    fn cube_a() -> IndexedMesh {
        Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_a build")
    }

    fn cube_b() -> IndexedMesh {
        Cube {
            origin: Point3r::new(-0.5, -0.5, -0.5),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_b build")
    }

    fn disk_a() -> IndexedMesh {
        Disk {
            center: Point3r::new(0.0, 0.0, 0.0),
            radius: 1.0,
            segments: 16,
        }
        .build()
        .expect("disk_a build")
    }

    fn disk_b() -> IndexedMesh {
        Disk {
            center: Point3r::new(0.5, 0.0, 0.0),
            radius: 1.0,
            segments: 16,
        }
        .build()
        .expect("disk_b build")
    }

    /// Assert a 3-D CSG result is watertight with a positive signed volume.
    fn assert_3d_watertight(mut mesh: IndexedMesh) {
        mesh.rebuild_edges();
        let report = check_watertight(&mesh.vertices, &mesh.faces, mesh.edges_ref().unwrap());
        assert!(
            report.is_watertight,
            "CSG result must be watertight: {} boundary edge(s), {} non-manifold edge(s)",
            report.boundary_edge_count, report.non_manifold_edge_count,
        );
        assert!(
            mesh.signed_volume() > 0.0,
            "CSG result must have positive signed volume (outward-oriented normals)",
        );
    }

    // ── sphere × cylinder (curved × curved — arrangement pipeline) ─────────────

    #[test]
    fn sphere_cylinder_union_is_watertight() {
        let result = csg_boolean_indexed(BooleanOp::Union, &sphere(), &cylinder())
            .expect("sphere ∪ cylinder");
        assert_3d_watertight(result);
    }

    #[test]
    fn sphere_cylinder_intersection_is_watertight() {
        let result = csg_boolean_indexed(BooleanOp::Intersection, &sphere(), &cylinder())
            .expect("sphere ∩ cylinder");
        assert_3d_watertight(result);
    }

    #[test]
    fn sphere_cylinder_difference_is_watertight() {
        let result = csg_boolean_indexed(BooleanOp::Difference, &sphere(), &cylinder())
            .expect("sphere \\ cylinder");
        assert_3d_watertight(result);
    }

    // ── cube × cube (flat faces — intersecting arrangement pipeline) ───────────

    #[test]
    fn cube_cube_union_is_watertight() {
        let result =
            csg_boolean_indexed(BooleanOp::Union, &cube_a(), &cube_b()).expect("cube ∪ cube");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cube_intersection_is_watertight() {
        let result = csg_boolean_indexed(BooleanOp::Intersection, &cube_a(), &cube_b())
            .expect("cube ∩ cube");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cube_difference_is_watertight() {
        let result =
            csg_boolean_indexed(BooleanOp::Difference, &cube_a(), &cube_b()).expect("cube \\ cube");
        assert_3d_watertight(result);
    }

    // ── disk × disk (coplanar — 2-D Sutherland-Hodgman pipeline) ───────────────
    // Disk operands are open surfaces; the coplanar path produces an open
    // surface result.  Only assert the operation completes without error.

    #[test]
    fn disk_disk_union_succeeds() {
        csg_boolean_indexed(BooleanOp::Union, &disk_a(), &disk_b())
            .expect("disk ∪ disk must not error");
    }

    #[test]
    fn disk_disk_intersection_succeeds() {
        csg_boolean_indexed(BooleanOp::Intersection, &disk_a(), &disk_b())
            .expect("disk ∩ disk must not error");
    }

    #[test]
    fn disk_disk_difference_succeeds() {
        csg_boolean_indexed(BooleanOp::Difference, &disk_a(), &disk_b())
            .expect("disk \\ disk must not error");
    }
}
