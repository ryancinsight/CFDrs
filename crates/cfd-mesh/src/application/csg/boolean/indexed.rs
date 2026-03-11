//! `IndexedMesh` wrapper API for CSG Boolean operations

use super::union_strategy::{concat_disjoint_meshes, spatial_union_order};
use crate::application::csg::boolean::BooleanOp;
use crate::application::csg::reconstruct;
use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::mesh::IndexedMesh;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// High-level boolean operation on two [`IndexedMesh`] objects.
///
/// Merges the vertex pools, runs the proven binary Boolean pipeline, and
/// reconstructs a fresh deduplicated `IndexedMesh`.
pub fn csg_boolean(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    csg_boolean_indexed(op, mesh_a, mesh_b)
}

/// Backward-compatible alias for the historical indexed binary Boolean entrypoint.
pub fn csg_boolean_indexed(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    let mut combined = VertexPool::for_csg();
    let (faces_a, faces_b) = remap_binary_face_soups(mesh_a, mesh_b, &mut combined);
    let is_coplanar = crate::application::csg::coplanar::detect_flat_plane(&faces_a, &combined)
        .is_some()
        && crate::application::csg::coplanar::detect_flat_plane(&faces_b, &combined).is_some();
    let result_faces =
        crate::application::csg::arrangement::boolean_csg::csg_boolean_unfinalized(
            op,
            &[faces_a.clone(), faces_b.clone()],
            &mut combined,
        )?;
    postprocess_boolean_mesh(result_faces, &combined, is_coplanar)
}

/// Best-effort indexed Boolean that preserves the historical tolerant path for
/// junction-heavy unions.
pub fn csg_boolean_indexed_tolerant(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    match csg_boolean_indexed(op, mesh_a, mesh_b) {
        Ok(mesh) => Ok(mesh),
        Err(MeshError::NotWatertight { .. }) => {
            let mut combined = VertexPool::for_csg();
            let (faces_a, faces_b) = remap_binary_face_soups(mesh_a, mesh_b, &mut combined);
            let is_coplanar = crate::application::csg::coplanar::detect_flat_plane(
                &faces_a,
                &combined,
            )
            .is_some()
                && crate::application::csg::coplanar::detect_flat_plane(&faces_b, &combined)
                    .is_some();
            let result_faces =
                crate::application::csg::arrangement::boolean_csg::csg_boolean_unfinalized(
                    op,
                    &[faces_a.clone(), faces_b.clone()],
                    &mut combined,
                )?;
            let mut mesh = reconstruct::reconstruct_mesh(&result_faces, &combined);
            mesh.recompute_normals();
            repair_boolean_mesh(&mut mesh, is_coplanar)?;
            Ok(mesh)
        }
        Err(error) => Err(error),
    }
}

/// Best-effort indexed Boolean used by higher-level assembly paths that can
/// accept a repaired-but-still-open intermediate mesh and enforce topology
/// policy at the final aggregate surface instead.
pub fn csg_boolean_indexed_best_effort(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    match csg_boolean_indexed(op, mesh_a, mesh_b) {
        Ok(mesh) => Ok(mesh),
        Err(MeshError::NotWatertight { .. }) => {
            let mut combined = VertexPool::for_csg();
            let (faces_a, faces_b) = remap_binary_face_soups(mesh_a, mesh_b, &mut combined);
            let is_coplanar = crate::application::csg::coplanar::detect_flat_plane(
                &faces_a,
                &combined,
            )
            .is_some()
                && crate::application::csg::coplanar::detect_flat_plane(&faces_b, &combined)
                    .is_some();
            let result_faces =
                crate::application::csg::arrangement::boolean_csg::csg_boolean_unfinalized(
                    op,
                    &[faces_a.clone(), faces_b.clone()],
                    &mut combined,
                )?;
            let mut mesh = reconstruct::reconstruct_mesh(&result_faces, &combined);
            mesh.recompute_normals();
            repair_boolean_mesh_with_policy(&mut mesh, is_coplanar, false)?;
            Ok(mesh)
        }
        Err(error) => Err(error),
    }
}

/// Compute an indexed Boolean union across an arbitrary number of meshes using
/// the canonical generalized arrangement engine.
pub fn csg_boolean_nary_union(meshes: &[IndexedMesh]) -> MeshResult<IndexedMesh> {
    if meshes.is_empty() {
        return Err(MeshError::EmptyBooleanResult {
            op: format!("{:?}", BooleanOp::Union),
        });
    }

    if meshes.len() == 1 {
        return Ok(meshes[0].clone());
    }

    match csg_boolean_nary(BooleanOp::Union, meshes) {
        Ok(mesh) => Ok(mesh),
        Err(MeshError::NotWatertight { .. }) => {
            let order = spatial_union_order(meshes);
            let mut accumulated = meshes[order[0]].clone();
            let mut accumulated_aabb = accumulated.bounding_box();
            for &idx in &order[1..] {
                let mesh = &meshes[idx];
                let mesh_aabb = mesh.bounding_box();
                accumulated = if !accumulated_aabb.intersects(&mesh_aabb) {
                    concat_disjoint_meshes(&accumulated, mesh)
                } else {
                    csg_boolean_indexed_tolerant(BooleanOp::Union, &accumulated, mesh)?
                };
                accumulated_aabb = accumulated.bounding_box();
            }
            Ok(accumulated)
        }
        Err(error) => Err(error),
    }
}

/// Compute an indexed Boolean across an arbitrary number of meshes using the
/// canonical generalized arrangement engine.
///
/// The first mesh is treated as the minuend for [`BooleanOp::Difference`], and
/// every later mesh is treated as a subtractive operand.
pub fn csg_boolean_nary(op: BooleanOp, meshes: &[IndexedMesh]) -> MeshResult<IndexedMesh> {
    if meshes.is_empty() {
        return Err(MeshError::EmptyBooleanResult {
            op: format!("{op:?}"),
        });
    }

    if meshes.len() == 1 {
        return Ok(meshes[0].clone());
    }

    let mut combined = VertexPool::for_csg();
    let face_soups = remap_nary_face_soups(meshes, &mut combined);
    let is_coplanar = face_soups.iter().all(|faces| {
        crate::application::csg::coplanar::detect_flat_plane(faces, &combined).is_some()
    });
    let result_faces = crate::application::csg::arrangement::boolean_csg::csg_boolean_unfinalized(
        op,
        &face_soups,
        &mut combined,
    )?;
    postprocess_boolean_mesh(result_faces, &combined, is_coplanar)
}

fn remap_binary_face_soups(
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
    combined: &mut VertexPool,
) -> (Vec<FaceData>, Vec<FaceData>) {
    use crate::domain::core::index::VertexId;
    use std::collections::HashMap;

    let mut remap_a: HashMap<VertexId, VertexId> = HashMap::with_capacity(mesh_a.vertices.len());
    for (old_id, _) in mesh_a.vertices.iter() {
        let pos = *mesh_a.vertices.position(old_id);
        let nrm = *mesh_a.vertices.normal(old_id);
        remap_a.insert(old_id, combined.insert_or_weld(pos, nrm));
    }

    let mut remap_b: HashMap<VertexId, VertexId> = HashMap::with_capacity(mesh_b.vertices.len());
    for (old_id, _) in mesh_b.vertices.iter() {
        let pos = *mesh_b.vertices.position(old_id);
        let nrm = *mesh_b.vertices.normal(old_id);
        remap_b.insert(old_id, combined.insert_or_weld(pos, nrm));
    }

    let faces_a: Vec<FaceData> = mesh_a
        .faces
        .iter()
        .map(|face| FaceData {
            vertices: face.vertices.map(|vertex_id| remap_a[&vertex_id]),
            region: face.region,
        })
        .collect();
    let faces_b: Vec<FaceData> = mesh_b
        .faces
        .iter()
        .map(|face| FaceData {
            vertices: face.vertices.map(|vertex_id| remap_b[&vertex_id]),
            region: face.region,
        })
        .collect();

    (faces_a, faces_b)
}

fn remap_nary_face_soups(meshes: &[IndexedMesh], combined: &mut VertexPool) -> Vec<Vec<FaceData>> {
    use crate::domain::core::index::VertexId;
    use std::collections::HashMap;

    let mut face_soups = Vec::with_capacity(meshes.len());
    for mesh in meshes {
        let mut remap: HashMap<VertexId, VertexId> = HashMap::with_capacity(mesh.vertices.len());
        for (old_id, _) in mesh.vertices.iter() {
            let pos = *mesh.vertices.position(old_id);
            let nrm = *mesh.vertices.normal(old_id);
            remap.insert(old_id, combined.insert_or_weld(pos, nrm));
        }

        let faces: Vec<FaceData> = mesh
            .faces
            .iter()
            .map(|face| FaceData {
                vertices: face.vertices.map(|vertex_id| remap[&vertex_id]),
                region: face.region,
            })
            .collect();
        face_soups.push(faces);
    }
    face_soups
}

fn postprocess_boolean_mesh(
    result_faces: Vec<FaceData>,
    combined: &VertexPool,
    is_coplanar: bool,
) -> MeshResult<IndexedMesh> {
    let mut mesh = reconstruct::reconstruct_mesh(&result_faces, combined);
    mesh.recompute_normals();
    repair_boolean_mesh(&mut mesh, is_coplanar)?;
    Ok(mesh)
}

fn repair_boolean_mesh(mesh: &mut IndexedMesh, is_coplanar: bool) -> MeshResult<()> {
    repair_boolean_mesh_with_policy(mesh, is_coplanar, true)
}

fn repair_boolean_mesh_with_policy(
    mesh: &mut IndexedMesh,
    is_coplanar: bool,
    require_watertight: bool,
) -> MeshResult<()> {
    if !is_coplanar {
        mesh.rebuild_edges();
        let mut report = crate::application::watertight::check::check_watertight(
            &mesh.vertices,
            &mesh.faces,
            mesh.edges_ref().unwrap(),
        );
        if !report.is_watertight {
            if report.is_closed && !report.orientation_consistent {
                mesh.orient_outward();
                mesh.rebuild_edges();
                report = crate::application::watertight::check::check_watertight(
                    &mesh.vertices,
                    &mesh.faces,
                    mesh.edges_ref().unwrap(),
                );
            }

            if !report.is_watertight && report.is_closed && !report.orientation_consistent {
                let edge_store =
                    crate::infrastructure::storage::edge_store::EdgeStore::from_face_store(
                        &mesh.faces,
                    );
                let _ = crate::domain::topology::orientation::fix_orientation(
                    &mut mesh.faces,
                    &edge_store,
                );
                mesh.rebuild_edges();
                mesh.orient_outward();
                mesh.rebuild_edges();
                report = crate::application::watertight::check::check_watertight(
                    &mesh.vertices,
                    &mesh.faces,
                    mesh.edges_ref().unwrap(),
                );
            }

            if !report.is_watertight
                && report.non_manifold_edge_count == 0
                && report.boundary_edge_count > 0
                && report.boundary_edge_count <= 512
            {
                let edge_store =
                    crate::infrastructure::storage::edge_store::EdgeStore::from_face_store(
                        &mesh.faces,
                    );
                let added = crate::application::watertight::seal::seal_boundary_loops(
                    &mut mesh.vertices,
                    &mut mesh.faces,
                    &edge_store,
                    crate::domain::core::index::RegionId::INVALID,
                );
                if added > 0 {
                    mesh.rebuild_edges();
                    mesh.orient_outward();
                    mesh.rebuild_edges();
                    report = crate::application::watertight::check::check_watertight(
                        &mesh.vertices,
                        &mesh.faces,
                        mesh.edges_ref().unwrap(),
                    );
                }
            }

            if !report.is_watertight && report.boundary_edge_count > 0 {
                let improved = crate::application::watertight::repair::MeshRepair::iterative_boundary_stitch(
                    &mut mesh.faces,
                    &mesh.vertices,
                    3,
                );
                if improved > 0 {
                    mesh.rebuild_edges();
                    mesh.orient_outward();
                    mesh.rebuild_edges();
                    report = crate::application::watertight::check::check_watertight(
                        &mesh.vertices,
                        &mesh.faces,
                        mesh.edges_ref().unwrap(),
                    );
                }
            }

            if require_watertight && !report.is_watertight {
                return Err(MeshError::NotWatertight {
                    count: report.boundary_edge_count + report.non_manifold_edge_count,
                });
            }
        }

        mesh.orient_outward();
        if mesh.signed_volume() < 0.0 {
            mesh.flip_faces();
            mesh.rebuild_edges();
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{Cube, Cylinder, Disk, PrimitiveMesh, UvSphere};

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

    fn component_count(mesh: &mut IndexedMesh) -> usize {
        use crate::domain::topology::connectivity::connected_components;
        use crate::domain::topology::AdjacencyGraph;

        mesh.rebuild_edges();
        let adjacency = AdjacencyGraph::build(&mesh.faces, mesh.edges_ref().unwrap());
        connected_components(&mesh.faces, &adjacency).len()
    }

    fn symmetric_parallel_cylinders(segments: usize) -> (IndexedMesh, IndexedMesh) {
        let radius = 0.6;
        let height = 3.0;
        let separation = radius;
        let cyl_a = Cylinder {
            base_center: Point3r::new(-separation / 2.0, -height / 2.0, 0.0),
            radius,
            height,
            segments,
        }
        .build()
        .expect("symmetric cyl_a build");
        let cyl_b = Cylinder {
            base_center: Point3r::new(separation / 2.0, -height / 2.0, 0.0),
            radius,
            height,
            segments,
        }
        .build()
        .expect("symmetric cyl_b build");
        (cyl_a, cyl_b)
    }

    fn planar_branch(angle_from_x: f64, radius: f64, height: f64, segments: usize) -> IndexedMesh {
        use crate::application::csg::CsgNode;
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius,
            height,
            segments,
        }
        .build()
        .expect("branch build");
        let rotation = UnitQuaternion::<f64>::from_axis_angle(
            &Vector3::z_axis(),
            angle_from_x - std::f64::consts::FRAC_PI_2,
        );
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(Box::new(raw))),
            iso: Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rotation),
        }
        .evaluate()
        .expect("branch transform")
    }

    fn planar_trunk(radius: f64, height: f64, extension: f64, segments: usize) -> IndexedMesh {
        use crate::application::csg::CsgNode;
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius,
            height: height + extension,
            segments,
        }
        .build()
        .expect("trunk build");
        let rotation = UnitQuaternion::<f64>::from_axis_angle(
            &Vector3::z_axis(),
            -std::f64::consts::FRAC_PI_2,
        );
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(Box::new(raw))),
            iso: Isometry3::from_parts(Translation3::new(-height, 0.0, 0.0), rotation),
        }
        .evaluate()
        .expect("trunk transform")
    }

    fn quadfurcation_meshes() -> Vec<IndexedMesh> {
        let radius = 0.5;
        let height = 3.0;
        let extension = radius * 0.10;
        let segments = 32;
        let mut meshes = vec![planar_trunk(radius, height, extension, segments)];
        for angle_deg in [60.0_f64, 20.0, -20.0, -60.0] {
            meshes.push(planar_branch(angle_deg.to_radians(), radius, height, segments));
        }
        meshes
    }

    fn trifurcation_meshes() -> Vec<IndexedMesh> {
        let radius = 0.5;
        let height = 3.0;
        let extension = radius * 0.10;
        let segments = 32;
        let mut meshes = vec![planar_trunk(radius, height, extension, segments)];
        for angle_deg in [45.0_f64, 90.0, -45.0] {
            meshes.push(planar_branch(angle_deg.to_radians(), radius, height, segments));
        }
        meshes
    }

    fn pentafurcation_meshes() -> Vec<IndexedMesh> {
        let radius = 0.5;
        let height = 3.0;
        let extension = radius * 0.10;
        let segments = 32;
        let mut meshes = vec![planar_trunk(radius, height, extension, segments)];
        for angle_deg in [60.0_f64, 30.0, 0.0, -30.0, -60.0] {
            meshes.push(planar_branch(angle_deg.to_radians(), radius, height, segments));
        }
        meshes
    }

    // ── sphere × cylinder (curved × curved — arrangement pipeline) ─────────────

    #[test]
    fn sphere_cylinder_union_is_watertight() {
        let result =
            csg_boolean(BooleanOp::Union, &sphere(), &cylinder()).expect("sphere ∪ cylinder");
        assert_3d_watertight(result);
    }

    #[test]
    fn sphere_cylinder_intersection_is_watertight() {
        let result = csg_boolean(BooleanOp::Intersection, &sphere(), &cylinder())
            .expect("sphere ∩ cylinder");
        assert_3d_watertight(result);
    }

    #[test]
    fn sphere_cylinder_difference_is_watertight() {
        let result =
            csg_boolean(BooleanOp::Difference, &sphere(), &cylinder()).expect("sphere \\ cylinder");
        assert_3d_watertight(result);
    }

    // ── cube × cube (flat faces — intersecting arrangement pipeline) ───────────

    #[test]
    fn cube_cube_union_is_watertight() {
        let result = csg_boolean(BooleanOp::Union, &cube_a(), &cube_b()).expect("cube ∪ cube");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cube_intersection_is_watertight() {
        let result =
            csg_boolean(BooleanOp::Intersection, &cube_a(), &cube_b()).expect("cube ∩ cube");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cube_difference_is_watertight() {
        let result =
            csg_boolean(BooleanOp::Difference, &cube_a(), &cube_b()).expect("cube \\ cube");
        assert_3d_watertight(result);
    }

    // ── cube × cylinder coplanar (caps flush with cube walls) ──────────────────

    fn cylinder_coplanar() -> IndexedMesh {
        Cylinder {
            base_center: Point3r::new(0.0, -1.0, 0.0),
            radius: 0.4,
            height: 2.0,
            segments: 16,
        }
        .build()
        .expect("cylinder_coplanar build")
    }

    /// Difference of cube minus a coplanar cylinder must be watertight.
    /// The cylinder end caps are coplanar with the cube's top and bottom walls.
    /// The 2-D coplanar pipeline must subtract circular discs from the square
    /// walls, producing annular rings (tunnel openings).
    #[test]
    fn cube_cylinder_coplanar_difference_is_watertight() {
        let result = csg_boolean(BooleanOp::Difference, &cube_a(), &cylinder_coplanar())
            .expect("cube \\\\ cylinder_coplanar");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cylinder_coplanar_union_is_watertight() {
        let result = csg_boolean(BooleanOp::Union, &cube_a(), &cylinder_coplanar())
            .expect("cube ∪ cylinder_coplanar");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cylinder_coplanar_intersection_is_watertight() {
        let result = csg_boolean(BooleanOp::Intersection, &cube_a(), &cylinder_coplanar())
            .expect("cube ∩ cylinder_coplanar");
        assert_3d_watertight(result);
    }

    // ── disk × disk (coplanar — 2-D Sutherland-Hodgman pipeline) ───────────────
    // Disk operands are open surfaces; the coplanar path produces an open
    // surface result.  Only assert the operation completes without error.

    #[test]
    fn disk_disk_union_succeeds() {
        csg_boolean(BooleanOp::Union, &disk_a(), &disk_b()).expect("disk ∪ disk must not error");
    }

    #[test]
    fn disk_disk_intersection_succeeds() {
        csg_boolean(BooleanOp::Intersection, &disk_a(), &disk_b())
            .expect("disk ∩ disk must not error");
    }

    #[test]
    fn disk_disk_difference_succeeds() {
        csg_boolean(BooleanOp::Difference, &disk_a(), &disk_b())
            .expect("disk \\ disk must not error");
    }

    #[test]
    fn symmetric_parallel_cylinder_intersection_is_single_watertight_component() {
        let (cyl_a, cyl_b) = symmetric_parallel_cylinders(64);
        let mut result =
            csg_boolean(BooleanOp::Intersection, &cyl_a, &cyl_b).expect("symmetric intersection");

        result.rebuild_edges();
        let report = check_watertight(&result.vertices, &result.faces, result.edges_ref().unwrap());
        assert!(
            report.is_watertight,
            "symmetric cylinder intersection must be watertight: boundary={}, non_manifold={}",
            report.boundary_edge_count,
            report.non_manifold_edge_count
        );
        assert_eq!(
            component_count(&mut result),
            1,
            "symmetric cylinder intersection must remain a single component",
        );

        let radius = 0.6;
        let height = 3.0;
        let theta = std::f64::consts::FRAC_PI_3;
        let overlap_area = 2.0 * radius * radius * (theta - theta.sin() * theta.cos());
        let expected = height * overlap_area;
        let relative_error = (result.signed_volume() - expected).abs() / expected;
        assert!(
            relative_error < 0.01,
            "symmetric cylinder intersection volume error {:.2}% exceeds 1%",
            relative_error * 100.0
        );
    }

    #[test]
    fn indexed_nary_quadfurcation_union_is_watertight_without_component_dropping() {
        let mut result =
            csg_boolean_nary(BooleanOp::Union, &quadfurcation_meshes()).expect("quadfurcation union");
        assert_eq!(
            component_count(&mut result),
            1,
            "quadfurcation union must be a single connected component",
        );
        assert_3d_watertight(result);
    }

    #[test]
    fn indexed_nary_trifurcation_union_is_watertight_without_component_dropping() {
        let mut result =
            csg_boolean_nary(BooleanOp::Union, &trifurcation_meshes()).expect("trifurcation union");
        assert_eq!(
            component_count(&mut result),
            1,
            "trifurcation union must be a single connected component",
        );
        assert_3d_watertight(result);
    }

    #[test]
    fn indexed_nary_pentafurcation_union_is_watertight_without_component_dropping() {
        let mut result =
            csg_boolean_nary(BooleanOp::Union, &pentafurcation_meshes()).expect("pentafurcation union");
        assert_eq!(
            component_count(&mut result),
            1,
            "pentafurcation union must be a single connected component",
        );
        assert_3d_watertight(result);
    }

    #[test]
    fn indexed_nary_union_is_permutation_invariant() {
        let forward = quadfurcation_meshes();
        let mut reversed = quadfurcation_meshes();
        reversed.reverse();

        let mut forward_union =
            csg_boolean_nary(BooleanOp::Union, &forward).expect("forward quadfurcation union");
        let mut reversed_union =
            csg_boolean_nary(BooleanOp::Union, &reversed).expect("reversed quadfurcation union");

        assert_3d_watertight(forward_union.clone());
        assert_3d_watertight(reversed_union.clone());
        assert_eq!(
            component_count(&mut forward_union),
            component_count(&mut reversed_union),
            "operand order must not change the number of connected components",
        );

        let forward_volume = forward_union.signed_volume();
        let reversed_volume = reversed_union.signed_volume();
        let relative_error =
            (forward_volume - reversed_volume).abs() / forward_volume.abs().max(1.0e-12);
        assert!(
            relative_error < 0.005,
            "operand order changed union volume by {:.2}%",
            relative_error * 100.0
        );
    }

    // ── Y-junction trunk difference (curved × curved, Difference) ──────────
    // Diagnostic: verify watertight trunk difference has outward-only normals.
    // The BFS seed is the extremal (max-X) face — by the Jordan-Brouwer theorem
    // its outward normal must have nx ≥ 0, so BFS correctly orients the mesh.
    #[test]
    fn cylinder_difference_normals_check() {
        use crate::application::csg::CsgNode;
        use crate::application::quality::normals::analyze_normals;
        use crate::domain::core::scalar::Point3r;
        use crate::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};
        use std::f64::consts::FRAC_PI_2;

        const R: f64 = 0.5;
        const H_TRUNK: f64 = 3.0;
        const H_BRANCH: f64 = 3.0;
        const EPS: f64 = R * 0.10;
        const SEGS: usize = 32;
        let theta = std::f64::consts::FRAC_PI_4;

        let trunk = {
            let raw = Cylinder {
                base_center: Point3r::new(0.0, 0.0, 0.0),
                radius: R,
                height: H_TRUNK + EPS,
                segments: SEGS,
            }
            .build()
            .unwrap();
            let rot = UnitQuaternion::<f64>::from_axis_angle(&Vector3::z_axis(), -FRAC_PI_2);
            let iso = Isometry3::from_parts(Translation3::new(-H_TRUNK, 0.0, 0.0), rot);
            CsgNode::Transform {
                node: Box::new(CsgNode::Leaf(Box::new(raw))),
                iso,
            }
            .evaluate()
            .unwrap()
        };
        let branch_up = {
            let raw = Cylinder {
                base_center: Point3r::new(0.0, 0.0, 0.0),
                radius: R,
                height: H_BRANCH,
                segments: SEGS,
            }
            .build()
            .unwrap();
            let rot = UnitQuaternion::<f64>::from_axis_angle(&Vector3::z_axis(), theta - FRAC_PI_2);
            let iso = Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rot);
            CsgNode::Transform {
                node: Box::new(CsgNode::Leaf(Box::new(raw))),
                iso,
            }
            .evaluate()
            .unwrap()
        };
        let branch_dn = {
            let raw = Cylinder {
                base_center: Point3r::new(0.0, 0.0, 0.0),
                radius: R,
                height: H_BRANCH,
                segments: SEGS,
            }
            .build()
            .unwrap();
            let rot =
                UnitQuaternion::<f64>::from_axis_angle(&Vector3::z_axis(), -theta - FRAC_PI_2);
            let iso = Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rot);
            CsgNode::Transform {
                node: Box::new(CsgNode::Leaf(Box::new(raw))),
                iso,
            }
            .evaluate()
            .unwrap()
        };
        let branches = csg_boolean(BooleanOp::Union, &branch_up, &branch_dn).unwrap();
        let mut result = csg_boolean(BooleanOp::Difference, &trunk, &branches).unwrap();
        let normals_before = analyze_normals(&result);
        eprintln!(
            "before orient_outward: outward={}, inward={}, degen={}",
            normals_before.outward_faces,
            normals_before.inward_faces,
            normals_before.degenerate_faces,
        );
        eprintln!("is_watertight_before={}", result.is_watertight());
        result.orient_outward();
        let normals_after = analyze_normals(&result);
        eprintln!(
            "after  orient_outward: outward={}, inward={}, degen={}",
            normals_after.outward_faces, normals_after.inward_faces, normals_after.degenerate_faces,
        );
        eprintln!("is_watertight_after={}", result.is_watertight());
        assert_eq!(
            normals_after.inward_faces, 0,
            "orient_outward must eliminate inward faces"
        );

        // Single connected component — retain_largest_component must have
        // stripped the 2 × 8-face phantom islands from the trunk difference.
        {
            use crate::domain::topology::connectivity::connected_components;
            use crate::domain::topology::AdjacencyGraph;
            result.rebuild_edges();
            let edges = result.edges_ref().unwrap();
            let adj = AdjacencyGraph::build(&result.faces, edges);
            let comps = connected_components(&result.faces, &adj);
            eprintln!("components={}", comps.len());
            assert_eq!(
                comps.len(),
                1,
                "trunk difference must be a single connected component; \
                 got {} (phantom islands not removed)",
                comps.len(),
            );
        }
        // Euler characteristic χ = 2 for a single genus-0 closed body.
        {
            use crate::application::watertight::check::check_watertight;
            result.rebuild_edges();
            let rpt =
                check_watertight(&result.vertices, &result.faces, result.edges_ref().unwrap());
            eprintln!("euler_characteristic={:?}", rpt.euler_characteristic);
            assert_eq!(
                rpt.euler_characteristic,
                Some(2),
                "trunk difference must have Euler χ = 2; got {:?}",
                rpt.euler_characteristic,
            );
        }
    }

}
