#![cfg(feature = "csg")]
use cfd_mesh::core::index::VertexId;
use cfd_mesh::core::scalar::Point3r;
use cfd_mesh::csg::{csg_boolean, csg_boolean_indexed, BooleanOp};
use cfd_mesh::geometry::primitives::{Cylinder, PrimitiveMesh, UvSphere};
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::storage::face_store::FaceData;
use cfd_mesh::storage::vertex_pool::VertexPool;
use cfd_mesh::watertight::check::check_watertight;
use std::collections::HashMap;

fn csg_boolean_indexed_unchecked(
    op: BooleanOp,
    mesh_a: &cfd_mesh::IndexedMesh,
    mesh_b: &cfd_mesh::IndexedMesh,
) -> cfd_mesh::IndexedMesh {
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

    let result_faces = csg_boolean(op, &faces_a, &faces_b, &mut combined).unwrap();
    cfd_mesh::csg::reconstruct::reconstruct_mesh(&result_faces, &combined)
}

#[test]
fn diag_cube_cube_union() {
    let cube_a = cfd_mesh::Cube {
        origin: Point3r::new(0.0, 0.0, 0.0),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()
    .unwrap();
    let cube_b = cfd_mesh::Cube {
        origin: Point3r::new(1.0, 0.0, 0.0),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()
    .unwrap();

    for (label, op) in [
        ("Union", BooleanOp::Union),
        ("Intersection", BooleanOp::Intersection),
        ("Difference", BooleanOp::Difference),
    ] {
        let result = csg_boolean_indexed_unchecked(op, &cube_a, &cube_b);
        let edges = EdgeStore::from_face_store(&result.faces);
        let report = check_watertight(&result.vertices, &result.faces, &edges);
        eprintln!("=== {label} ===");
        eprintln!("Vertices: {}", result.vertex_count());
        eprintln!("Faces   : {}", result.face_count());
        eprintln!("Edges   : {}", edges.len());
        eprintln!("Watertight: {}", report.is_watertight);
        eprintln!("Closed    : {}", report.is_closed);
        eprintln!("Boundary  : {}", report.boundary_edge_count);
        eprintln!("NonManifld: {}", report.non_manifold_edge_count);
        eprintln!("Orient OK : {}", report.orientation_consistent);
        eprintln!(
            "Euler     : {:?} (expected {})",
            report.euler_characteristic, report.euler_expected
        );
        eprintln!("Volume    : {}", report.signed_volume);

        let mut plane_counts: std::collections::BTreeMap<String, usize> =
            std::collections::BTreeMap::new();
        for face in result.faces.iter() {
            let a = result.vertices.position(face.vertices[0]);
            let b = result.vertices.position(face.vertices[1]);
            let c = result.vertices.position(face.vertices[2]);
            let eps = 1e-9;
            let key = if (a.x - b.x).abs() < eps && (a.x - c.x).abs() < eps {
                format!("x={:.3}", a.x)
            } else if (a.y - b.y).abs() < eps && (a.y - c.y).abs() < eps {
                format!("y={:.3}", a.y)
            } else if (a.z - b.z).abs() < eps && (a.z - c.z).abs() < eps {
                format!("z={:.3}", a.z)
            } else {
                "other".to_string()
            };
            *plane_counts.entry(key).or_default() += 1;
        }
        eprintln!("Planes    : {:?}", plane_counts);

        if label != "Union" {
            for plane in ["x=1.000", "x=2.000", "y=2.000", "z=2.000"] {
                eprintln!("  Faces on {plane}:");
                for face in result.faces.iter() {
                    let a = result.vertices.position(face.vertices[0]);
                    let b = result.vertices.position(face.vertices[1]);
                    let c = result.vertices.position(face.vertices[2]);
                    let eps = 1e-9;
                    let key = if (a.x - b.x).abs() < eps && (a.x - c.x).abs() < eps {
                        format!("x={:.3}", a.x)
                    } else if (a.y - b.y).abs() < eps && (a.y - c.y).abs() < eps {
                        format!("y={:.3}", a.y)
                    } else if (a.z - b.z).abs() < eps && (a.z - c.z).abs() < eps {
                        format!("z={:.3}", a.z)
                    } else {
                        "other".to_string()
                    };
                    if key == plane {
                        eprintln!(
                            "    ({:.1},{:.1},{:.1}) ({:.1},{:.1},{:.1}) ({:.1},{:.1},{:.1})",
                            a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z
                        );
                    }
                }
            }
        }

        // Print all non-manifold edges with coordinates
        let mut bad_count = 0usize;
        for edge_data in edges.iter() {
            if edge_data.faces.len() != 2 {
                bad_count += 1;
                if bad_count <= 20 {
                    let (vi, vj) = edge_data.vertices;
                    let pi = result.vertices.position(vi);
                    let pj = result.vertices.position(vj);
                    eprintln!(
                        "  bad edge: ({:.4},{:.4},{:.4})-({:.4},{:.4},{:.4}) face_count={}",
                        pi.x,
                        pi.y,
                        pi.z,
                        pj.x,
                        pj.y,
                        pj.z,
                        edge_data.faces.len()
                    );
                }
            }
        }
        if bad_count > 20 {
            eprintln!("  ... {} more bad edges", bad_count - 20);
        }
        eprintln!();
    }
}

#[test]
fn diag_sphere_cylinder() {
    let sphere = UvSphere {
        radius: 1.0,
        center: Point3r::new(0.0, 0.0, 0.0),
        stacks: 64,
        segments: 32,
    }
    .build()
    .unwrap();
    let cylinder = Cylinder {
        base_center: Point3r::new(0.0, -1.5, 0.0),
        radius: 0.4,
        height: 3.0,
        segments: 64,
    }
    .build()
    .unwrap();

    for (label, op) in [
        ("Union", BooleanOp::Union),
        ("Intersection", BooleanOp::Intersection),
        ("Difference", BooleanOp::Difference),
    ] {
        let result = csg_boolean_indexed(op, &sphere, &cylinder).unwrap();
        let edges = EdgeStore::from_face_store(&result.faces);
        let report = check_watertight(&result.vertices, &result.faces, &edges);
        eprintln!("=== Sphere/Cylinder {label} ===");
        eprintln!("Vertices: {}", result.vertex_count());
        eprintln!("Faces   : {}", result.face_count());
        eprintln!("Boundary: {}", report.boundary_edge_count);
        eprintln!("NonMan  : {}", report.non_manifold_edge_count);
        eprintln!("OrientOK: {}", report.orientation_consistent);
        eprintln!("Volume  : {}", report.signed_volume);
    }
}
