use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::domain::core::scalar::Point3r;
use cfd_mesh::domain::geometry::primitives::{Disk, PrimitiveMesh};
use cfd_mesh::domain::topology::manifold;
use cfd_mesh::infrastructure::storage::edge_store::EdgeStore;

fn disk(cx: f64, r: f64, n: usize) -> cfd_mesh::IndexedMesh {
    Disk {
        center: Point3r::new(cx, 0., 0.),
        radius: r,
        segments: n,
    }
    .build()
    .unwrap()
}

fn check_mesh(name: &str, mesh: &cfd_mesh::IndexedMesh) {
    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = manifold::check_manifold(&edges);
    println!(
        "{}: faces={} boundary_edges={} non_manifold_edges={}",
        name,
        mesh.faces.len(),
        report.boundary_edges,
        report.non_manifold_edges
    );

    // Check for duplicate edges or unexpected boundary edges internal to the mesh.
    // In a simple 2D polygon like union of two overlapping disks, non_manifold should be 0.
    if report.non_manifold_edges > 0 {
        println!(
            "  WARNING: {} has {} non-manifold edges (duplicate/overlap?)",
            name, report.non_manifold_edges
        );
    }
}

fn main() {
    let (r, d) = (1.0_f64, 1.0_f64);
    let a = disk(-d / 2., r, 128);
    let b = disk(d / 2., r, 128);

    let u = csg_boolean_indexed(BooleanOp::Union, &a, &b).unwrap();
    check_mesh("Union", &u);

    let i = csg_boolean_indexed(BooleanOp::Intersection, &a, &b).unwrap();
    check_mesh("Intersection", &i);

    let diff = csg_boolean_indexed(BooleanOp::Difference, &a, &b).unwrap();
    check_mesh("Difference", &diff);
}
