//! Export: Disk → `outputs/primitives/disk.stl`
//!
//! Builds a flat circular disk (r = 1 mm, 64 segments) using
//! [`cfd_mesh::Disk`], validates mesh properties, and exports binary STL.
//!
//! ## Open-surface semantics
//!
//! Unlike closed primitives (sphere, cylinder, cube), a [`Disk`] is an
//! **open surface**: it has a boundary rim of `segments` edges and zero
//! signed volume.  The mesh is therefore intentionally **not watertight** —
//! `is_watertight` will be `false` and `boundary_edge_count == segments`.
//!
//! ## Geometry
//!
//! ```text
//!  r = 1 mm,  segments = 64
//!  A_N = (N/2) r² sin(2π/N)  →  A_64 ≈ π · 1²  (< 0.04% error)
//! ```
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_disk
//! ```

use std::f64::consts::PI;
use std::fs;
use std::io::BufWriter;

use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::infrastructure::storage::edge_store::EdgeStore;
use cfd_mesh::Disk;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Disk r=1 mm, 64 segments  (via Disk builder)");
    println!("  (Open surface — has boundary rim, zero volume by design)");
    println!("=================================================================");

    let r: f64 = 1.0;
    let n: usize = 64;

    let disk = Disk {
        center: cfd_mesh::domain::core::scalar::Point3r::origin(),
        radius: r,
        segments: n,
    };
    let mesh = disk.build()?;

    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);

    // Theoretical area of a regular N-gon inscribed in a circle of radius r:
    //   A_N = (N/2) · r² · sin(2π/N)
    let a_theoretical = (n as f64 / 2.0) * r * r * (2.0 * PI / n as f64).sin();
    let a_expected = PI * r * r;
    let a_err = (a_theoretical - a_expected).abs() / a_expected * 100.0;

    // Compute actual mesh area from face cross products.
    let a_mesh: f64 = mesh
        .faces
        .iter()
        .map(|f| {
            let a = mesh.vertices.position(f.vertices[0]);
            let b = mesh.vertices.position(f.vertices[1]);
            let c = mesh.vertices.position(f.vertices[2]);
            let ab = b - a;
            let ac = c - a;
            ab.cross(&ac).norm() * 0.5
        })
        .sum();
    let a_mesh_err = (a_mesh - a_expected).abs() / a_expected * 100.0;

    println!("  Vertices           : {}", mesh.vertices.len());
    println!(
        "  Faces              : {} (expected {})",
        mesh.faces.len(),
        n
    );
    println!();
    println!(
        "  Area (theoretical) : {:.6} mm²  (N-gon formula, err {:.4}%)",
        a_theoretical, a_err
    );
    println!(
        "  Area (mesh)        : {:.6} mm²  (from face cross products, err {:.4}%)",
        a_mesh, a_mesh_err
    );
    println!("  Expected π·r²      : {:.6} mm²", a_expected);
    println!();
    println!("  Open-surface properties:");
    println!(
        "  Watertight         : {}  (expected false — open surface)",
        report.is_watertight
    );
    println!(
        "  Boundary edges     : {:?}  (expected {})",
        report.boundary_edge_count, n
    );
    println!(
        "  Signed volume      : {:.6} mm³  (expected 0 — flat surface)",
        report.signed_volume
    );
    println!(
        "  Euler χ            : {:?}  (expected 1 for disk — open manifold)",
        report.euler_characteristic
    );

    // Sanity checks
    assert_eq!(
        mesh.face_count(),
        n,
        "Disk should have exactly {} triangles",
        n
    );
    assert_eq!(
        mesh.vertices.len(),
        n + 1,
        "Disk should have {} + 1 = {} vertices",
        n,
        n + 1
    );
    assert!(
        !report.is_watertight,
        "Disk is an open surface — not watertight"
    );
    // At N=64, A_64 = (64/2)·r²·sin(2π/64) ≈ 0.161% below πr² — the inscribed-polygon
    // discretisation error.  We assert < 0.2% to allow a comfortable margin at this resolution.
    assert!(
        a_mesh_err < 0.2,
        "Area error {:.4}% exceeds 0.2% at N={}",
        a_mesh_err,
        n
    );
    println!();
    println!("  All assertions passed.");

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;
    let stl_path = out_dir.join("disk.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL       : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
