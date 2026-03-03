//! CSG (Constructive Solid Geometry) operations example
//!
//! Demonstrates union, intersection, and difference on representative
//! primitive pairs using the Mesh Arrangement pipeline.  Each result is
//! validated for volume accuracy, watertightness, Euler characteristic,
//! connected component count, and normal orientation.  Geometry issues
//! are flagged clearly so they can be investigated.
//!
//! ## Run
//!
//! ```sh
//! cargo run --example csg_operations --features csg
//! # Enable trace output for seam-repair diagnostics:
//! CSG_TRACE=1 cargo run --example csg_operations --features csg
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

#![allow(missing_docs)]

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh, UvSphere};
use cfd_mesh::domain::topology::connectivity::connected_components;
use cfd_mesh::domain::topology::AdjacencyGraph;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::{analyze_normals, Cube, IndexedMesh};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Operations: Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — volume + topology validation)");
    println!("=================================================================");
    println!();

    let out_dir = std::path::Path::new("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    // ── Pair 1: Cube × Cylinder ───────────────────────────────────────────────
    {
        println!("── Pair 1: Cube × Cylinder ──────────────────────────────────────");
        let r: f64 = 0.6;
        let h_cyl: f64 = 3.0;
        let h_cube: f64 = 2.0;
        let v_cube = h_cube.powi(3);
        let v_cyl = std::f64::consts::PI * r * r * h_cyl;
        let v_overlap = std::f64::consts::PI * r * r * h_cube;

        let cube = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: h_cube,
            height: h_cube,
            depth: h_cube,
        }
        .build()?;
        let cylinder = Cylinder {
            base_center: Point3r::new(1.0, -0.5, 1.0),
            radius: r,
            height: h_cyl,
            segments: 64,
        }
        .build()?;

        println!(
            "  Cube [0,2]³ mm  ({} faces)  +  Cylinder r={r}, h={h_cyl}  ({} faces)",
            cube.face_count(),
            cylinder.face_count()
        );
        println!();

        let t0 = Instant::now();
        let mut res = csg_boolean_indexed(BooleanOp::Union, &cube, &cylinder)?;
        report("Union (A ∪ B)", &mut res, v_cube + v_cyl - v_overlap, 0.05, t0.elapsed().as_millis());
        write_stl(&res, &out_dir.join("csg_ops_cube_cylinder_union.stl"))?;

        let t0 = Instant::now();
        let mut res = csg_boolean_indexed(BooleanOp::Intersection, &cube, &cylinder)?;
        report("Intersection (A ∩ B)", &mut res, v_overlap, 0.05, t0.elapsed().as_millis());
        write_stl(&res, &out_dir.join("csg_ops_cube_cylinder_intersection.stl"))?;

        let t0 = Instant::now();
        let mut res = csg_boolean_indexed(BooleanOp::Difference, &cube, &cylinder)?;
        report("Difference (A \\ B)", &mut res, v_cube - v_overlap, 0.05, t0.elapsed().as_millis());
        write_stl(&res, &out_dir.join("csg_ops_cube_cylinder_difference.stl"))?;
        println!();
    }

    // ── Pair 2: Sphere × Cylinder ─────────────────────────────────────────────
    {
        println!("── Pair 2: Sphere × Cylinder ────────────────────────────────────");
        let r_s: f64 = 1.0;
        let r_c: f64 = 0.4;
        let h_c: f64 = 3.0;
        let v_sphere = (4.0 / 3.0) * std::f64::consts::PI * r_s.powi(3);
        let v_cyl = std::f64::consts::PI * r_c * r_c * h_c;
        // Approximate overlap: cylinder clipped to sphere height 2*r_s
        let v_overlap_approx = std::f64::consts::PI * r_c * r_c * 2.0 * r_s;

        let sphere = UvSphere {
            radius: r_s,
            center: Point3r::origin(),
            segments: 32,
            stacks: 16,
        }
        .build()?;
        let cyl = Cylinder {
            base_center: Point3r::new(0.0, -r_s, 0.0),
            radius: r_c,
            height: h_c,
            segments: 32,
        }
        .build()?;

        println!(
            "  Sphere r={r_s}  ({} faces)  +  Cylinder r={r_c}, h={h_c}  ({} faces)",
            sphere.face_count(),
            cyl.face_count()
        );
        println!();

        let t0 = Instant::now();
        let mut res = csg_boolean_indexed(BooleanOp::Union, &sphere, &cyl)?;
        report(
            "Union (A ∪ B)",
            &mut res,
            v_sphere + v_cyl - v_overlap_approx,
            0.10, // 10% tolerance — approximate overlap
            t0.elapsed().as_millis(),
        );
        write_stl(&res, &out_dir.join("csg_ops_sphere_cylinder_union.stl"))?;

        let t0 = Instant::now();
        let mut res = csg_boolean_indexed(BooleanOp::Difference, &sphere, &cyl)?;
        report(
            "Difference (A \\ B)",
            &mut res,
            v_sphere - v_overlap_approx,
            0.10,
            t0.elapsed().as_millis(),
        );
        write_stl(&res, &out_dir.join("csg_ops_sphere_cylinder_difference.stl"))?;
        println!();
    }

    println!("  STL outputs written to outputs/csg/csg_ops_*.stl");
    println!("=================================================================");
    Ok(())
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn report(label: &str, mesh: &mut IndexedMesh, expected: f64, tol: f64, ms: u128) {
    let vol = mesh.signed_volume();
    let n = analyze_normals(mesh);
    let err = (vol - expected).abs() / expected.abs().max(1e-12);
    let vol_status = if err <= tol { "PASS" } else { "FAIL" };

    mesh.rebuild_edges();
    let wt = check_watertight(&mesh.vertices, &mesh.faces, mesh.edges_ref().unwrap());
    let adj = AdjacencyGraph::build(&mesh.faces, mesh.edges_ref().unwrap());
    let n_comps = connected_components(&mesh.faces, &adj).len();

    let chi_ok = wt.euler_characteristic == Some(2);
    let comps_ok = n_comps == 1;
    let norm_ok = n.inward_faces == 0;
    let any_issue = !wt.is_watertight || !chi_ok || !comps_ok || !norm_ok;

    println!("  ── {label} ──");
    println!("    Faces      : {}", mesh.face_count());
    println!("    Volume     : {vol:.4} mm³  (expected {expected:.4})");
    println!("    Vol error  : {:.2}%  [{vol_status}]", err * 100.0);
    println!(
        "    Watertight : {}  (boundary={}, non-manifold={})",
        wt.is_watertight, wt.boundary_edge_count, wt.non_manifold_edge_count
    );
    println!(
        "    Euler χ    : {:?}  (expected 2)  [{}]",
        wt.euler_characteristic,
        if chi_ok { "PASS" } else { "WARN" }
    );
    println!(
        "    Components : {n_comps}  [{}]",
        if comps_ok { "PASS" } else { "WARN phantom islands" }
    );
    println!(
        "    Normals    : outward={}, inward={} ({:.1}%), degen={}  [{}]",
        n.outward_faces,
        n.inward_faces,
        if mesh.face_count() > 0 {
            n.inward_faces as Real / mesh.face_count() as Real * 100.0
        } else {
            0.0
        },
        n.degenerate_faces,
        if norm_ok { "PASS" } else { "WARN" }
    );
    println!(
        "    Alignment  : mean={:.4}  min={:.4}",
        n.face_vertex_alignment_mean, n.face_vertex_alignment_min
    );
    println!("    Elapsed    : {ms} ms");

    if any_issue {
        println!("    *** GEOMETRY ISSUES DETECTED ***");
        if !wt.is_watertight {
            println!(
                "       - Not watertight: {} boundary + {} non-manifold edge(s)",
                wt.boundary_edge_count, wt.non_manifold_edge_count
            );
        }
        if !chi_ok {
            println!(
                "       - Euler χ = {:?} (expected 2): phantom islands or non-manifold topology",
                wt.euler_characteristic
            );
        }
        if !comps_ok {
            println!(
                "       - {} connected component(s): {} phantom island(s) present",
                n_comps,
                n_comps.saturating_sub(1)
            );
        }
        if !norm_ok {
            println!(
                "       - {}/{} face(s) with inward normals: winding order errors",
                n.inward_faces,
                mesh.face_count()
            );
        }
    }
}

fn write_stl(mesh: &IndexedMesh, path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = fs::File::create(path)?;
    let mut w = BufWriter::new(file);
    stl::write_binary_stl(&mut w, &mesh.vertices, &mesh.faces)?;
    Ok(())
}
