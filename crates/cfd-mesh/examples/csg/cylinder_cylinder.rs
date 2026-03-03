//! CSG Cylinder–Cylinder (alias for symmetric variant).
//!
//! This file is kept for backward compatibility.  Prefer running
//! `csg_cylinder_cylinder_symmetric` or `csg_cylinder_cylinder_asymmetric`.
//!
//! Two equal-radius, equal-height Y-axis cylinders (r = 0.6 mm, h = 3 mm)
//! whose axes are parallel but offset by exactly one radius (Δx = 0.6 mm = r).
//! The cylinders partially overlap along their full height, forming a
//! **cylindrical lens** cross-section.  Both operands are curved and their cap
//! faces are **exactly coplanar** — this exercises the coplanar-face
//! classification path in the Mesh Arrangement pipeline.
//!
//! ## Geometry
//!
//! ```text
//! Cylinder A : base (−0.3, −1.5, 0), r = 0.6, h = 3 mm  →  Y ∈ [−1.5, 1.5]
//! Cylinder B : base (+0.3, −1.5, 0), r = 0.6, h = 3 mm  →  Y ∈ [−1.5, 1.5]
//!   axis separation d = 0.6 mm = r  →  θ = arccos(d / 2r) = arccos(0.5) = π/3
//!   cap planes exactly coplanar at y = ±1.5
//! ```
//!
//! The intersection cross-section is a symmetric lens of two equal circular
//! segments.  For d = r:
//!
//! ```text
//! θ         = π / 3
//! A_segment = r² (θ − sin θ cos θ)  =  r² (π/3 − √3/4)
//! A_lens    = 2 · A_segment
//! V_∩       = A_lens · h  =  2 h r² (π/3 − √3/4)
//! ```
//!
//! With r = 0.6, h = 3:
//!   A_lens ≈ 2 · 0.36 · (1.0472 − 0.4330) ≈ 0.4422 mm²
//!   V_∩    ≈ 0.4422 · 3 ≈ 1.3267 mm³
//!
//! | Operation | Expected (mm³)          | Pipeline    |
//! |-----------|------------------------|-------------|
//! | A ∪ B     | 2·V_cyl − V_∩          | Arrangement |
//! | A ∩ B     | V_∩ ≈ 1.3267           | Arrangement |
//! | A \ B     | V_cyl − V_∩ ≈ 2.0662   | Arrangement |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::{analyze_normals, IndexedMesh};
use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::topology::connectivity::connected_components;
use cfd_mesh::domain::topology::AdjacencyGraph;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cylinder–Cylinder: Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — both operands curved, coplanar caps)");
    println!("=================================================================");

    let r: f64 = 0.6; // both cylinder radii
    let h: f64 = 3.0; // both cylinder heights (equal → coplanar cap planes)
    let d: f64 = r; // axis separation = r  →  θ = π/3

    let v_cyl = std::f64::consts::PI * r * r * h;

    // Cylindrical-lens cross-section (d = r  →  θ = π/3):
    //   A_segment = r²(θ − sin θ · cos θ)
    //   V_∩ = 2 · h · A_segment
    let theta = (d / (2.0 * r)).acos(); // π/3
    let a_seg = r * r * (theta - theta.sin() * theta.cos());
    let v_intersect = 2.0 * h * a_seg;

    println!("  Cylinder A : base (−0.3,−1.5,0), r={r}, h={h}  V = {v_cyl:.4} mm³");
    println!("  Cylinder B : base (+0.3,−1.5,0), r={r}, h={h}  V = {v_cyl:.4} mm³");
    println!("  Axis sep d={d} mm (=r → θ=π/3); caps coplanar at y=±1.5");
    println!("  Overlap    : cylindrical lens cross-section      V = {v_intersect:.4} mm³");
    println!();
    println!("  Expected volumes:");
    println!("    Union        : {:.4} mm³", 2.0 * v_cyl - v_intersect);
    println!("    Intersection : {v_intersect:.4} mm³");
    println!("    Difference   : {:.4} mm³", v_cyl - v_intersect);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();
    let cyl_a = Cylinder {
        base_center: Point3r::new(-d / 2.0, -h / 2.0, 0.0),
        radius: r,
        height: h,
        segments: 64,
    }
    .build()?;
    let cyl_b = Cylinder {
        base_center: Point3r::new(d / 2.0, -h / 2.0, 0.0),
        radius: r,
        height: h,
        segments: 64,
    }
    .build()?;
    println!(
        "  Mesh built: {} + {} cylinder faces  ({} ms)",
        cyl_a.face_count(),
        cyl_b.face_count(),
        t_build.elapsed().as_millis()
    );
    println!();

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Union, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report(
            "Union (A ∪ B)",
            &mut result,
            2.0 * v_cyl - v_intersect,
            0.05,
            ms,
        );
        write_stl(&result, &out_dir.join("cylinder_cylinder_union.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Intersection, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &mut result, v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_intersection.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Difference, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report(
            "Difference (A \\ B)",
            &mut result,
            v_cyl - v_intersect,
            0.05,
            ms,
        );
        write_stl(&result, &out_dir.join("cylinder_cylinder_difference.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_difference.stl");
        println!();
    }

    println!("=================================================================");
    Ok(())
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn report(label: &str, mesh: &mut IndexedMesh, expected: f64, tol: f64, ms: u128) {
    let vol = mesh.signed_volume();
    let n = analyze_normals(mesh);
    let err = (vol - expected).abs() / expected.abs().max(1e-12);
    let vol_status = if err <= tol { "PASS" } else { "FAIL" };

    // Topology checks: Euler χ, connected components, boundary/non-manifold edges.
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
