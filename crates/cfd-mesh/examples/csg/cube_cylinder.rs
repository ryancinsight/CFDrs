//! CSG Cube–Cylinder: union, intersection, and difference
//!
//! A 2×2×2 mm axis-aligned cube and a cylinder (r = 0.6 mm, h = 3 mm) aligned
//! with the Y axis and centred on the cube's XZ face at `(1, −0.5, 1)`.
//! The cylinder punches straight through the cube from below to above.
//! Because one operand is curved the **Mesh Arrangement pipeline** is used.
//!
//! ## Geometry
//!
//! ```text
//! Cube      A : [0, 2]³ mm           V_A = 8 mm³
//! Cylinder  B : base (1, −0.5, 1),   V_B = π·r²·h = π·0.36·3 ≈ 3.3929 mm³
//!               r = 0.6, h = 3 (−Y to +Y through cube)
//! Overlap     : cylinder clipped to cube height [0, 2]
//!               V_∩ = π·r²·2 = π·0.36·2 ≈ 2.2619 mm³
//! ```
//!
//! | Operation | Expected (mm³)          | Pipeline    |
//! |-----------|------------------------|-------------|
//! | A ∪ B     | V_A + V_B − V_∩        | Arrangement |
//! | A ∩ B     | V_∩ ≈ 2.2619           | Arrangement |
//! | A \ B     | V_A − V_∩ ≈ 5.7381     | Arrangement |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cube_cylinder
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::application::csg::boolean::{csg_boolean, BooleanOp};
use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
use cfd_mesh::domain::topology::connectivity::connected_components;
use cfd_mesh::domain::topology::AdjacencyGraph;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::{analyze_normals, Cube, IndexedMesh};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cube–Cylinder: Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — curved surface)");
    println!("=================================================================");

    let r: f64 = 0.6;
    let h_cyl: f64 = 3.0; // cylinder height (extends ±0.5 mm beyond cube faces)
    let h_cube: f64 = 2.0; // cube height along Y

    let v_cube = h_cube.powi(3);
    let v_cyl = std::f64::consts::PI * r * r * h_cyl;
    // Overlap: cylinder clipped to the cube's Y extent [0, 2]
    let v_overlap = std::f64::consts::PI * r * r * h_cube;

    println!("  Cube     A : [0,2]³ mm                              V = {v_cube:.4} mm³");
    println!("  Cylinder B : base (1,−0.5,1), r={r}, h={h_cyl}  V = {v_cyl:.4} mm³");
    println!("  Overlap    : cylinder ∩ cube [Y=0..2]               V = {v_overlap:.4} mm³");
    println!();
    println!("  Expected volumes:");
    println!("    Union        : {:.4} mm³", v_cube + v_cyl - v_overlap);
    println!("    Intersection : {v_overlap:.4} mm³");
    println!("    Difference   : {:.4} mm³", v_cube - v_overlap);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();
    let cube = Cube {
        origin: Point3r::new(0.0, 0.0, 0.0),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()?;
    // Cylinder aligned with +Y, base at (1, −0.5, 1) so it extends from Y=−0.5 to Y=2.5,
    // passing cleanly through the full height of the cube.
    let cylinder = Cylinder {
        base_center: Point3r::new(1.0, -0.5, 1.0),
        radius: r,
        height: h_cyl,
        segments: 64,
    }
    .build()?;
    println!(
        "  Mesh built: {} cube + {} cylinder faces  ({} ms)",
        cube.face_count(),
        cylinder.face_count(),
        t_build.elapsed().as_millis()
    );
    println!();

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean(BooleanOp::Union, &cube, &cylinder)?;
        let ms = t0.elapsed().as_millis();
        report(
            "Union (A ∪ B)",
            &mut result,
            v_cube + v_cyl - v_overlap,
            0.05,
            ms,
            0, // genus 0 — no holes
        );
        write_stl(&result, &out_dir.join("cube_cylinder_union.stl"))?;
        println!("  STL: outputs/csg/cube_cylinder_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean(BooleanOp::Intersection, &cube, &cylinder)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &mut result, v_overlap, 0.05, ms, 0);
        write_stl(&result, &out_dir.join("cube_cylinder_intersection.stl"))?;
        println!("  STL: outputs/csg/cube_cylinder_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean(BooleanOp::Difference, &cube, &cylinder)?;
        let ms = t0.elapsed().as_millis();
        report(
            "Difference (A \\ B)",
            &mut result,
            v_cube - v_overlap,
            0.05,
            ms,
            1, // genus 1 — cylinder punches a through-hole (torus topology)
        );
        write_stl(&result, &out_dir.join("cube_cylinder_difference.stl"))?;
        println!("  STL: outputs/csg/cube_cylinder_difference.stl");
        println!();
    }

    println!("=================================================================");
    Ok(())
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn report(label: &str, mesh: &mut IndexedMesh, expected: f64, tol: f64, ms: u128, genus: i64) {
    let vol = mesh.signed_volume();
    let n = analyze_normals(mesh);
    let err = (vol - expected).abs() / expected.abs().max(1e-12);
    let vol_status = if err <= tol { "PASS" } else { "FAIL" };

    // Topology checks: Euler χ, connected components, boundary/non-manifold edges.
    mesh.rebuild_edges();
    let wt = check_watertight(&mesh.vertices, &mesh.faces, mesh.edges_ref().unwrap());
    let adj = AdjacencyGraph::build(&mesh.faces, mesh.edges_ref().unwrap());
    let n_comps = connected_components(&mesh.faces, &adj).len();

    let expected_chi = 2 * (1 - genus);
    let chi_ok = wt.euler_characteristic == Some(expected_chi);
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
        "    Euler χ    : {:?}  (expected {}, genus {})  [{}]",
        wt.euler_characteristic,
        expected_chi,
        genus,
        if chi_ok { "PASS" } else { "WARN" }
    );
    println!(
        "    Components : {n_comps}  [{}]",
        if comps_ok {
            "PASS"
        } else {
            "WARN phantom islands"
        }
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
                "       - Euler χ = {:?} (expected {}, genus {}): topology mismatch",
                wt.euler_characteristic, expected_chi, genus
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
