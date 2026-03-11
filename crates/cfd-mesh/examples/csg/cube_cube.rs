//! CSG Cube–Cube: union, intersection, and difference
//!
//! All three Boolean operations between two overlapping 2×2×2 mm axis-aligned
//! cubes.  Because both operands are flat-faced the **BSP pipeline** is used
//! (zero discretisation error for axis-aligned planes).
//!
//! ## Geometry
//!
//! ```text
//! Cube A : [0, 2]³  mm        V_A = 8 mm³
//! Cube B : [1, 3]×[0,2]×[0,2]   V_B = 8 mm³
//! Overlap: [1, 2]×[0,2]×[0,2]   V_∩ = 4 mm³   (1×2×2 slab)
//! ```
//!
//! | Operation | Expected (mm³) | Pipeline |
//! |-----------|---------------|----------|
//! | A ∪ B     | 8 + 8 − 4 = 12 | BSP (flat-face) |
//! | A ∩ B     | 4              | BSP |
//! | A \ B     | 8 − 4 = 4      | BSP |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cube_cube
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::application::csg::boolean::{csg_boolean, BooleanOp};
use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::domain::topology::connectivity::connected_components;
use cfd_mesh::domain::topology::AdjacencyGraph;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::{analyze_normals, Cube, IndexedMesh};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cube–Cube: Union | Intersection | Difference");
    println!("  (BSP pipeline — exact for flat-face geometry)");
    println!("=================================================================");

    // Cube A: [0,2]³    Cube B: [1,3]×[0,2]×[0,2]
    // Overlap slab: [1,2]×[0,2]×[0,2] = 1×2×2 = 4 mm³
    let v_a = 8.0_f64;
    let v_b = 8.0_f64;
    let v_overlap = 4.0_f64;

    println!("  Cube A : [0,2]³ mm              V = {v_a:.4} mm³");
    println!("  Cube B : [1,3]×[0,2]×[0,2] mm  V = {v_b:.4} mm³");
    println!("  Overlap: [1,2]×[0,2]×[0,2] mm  V = {v_overlap:.4} mm³");
    println!();
    println!("  Expected volumes:");
    println!("    Union        : {:.4} mm³", v_a + v_b - v_overlap);
    println!("    Intersection : {v_overlap:.4} mm³");
    println!("    Difference   : {:.4} mm³", v_a - v_overlap);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();
    let cube_a = Cube {
        origin: Point3r::new(0.0, 0.0, 0.0),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()?;
    let cube_b = Cube {
        origin: Point3r::new(1.0, 0.0, 0.0),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()?;
    println!(
        "  Mesh built: {} + {} faces  ({} ms)",
        cube_a.face_count(),
        cube_b.face_count(),
        t_build.elapsed().as_millis()
    );
    println!();

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean(BooleanOp::Union, &cube_a, &cube_b)?;
        let ms = t0.elapsed().as_millis();
        report(
            "Union (A ∪ B)",
            &mut result,
            v_a + v_b - v_overlap,
            0.01,
            ms,
        );
        write_stl(&result, &out_dir.join("cube_cube_union.stl"))?;
        println!("  STL: outputs/csg/cube_cube_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean(BooleanOp::Intersection, &cube_a, &cube_b)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &mut result, v_overlap, 0.01, ms);
        write_stl(&result, &out_dir.join("cube_cube_intersection.stl"))?;
        println!("  STL: outputs/csg/cube_cube_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean(BooleanOp::Difference, &cube_a, &cube_b)?;
        let ms = t0.elapsed().as_millis();
        report(
            "Difference (A \\ B)",
            &mut result,
            v_a - v_overlap,
            0.01,
            ms,
        );
        write_stl(&result, &out_dir.join("cube_cube_difference.stl"))?;
        println!("  STL: outputs/csg/cube_cube_difference.stl");
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
