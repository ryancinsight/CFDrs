//! CSG primitives demonstration
//!
//! Builds each supported primitive, validates volume, watertightness, Euler
//! characteristic, and normal orientation, then flags any geometry issues.
//! No boolean operations are performed — this is a baseline quality check
//! for the tessellation of each primitive type.
//!
//! ## Run
//!
//! ```sh
//! cargo run --example csg_primitives_demo --features csg
//! ```

#![allow(missing_docs)]

use std::f64::consts::PI;
use std::time::Instant;

use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::geometry::primitives::{Cone, Cylinder, PrimitiveMesh, UvSphere};
use cfd_mesh::domain::topology::connectivity::connected_components;
use cfd_mesh::domain::topology::AdjacencyGraph;
use cfd_mesh::{analyze_normals, Cube, IndexedMesh};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Primitives — volume, watertightness & topology validation");
    println!("=================================================================");
    println!();

    // ── Cube ──────────────────────────────────────────────────────────────────
    {
        let side = 2.0_f64;
        let expected_vol = side.powi(3);
        let t0 = Instant::now();
        let mut mesh = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: side,
            height: side,
            depth: side,
        }
        .build()?;
        report_primitive("Cube", &mut mesh, expected_vol, t0.elapsed().as_millis());
    }

    // ── UV Sphere ─────────────────────────────────────────────────────────────
    {
        let r = 1.0_f64;
        let expected_vol = (4.0 / 3.0) * PI * r.powi(3);
        let t0 = Instant::now();
        let mut mesh = UvSphere {
            radius: r,
            center: Point3r::origin(),
            segments: 32,
            stacks: 16,
        }
        .build()?;
        report_primitive(
            "UV Sphere (r=1, seg=32, stk=16)",
            &mut mesh,
            expected_vol,
            t0.elapsed().as_millis(),
        );
    }

    // ── Cylinder ──────────────────────────────────────────────────────────────
    {
        let r = 0.5_f64;
        let h = 2.0_f64;
        let expected_vol = PI * r * r * h;
        let t0 = Instant::now();
        let mut mesh = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: r,
            height: h,
            segments: 64,
        }
        .build()?;
        report_primitive(
            "Cylinder (r=0.5, h=2, seg=64)",
            &mut mesh,
            expected_vol,
            t0.elapsed().as_millis(),
        );
    }

    // ── Cone ──────────────────────────────────────────────────────────────────
    {
        let r = 1.0_f64;
        let h = 2.0_f64;
        let expected_vol = PI * r * r * h / 3.0;
        let t0 = Instant::now();
        let mut mesh = Cone {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: r,
            height: h,
            segments: 64,
        }
        .build()?;
        report_primitive(
            "Cone (r=1, h=2, seg=64)",
            &mut mesh,
            expected_vol,
            t0.elapsed().as_millis(),
        );
    }

    println!("=================================================================");
    Ok(())
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn report_primitive(label: &str, mesh: &mut IndexedMesh, expected_vol: f64, ms: u128) {
    let tol = 0.02; // 2% volume tolerance for tessellation approximation
    let vol = mesh.signed_volume();
    let n = analyze_normals(mesh);
    let err = (vol - expected_vol).abs() / expected_vol.abs().max(1e-12);
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
    println!("    Volume     : {vol:.4} mm³  (expected {expected_vol:.4})");
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
        if comps_ok { "PASS" } else { "WARN" }
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
    println!("    Built in   : {ms} ms");

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
                "       - Euler χ = {:?} (expected 2): topology defect in primitive tessellation",
                wt.euler_characteristic
            );
        }
        if !comps_ok {
            println!(
                "       - {} connected component(s): unexpected disconnected patches",
                n_comps
            );
        }
        if !norm_ok {
            println!(
                "       - {}/{} face(s) with inward normals: winding error in primitive",
                n.inward_faces,
                mesh.face_count()
            );
        }
    }
    println!();
}
