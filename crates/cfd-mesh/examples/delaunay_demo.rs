//! Delaunay triangulation demonstration.
//!
//! Demonstrates the core `cfd-mesh` Delaunay / CDT / Ruppert refinement pipeline:
//!
//! 1. **Delaunay triangulation** — Bowyer-Watson incremental insertion
//! 2. **Constrained Delaunay (CDT)** — PSLG constraint enforcement
//! 3. **Ruppert refinement** — Quality mesh generation with guarantees
//! 4. **Export to IndexedMesh** — Integration with the cfd-mesh ecosystem
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example delaunay_demo
//! ```

use std::f64::consts::PI;
use std::fs;
use std::io::BufWriter;

use cfd_mesh::application::delaunay::convert::indexed_mesh::to_indexed_mesh;
use cfd_mesh::application::delaunay::{Cdt, Pslg, RuppertRefiner, TriangleQuality};
use cfd_mesh::infrastructure::io::stl;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Delaunay Triangulation Demo (Bowyer-Watson + CDT + Ruppert)");
    println!("=================================================================");
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("delaunay");
    fs::create_dir_all(&out_dir)?;

    // ── 1. Simple Delaunay: Random Points ───────────────────────────────────
    println!("1. Simple Delaunay Triangulation (100 random points)");
    println!("   ----------------------------------------------------");

    let points: Vec<(f64, f64)> = (0..100)
        .map(|i| {
            // Deterministic pseudo-random via LCG
            let rng = 42_u64.wrapping_mul((i + 1) as u64);
            let x = ((rng >> 33) as f64 / (1u64 << 31) as f64) * 10.0;
            let y =
                ((rng.wrapping_mul(6364136223846793005) >> 33) as f64 / (1u64 << 31) as f64) * 10.0;
            (x, y)
        })
        .collect();

    // Build Delaunay triangulation directly from points
    let dt =
        cfd_mesh::application::delaunay::triangulation::DelaunayTriangulation::from_points(&points);

    println!("   Vertices : {}", dt.vertex_count());
    println!("   Triangles: {}", dt.triangle_count());
    println!("   Delaunay : {}", dt.is_delaunay());
    println!();

    // ── 2. CDT: Square with Hole ─────────────────────────────────────────────
    println!("2. CDT: Square Domain with Rectangular Hole");
    println!("   -----------------------------------------");

    let mut pslg = Pslg::new();

    // Outer square boundary (CCW)
    let v0 = pslg.add_vertex(0.0, 0.0);
    let v1 = pslg.add_vertex(10.0, 0.0);
    let v2 = pslg.add_vertex(10.0, 10.0);
    let v3 = pslg.add_vertex(0.0, 10.0);

    pslg.add_segment(v0, v1);
    pslg.add_segment(v1, v2);
    pslg.add_segment(v2, v3);
    pslg.add_segment(v3, v0);

    // Inner square hole (CW for hole)
    let h0 = pslg.add_vertex(3.0, 3.0);
    let h1 = pslg.add_vertex(7.0, 3.0);
    let h2 = pslg.add_vertex(7.0, 7.0);
    let h3 = pslg.add_vertex(3.0, 7.0);

    pslg.add_segment(h0, h3); // CW order for hole
    pslg.add_segment(h3, h2);
    pslg.add_segment(h2, h1);
    pslg.add_segment(h1, h0);

    // Mark the hole interior
    pslg.add_hole(5.0, 5.0);

    let cdt = Cdt::from_pslg(&pslg);

    println!("   PSLG vertices  : {}", pslg.vertices().len());
    println!("   PSLG segments  : {}", pslg.segments().len());
    println!(
        "   CDT triangles  : {}",
        cdt.triangulation().triangle_count()
    );
    println!("   CDT is Delaunay: {}", cdt.triangulation().is_delaunay());
    println!();

    // ── 3. CDT: L-Shaped Domain ─────────────────────────────────────────────
    println!("3. CDT: L-Shaped Domain");
    println!("   ---------------------");

    let mut pslg_l = Pslg::new();

    // L-shaped polygon (re-entrant corner)
    let l0 = pslg_l.add_vertex(0.0, 0.0);
    let l1 = pslg_l.add_vertex(2.0, 0.0);
    let l2 = pslg_l.add_vertex(2.0, 1.0);
    let l3 = pslg_l.add_vertex(1.0, 1.0);
    let l4 = pslg_l.add_vertex(1.0, 2.0);
    let l5 = pslg_l.add_vertex(0.0, 2.0);

    pslg_l.add_segment(l0, l1);
    pslg_l.add_segment(l1, l2);
    pslg_l.add_segment(l2, l3);
    pslg_l.add_segment(l3, l4);
    pslg_l.add_segment(l4, l5);
    pslg_l.add_segment(l5, l0);

    let cdt_l = Cdt::from_pslg(&pslg_l);

    println!("   PSLG vertices  : {}", pslg_l.vertices().len());
    println!(
        "   CDT triangles  : {}",
        cdt_l.triangulation().triangle_count()
    );
    println!(
        "   CDT is Delaunay: {}",
        cdt_l.triangulation().is_delaunay()
    );
    println!();

    // ── 4. Ruppert Refinement: Quality Mesh ─────────────────────────────────
    println!("4. Ruppert Refinement: Unit Square with Quality Guarantee");
    println!("   --------------------------------------------------------");

    let mut pslg_sq = Pslg::new();

    let s0 = pslg_sq.add_vertex(0.0, 0.0);
    let s1 = pslg_sq.add_vertex(1.0, 0.0);
    let s2 = pslg_sq.add_vertex(1.0, 1.0);
    let s3 = pslg_sq.add_vertex(0.0, 1.0);

    pslg_sq.add_segment(s0, s1);
    pslg_sq.add_segment(s1, s2);
    pslg_sq.add_segment(s2, s3);
    pslg_sq.add_segment(s3, s0);

    let cdt_sq = Cdt::from_pslg(&pslg_sq);

    // Count skinny triangles before refinement
    let skinny_before = count_skinny_triangles(&cdt_sq, 1.5);

    let mut refiner = RuppertRefiner::new(cdt_sq);
    refiner.set_max_ratio(1.414); // √2 — guarantees termination
    refiner.set_max_steiner(1000);
    let steiner_count = refiner.refine();

    let refined = refiner.into_cdt();

    // Count skinny triangles after refinement
    let skinny_after = count_skinny_triangles(&refined, 1.5);

    println!("   Steiner points added: {}", steiner_count);
    println!(
        "   Final triangles     : {}",
        refined.triangulation().triangle_count()
    );
    println!(
        "   Skinny triangles before (ratio > 1.5): {}",
        skinny_before
    );
    println!("   Skinny triangles after  (ratio > 1.5): {}", skinny_after);
    println!(
        "   Delaunay property   : {}",
        refined.triangulation().is_delaunay()
    );
    println!();

    // ── 5. Export to STL via IndexedMesh ────────────────────────────────────
    println!("5. Export: Refined Square → STL");
    println!("   -----------------------------");

    let indexed = to_indexed_mesh(refined.triangulation());

    let stl_path = out_dir.join("refined_square.stl");
    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &indexed.vertices, &indexed.faces)?;

    println!("   Output vertices: {}", indexed.vertices.len());
    println!("   Output faces   : {}", indexed.faces.len());
    println!("   STL written    : {}", stl_path.display());
    println!();

    // ── 6. Quality Statistics ───────────────────────────────────────────────
    println!("6. Quality Statistics for Refined Mesh");
    println!("   ------------------------------------");

    let stats = compute_quality_stats(refined.triangulation());

    println!("   Min angle       : {:.2}°", stats.min_angle_deg);
    println!("   Max angle       : {:.2}°", stats.max_angle_deg);
    println!("   Min radius/edge : {:.4}", stats.min_radius_edge);
    println!("   Max radius/edge : {:.4}", stats.max_radius_edge);
    println!("   Mean area       : {:.6}", stats.mean_area);
    println!();

    println!("=================================================================");
    println!("  Delaunay triangulation demo complete.");
    println!("=================================================================");

    Ok(())
}

/// Count triangles with radius-edge ratio exceeding the threshold.
fn count_skinny_triangles(cdt: &Cdt, max_ratio: f64) -> usize {
    let dt = cdt.triangulation();
    let mut count = 0;

    for (_, tri) in dt.interior_triangles() {
        let a = dt.vertex(tri.vertices[0]);
        let b = dt.vertex(tri.vertices[1]);
        let c = dt.vertex(tri.vertices[2]);

        let q = TriangleQuality::compute(a, b, c);
        if q.radius_edge_ratio > max_ratio {
            count += 1;
        }
    }

    count
}

/// Quality statistics for a triangulation.
struct QualityStats {
    min_angle_deg: f64,
    max_angle_deg: f64,
    min_radius_edge: f64,
    max_radius_edge: f64,
    mean_area: f64,
}

/// Compute quality statistics for a triangulation.
fn compute_quality_stats(
    dt: &cfd_mesh::application::delaunay::triangulation::DelaunayTriangulation,
) -> QualityStats {
    let mut min_angle: f64 = f64::MAX;
    let mut min_re: f64 = f64::MAX;
    let mut max_re: f64 = 0.0;
    let mut total_area: f64 = 0.0;
    let mut count = 0;

    for (_, tri) in dt.interior_triangles() {
        let a = dt.vertex(tri.vertices[0]);
        let b = dt.vertex(tri.vertices[1]);
        let c = dt.vertex(tri.vertices[2]);

        let q = TriangleQuality::compute(a, b, c);

        min_angle = min_angle.min(q.min_angle_deg());
        min_re = min_re.min(q.radius_edge_ratio);
        max_re = max_re.max(q.radius_edge_ratio);
        total_area += q.area;
        count += 1;
    }

    QualityStats {
        min_angle_deg: min_angle,
        max_angle_deg: 180.0 - 2.0 * min_angle, // Approximate from min_angle
        min_radius_edge: min_re,
        max_radius_edge: max_re,
        mean_area: if count > 0 {
            total_area / count as f64
        } else {
            0.0
        },
    }
}
