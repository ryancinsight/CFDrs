//! STL mesh inspection and watertightness analysis tool.
//!
//! ## Usage
//!
//! ```sh
//! inspect_stl <file.stl>
//! ```
//!
//! ## Output
//!
//! Reports mesh statistics, topology health (open/non-manifold edges),
//! volume and surface area, bounding box, degenerate face count, and
//! a normal orientation heuristic (outward vs. inward-facing faces).
//!
//! Exits with code 0 if the mesh is watertight and manifold, 1 otherwise.
//!
//! ## Algorithms
//!
//! - **Watertightness**: edge manifold test — each edge must appear exactly
//!   twice (once per adjacent face). Open edges indicate holes; ≥3 indicates
//!   non-manifold topology.
//! - **Volume**: signed tetrahedron sum from origin (Gauss divergence theorem):
//!   `V = (1/6) Σ p0·(p1 × p2)`
//! - **Normal orientation**: a face's geometric normal `(p1-p0)×(p2-p0)` should
//!   point outward (away from the mesh centroid). Inward-facing faces indicate
//!   a winding error.
//! - **Aspect ratio**: `longest_edge / shortest_edge` per triangle — high
//!   aspect ratios (> 10) indicate poor mesh quality.

#![allow(clippy::items_after_statements)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::field_reassign_with_default)]

use cfd_mesh::domain::core::scalar::{Point3r, Real, Vector3r};
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::BufReader;

/// Complete mesh analysis report.
#[derive(Debug, Default)]
struct StlStats {
    triangle_count: usize,
    vertex_count: usize,
    degenerate_faces: usize,
    open_edges: usize,
    non_manifold_edges: usize,
    inward_facing_faces: usize,
    volume: Real,
    area: Real,
    max_aspect_ratio: Real,
    mean_aspect_ratio: Real,
    min_corner: Point3r,
    max_corner: Point3r,
    centroid: Point3r,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: inspect_stl <file.stl>");
        eprintln!();
        eprintln!("  Analyzes an STL mesh file and reports:");
        eprintln!("    - Watertightness (open/non-manifold edges)");
        eprintln!("    - Volume and surface area");
        eprintln!("    - Bounding box and centroid");
        eprintln!("    - Degenerate face count");
        eprintln!("    - Normal orientation (outward vs inward)");
        eprintln!("    - Triangle quality (aspect ratio)");
        eprintln!();
        eprintln!("  Exit code: 0 = watertight manifold, 1 = issues detected");
        std::process::exit(1);
    }

    let path = &args[1];
    println!("═══════════════════════════════════════════════");
    println!("  STL Inspector: {path}");
    println!("═══════════════════════════════════════════════");

    let mut pool =
        cfd_mesh::infrastructure::storage::vertex_pool::VertexPool::default_millifluidic();
    let mut faces = cfd_mesh::infrastructure::storage::face_store::FaceStore::new();
    let region = cfd_mesh::domain::core::index::RegionId::from_usize(0);

    // Detect ASCII vs binary STL.
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut start = [0u8; 5];
    use std::io::Read;
    reader.read_exact(&mut start)?;

    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let is_ascii = &start == b"solid";

    if is_ascii {
        cfd_mesh::infrastructure::io::stl::read_ascii_stl(reader, &mut pool, &mut faces, region)?;
    } else {
        read_binary_stl(&mut reader, &mut pool, &mut faces, region)?;
    }

    println!("  Format: {}", if is_ascii { "ASCII" } else { "Binary" });
    println!();

    let stats = analyze_mesh(&pool, &faces);
    print_report(&stats);

    let has_issues =
        stats.open_edges > 0 || stats.non_manifold_edges > 0 || stats.degenerate_faces > 0;

    if has_issues {
        println!("  ✗  Mesh has issues — see report above.");
        std::process::exit(1);
    } else {
        println!("  ✓  Mesh is watertight and manifold.");
        std::process::exit(0);
    }
}

// ── Binary STL reader ─────────────────────────────────────────────────────────

fn read_binary_stl<R: std::io::Read>(
    reader: &mut R,
    pool: &mut cfd_mesh::infrastructure::storage::vertex_pool::VertexPool,
    faces: &mut cfd_mesh::infrastructure::storage::face_store::FaceStore,
    region: cfd_mesh::domain::core::index::RegionId,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut header = [0u8; 80];
    reader.read_exact(&mut header)?;

    let mut count_buf = [0u8; 4];
    reader.read_exact(&mut count_buf)?;
    let count = u32::from_le_bytes(count_buf);

    for _ in 0..count {
        // 12 bytes normal + 36 bytes vertices (3×12) + 2 bytes attribute
        let mut buf = [0u8; 50];
        reader.read_exact(&mut buf)?;

        let nx = f32::from_le_bytes(buf[0..4].try_into()?);
        let ny = f32::from_le_bytes(buf[4..8].try_into()?);
        let nz = f32::from_le_bytes(buf[8..12].try_into()?);
        let normal = Vector3r::new(Real::from(nx), Real::from(ny), Real::from(nz));

        let mut vs = [cfd_mesh::domain::core::index::VertexId::new(0); 3];
        for i in 0..3 {
            let o = 12 + i * 12;
            let vx = f32::from_le_bytes(buf[o..o + 4].try_into()?);
            let vy = f32::from_le_bytes(buf[o + 4..o + 8].try_into()?);
            let vz = f32::from_le_bytes(buf[o + 8..o + 12].try_into()?);
            let p = Point3r::new(Real::from(vx), Real::from(vy), Real::from(vz));
            vs[i] = pool.insert_or_weld(p, normal);
        }

        faces.push(cfd_mesh::infrastructure::storage::face_store::FaceData {
            vertices: vs,
            region,
        });
    }
    Ok(())
}

// ── Analysis ──────────────────────────────────────────────────────────────────

fn analyze_mesh(
    pool: &cfd_mesh::infrastructure::storage::vertex_pool::VertexPool,
    faces: &cfd_mesh::infrastructure::storage::face_store::FaceStore,
) -> StlStats {
    let mut stats = StlStats {
        min_corner: Point3r::new(Real::MAX, Real::MAX, Real::MAX),
        max_corner: Point3r::new(Real::MIN, Real::MIN, Real::MIN),
        ..Default::default()
    };
    stats.triangle_count = faces.len();
    stats.vertex_count = pool.len();

    let mut edge_counts: HashMap<(u32, u32), usize> = HashMap::new();
    let mut centroid_sum = Vector3r::zeros();
    let mut aspect_sum: Real = 0.0;

    // First pass: geometry stats and edge census.
    for (_, face) in faces.iter_enumerated() {
        let v = face.vertices;
        let p0 = *pool.position(v[0]);
        let p1 = *pool.position(v[1]);
        let p2 = *pool.position(v[2]);

        // Bounds
        for p in [p0, p1, p2] {
            stats.min_corner.x = stats.min_corner.x.min(p.x);
            stats.min_corner.y = stats.min_corner.y.min(p.y);
            stats.min_corner.z = stats.min_corner.z.min(p.z);
            stats.max_corner.x = stats.max_corner.x.max(p.x);
            stats.max_corner.y = stats.max_corner.y.max(p.y);
            stats.max_corner.z = stats.max_corner.z.max(p.z);
        }

        let e1 = p1 - p0;
        let e2 = p2 - p0;
        let e3 = p2 - p1;
        let cross = e1.cross(&e2);
        let area = cross.norm() * 0.5;
        stats.area += area;

        if area < 1e-9 {
            stats.degenerate_faces += 1;
        }

        // Signed volume contribution (Gauss divergence theorem).
        stats.volume += p0.coords.dot(&p1.coords.cross(&p2.coords)) / 6.0;

        // Centroid accumulation (weighted by area).
        let fc = (p0.coords + p1.coords + p2.coords) / 3.0;
        centroid_sum += fc * area;

        // Aspect ratio: longest / shortest edge.
        let l1 = e1.norm();
        let l2 = e2.norm();
        let l3 = e3.norm();
        let longest = l1.max(l2).max(l3);
        let shortest = l1.min(l2).min(l3);
        let ar = if shortest > 1e-20 {
            longest / shortest
        } else {
            Real::MAX
        };
        if ar.is_finite() && ar > stats.max_aspect_ratio {
            stats.max_aspect_ratio = ar;
        }
        aspect_sum += if ar.is_finite() { ar } else { 0.0 };

        // Edge manifold census.
        for i in 0..3 {
            let va = v[i].raw();
            let vb = v[(i + 1) % 3].raw();
            let key = (va.min(vb), va.max(vb));
            *edge_counts.entry(key).or_insert(0) += 1;
        }
    }

    // Centroid
    if stats.area > 1e-20 {
        let c = centroid_sum / stats.area;
        stats.centroid = Point3r::new(c.x, c.y, c.z);
    }
    stats.mean_aspect_ratio = if stats.triangle_count > 0 {
        aspect_sum / stats.triangle_count as Real
    } else {
        0.0
    };

    // Edge topology review.
    for &count in edge_counts.values() {
        if count == 1 {
            stats.open_edges += 1;
        } else if count > 2 {
            stats.non_manifold_edges += 1;
        }
    }

    // Second pass: normal orientation check.
    // A face is considered inward-facing if its geometric normal points
    // toward the mesh centroid rather than away from it.
    // This heuristic is reliable for convex and near-convex meshes.
    for (_, face) in faces.iter_enumerated() {
        let v = face.vertices;
        let p0 = *pool.position(v[0]);
        let p1 = *pool.position(v[1]);
        let p2 = *pool.position(v[2]);
        let geom_normal = (p1 - p0).cross(&(p2 - p0));
        let face_centroid = (p0.coords + p1.coords + p2.coords) / 3.0;
        let to_centroid = stats.centroid.coords - face_centroid;
        // Inward facing: normal points toward mesh centroid.
        if geom_normal.dot(&to_centroid) > 0.0 {
            stats.inward_facing_faces += 1;
        }
    }

    stats
}

// ── Report ────────────────────────────────────────────────────────────────────

fn print_report(stats: &StlStats) {
    let watertight = stats.open_edges == 0 && stats.non_manifold_edges == 0;
    let vol_sign = if stats.volume < 0.0 {
        " ⚠ INVERTED"
    } else {
        ""
    };
    let n_total = stats.triangle_count;
    let inward_pct = if n_total > 0 {
        100.0 * stats.inward_facing_faces as Real / n_total as Real
    } else {
        0.0
    };

    println!("  ── Geometry ─────────────────────────────────");
    println!("  Triangles         : {}", stats.triangle_count);
    println!("  Vertices          : {}", stats.vertex_count);
    println!(
        "  Bounds (min)      : [{:.4}, {:.4}, {:.4}]",
        stats.min_corner.x, stats.min_corner.y, stats.min_corner.z
    );
    println!(
        "  Bounds (max)      : [{:.4}, {:.4}, {:.4}]",
        stats.max_corner.x, stats.max_corner.y, stats.max_corner.z
    );
    println!(
        "  Centroid          : [{:.4}, {:.4}, {:.4}]",
        stats.centroid.x, stats.centroid.y, stats.centroid.z
    );
    println!("  Surface area      : {:.6}", stats.area);
    println!("  Volume            : {:.6}{vol_sign}", stats.volume);
    println!();

    println!("  ── Topology ─────────────────────────────────");
    println!(
        "  Watertight        : {}",
        if watertight { "YES ✓" } else { "NO  ✗" }
    );
    println!("  Open edges        : {}", stats.open_edges);
    println!("  Non-manifold edges: {}", stats.non_manifold_edges);
    println!();

    println!("  ── Quality ──────────────────────────────────");
    println!("  Degenerate faces  : {}", stats.degenerate_faces);
    println!(
        "  Inward normals    : {}  ({:.1}%)",
        stats.inward_facing_faces, inward_pct
    );
    println!("  Aspect ratio (max): {:.2}", stats.max_aspect_ratio);
    println!("  Aspect ratio (avg): {:.2}", stats.mean_aspect_ratio);
    println!("═══════════════════════════════════════════════");
    println!();
}
