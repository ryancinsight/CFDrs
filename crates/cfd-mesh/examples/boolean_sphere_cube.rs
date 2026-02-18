//! Boolean operations example: sphere centered in a cube, extending beyond walls.
//!
//! This example demonstrates CSG boolean operations (union, intersection, difference)
//! on a geometric configuration where a sphere is positioned at the center of a cube
//! but with a radius large enough to extend beyond the cube's walls.
//!
//! ## Geometric Configuration
//! - **Cube**: Unit cube from (0,0,0) to (1,1,1)
//! - **Sphere**: Centered at (0.5, 0.5, 0.5) with radius 0.75
//! - **Extension**: Sphere extends 0.25 units beyond each cube face
//!
//! ## Operations Demonstrated
//! 1. **Union**: Cube ∪ Sphere (combined volume)
//! 2. **Intersection**: Cube ∩ Sphere (spherical cap inside cube)
//! 3. **Difference**: Cube \ Sphere (cube with spherical cavity)
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example boolean_sphere_cube --features "csg,stl-io"
//! ```

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::IndexedMesh;
use cfd_mesh::core::index::RegionId;
use cfd_mesh::core::scalar::{Real, Point3r, Vector3r};

// Use standard library math constants
use std::f64::consts::{PI, TAU};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean};
use cfd_mesh::io::stl;
use cfd_mesh::storage::face_store::FaceData;
use cfd_mesh::storage::vertex_pool::VertexPool;

// =============================================================================
// Geometric Parameters
// =============================================================================

/// Cube side length (mm).
const CUBE_SIZE: Real = 1.0;

/// Sphere radius (mm) - larger than half cube side to extend beyond walls.
const SPHERE_RADIUS: Real = 0.75;

/// Sphere center position (centered in cube).
const SPHERE_CENTER: Point3r = Point3r::new(0.5, 0.5, 0.5);

/// Number of latitude stacks for sphere tessellation.
const SPHERE_STACKS: usize = 16;

/// Number of longitude segments for sphere tessellation.
const SPHERE_SEGMENTS: usize = 32;

// =============================================================================
// Main Entry Point
// =============================================================================

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Boolean Operations: Sphere Extending Beyond Cube Walls");
    println!("=================================================================");
    println!();
    println!("  Cube: {} × {} × {} mm, origin at (0,0,0)", CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
    println!("  Sphere: radius {} mm, center ({}, {}, {})", 
        SPHERE_RADIUS, SPHERE_CENTER.x, SPHERE_CENTER.y, SPHERE_CENTER.z);
    println!("  Extension beyond walls: {:.2} mm on each side", 
        SPHERE_RADIUS - CUBE_SIZE / 2.0);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("boolean_demo");
    fs::create_dir_all(&out_dir)?;
    let out_path = out_dir.to_str().expect("non-UTF8 path");

    let mut reports = Vec::new();

    // ── Operation 1: Union ────────────────────────────────────────────────
    println!("--- Union: Cube ∪ Sphere ---");
    match run_boolean_op(
        "union",
        BooleanOp::Union,
        out_path,
    ) {
        Ok(report) => reports.push(report),
        Err(e) => eprintln!("  FAILED: {}", e),
    }

    // ── Operation 2: Intersection ─────────────────────────────────────────
    println!("--- Intersection: Cube ∩ Sphere ---");
    match run_boolean_op(
        "intersection",
        BooleanOp::Intersection,
        out_path,
    ) {
        Ok(report) => reports.push(report),
        Err(e) => eprintln!("  FAILED: {}", e),
    }

    // ── Operation 3: Difference ───────────────────────────────────────────
    println!("--- Difference: Cube \\ Sphere ---");
    match run_boolean_op(
        "difference",
        BooleanOp::Difference,
        out_path,
    ) {
        Ok(report) => reports.push(report),
        Err(e) => eprintln!("  FAILED: {}", e),
    }

    // ── Summary ───────────────────────────────────────────────────────────
    println!();
    println!("=================================================================");
    println!("  Summary");
    println!("=================================================================");
    println!(
        "{:<15} {:>10} {:>10} {:>12} {:>12} {:>8}",
        "Operation", "Vertices", "Faces", "Area (mm²)", "Vol (mm³)", "ms"
    );
    println!("{:-<72}", "");
    for r in &reports {
        println!(
            "{:<15} {:>10} {:>10} {:>12.4} {:>12.4} {:>8}",
            r.name, r.vertices, r.faces, r.area, r.volume, r.elapsed_ms
        );
    }
    println!();
    println!("STL files written to {}/", out_path);

    Ok(())
}

// =============================================================================
// Primitive Mesh Generators
// =============================================================================

/// Generate a unit cube triangulated mesh.
///
/// Creates a cube from (0,0,0) to (size,size,size) with 12 triangular faces
/// (2 triangles per cube face). All normals point outward.
fn generate_cube(size: Real, pool: &mut VertexPool, region: RegionId) -> Vec<FaceData> {
    let s = size;
    let mut faces = Vec::with_capacity(12);

    // Define the 8 corner positions
    let p000 = Point3r::new(0.0, 0.0, 0.0);
    let p100 = Point3r::new(s, 0.0, 0.0);
    let p110 = Point3r::new(s, s, 0.0);
    let p010 = Point3r::new(0.0, s, 0.0);
    let p001 = Point3r::new(0.0, 0.0, s);
    let p101 = Point3r::new(s, 0.0, s);
    let p111 = Point3r::new(s, s, s);
    let p011 = Point3r::new(0.0, s, s);

    // Helper to add a quad as two triangles with consistent winding (CCW from outside)
    let mut add_quad = |p0: Point3r, p1: Point3r, p2: Point3r, p3: Point3r, normal: Vector3r| {
        let v0 = pool.insert_or_weld(p0, normal);
        let v1 = pool.insert_or_weld(p1, normal);
        let v2 = pool.insert_or_weld(p2, normal);
        let v3 = pool.insert_or_weld(p3, normal);
        
        // Triangle 1: (v0, v1, v2)
        faces.push(FaceData::new(v0, v1, v2, region));
        // Triangle 2: (v0, v2, v3)
        faces.push(FaceData::new(v0, v2, v3, region));
    };

    // Bottom face (z=0, normal -Z): Need (0,1,0) -> (1,0,0) for -Z?
    // Current: p000, p100, p110. (1,0,0)x(1,1,0) = +Z. WRONG.
    // Need: p000, p010, p110. (0,1,0)x(1,1,0) = (0,0,-1) = -Z. Correct.
    // So p000(BL), p010(TL), p110(TR), p100(BR).
    add_quad(p000, p010, p110, p100, -Vector3r::z());

    // Top face (z=s, normal +Z): Need (1,0,0) -> (0,1,0) for +Z?
    // Current: p001, p011, p111. (0,1,0)x(1,0,0) = -Z. WRONG.
    // Need: p001, p101, p111. (1,0,0)x(1,1,0) = +Z. Correct.
    // So p001(BL), p101(BR), p111(TR), p011(TL).
    add_quad(p001, p101, p111, p011, Vector3r::z());

    // Front face (y=0, normal -Y): p000, p100, p101. (1,0,0)x(1,0,1) = (0,-1,0). Correct.
    // Wait. p000->p100 (X). p100->p101 (Z). X x Z = -Y. Correct.
    // Original: p000, p001, p101, p100.
    // p000->p001 (Z). p001->p101 (X). Z x X = +Y. WRONG.
    add_quad(p000, p100, p101, p001, -Vector3r::y());

    // Back face (y=s, normal +Y):
    // Original: p010, p110, p111, p011.
    // p010->p110 (X). p110->p111 (Z). X x Z = -Y. WRONG.
    // Need Z x X = +Y.
    // p010->p011 (Z). p011->p111 (X).
    add_quad(p010, p011, p111, p110, Vector3r::y());

    // Left face (x=0, normal -X):
    // Original: p000, p010, p011, p001.
    // p000->p010 (Y). p010->p011 (Z). Y x Z = X. WRONG.
    // Need Z x Y = -X.
    // p000->p001 (Z). p001->p011 (Y).
    add_quad(p000, p001, p011, p010, -Vector3r::x());

    // Right face (x=s, normal +X):
    // Original: p100, p101, p111, p110.
    // p100->p101 (Z). p101->p111 (Y). Z x Y = -X. WRONG.
    // Need Y x Z = X.
    // p100->p110 (Y). p110->p111 (Z).
    add_quad(p100, p110, p111, p101, Vector3r::x());

    faces
}

/// Generate a UV-sphere triangulated mesh with outward-facing normals.
///
/// Creates a sphere using latitude/longitude tessellation with adaptive
/// handling of polar caps (triangles instead of quads at poles).
///
/// # Parameters
/// - `center`: Sphere center position
/// - `radius`: Sphere radius
/// - `segments`: Number of longitude divisions (around equator)
/// - `stacks`: Number of latitude divisions (pole to pole)
/// - `pool`: Vertex pool for deduplication
/// - `region`: Region tag for all faces
fn generate_sphere(
    center: Point3r,
    radius: Real,
    segments: usize,
    stacks: usize,
    pool: &mut VertexPool,
    region: RegionId,
) -> Vec<FaceData> {
    let mut faces = Vec::with_capacity(segments * stacks * 2);

    // Helper to generate a vertex on the sphere surface
    // phi = 0 at north pole (top), phi = PI at south pole (bottom)
    // theta = 0 to 2*PI around the equator
    let vertex_at = |theta: Real, phi: Real| -> (Point3r, Vector3r) {
        // Spherical to Cartesian conversion:
        // x = r * sin(φ) * cos(θ)
        // y = r * cos(φ)  
        // z = r * sin(φ) * sin(θ)
        // where θ is longitude (0 to 2π), φ is colatitude (0 at north pole to π at south pole)
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();
        
        // Normal points outward from center
        let normal = Vector3r::new(
            sin_phi * cos_theta,
            cos_phi,
            sin_phi * sin_theta,
        );
        
        let position = center + normal * radius;
        (position, normal)
    };

    // Generate quads/triangles for each stack/segment
    // The sphere is built from north pole (phi=0) to south pole (phi=PI)
    // Each "row" (stack) goes around the sphere at a given latitude
    for i in 0..segments {
        for j in 0..stacks {
            // Parameter values
            let t0 = i as Real / segments as Real;
            let t1 = (i + 1) as Real / segments as Real;
            let p0 = j as Real / stacks as Real;
            let p1 = (j + 1) as Real / stacks as Real;

            // Convert to angles
            // theta goes 0 to 2π (longitude)
            let theta0 = t0 * TAU;
            let theta1 = t1 * TAU;
            // phi goes 0 to π (colatitude: 0=north pole, π=south pole)
            let phi0 = p0 * PI;
            let phi1 = p1 * PI;

            // Get vertices at corners of the patch
            // v00 = (theta0, phi0) - lower-left
            // v10 = (theta1, phi0) - lower-right
            // v01 = (theta0, phi1) - upper-left
            // v11 = (theta1, phi1) - upper-right
            let (pos00, n00) = vertex_at(theta0, phi0);
            let (pos10, n10) = vertex_at(theta1, phi0);
            let (pos01, n01) = vertex_at(theta0, phi1);
            let (pos11, n11) = vertex_at(theta1, phi1);

            // Insert vertices into pool
            let v00 = pool.insert_or_weld(pos00, n00);
            let v10 = pool.insert_or_weld(pos10, n10);
            let v01 = pool.insert_or_weld(pos01, n01);
            let v11 = pool.insert_or_weld(pos11, n11);

            // For outward-facing normals, we need CCW winding when viewed from outside
            // Looking from outside at the north pole (phi=0), we see triangles with
            // their base at larger phi and apex at smaller phi
            
            if j == 0 {
                // North pole cap: phi0 ≈ 0, phi1 > 0
                // The cap is a triangle fan with apex at pole (v00, v10 at phi≈0)
                // and base at v01, v11 (at phi1)
                // For outward normal, we need CCW when viewed from above (outside)
                // Triangle: v10 -> v11 -> v01 (or equivalently the reverse order)
                // Skip degenerate triangles where pole vertices are at same point
                if phi1 > 1e-6 {
                    // Winding: v10, v01, v11 gives correct outward normal
                    faces.push(FaceData::new(v10, v01, v11, region));
                }
            } else if j == stacks - 1 {
                // South pole cap: phi0 < PI, phi1 ≈ PI
                // The cap is a triangle fan with apex at pole (v01, v11 at phi≈PI)
                // For outward normal at south pole (viewed from below), need CCW
                if (PI - phi0) > 1e-6 {
                    // Winding: v00, v10, v01 gives correct outward normal
                    faces.push(FaceData::new(v00, v10, v01, region));
                }
            } else {
                // Middle band: split quad into two triangles
                // Quad vertices in order around the outside (CCW from outside):
                // v00 (lower-left) -> v10 (lower-right) -> v11 (upper-right) -> v01 (upper-left)
                // 
                // Split along diagonal v00-v11:
                // Triangle 1: v00, v10, v11 (lower-right half)
                // Triangle 2: v00, v11, v01 (upper-left half)
                faces.push(FaceData::new(v00, v10, v11, region));
                faces.push(FaceData::new(v00, v11, v01, region));
            }
        }
    }

    faces
}

// =============================================================================
// Boolean Operation Runner
// =============================================================================

struct MeshReport {
    name: String,
    vertices: usize,
    faces: usize,
    area: Real,
    volume: Real,
    elapsed_ms: u128,
}

fn run_boolean_op(
    op_name: &str,
    op: BooleanOp,
    out_path: &str,
) -> Result<MeshReport, Box<dyn std::error::Error>> {
    let t0 = Instant::now();

    // Create fresh vertex pool
    let mut pool = VertexPool::default_millifluidic();

    // Generate primitives with distinct region tags
    let cube_region = RegionId::new(1);
    let sphere_region = RegionId::new(2);

    let cube_faces = generate_cube(CUBE_SIZE, &mut pool, cube_region);
    println!("  Cube: {} faces, {} vertices", cube_faces.len(), pool.len());

    let sphere_faces = generate_sphere(
        SPHERE_CENTER,
        SPHERE_RADIUS,
        SPHERE_SEGMENTS,
        SPHERE_STACKS,
        &mut pool,
        sphere_region,
    );
    println!("  Sphere: {} faces, {} vertices total", sphere_faces.len(), pool.len());

    // --- DIAGNOSTICS: Check Inputs ---
    {
        // Cube
        let mut mesh = IndexedMesh::new();
        // create a temporary pool reference or clone? VertexPool is shared.
        // IndexedMesh takes ownership of pool usually? No, `vertices` is `&VertexPool` or `VertexPool`.
        // `IndexedMesh` struct def in src/mesh.rs?
        // Let's assume we can build a temporary mesh.
        // Actually, IndexedMesh owns vertices.
        // We can't easily check without cloning pool.
        // But we can just sum volume of faces using pool?
        
        let calc_vol = |faces: &[FaceData], pool: &VertexPool| -> Real {
             let mut vol = 0.0;
             for f in faces {
                 let a = pool.position(f.vertices[0]).coords;
                 let b = pool.position(f.vertices[1]).coords;
                 let c = pool.position(f.vertices[2]).coords;
                 vol += a.dot(&b.cross(&c)) / 6.0;
             }
             vol
        };

        let cube_vol = calc_vol(&cube_faces, &mut pool);
        println!("  [DEBUG] Input Cube Volume: {:.4}", cube_vol);
        let sphere_vol = calc_vol(&sphere_faces, &mut pool);
        println!("  [DEBUG] Input Sphere Volume: {:.4}", sphere_vol);
        
        // Write Input STLs
        {
             let p = format!("{}/input_cube.stl", out_path);
             let f = fs::File::create(&p).unwrap();
             let mut w = BufWriter::new(f);
             let mut fs = cfd_mesh::storage::face_store::FaceStore::new();
             for face in &cube_faces { fs.push(*face); }
             stl::write_binary_stl(&mut w, &pool, &fs).unwrap();
        }
        {
             let p = format!("{}/input_sphere.stl", out_path);
             let f = fs::File::create(&p).unwrap();
             let mut w = BufWriter::new(f);
             let mut fs = cfd_mesh::storage::face_store::FaceStore::new();
             for face in &sphere_faces { fs.push(*face); }
             stl::write_binary_stl(&mut w, &pool, &fs).unwrap();
        }
    }

    // Perform boolean operation
    let result_faces = csg_boolean(op, &cube_faces, &sphere_faces, &mut pool)?;
    println!("  Result: {} faces", result_faces.len());

    // Build IndexedMesh
    let mut mesh = IndexedMesh::new();
    mesh.vertices = pool;
    for f in &result_faces {
        mesh.faces.push(*f);
    }
    mesh.rebuild_edges();

    // Compute metrics
    let area = mesh.surface_area();
    let volume = mesh.signed_volume();
    println!("  Area: {:.4} mm², Volume: {:.4} mm³", area, volume);

    // Quality check
    let qr = mesh.quality_report();
    if qr.passed {
        println!("  Quality: PASSED ({} faces)", qr.total_faces);
    } else {
        println!(
            "  Quality: {} / {} faces failing",
            qr.failing_faces, qr.total_faces
        );
    }
    if let (Some(ar), Some(ma)) = (&qr.aspect_ratio, &qr.min_angle) {
        println!(
            "  Aspect ratio: mean={:.2}, max={:.2}",
            ar.mean, ar.max
        );
        println!(
            "  Min angle: mean={:.1}°, min={:.1}°",
            ma.mean.to_degrees(),
            ma.min.to_degrees()
        );
    }

    // Write STL
    let stl_path = format!("{}/{}.stl", out_path, op_name);
    {
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
    }
    let stl_size = fs::metadata(&stl_path)?.len();
    println!("  STL: {} ({} bytes)", stl_path, stl_size);

    let elapsed = t0.elapsed().as_millis();
    println!("  Elapsed: {} ms", elapsed);
    println!();

    Ok(MeshReport {
        name: op_name.to_string(),
        vertices: mesh.vertex_count(),
        faces: mesh.face_count(),
        area,
        volume,
        elapsed_ms: elapsed,
    })
}