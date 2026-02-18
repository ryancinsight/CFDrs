//! Provable Mesh Diagnostics
//!
//! This example generates meshes with specific, known defects and verifies that
//! the `cfd-mesh` diagnostic tools correctly identify them.
//!
//! ## Test Cases
//! 1. **Manifold Cube**: Control group (should pass all checks).
//! 2. **Punctured Cube**: Missing one face (should fail closedness).
//! 3. **Non-Manifold Edge**: "Book" structure (3 faces sharing an edge).
//! 4. **Inverted Face**: One face has wrong winding (should fail orientation).
//! 5. **Degenerate Face**: Zero-area triangle (should fail quality).
//! 6. **Bow Tie**: Self-intersecting geometry (single vertex singularity).

use std::fs;
use std::io::BufWriter;
use std::path::Path;

use cfd_mesh::IndexedMesh;
use cfd_mesh::core::index::{FaceId, RegionId};
use cfd_mesh::core::scalar::{Point3r, Vector3r};
use cfd_mesh::storage::face_store::FaceData;
use cfd_mesh::storage::vertex_pool::VertexPool;
use cfd_mesh::io::stl;
use cfd_mesh::quality::validation::MeshValidator;
use cfd_mesh::watertight::check::check_watertight;

// =============================================================================
// Main
// =============================================================================

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Provable Mesh Diagnostics");
    println!("=================================================================");
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("provable_diagnostics");
    fs::create_dir_all(&out_dir)?;

    run_test("Manifold Cube", generate_cube, &out_dir, true)?;
    run_test("Punctured Cube", generate_punctured_cube, &out_dir, false)?;
    run_test("Non-Manifold Edge", generate_book, &out_dir, false)?;
    run_test("Inverted Face", generate_inverted_face_cube, &out_dir, false)?;
    run_test("Degenerate Face", generate_degenerate_face, &out_dir, false)?;
    // run_test("Bow Tie", generate_bow_tie, &out_dir, false)?; // Optional extra

    println!("=================================================================");
    println!("  All diagnostics completed.");
    println!("=================================================================");

    Ok(())
}

fn run_test<F>(
    name: &str,
    generator: F,
    out_dir: &Path,
    expect_success: bool,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: Fn() -> IndexedMesh,
{
    println!("--- Test Case: {} ---", name);
    let mut mesh = generator();
    
    // Ensure topology is built
    mesh.rebuild_edges();
    mesh.recompute_normals(); // For visualization

    // 1. Watertightness Check
    // Use edges_ref() to avoid mutable borrow conflict
    let edges = mesh.edges_ref().expect("Edges should be built");
    let wt = check_watertight(&mesh.vertices, &mesh.faces, edges);
    
    println!("  [Watertight] Total Edges: {}", edges.len());
    println!("  [Watertight] Boundary Edges: {}", wt.boundary_edge_count);
    println!("  [Watertight] Non-Manifold Edges: {}", wt.non_manifold_edge_count);
    println!("  [Watertight] Closed Manifold: {}", wt.is_closed);
    println!("  [Watertight] Consistent Orientation: {}", wt.orientation_consistent);
    println!("  [Watertight] Signed Volume: {:.4}", wt.signed_volume);
    
    // 2. Quality Check
    let validator = MeshValidator::default();
    let quality = validator.validate(&mesh.faces, &mesh.vertices);
    println!("  [Quality]    Failing Faces: {} / {}", quality.failing_faces, quality.total_faces);

    // 3. Topology (Euler Characteristic)
    // Chi = V - E + F. For a closed sphere-like surface (genus 0), Chi should be 2.
    let n_v = mesh.vertices.len();
    let n_e = edges.len();
    let n_f = mesh.faces.len();
    let chi = (n_v as i64) - (n_e as i64) + (n_f as i64);
    println!("  [Topology]   V={}, E={}, F={}", n_v, n_e, n_f);
    println!("  [Topology]   Euler Characteristic (Chi): {}", chi);

    // 4. Overall Verification
    let passed = wt.is_watertight && quality.passed && chi == 2;
    let status = if passed { "PASSED" } else { "FAILED" };
    
    println!("  [Result]     Status: {}", status);
    
    let mut repair_attempted = false;
    let mut repair_success = false;

    if !passed {
        println!("\n  [Repair]     Attempting repairs...");
        
        let mut changes_made = false;

        // Repair 1: Degenerate Faces
        if quality.failing_faces > 0 {
            use cfd_mesh::watertight::repair::MeshRepair;
            let removed = MeshRepair::remove_degenerate_faces(&mesh.faces, &mesh.vertices);
            if !removed.is_empty() {
                println!("  [Repair]     Removed {} degenerate faces.", removed.len());
                // In a real app we'd need to remove these from face_store.
                // Since this example doesn't have a mutable remove API readily exposed on IndexedMesh for FaceId list,
                // we simulates it or needs to extend API. 
                // Wait, we don't have easy removal by ID list in the current API exposed to example.
                // Let's check if we can filter the faces.
                // Actually, let's just use the fact that we can rebuild the mesh or use internal features if public.
                // For this example, let's skip complex removal if API is missing and just note it, 
                // OR we can implement a simple removal loop if FaceStore supports swap_remove or similar.
                // Checking FaceStore - it has `faces` Vec.
                // We will implement a naive reconstruction for this example.
                
                let mut new_faces = std::collections::HashSet::new();
                for f in mesh.faces.as_slice() {
                    new_faces.insert(f.clone());
                }
                
                let mut valid_faces = Vec::new();
                let removed_set: std::collections::HashSet<_> = removed.iter().collect();
                for (i, face) in mesh.faces.iter().enumerate() {
                    let fid = FaceId::from_usize(i);
                    if !removed_set.contains(&fid) {
                        valid_faces.push(face.clone());
                    }
                }
                mesh.faces.clear();
                for f in valid_faces {
                    mesh.faces.push(f);
                }
                changes_made = true;
            }
        }

        // Repair 2: Closure (Sealing)
        if !wt.is_closed {
            use cfd_mesh::watertight::seal::seal_boundary_loops;
            // Edges need to be fresh
            mesh.rebuild_edges();
            // We need to clone the edges to break the borrow on mesh.
            // EdgeStore implements Clone.
            let edges_clone = mesh.edges_ref().unwrap().clone();
            
            let added = seal_boundary_loops(
                &mut mesh.vertices,
                &mut mesh.faces,
                &edges_clone,
                RegionId::new(0)
            );
            if added > 0 {
                println!("  [Repair]     Sealed boundary loops: {} new faces.", added);
                changes_made = true;
                // Edges are now stale, force rebuild next time if needed
            }
        }

        // Repair 3: Orientation
        // Note: Orientation check might be invalid if mesh is not closed/manifold.
        // We re-check topology first.
        mesh.rebuild_edges();
        let edges_for_check = mesh.edges_ref().unwrap();
        let wt_after_seal = check_watertight(&mesh.vertices, &mesh.faces, edges_for_check);
        
        if wt_after_seal.is_closed && !wt_after_seal.orientation_consistent {
             use cfd_mesh::watertight::repair::MeshRepair;
             // We need mutable access to faces and existing edges.
             // Edges are safe to clone as they are just indices.
             let edges_clone = mesh.edges_ref().unwrap().clone();
             let flipped = MeshRepair::fix_orientations(&mut mesh.faces, &edges_clone);
             
             if flipped > 0 {
                 println!("  [Repair]     Flipped {} faces for consistency.", flipped);
                 changes_made = true;
             }
             
             // Check volume ensuring outward normals
             // fix_orientations ensures consistency, but might result in all-inward normals (negative volume).
             // If so, we flip everything.
             let current_volume = compute_signed_volume(&mesh.vertices, &mesh.faces);
             if current_volume < 0.0 {
                 println!("  [Repair]     Volume is negative ({:.4}). Flipping all faces to ensure outward normals.", current_volume);
                 for face in mesh.faces.iter_mut() {
                     face.flip();
                 }
                 changes_made = true;
             }
        }

        if changes_made {
             repair_attempted = true;
             // Re-verify
             mesh.rebuild_edges();
             mesh.recompute_normals();
             
             let wt2 = check_watertight(&mesh.vertices, &mesh.faces, mesh.edges_ref().unwrap());
             let q2 = validator.validate(&mesh.faces, &mesh.vertices);
             
             let n_v2 = mesh.vertices.len();
             let n_e2 = mesh.edges_ref().unwrap().len();
             let n_f2 = mesh.faces.len();
             let chi2 = (n_v2 as i64) - (n_e2 as i64) + (n_f2 as i64);
             
             repair_success = wt2.is_watertight && q2.passed && chi2 == 2;
             let new_status = if repair_success { "REPAIRED" } else { "STILL FAILED" };
             println!("  [Repair]     Result: {}", new_status);
             
             if repair_success {
                 // Write repaired STL
                 let filename = name.replace(" ", "_").to_lowercase();
                 let stl_path = out_dir.join(format!("{}_repaired.stl", filename));
                 let file = fs::File::create(&stl_path)?;
                 let mut writer = BufWriter::new(file);
                 stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
                 println!("  [Output]     Wrote {}", stl_path.display());
             }
        } else {
             println!("  [Repair]     No repair strategy available or applicable.");
        }
    }

    if expect_success && !passed {
        println!("  [WARNING]    Expected PASS but got FAIL.");
    } else if !expect_success && passed {
        println!("  [WARNING]    Expected FAIL but got PASS.");
    } else if !expect_success && repair_attempted && !repair_success {
         println!("  [INFO]       Repair attempted but failed (structure may be unfixable automatically).");
    } else if !expect_success && repair_success {
         println!("  [SUCCESS]    Defect was automatically repaired.");
    } else {
        println!("  [OK]         Outcome matches expectation.");
    }

    println!(); // Spacing

    Ok(())
}

fn compute_signed_volume(vertices: &VertexPool, faces: &cfd_mesh::storage::face_store::FaceStore) -> f64 {
    let mut vol = 0.0;
    for face in faces.iter() {
        let p0 = vertices.position(face.vertices[0]);
        let p1 = vertices.position(face.vertices[1]);
        let p2 = vertices.position(face.vertices[2]);
        // V = (p0 . (p1 x p2)) / 6
        // Need to use coords for vector operations on Points
        vol += p0.coords.dot(&p1.coords.cross(&p2.coords)) / 6.0;
    }
    vol as f64
}

// =============================================================================
// Generators
// =============================================================================

/// 1. Manifold Cube: Standard 12-triangle cube.
fn generate_cube() -> IndexedMesh {
    let mut pool = VertexPool::default();
    let region = RegionId::new(1);
    let s = 1.0;
    
    // Vertices
    let p000 = Point3r::new(0.0, 0.0, 0.0);
    let p100 = Point3r::new(s, 0.0, 0.0);
    let p110 = Point3r::new(s, s, 0.0);
    let p010 = Point3r::new(0.0, s, 0.0);
    let p001 = Point3r::new(0.0, 0.0, s);
    let p101 = Point3r::new(s, 0.0, s);
    let p111 = Point3r::new(s, s, s);
    let p011 = Point3r::new(0.0, s, s);

    let mut faces = Vec::new();
    // Helper
    let mut add_quad = |p0, p1, p2, p3, n| {
        let v0 = pool.insert_or_weld(p0, n);
        let v1 = pool.insert_or_weld(p1, n);
        let v2 = pool.insert_or_weld(p2, n);
        let v3 = pool.insert_or_weld(p3, n);
        faces.push(FaceData::new(v0, v1, v2, region));
        faces.push(FaceData::new(v0, v2, v3, region));
    };

    add_quad(p000, p010, p110, p100, -Vector3r::z()); // Bottom
    add_quad(p001, p101, p111, p011, Vector3r::z());  // Top
    add_quad(p000, p100, p101, p001, -Vector3r::y()); // Front
    add_quad(p010, p011, p111, p110, Vector3r::y());  // Back
    add_quad(p000, p001, p011, p010, -Vector3r::x()); // Left
    add_quad(p100, p110, p111, p101, Vector3r::x());  // Right

    let mut mesh = IndexedMesh::new();
    mesh.vertices = pool;
    for f in faces { mesh.faces.push(f); }
    mesh
}

/// 2. Punctured Cube: Manifold cube minus the top face (2 triangles).
fn generate_punctured_cube() -> IndexedMesh {
    let mut mesh = generate_cube();
    // Remove the last 2 faces (Top face)
    let _ = mesh.faces.pop();
    let _ = mesh.faces.pop();
    mesh
}

/// 3. Non-Manifold Edge: Two triangles sharing an edge, plus a third "fin" triangle sharing the same edge.
/// Like a book binding with 3 pages.
fn generate_book() -> IndexedMesh {
    let mut pool = VertexPool::default();
    let region = RegionId::new(1);
    
    // Shared edge vertices
    let v_a = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), Vector3r::z());
    let v_b = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), Vector3r::z());

    // Page vertices
    let v_p1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), Vector3r::z());
    let v_p2 = pool.insert_or_weld(Point3r::new(-1.0, 0.0, 0.0), Vector3r::z());
    let v_p3 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 1.0), Vector3r::x());

    let mut mesh = IndexedMesh::new();
    mesh.vertices = pool;
    mesh.faces.push(FaceData::new(v_a, v_p1, v_b, region));
    mesh.faces.push(FaceData::new(v_a, v_b, v_p2, region)); // Note winding might be inconsistent, but edge is shared 3x
    mesh.faces.push(FaceData::new(v_a, v_b, v_p3, region));
    
    mesh
}

/// 4. Inverted Face: Cube where one face has flipped winding.
fn generate_inverted_face_cube() -> IndexedMesh {
    let mut mesh = generate_cube();
    // Flip the first face
    let f = mesh.faces.get_mut(FaceId::from_usize(0));
    let tmp = f.vertices[1];
    f.vertices[1] = f.vertices[2];
    f.vertices[2] = tmp;
    mesh
}

/// 5. Degenerate Face: A standard cube plus one extra degenerate triangle (using existing vertices).
fn generate_degenerate_face() -> IndexedMesh {
    let mut mesh = generate_cube();
    
    // Add a degenerate triangle using existing vertices (v0, v0, v1)
    // This ensures that when removed, we don't leave unused vertices that mess up Chi.
    let region = RegionId::new(1);
    let v0 = FaceId::from_usize(0); // arbitrary existing
    // We need vertex IDs, not face IDs.
    // generate_cube creates 8 vertices.
    let v_a = cfd_mesh::core::index::VertexId::from_usize(0);
    let v_b = cfd_mesh::core::index::VertexId::from_usize(1);
    
    // Degenerate face: v_a, v_a, v_b (two vertices are same)
    mesh.faces.push(FaceData::new(v_a, v_a, v_b, region));
    mesh
}
