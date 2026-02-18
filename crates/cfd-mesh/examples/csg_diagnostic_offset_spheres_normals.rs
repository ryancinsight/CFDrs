//! CSG Diagnostic: Offset overlapping spheres with normal-orientation analysis.
//!
//! Sphere A: radius 1.0 mm, centered at (0,0,0)
//! Sphere B: radius 1.0 mm, centered at (0.5,0.5,0.5)
//! Center distance: √0.75 ≈ 0.866 mm
//!
//! Expected volumes (analytical, two equal spheres, r=1, d=√0.75):
//!   h = r - d/2 ≈ 0.5670
//!   V_cap = (π h² / 3)(3r - h) ≈ 0.8191
//! - Union:        ≈ 6.739 mm³
//! - Intersection: ≈ 1.638 mm³
//! - Difference:   ≈ 2.551 mm³
//!
//! Run with:
//! cargo run --manifest-path CFDrs/Cargo.toml -p cfd-mesh --example csg_diagnostic_offset_spheres_normals --features "csg,stl-io"

use std::f64::consts::{PI, TAU};
use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::index::{RegionId, VertexId};
use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{csg_boolean, BooleanOp};
use cfd_mesh::geometry::normal::triangle_normal;
use cfd_mesh::io::stl;
use cfd_mesh::storage::face_store::FaceData;
use cfd_mesh::storage::vertex_pool::VertexPool;
use cfd_mesh::watertight::check::check_watertight;
use cfd_mesh::IndexedMesh;

// ─── Analytical expected values ───────────────────────────────────────────────
// Two unit spheres, d = √0.75, h = 1 - d/2
const EXPECTED_UNION_VOLUME: Real = 6.739;
const EXPECTED_INTERSECTION_VOLUME: Real = 1.638;
const EXPECTED_DIFFERENCE_VOLUME: Real = 2.551;

// ─── Pass / fail thresholds ───────────────────────────────────────────────────
/// Tessellated spheres carry approximation error; allow generous volume tolerance.
const MAX_VOLUME_ERROR: Real = 0.25;
const MAX_INWARD_FACE_FRACTION: Real = 0.05;
const MIN_FACE_VERTEX_ALIGNMENT: Real = 0.50;
const MIN_VERTEX_NORMAL_LEN_MEAN: Real = 0.70;

/// Longitude divisions for sphere tessellation.
const SPHERE_SEGMENTS: usize = 32;
/// Latitude divisions for sphere tessellation.
const SPHERE_STACKS: usize = 16;

#[derive(Debug)]
struct NormalAnalysis {
    outward_faces: usize,
    inward_faces: usize,
    degenerate_faces: usize,
    face_vertex_alignment_mean: Real,
    face_vertex_alignment_min: Real,
    vertex_normal_len_mean: Real,
    vertex_normal_len_min: Real,
    vertex_normal_len_max: Real,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Diagnostic: Offset Spheres + Normal Analysis");
    println!("=================================================================");
    println!("  Sphere A: r=1.0 mm @ (0,0,0)");
    println!("  Sphere B: r=1.0 mm @ (0.5,0.5,0.5)");
    println!("  Center distance: {:.4} mm", (0.75_f64).sqrt());
    println!("  Tessellation: {} segments × {} stacks", SPHERE_SEGMENTS, SPHERE_STACKS);

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir
        .join("outputs")
        .join("csg_diagnostic_offset_spheres_normals");
    fs::create_dir_all(&out_dir)?;
    let out_path = out_dir.to_str().expect("non-UTF8 path");

    run_case("union", BooleanOp::Union, EXPECTED_UNION_VOLUME, out_path)?;
    run_case(
        "intersection",
        BooleanOp::Intersection,
        EXPECTED_INTERSECTION_VOLUME,
        out_path,
    )?;
    run_case(
        "difference",
        BooleanOp::Difference,
        EXPECTED_DIFFERENCE_VOLUME,
        out_path,
    )?;

    println!("=================================================================");
    println!("  Complete");
    println!("=================================================================");

    Ok(())
}

/// Boundary-only T-junction seam repair.
///
/// After BSP boolean operations on curved surfaces, the intersection seam
/// alternates between cut vertices from mesh A (on sphere-B face planes) and
/// cut vertices from mesh B (on sphere-A face planes).  Adjacent pairs are
/// separated by ~0.089 mm (half a tessellation step at the intersection
/// latitude) — far larger than the 1e-4 mm BSP weld tolerance.
///
/// This function:
///  1. Identifies boundary vertices (vertices on edges used by exactly 1 face).
///  2. Uses union-find to merge boundary vertex pairs within `tolerance`.
///  3. Interior vertices (including thin pole-adjacent edges) are untouched.
///
/// Safety: min edge on sphere surface at seam latitude ≈ 0.177 mm >> 0.1 mm
/// tolerance, so no legitimate edges are collapsed.
fn repair_tj_seam(faces: &[FaceData], pool: &VertexPool, tolerance: Real) -> Vec<FaceData> {
    use std::collections::{HashMap, HashSet};

    // Build edge → face-use count.
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for face in faces {
        let verts = face.vertices;
        for &(u, v) in &[
            (verts[0].as_usize(), verts[1].as_usize()),
            (verts[1].as_usize(), verts[2].as_usize()),
            (verts[2].as_usize(), verts[0].as_usize()),
        ] {
            let key = if u < v { (u, v) } else { (v, u) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }

    // Collect boundary vertex indices.
    let mut bvid_set: HashSet<usize> = HashSet::new();
    for (&(u, v), &cnt) in &edge_count {
        if cnt == 1 {
            bvid_set.insert(u);
            bvid_set.insert(v);
        }
    }
    let mut bvids: Vec<usize> = bvid_set.into_iter().collect();
    bvids.sort_unstable();
    let n = bvids.len();
    if n == 0 {
        return faces.to_vec();
    }

    // Union-find: merge boundary vertices within `tolerance`.
    let mut parent: Vec<usize> = (0..n).collect();
    let find_root = |parent: &Vec<usize>, mut x: usize| -> usize {
        while parent[x] != x {
            x = parent[x];
        }
        x
    };
    let tol_sq = tolerance * tolerance;
    for i in 0..n {
        let pi = pool.position(VertexId::new(bvids[i] as u32));
        for j in (i + 1)..n {
            let pj = pool.position(VertexId::new(bvids[j] as u32));
            let dx = pi.x - pj.x;
            let dy = pi.y - pj.y;
            let dz = pi.z - pj.z;
            if dx * dx + dy * dy + dz * dz <= tol_sq {
                let ri = find_root(&parent, i);
                let rj = find_root(&parent, j);
                if ri != rj {
                    parent[rj] = ri;
                }
            }
        }
    }

    // Build merge map: old index → representative index (only for changed vertices).
    let mut merge: HashMap<usize, usize> = HashMap::new();
    for i in 0..n {
        let root = find_root(&parent, i);
        if root != i {
            merge.insert(bvids[i], bvids[root]);
        }
    }

    // Apply merge map to face vertex IDs; drop collapsed triangles.
    faces
        .iter()
        .filter_map(|face| {
            let verts: [VertexId; 3] = face.vertices.map(|vid| {
                let idx = vid.as_usize();
                VertexId::new(*merge.get(&idx).unwrap_or(&idx) as u32)
            });
            if verts[0] == verts[1] || verts[1] == verts[2] || verts[0] == verts[2] {
                None
            } else {
                Some(FaceData { vertices: verts, region: face.region })
            }
        })
        .collect()
}

/// Compute signed volume directly from a face list + pool (no IndexedMesh needed).
fn signed_volume_from_faces(faces: &[FaceData], pool: &VertexPool) -> Real {
    faces
        .iter()
        .map(|f| {
            let a = pool.position(f.vertices[0]).coords;
            let b = pool.position(f.vertices[1]).coords;
            let c = pool.position(f.vertices[2]).coords;
            a.dot(&b.cross(&c)) / 6.0
        })
        .sum()
}

fn run_case(
    name: &str,
    op: BooleanOp,
    expected_volume: Real,
    out_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n--- {} ---", name.to_uppercase());

    let t0 = Instant::now();
    let mut pool = VertexPool::default_millifluidic();

    let sphere_a = generate_sphere(
        Point3r::origin(),
        1.0,
        SPHERE_SEGMENTS,
        SPHERE_STACKS,
        &mut pool,
        RegionId::new(1),
    );
    let sphere_b = generate_sphere(
        Point3r::new(0.5, 0.5, 0.5),
        1.0,
        SPHERE_SEGMENTS,
        SPHERE_STACKS,
        &mut pool,
        RegionId::new(2),
    );

    // ── Input mesh validation ──────────────────────────────────────────────
    let analytical_sphere_vol: Real = 4.0 / 3.0 * PI;
    let vol_a = signed_volume_from_faces(&sphere_a, &pool);
    let vol_b = signed_volume_from_faces(&sphere_b, &pool);
    println!(
        "  Input A volume: {:.4} mm³ (analytic {:.4}, err {:.4})",
        vol_a,
        analytical_sphere_vol,
        (vol_a - analytical_sphere_vol).abs()
    );
    println!(
        "  Input B volume: {:.4} mm³ (analytic {:.4}, err {:.4})",
        vol_b,
        analytical_sphere_vol,
        (vol_b - analytical_sphere_vol).abs()
    );
    if vol_a <= 0.0 {
        println!("  WARNING: Sphere A has non-positive volume — winding is inverted!");
    }
    if vol_b <= 0.0 {
        println!("  WARNING: Sphere B has non-positive volume — winding is inverted!");
    }

    let result_faces = csg_boolean(op, &sphere_a, &sphere_b, &mut pool)?;

    // ── Phase 1: raw BSP output ────────────────────────────────────────────
    let mut mesh = IndexedMesh::new();
    mesh.vertices = pool;
    for f in &result_faces {
        mesh.faces.push(*f);
    }
    mesh.rebuild_edges();
    {
        let wt_raw = {
            let edges = mesh.edges_ref().unwrap();
            check_watertight(&mesh.vertices, &mesh.faces, edges)
        };
        println!(
            "  [BSP raw] faces={}, boundary_edges={}, non_manifold_edges={}",
            mesh.face_count(), wt_raw.boundary_edge_count, wt_raw.non_manifold_edge_count
        );
    }

    // ── Phase 2: boundary-only T-junction seam repair ──────────────────────
    // Merge boundary vertex pairs within 0.1 mm.
    // Gap between matched seam pairs ≈ 0.089 mm; min sphere edge at seam ≈ 0.177 mm.
    const TJ_SEAM_TOLERANCE: Real = 0.1;
    let rw_faces = repair_tj_seam(mesh.faces.as_slice(), &mesh.vertices, TJ_SEAM_TOLERANCE);
    // Replace the face list in-place (vertex pool stays the same).
    mesh.faces = cfd_mesh::storage::face_store::FaceStore::new();
    for f in &rw_faces {
        mesh.faces.push(*f);
    }
    mesh.rebuild_edges();
    mesh.recompute_normals();

    let area = mesh.surface_area();
    let volume = mesh.signed_volume();
    let volume_err = (volume - expected_volume).abs();
    let normals = analyze_normals(&mesh);

    // ── Detailed watertight report ─────────────────────────────────────────
    let wt = {
        let edges = mesh.edges_ref().expect("edges just rebuilt");
        check_watertight(&mesh.vertices, &mesh.faces, edges)
    };

    println!("  Faces:        {}", mesh.face_count());
    println!("  Vertices:     {}", mesh.vertex_count());
    println!("  Surface Area: {:.4} mm²", area);
    println!(
        "  Volume:       {:.4} mm³ (expected {:.3}, err {:.4})",
        volume, expected_volume, volume_err
    );
    println!(
        "  Watertight:   {} (boundary_edges={}, non_manifold_edges={})",
        wt.is_watertight, wt.boundary_edge_count, wt.non_manifold_edge_count
    );

    println!("  Normal analysis:");
    println!(
        "    Face orientation estimate: outward={}, inward={}, degenerate={}",
        normals.outward_faces, normals.inward_faces, normals.degenerate_faces
    );
    println!(
        "    Face↔vertex alignment: mean={:.4}, min={:.4}",
        normals.face_vertex_alignment_mean, normals.face_vertex_alignment_min
    );
    println!(
        "    Vertex |n| stats: mean={:.4}, min={:.4}, max={:.4}",
        normals.vertex_normal_len_mean,
        normals.vertex_normal_len_min,
        normals.vertex_normal_len_max
    );

    let total_faces = mesh.face_count();
    let inward_fraction = if total_faces > 0 {
        normals.inward_faces as Real / total_faces as Real
    } else {
        1.0
    };

    let volume_ok = volume_err <= MAX_VOLUME_ERROR;
    let watertight_ok = wt.is_watertight;
    let inward_ok = inward_fraction <= MAX_INWARD_FACE_FRACTION;
    let alignment_ok = normals.face_vertex_alignment_min >= MIN_FACE_VERTEX_ALIGNMENT;
    let normal_len_ok = normals.vertex_normal_len_mean >= MIN_VERTEX_NORMAL_LEN_MEAN;
    let status = if volume_ok && watertight_ok && inward_ok && alignment_ok && normal_len_ok {
        "PASS"
    } else {
        "FAIL"
    };
    let mut why_fail = Vec::new();
    if !volume_ok       { why_fail.push(format!("vol_err={:.4}", volume_err)); }
    if !watertight_ok   { why_fail.push(format!("boundary_edges={}", wt.boundary_edge_count)); }
    if !inward_ok       { why_fail.push(format!("inward={:.3}", inward_fraction)); }
    if !alignment_ok    { why_fail.push(format!("align_min={:.4}", normals.face_vertex_alignment_min)); }
    if !normal_len_ok   { why_fail.push(format!("|n|_mean={:.4}", normals.vertex_normal_len_mean)); }
    println!(
        "  Status:       {} ({})",
        status,
        if why_fail.is_empty() { "all checks passed".to_string() } else { why_fail.join(", ") }
    );

    let stl_path = format!("{}/{}.stl", out_path, name);
    {
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
    }
    println!("  STL: {}", stl_path);
    println!("  Elapsed: {} ms", t0.elapsed().as_millis());

    Ok(())
}

fn analyze_normals(mesh: &IndexedMesh) -> NormalAnalysis {
    let mut centroid_sum = Vector3r::zeros();
    let mut vertex_count = 0usize;
    for (_, v) in mesh.vertices.iter() {
        centroid_sum += v.position.coords;
        vertex_count += 1;
    }
    let mesh_center = if vertex_count > 0 {
        Point3r::from(centroid_sum / vertex_count as Real)
    } else {
        Point3r::origin()
    };

    let mut outward_faces = 0usize;
    let mut inward_faces = 0usize;
    let mut degenerate_faces = 0usize;

    let mut align_sum = 0.0;
    let mut align_count = 0usize;
    let mut align_min: Real = 1.0;

    for face in mesh.faces.iter() {
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);

        let Some(face_n) = triangle_normal(a, b, c) else {
            degenerate_faces += 1;
            continue;
        };

        let face_center = Point3r::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        let to_face = face_center - mesh_center;
        if to_face.norm() > 1e-12 {
            if face_n.dot(&to_face.normalize()) >= 0.0 {
                outward_faces += 1;
            } else {
                inward_faces += 1;
            }
        }

        let avg_vertex_n = (*mesh.vertices.normal(face.vertices[0])
            + *mesh.vertices.normal(face.vertices[1])
            + *mesh.vertices.normal(face.vertices[2]))
            / 3.0;
        let avg_len = avg_vertex_n.norm();
        if avg_len > 1e-12 {
            let alignment = face_n.dot(&(avg_vertex_n / avg_len));
            align_sum += alignment;
            align_count += 1;
            align_min = align_min.min(alignment);
        }
    }

    let mut vlen_sum: Real = 0.0;
    let mut vlen_count = 0usize;
    let mut vlen_min = Real::MAX;
    let mut vlen_max: Real = 0.0;
    for (_, v) in mesh.vertices.iter() {
        let l = v.normal.norm();
        vlen_sum += l;
        vlen_count += 1;
        vlen_min = vlen_min.min(l);
        vlen_max = vlen_max.max(l);
    }

    NormalAnalysis {
        outward_faces,
        inward_faces,
        degenerate_faces,
        face_vertex_alignment_mean: if align_count > 0 {
            align_sum / align_count as Real
        } else {
            0.0
        },
        face_vertex_alignment_min: if align_count > 0 { align_min } else { 0.0 },
        vertex_normal_len_mean: if vlen_count > 0 {
            vlen_sum / vlen_count as Real
        } else {
            0.0
        },
        vertex_normal_len_min: if vlen_count > 0 { vlen_min } else { 0.0 },
        vertex_normal_len_max: if vlen_count > 0 { vlen_max } else { 0.0 },
    }
}

/// Generate a UV-sphere mesh with outward-facing normals.
///
/// Uses latitude/longitude (stacked ring) tessellation. Poles are capped
/// with triangles to avoid degenerate quads.
fn generate_sphere(
    center: Point3r,
    radius: Real,
    segments: usize,
    stacks: usize,
    pool: &mut VertexPool,
    region: RegionId,
) -> Vec<FaceData> {
    let mut faces = Vec::with_capacity(segments * stacks * 2);

    // Spherical → Cartesian. phi = colatitude (0 at north pole, π at south pole).
    let vertex_at = |theta: Real, phi: Real| -> (Point3r, Vector3r) {
        let normal = Vector3r::new(
            phi.sin() * theta.cos(),
            phi.cos(),
            phi.sin() * theta.sin(),
        );
        (center + normal * radius, normal)
    };

    for i in 0..segments {
        for j in 0..stacks {
            let theta0 = (i as Real / segments as Real) * TAU;
            let theta1 = ((i + 1) as Real / segments as Real) * TAU;
            let phi0 = (j as Real / stacks as Real) * PI;
            let phi1 = ((j + 1) as Real / stacks as Real) * PI;

            let (pos00, n00) = vertex_at(theta0, phi0);
            let (pos10, n10) = vertex_at(theta1, phi0);
            let (pos01, n01) = vertex_at(theta0, phi1);
            let (pos11, n11) = vertex_at(theta1, phi1);

            let v00 = pool.insert_or_weld(pos00, n00);
            let v10 = pool.insert_or_weld(pos10, n10);
            let v01 = pool.insert_or_weld(pos01, n01);
            let v11 = pool.insert_or_weld(pos11, n11);

            if j == 0 {
                // North pole cap — one triangle per segment.
                // Winding: (pole, v11, v01) gives outward normal (+Y at north pole).
                // (pole, v01, v11) would be inverted — a common mistake.
                faces.push(FaceData::new(v10, v11, v01, region));
            } else if j == stacks - 1 {
                // South pole cap — one triangle per segment.
                faces.push(FaceData::new(v00, v10, v01, region));
            } else {
                // Middle band — split quad into two triangles.
                faces.push(FaceData::new(v00, v10, v11, region));
                faces.push(FaceData::new(v00, v11, v01, region));
            }
        }
    }

    faces
}
