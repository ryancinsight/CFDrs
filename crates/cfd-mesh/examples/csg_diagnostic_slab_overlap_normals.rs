//! CSG Diagnostic: Slab overlap cubes with normal-orientation analysis.
//!
//! Cube A: 2×2×2 at (0,0,0)
//! Cube B: 2×2×2 at (1.5,0,0)
//! Overlap slab: 0.5×2×2 = 2.0 mm³
//!
//! Expected volumes:
//! - Union:        14.0 mm³
//! - Intersection:  2.0 mm³
//! - Difference:    6.0 mm³
//!
//! Run with:
//! cargo run --manifest-path CFDrs/Cargo.toml -p cfd-mesh --example csg_diagnostic_slab_overlap_normals --features "csg,stl-io"

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::index::RegionId;
use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{csg_boolean, BooleanOp};
use cfd_mesh::geometry::normal::triangle_normal;
use cfd_mesh::io::stl;
use cfd_mesh::storage::face_store::FaceData;
use cfd_mesh::storage::vertex_pool::VertexPool;
use cfd_mesh::IndexedMesh;

const EXPECTED_UNION_VOLUME: Real = 14.0;
const EXPECTED_INTERSECTION_VOLUME: Real = 2.0;
const EXPECTED_DIFFERENCE_VOLUME: Real = 6.0;

const MAX_VOLUME_ERROR: Real = 0.1;
const MAX_INWARD_FACE_FRACTION: Real = 0.05;
const MIN_FACE_VERTEX_ALIGNMENT: Real = 0.50;
const MIN_VERTEX_NORMAL_LEN_MEAN: Real = 0.70;

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
    println!("  CSG Diagnostic: Slab Overlap Cubes + Normal Analysis");
    println!("=================================================================");
    println!("  Cube A: 2×2×2 mm @ (0,0,0)");
    println!("  Cube B: 2×2×2 mm @ (1.5,0,0)");
    println!("  Overlap volume: 2.0 mm³");

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir
        .join("outputs")
        .join("csg_diagnostic_slab_overlap_normals");
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

fn run_case(
    name: &str,
    op: BooleanOp,
    expected_volume: Real,
    out_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n--- {} ---", name.to_uppercase());

    let t0 = Instant::now();
    let mut pool = VertexPool::default_millifluidic();

    let cube_a = generate_cube(2.0, Point3r::origin(), &mut pool, RegionId::new(1));
    let cube_b = generate_cube(
        2.0,
        Point3r::new(1.5, 0.0, 0.0),
        &mut pool,
        RegionId::new(2),
    );

    let result_faces = csg_boolean(op, &cube_a, &cube_b, &mut pool)?;

    let mut mesh = IndexedMesh::new();
    mesh.vertices = pool;
    for f in &result_faces {
        mesh.faces.push(*f);
    }
    mesh.rebuild_edges();
    mesh.recompute_normals();

    let area = mesh.surface_area();
    let volume = mesh.signed_volume();
    let volume_err = (volume - expected_volume).abs();
    let normals = analyze_normals(&mesh);

    println!("  Faces:        {}", mesh.face_count());
    println!("  Vertices:     {}", mesh.vertex_count());
    println!("  Surface Area: {:.4} mm²", area);
    println!(
        "  Volume:       {:.4} mm³ (expected {:.3}, err {:.4})",
        volume, expected_volume, volume_err
    );
    println!("  Watertight:   {}", mesh.is_watertight());

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
    let inward_ok = inward_fraction <= MAX_INWARD_FACE_FRACTION;
    let alignment_ok = normals.face_vertex_alignment_min >= MIN_FACE_VERTEX_ALIGNMENT;
    let normal_len_ok = normals.vertex_normal_len_mean >= MIN_VERTEX_NORMAL_LEN_MEAN;
    let status = if volume_ok && inward_ok && alignment_ok && normal_len_ok {
        "PASS"
    } else {
        "FAIL"
    };
    println!(
        "  Status:       {} (vol_err={:.4}, inward={:.3}, align_min={:.4}, |n|_mean={:.4})",
        status,
        volume_err,
        inward_fraction,
        normals.face_vertex_alignment_min,
        normals.vertex_normal_len_mean
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

fn generate_cube(
    size: Real,
    origin: Point3r,
    pool: &mut VertexPool,
    region: RegionId,
) -> Vec<FaceData> {
    let s = size;
    let o = origin;
    let mut faces = Vec::with_capacity(12);

    let p000 = Point3r::new(o.x, o.y, o.z);
    let p100 = Point3r::new(o.x + s, o.y, o.z);
    let p110 = Point3r::new(o.x + s, o.y + s, o.z);
    let p010 = Point3r::new(o.x, o.y + s, o.z);
    let p001 = Point3r::new(o.x, o.y, o.z + s);
    let p101 = Point3r::new(o.x + s, o.y, o.z + s);
    let p111 = Point3r::new(o.x + s, o.y + s, o.z + s);
    let p011 = Point3r::new(o.x, o.y + s, o.z + s);

    let mut add_quad = |p0: Point3r, p1: Point3r, p2: Point3r, p3: Point3r, n: Vector3r| {
        let v0 = pool.insert_or_weld(p0, n);
        let v1 = pool.insert_or_weld(p1, n);
        let v2 = pool.insert_or_weld(p2, n);
        let v3 = pool.insert_or_weld(p3, n);
        faces.push(FaceData::new(v0, v1, v2, region));
        faces.push(FaceData::new(v0, v2, v3, region));
    };

    add_quad(p000, p010, p110, p100, -Vector3r::z());
    add_quad(p001, p101, p111, p011, Vector3r::z());
    add_quad(p000, p100, p101, p001, -Vector3r::y());
    add_quad(p010, p011, p111, p110, Vector3r::y());
    add_quad(p000, p001, p011, p010, -Vector3r::x());
    add_quad(p100, p110, p111, p101, Vector3r::x());

    faces
}
