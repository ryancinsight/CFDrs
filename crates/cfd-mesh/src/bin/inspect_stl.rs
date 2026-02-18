use std::env;
use std::fs::File;
use std::io::BufReader;
use std::collections::{HashMap, HashSet};
use cfd_mesh::core::scalar::{Real, Point3r, Vector3r};

#[derive(Debug, Default)]
struct StlStats {
    triangle_count: usize,
    vertex_count: usize,
    degenerate_faces: usize,
    open_edges: usize,
    non_manifold_edges: usize,
    inverted_normals: usize, // Heuristic check
    volume: Real,
    area: Real,
    min_corner: Point3r,
    max_corner: Point3r,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: inspect_stl <file.stl>");
        std::process::exit(1);
    }

    let path = &args[1];
    println!("Inspecting: {}", path);

    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    // We use a simplified STL parser here or crate support
    // For now, let's use the cfd-mesh IO if available, or write a quick parser.
    // Since `cfd-mesh` has `read_ascii_stl`, we might want to expose a generic read.
    // But `cfd-mesh` parser requires VertexPool and FaceStore.
    
    // Let's use `cfd_mesh::io::stl` logic but adapted, or rely on `cfd_mesh`.
    let mut pool = cfd_mesh::storage::vertex_pool::VertexPool::default_millifluidic();
    let mut faces = cfd_mesh::storage::face_store::FaceStore::new();
    let region = cfd_mesh::core::index::RegionId::new(0);

    // Try reading (assuming ASCII for now based on example output, but example used binary write!)
    // Wait, example output said "STL files written to ...".
    // boolean_sphere_cube writes binary STL.
    // `io::stl` has `write_binary_stl` but no `read_binary_stl` in the file I viewed?
    // Let's check `src/io/stl.rs` again. It had `read_ascii_stl`. Does it have binary read?
    // I only viewed 150 lines. It might be missing.

    // If no binary read, I'll implement a quick binary reader here.
    
    // Check if file is binary (starts with 80 bytes header, then u32 count).
    // ASCII usually starts with "solid".
    
    let file = File::open(path)?;
    let meta = file.metadata()?;
    let len = meta.len();
    
    let mut reader = BufReader::new(file);
    // Simple heuristic: read first 5 bytes.
    let mut start = [0u8; 5];
    use std::io::Read;
    reader.read_exact(&mut start)?;
    
    // Reset
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    let is_ascii = &start == b"solid";

    if is_ascii {
        cfd_mesh::io::stl::read_ascii_stl(reader, &mut pool, &mut faces, region)?;
    } else {
        // Binary reader impl
        read_binary_stl(&mut reader, &mut pool, &mut faces, region)?;
    }

    let stats = analyze_mesh(&pool, &faces);
    print_report(&stats);

    if stats.open_edges > 0 || stats.non_manifold_edges > 0 {
        std::process::exit(1);
    }

    Ok(())
}

fn read_binary_stl<R: std::io::Read>(
    reader: &mut R,
    pool: &mut cfd_mesh::storage::vertex_pool::VertexPool,
    faces: &mut cfd_mesh::storage::face_store::FaceStore,
    region: cfd_mesh::core::index::RegionId,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut header = [0u8; 80];
    reader.read_exact(&mut header)?;
    
    let mut count_buf = [0u8; 4];
    reader.read_exact(&mut count_buf)?;
    let count = u32::from_le_bytes(count_buf);

    for _ in 0..count {
        let mut buf = [0u8; 50]; // 12 (normal) + 36 (verts) + 2 (attr)
        reader.read_exact(&mut buf)?;
        
        let nx = f32::from_le_bytes(buf[0..4].try_into()?);
        let ny = f32::from_le_bytes(buf[4..8].try_into()?);
        let nz = f32::from_le_bytes(buf[8..12].try_into()?);
        let normal = Vector3r::new(nx as Real, ny as Real, nz as Real);

        let mut vs = [cfd_mesh::core::index::VertexId::new(0); 3];
        for i in 0..3 {
            let offset = 12 + i * 12;
            let vx = f32::from_le_bytes(buf[offset..offset+4].try_into()?);
            let vy = f32::from_le_bytes(buf[offset+4..offset+8].try_into()?);
            let vz = f32::from_le_bytes(buf[offset+8..offset+12].try_into()?);
            let p = Point3r::new(vx as Real, vy as Real, vz as Real);
            
            // We use the normal from STL? Or recalc?
            // STL normal is often garbage. Let's store it but rely on geom.
            vs[i] = pool.insert_or_weld(p, normal);
        }
        
        faces.push(cfd_mesh::storage::face_store::FaceData {
            vertices: vs,
            region,
        });
    }
    Ok(())
}

fn analyze_mesh(
    pool: &cfd_mesh::storage::vertex_pool::VertexPool,
    faces: &cfd_mesh::storage::face_store::FaceStore,
) -> StlStats {
    let mut stats = StlStats::default();
    stats.triangle_count = faces.len();
    stats.vertex_count = pool.len();
    stats.min_corner = Point3r::new(Real::MAX, Real::MAX, Real::MAX);
    stats.max_corner = Point3r::new(Real::MIN, Real::MIN, Real::MIN);

    let mut edge_counts = HashMap::new();

    for (_, face) in faces.iter_enumerated() {
        let v = face.vertices;
        let p0 = pool.position(v[0]);
        let p1 = pool.position(v[1]);
        let p2 = pool.position(v[2]);

        // Bounds
        for p in [p0, p1, p2] {
            stats.min_corner.x = stats.min_corner.x.min(p.x);
            stats.min_corner.y = stats.min_corner.y.min(p.y);
            stats.min_corner.z = stats.min_corner.z.min(p.z);
            stats.max_corner.x = stats.max_corner.x.max(p.x);
            stats.max_corner.y = stats.max_corner.y.max(p.y);
            stats.max_corner.z = stats.max_corner.z.max(p.z);
        }

        // Area & Volume
        let edge1 = p1 - p0;
        let edge2 = p2 - p0;
        let cross = edge1.cross(&edge2);
        let area = cross.norm() * 0.5;
        stats.area += area;
        
        if area < 1e-9 {
            stats.degenerate_faces += 1;
        }

        // Signed volume (tetrahedron from origin)
        stats.volume += p0.coords.dot(&p1.coords.cross(&p2.coords)) / 6.0;

        // Edges (sorted vertex pairs)
        for i in 0..3 {
            let v_start = v[i].min(v[(i+1)%3]);
            let v_end = v[i].max(v[(i+1)%3]);
            *edge_counts.entry((v_start, v_end)).or_insert(0) += 1;
        }
    }

    for &count in edge_counts.values() {
        if count == 1 {
            stats.open_edges += 1;
        } else if count > 2 {
            stats.non_manifold_edges += 1;
        }
    }

    stats
}

fn print_report(stats: &StlStats) {
    println!("--- Mesh Report ---");
    println!("Triangles: {}", stats.triangle_count);
    println!("Vertices:  {}", stats.vertex_count);
    println!("Bounds:    [{:.2}, {:.2}, {:.2}] to [{:.2}, {:.2}, {:.2}]", 
        stats.min_corner.x, stats.min_corner.y, stats.min_corner.z,
        stats.max_corner.x, stats.max_corner.y, stats.max_corner.z);
    println!("Volume:    {:.4}", stats.volume);
    println!("Area:      {:.4}", stats.area);
    println!("Legacy Qs:");
    println!("  Degenerate faces:   {}", stats.degenerate_faces);
    println!("  Open edges (holes): {}", stats.open_edges);
    println!("  Non-manifold edges: {}", stats.non_manifold_edges);
    
    if stats.volume < 0.0 {
        println!("WARNING: Negative volume! Mesh likely inverted (inside-out).");
    }
}
