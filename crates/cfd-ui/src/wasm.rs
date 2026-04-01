//! WASM entry points for the CFD mesh visualizer.
//!
//! Exports mesh geometry (vertices + triangle indices) as JSON to JavaScript,
//! where Three.js handles all WebGL rendering. This decouples the rendering
//! backend from `wgpu` version constraints.

use wasm_bindgen::prelude::*;
use nalgebra::Point3;
use serde::Serialize;

use cfd_mesh::application::delaunay::dim3::sdf::FiniteCylinderSdf;
use cfd_mesh::application::delaunay::dim3::SdfMesher;

fn log(msg: &str) {
    web_sys::console::log_1(&wasm_bindgen::JsValue::from_str(msg));
}

/// Serializable mesh representation for JavaScript consumption.
#[derive(Serialize)]
struct JsMeshData {
    /// Flat array of vertex positions: [x0, y0, z0, x1, y1, z1, ...]
    vertices: Vec<f32>,
    /// Flat array of triangle indices: [i0, i1, i2, ...]
    indices: Vec<u32>,
    /// Flat array of vertex normals: [nx0, ny0, nz0, ...]
    normals: Vec<f32>,
    /// Number of vertices
    vertex_count: usize,
    /// Number of triangles
    triangle_count: usize,
}

#[wasm_bindgen(start)]
pub fn main_js() -> Result<(), JsValue> {
    console_error_panic_hook::set_once();
    log("main_js: panic hook set");
    Ok(())
}

/// Minimal test function to verify WASM loads and executes.
#[wasm_bindgen]
pub fn ping() -> String {
    log("ping: called");
    "pong".to_string()
}

/// Generate a Delaunay tetrahedralized cylinder mesh and return it as JSON.
///
/// The mesh is generated using the `BowyerWatson3D` kernel with a
/// `FiniteCylinderSdf` domain. The output contains the surface triangle
/// mesh extracted from the tetrahedral boundary.
#[wasm_bindgen]
pub fn generate_delaunay_cylinder(resolution: f64) -> String {
    log("generate_delaunay_cylinder: entry");
    let sdf = FiniteCylinderSdf::new(
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        0.5,
    );
    // Clamp to WASM-safe range: minimum 0.25 to keep seed count tractable.
    let res = resolution.clamp(0.25, 2.0);
    log(&format!("generate_delaunay_cylinder: creating mesher with res={res}"));
    let mut mesher = SdfMesher::<f64>::new(res);
    mesher.snap_iterations = 0;
    log("generate_delaunay_cylinder: calling build_volume");
    let mesh = mesher.build_volume(&sdf);
    log(&format!("generate_delaunay_cylinder: mesh has {} verts, {} faces",
        mesh.vertex_count(), mesh.face_count()));

    indexed_mesh_to_json(&mesh)
}

/// Generate a Delaunay tetrahedralized sphere mesh and return it as JSON.
///
/// Uses an offset center to avoid axis-aligned co-spherical degeneracies
/// in the `BowyerWatson3D` kernel.
#[wasm_bindgen]
pub fn generate_delaunay_sphere(resolution: f64) -> String {
    use cfd_mesh::application::delaunay::dim3::sdf::SphereSdf;
    log("generate_delaunay_sphere: entry");
    // Offset center by a small irrational displacement to break grid alignment symmetries
    let sdf = SphereSdf {
        center: Point3::new(0.037, 0.019, 0.013),
        radius: 1.0,
    };
    // Clamp resolution to WASM-safe range.
    // Minimum 0.25: at this resolution a unit sphere BCC lattice produces ~500 seed
    // points which fits within WASM linear memory constraints.
    // Maximum 2.0: coarser than this produces degenerate geometry.
    let res = resolution.clamp(0.25, 2.0);
    log(&format!("generate_delaunay_sphere: creating mesher with res={res}"));
    let mut mesher = SdfMesher::<f64>::new(res);
    mesher.snap_iterations = 0;
    log("generate_delaunay_sphere: calling build_volume");
    let mesh = mesher.build_volume(&sdf);
    log(&format!("generate_delaunay_sphere: mesh has {} verts, {} faces",
        mesh.vertex_count(), mesh.face_count()));

    indexed_mesh_to_json(&mesh)
}

/// Convert an `IndexedMesh` to a JSON string with vertices, indices, normals.
fn indexed_mesh_to_json(mesh: &cfd_mesh::domain::mesh::IndexedMesh) -> String {
    let n_verts = mesh.vertex_count();
    let mut vertices = Vec::with_capacity(n_verts * 3);

    for (_, vdata) in mesh.vertices.iter() {
        vertices.push(vdata.position.x as f32);
        vertices.push(vdata.position.y as f32);
        vertices.push(vdata.position.z as f32);
    }

    let n_faces = mesh.face_count();
    let mut indices = Vec::with_capacity(n_faces * 3);
    for (_, face) in mesh.faces.iter_enumerated() {
        indices.push(face.vertices[0].as_usize() as u32);
        indices.push(face.vertices[1].as_usize() as u32);
        indices.push(face.vertices[2].as_usize() as u32);
    }

    // Recompute per-vertex normals from face geometry for fidelity
    let mut computed_normals = vec![0.0f32; n_verts * 3];
    for (_, face) in mesh.faces.iter_enumerated() {
        let i0 = face.vertices[0].as_usize();
        let i1 = face.vertices[1].as_usize();
        let i2 = face.vertices[2].as_usize();
        let v0 = mesh.vertices.position(face.vertices[0]);
        let v1 = mesh.vertices.position(face.vertices[1]);
        let v2 = mesh.vertices.position(face.vertices[2]);
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let n = e1.cross(&e2);
        for idx in [i0, i1, i2] {
            computed_normals[idx * 3] += n.x as f32;
            computed_normals[idx * 3 + 1] += n.y as f32;
            computed_normals[idx * 3 + 2] += n.z as f32;
        }
    }
    for chunk in computed_normals.chunks_exact_mut(3) {
        let len = (chunk[0] * chunk[0] + chunk[1] * chunk[1] + chunk[2] * chunk[2]).sqrt();
        if len > 1e-12 {
            chunk[0] /= len;
            chunk[1] /= len;
            chunk[2] /= len;
        }
    }

    let data = JsMeshData {
        vertices,
        indices,
        normals: computed_normals,
        vertex_count: n_verts,
        triangle_count: n_faces,
    };
    serde_json::to_string(&data).expect("mesh serialization failed")
}
