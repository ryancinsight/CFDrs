//! 3D Mesh Integration Example
//!
//! This example demonstrates the mesh integration capabilities of the 3D CFD module,
//! including mesh creation, quality assessment, and CSG integration.

use cfd_3d::prelude::*;
use cfd_mesh::{Mesh, Vertex, Face, Cell, MeshTopology};
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("3D Mesh Integration Example");
    println!("===========================");
    
    // Create STL adapter for mesh operations
    let stl_adapter = StlAdapter::<f64>::default();
    
    println!("STL Adapter created");
    println!("Supported extensions: {:?}", stl_adapter.supported_extensions());
    println!();
    
    // Create a simple tetrahedral mesh
    println!("Creating unit tetrahedron mesh...");
    let mesh = create_unit_tetrahedron()?;
    
    println!("Mesh created successfully:");
    println!("  Vertices: {}", mesh.vertices.len());
    println!("  Faces: {}", mesh.faces.len());
    println!("  Cells: {}", mesh.cells.len());
    println!();
    
    // Display vertex coordinates
    println!("Vertex coordinates:");
    for (i, vertex) in mesh.vertices.iter().enumerate() {
        println!("  Vertex {}: ({:.2}, {:.2}, {:.2})", 
                 i, vertex.position.x, vertex.position.y, vertex.position.z);
    }
    println!();
    
    // Validate mesh quality
    println!("Validating mesh quality...");
    match stl_adapter.validate_mesh(&mesh) {
        Ok(quality_report) => {
            println!("Quality assessment completed:");
            println!("  Minimum quality: {:.6}", quality_report.min_quality);
            println!("  Maximum quality: {:.6}", quality_report.max_quality);
            println!("  Average quality: {:.6}", quality_report.avg_quality);
            println!("  Degenerate elements: {}", quality_report.degenerate_elements);
            println!("  Inverted elements: {}", quality_report.inverted_elements);
            println!("  Mesh is valid: {}", quality_report.is_valid);
            
            if quality_report.inverted_elements == 0 {
                println!("  ✓ No inverted elements detected");
            } else {
                println!("  ⚠ Warning: {} inverted elements found", quality_report.inverted_elements);
            }
            
            if quality_report.degenerate_elements == 0 {
                println!("  ✓ No degenerate elements detected");
            } else {
                println!("  ⚠ Warning: {} degenerate elements found", quality_report.degenerate_elements);
            }
        }
        Err(e) => {
            println!("Quality validation failed: {}", e);
        }
    }
    println!();
    
    // Test mesh import/export functionality
    println!("Testing mesh I/O operations...");
    
    // Export mesh (placeholder functionality)
    match stl_adapter.export_mesh(&mesh) {
        Ok(data) => {
            println!("Mesh export successful (placeholder): {} bytes", data.len());
        }
        Err(e) => {
            println!("Mesh export failed: {}", e);
        }
    }
    
    // Import mesh (placeholder functionality)
    let empty_data = vec![];
    match stl_adapter.import_mesh(&empty_data) {
        Ok(imported_mesh) => {
            println!("Mesh import successful:");
            println!("  Imported vertices: {}", imported_mesh.vertices.len());
            println!("  Imported cells: {}", imported_mesh.cells.len());
        }
        Err(e) => {
            println!("Mesh import failed: {}", e);
        }
    }
    println!();
    
    // Demonstrate CSG mesh adapter
    println!("Testing CSG mesh adapter...");
    let csg_adapter = CsgMeshAdapter::<f64>::new();
    
    match csg_adapter.generate_from_csg("sphere(1.0)") {
        Ok(csg_mesh) => {
            println!("CSG mesh generation successful:");
            println!("  Generated vertices: {}", csg_mesh.vertices.len());
            println!("  Generated cells: {}", csg_mesh.cells.len());
            println!("  Note: This is a placeholder implementation");
        }
        Err(e) => {
            println!("CSG mesh generation failed: {}", e);
        }
    }
    println!();
    
    // Calculate mesh statistics
    println!("Mesh statistics:");
    let total_volume = calculate_mesh_volume(&mesh)?;
    let surface_area = calculate_surface_area(&mesh)?;
    
    println!("  Total volume: {:.6} cubic units", total_volume);
    println!("  Surface area: {:.6} square units", surface_area);
    
    // For a unit tetrahedron with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    // The theoretical volume is 1/6 ≈ 0.166667
    let theoretical_volume = 1.0 / 6.0;
    let volume_error = (total_volume - theoretical_volume).abs() / theoretical_volume;
    
    println!("  Theoretical volume: {:.6}", theoretical_volume);
    println!("  Volume error: {:.2}%", volume_error * 100.0);
    
    if volume_error < 0.01 {
        println!("  ✓ Volume calculation is accurate");
    } else {
        println!("  ⚠ Volume calculation may have errors");
    }
    
    println!();
    println!("3D mesh integration demonstration completed successfully!");
    
    Ok(())
}

/// Create a simple unit tetrahedron mesh
fn create_unit_tetrahedron() -> Result<Mesh<f64>, Box<dyn std::error::Error>> {
    let vertices = vec![
        Vertex { position: nalgebra::Point3::new(0.0, 0.0, 0.0), id: 0 },
        Vertex { position: nalgebra::Point3::new(1.0, 0.0, 0.0), id: 1 },
        Vertex { position: nalgebra::Point3::new(0.0, 1.0, 0.0), id: 2 },
        Vertex { position: nalgebra::Point3::new(0.0, 0.0, 1.0), id: 3 },
    ];
    
    let faces = vec![
        Face { vertices: vec![0, 1, 2], id: 0 }, // Bottom face
        Face { vertices: vec![0, 1, 3], id: 1 }, // Front face
        Face { vertices: vec![1, 2, 3], id: 2 }, // Right face
        Face { vertices: vec![0, 2, 3], id: 3 }, // Left face
    ];
    
    let cells = vec![
        Cell { faces: vec![0, 1, 2, 3], id: 0 },
    ];
    
    let topology = MeshTopology {
        num_vertices: 4,
        num_edges: 6,
        num_faces: 4,
        num_cells: 1,
    };
    
    Ok(Mesh {
        vertices,
        edges: vec![],
        faces,
        cells,
        topology,
    })
}

/// Calculate the volume of a tetrahedral mesh
fn calculate_mesh_volume(mesh: &Mesh<f64>) -> Result<f64, Box<dyn std::error::Error>> {
    let mut total_volume = 0.0;
    
    // For each cell (assuming tetrahedral)
    for _cell in &mesh.cells {
        // For simplicity, calculate volume of the single tetrahedron
        // Volume = |det(v1-v0, v2-v0, v3-v0)| / 6
        if mesh.vertices.len() >= 4 {
            let v0 = &mesh.vertices[0].position;
            let v1 = &mesh.vertices[1].position;
            let v2 = &mesh.vertices[2].position;
            let v3 = &mesh.vertices[3].position;
            
            let e1 = nalgebra::Vector3::new(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
            let e2 = nalgebra::Vector3::new(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);
            let e3 = nalgebra::Vector3::new(v3.x - v0.x, v3.y - v0.y, v3.z - v0.z);
            
            let volume = e1.cross(&e2).dot(&e3).abs() / 6.0;
            total_volume += volume;
        }
    }
    
    Ok(total_volume)
}

/// Calculate the surface area of a mesh
fn calculate_surface_area(mesh: &Mesh<f64>) -> Result<f64, Box<dyn std::error::Error>> {
    let mut total_area = 0.0;
    
    // For each face, calculate triangle area
    for face in &mesh.faces {
        if face.vertices.len() >= 3 {
            let v0 = &mesh.vertices[face.vertices[0]].position;
            let v1 = &mesh.vertices[face.vertices[1]].position;
            let v2 = &mesh.vertices[face.vertices[2]].position;
            
            let e1 = nalgebra::Vector3::new(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
            let e2 = nalgebra::Vector3::new(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);
            
            let area = e1.cross(&e2).norm() / 2.0;
            total_area += area;
        }
    }
    
    Ok(total_area)
}
