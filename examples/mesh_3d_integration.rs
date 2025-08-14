//! 3D Mesh Integration Example
//!
//! This example demonstrates the basic mesh generation capabilities of the 3D CFD module,
//! including mesh creation and validation. CSG boolean operations are not currently implemented.

use cfd_mesh::{Mesh, Vertex, Face, Cell, MeshTopology, csg::CsgOperator};
use nalgebra::Point3;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("3D Mesh Integration Example");
    println!("===========================");
    
    // Create CSG operator for mesh generation (basic primitives only)
    let csg_operator = CsgOperator::<f64>::new();
    
    println!("CSG Operator created for basic mesh generation");
    println!();
    
    // Create a sphere mesh using the CSG operator
    println!("Creating sphere mesh...");
    match csg_operator.create_sphere(1.0, 8) {
        Ok(sphere_mesh) => {
            println!("Sphere mesh created successfully:");
            println!("  Vertices: {}", sphere_mesh.vertices.len());
            println!("  Faces: {}", sphere_mesh.faces.len());
            println!("  Cells: {}", sphere_mesh.cells.len());
            
            // Validate the mesh
            match csg_operator.validate_mesh(&sphere_mesh) {
                Ok(is_valid) => {
                    println!("  Mesh validation: {}", if is_valid { "PASSED" } else { "FAILED" });
                }
                Err(e) => {
                    println!("  Mesh validation error: {}", e);
                }
            }
        }
        Err(e) => {
            println!("Sphere mesh creation failed: {}", e);
        }
    }
    println!();
    
    // Create a box mesh using the CSG operator  
    println!("Creating box mesh...");
    match csg_operator.create_box(2.0, 1.0, 1.5) {
        Ok(box_mesh) => {
            println!("Box mesh created successfully:");
            println!("  Vertices: {}", box_mesh.vertices.len());
            println!("  Faces: {}", box_mesh.faces.len());
            println!("  Cells: {}", box_mesh.cells.len());
            
            // Validate the mesh
            match csg_operator.validate_mesh(&box_mesh) {
                Ok(is_valid) => {
                    println!("  Mesh validation: {}", if is_valid { "PASSED" } else { "FAILED" });
                }
                Err(e) => {
                    println!("  Mesh validation error: {}", e);
                }
            }
        }
        Err(e) => {
            println!("Box mesh creation failed: {}", e);
        }
    }
    println!();
    
    // Create a manual tetrahedral mesh
    println!("Creating manual tetrahedral mesh...");
    let tet_mesh = create_unit_tetrahedron()?;
    
    println!("Tetrahedron mesh created successfully:");
    println!("  Vertices: {}", tet_mesh.vertices.len());
    println!("  Faces: {}", tet_mesh.faces.len());
    println!("  Cells: {}", tet_mesh.cells.len());
    
    // Display vertex coordinates
    println!("Vertex coordinates:");
    for (i, vertex) in tet_mesh.vertices.iter().enumerate() {
        println!("  Vertex {}: ({:.2}, {:.2}, {:.2})", 
                 i, vertex.position.x, vertex.position.y, vertex.position.z);
    }
    
    // Validate the manual mesh
    match csg_operator.validate_mesh(&tet_mesh) {
        Ok(is_valid) => {
            println!("  Mesh validation: {}", if is_valid { "PASSED" } else { "FAILED" });
        }
        Err(e) => {
            println!("  Mesh validation error: {}", e);
        }
    }
    println!();
    
    println!("NOTE: CSG boolean operations (union, intersection, difference) are not currently implemented.");
    println!("This is a documented limitation. Only basic primitive generation and validation are available.");
    println!();
    
    println!("Example completed successfully!");
    
    Ok(())
}

/// Create a unit tetrahedron mesh
fn create_unit_tetrahedron() -> std::result::Result<Mesh<f64>, Box<dyn std::error::Error>> {
    let mut mesh = Mesh::new();
    
    // Define vertices of a unit tetrahedron
    mesh.vertices = vec![
        Vertex { id: 0, position: Point3::new(0.0, 0.0, 0.0) },
        Vertex { id: 1, position: Point3::new(1.0, 0.0, 0.0) },
        Vertex { id: 2, position: Point3::new(0.5, 0.866, 0.0) },
        Vertex { id: 3, position: Point3::new(0.5, 0.289, 0.816) },
    ];
    
    // Define faces (triangles)
    mesh.faces = vec![
        Face { id: 0, vertices: vec![0, 1, 2] }, // Base
        Face { id: 1, vertices: vec![0, 1, 3] }, // Side 1
        Face { id: 2, vertices: vec![1, 2, 3] }, // Side 2
        Face { id: 3, vertices: vec![2, 0, 3] }, // Side 3
    ];
    
    // Define the single tetrahedral cell
    mesh.cells = vec![
        Cell { id: 0, faces: vec![0, 1, 2, 3] },
    ];
    
    // Update topology
    mesh.topology = MeshTopology {
        num_vertices: mesh.vertices.len(),
        num_edges: 6,  // A tetrahedron has 6 edges
        num_faces: mesh.faces.len(),
        num_cells: mesh.cells.len(),
    };
    
    Ok(mesh)
}


