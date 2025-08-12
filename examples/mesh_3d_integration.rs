//! 3D Mesh Integration Example
//!
//! This example demonstrates the CSG mesh integration capabilities of the 3D CFD module,
//! including mesh creation, CSG operations, and mesh quality assessment.

use cfd_suite::prelude::*;
use cfd_mesh::{Mesh, Vertex, Face, Cell, MeshTopology, csg::CsgMeshAdapter};
use nalgebra::Point3;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("3D Mesh Integration Example");
    println!("===========================");
    
    // Create CSG adapter for mesh operations
    let csg_adapter = CsgMeshAdapter::<f64>::new();
    
    println!("CSG Mesh Adapter created");
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
    
    // Create a second mesh for CSG operations
    println!("Creating a second tetrahedron for CSG operations...");
    let mesh2 = create_offset_tetrahedron(0.5)?;
    
    // Perform CSG union
    println!("Performing CSG union operation...");
    match csg_adapter.union(&mesh, &mesh2) {
        Ok(union_mesh) => {
            println!("Union mesh created:");
            println!("  Vertices: {}", union_mesh.vertices.len());
            println!("  Faces: {}", union_mesh.faces.len());
        }
        Err(e) => {
            println!("Union operation failed: {}", e);
        }
    }
    println!();
    
    // Perform CSG intersection
    println!("Performing CSG intersection operation...");
    match csg_adapter.intersection(&mesh, &mesh2) {
        Ok(intersection_mesh) => {
            println!("Intersection mesh created:");
            println!("  Vertices: {}", intersection_mesh.vertices.len());
            println!("  Faces: {}", intersection_mesh.faces.len());
        }
        Err(e) => {
            println!("Intersection operation failed: {}", e);
        }
    }
    println!();
    
    // Perform CSG difference
    println!("Performing CSG difference operation...");
    match csg_adapter.difference(&mesh, &mesh2) {
        Ok(difference_mesh) => {
            println!("Difference mesh created:");
            println!("  Vertices: {}", difference_mesh.vertices.len());
            println!("  Faces: {}", difference_mesh.faces.len());
        }
        Err(e) => {
            println!("Difference operation failed: {}", e);
        }
    }
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

/// Create an offset tetrahedron mesh
fn create_offset_tetrahedron(offset: f64) -> std::result::Result<Mesh<f64>, Box<dyn std::error::Error>> {
    let mut mesh = Mesh::new();
    
    // Define vertices of an offset tetrahedron
    mesh.vertices = vec![
        Vertex { id: 0, position: Point3::new(offset, offset, offset) },
        Vertex { id: 1, position: Point3::new(1.0 + offset, offset, offset) },
        Vertex { id: 2, position: Point3::new(0.5 + offset, 0.866 + offset, offset) },
        Vertex { id: 3, position: Point3::new(0.5 + offset, 0.289 + offset, 0.816 + offset) },
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
