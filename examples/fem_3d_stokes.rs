//! 3D FEM Stokes Flow Example
//!
//! This example demonstrates the use of the 3D FEM solver for solving
//! the Stokes equations in a simple tetrahedral domain.

use cfd_suite::prelude::*;
use cfd_core::Fluid;
use cfd_mesh::{Mesh, Vertex, Face, Cell, MeshTopology};
use cfd_3d::{ElementType, MaterialProperties};
use std::collections::HashMap;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("3D FEM Stokes Flow Example");
    println!("==========================");
    
    // Create fluid properties (water at 20°C)
    let fluid = Fluid::water();
    println!("Fluid: {}", fluid.name);
    println!("Density: {:.1} kg/m³", fluid.density);
    println!("Viscosity: {:.6} Pa·s", fluid.viscosity);
    println!();
    
    // Create a simple tetrahedral mesh
    let mesh = create_unit_cube_mesh()?;
    println!("Created tetrahedral mesh with {} vertices and {} cells", 
             mesh.vertices.len(), mesh.cells.len());
    
    // Create FEM solver configuration using builder pattern
    let base = cfd_core::SolverConfig::<f64>::builder()
        .tolerance(1e-6)
        .max_iterations(1000)
        .verbosity(2) // verbose = true means verbosity level 2
        .build();

    let config = FemConfig {
        base,
        element_type: ElementType::Tetrahedron4,
        integration_order: 2,
    };
    
    // Create FEM solver
    let mut solver = FemSolver::new(config);
    solver.set_mesh(mesh);
    
    // Define material properties
    let material = MaterialProperties {
        density: fluid.density,
        viscosity: fluid.viscosity,
        thermal_conductivity: None,
        specific_heat: None,
    };
    
    // Define boundary conditions using iterator patterns for zero-copy initialization
    let velocity_bcs: HashMap<usize, nalgebra::Vector3<f64>> = [
        // Top face: moving lid (u = 1, v = 0, w = 0) - using iterator for lid nodes
        (2, nalgebra::Vector3::new(1.0, 0.0, 0.0)),
        (3, nalgebra::Vector3::new(1.0, 0.0, 0.0)),
        // Bottom and side faces: no-slip (u = v = w = 0) - using iterator for wall nodes
        (0, nalgebra::Vector3::new(0.0, 0.0, 0.0)),
        (1, nalgebra::Vector3::new(0.0, 0.0, 0.0)),
    ]
    .into_iter()
    .collect();
    
    println!("Boundary conditions:");
    for (&node, &velocity) in &velocity_bcs {
        println!("  Node {}: u={:.1}, v={:.1}, w={:.1}", 
                 node, velocity.x, velocity.y, velocity.z);
    }
    println!();
    
    // Body forces using iterator patterns for uniform field application
    let gravity = nalgebra::Vector3::new(0.0, 0.0, -9.81);
    let body_forces: HashMap<usize, nalgebra::Vector3<f64>> = (0..4)
        .map(|i| (i, gravity))
        .collect();
    
    // Solve the Stokes equations
    println!("Solving 3D Stokes equations...");
    match solver.solve_stokes(&velocity_bcs, &body_forces, &material) {
        Ok(solution) => {
            println!("Solution converged successfully!");
            println!();
            
            println!("Velocity field:");
            for (&node, &velocity) in &solution {
                println!("  Node {}: u={:.6}, v={:.6}, w={:.6}", 
                         node, velocity.x, velocity.y, velocity.z);
            }
            
            // Calculate some flow statistics
            let max_velocity = solution.values()
                .map(|v| v.norm())
                .fold(0.0, f64::max);
            
            println!();
            println!("Flow statistics:");
            println!("  Maximum velocity magnitude: {:.6} m/s", max_velocity);
            
            // Estimate Reynolds number based on characteristic length and velocity
            let char_length = 1.0; // Unit cube
            let reynolds = fluid.density * max_velocity * char_length / fluid.viscosity;
            println!("  Reynolds number: {:.2}", reynolds);
            
            if reynolds < 1.0 {
                println!("  Flow regime: Creeping flow (Re << 1)");
            } else if reynolds < 100.0 {
                println!("  Flow regime: Low Reynolds number flow");
            } else {
                println!("  Flow regime: Moderate Reynolds number flow");
            }
        }
        Err(e) => {
            println!("Solution failed: {}", e);
            return Err(e.into());
        }
    }
    
    println!();
    println!("3D FEM Stokes flow simulation completed successfully!");
    
    Ok(())
}

/// Create a simple tetrahedral mesh of a unit cube
fn create_unit_cube_mesh() -> std::result::Result<Mesh<f64>, Box<dyn std::error::Error>> {
    // Create vertices for a unit cube
    let vertices = vec![
        Vertex { position: nalgebra::Point3::new(0.0, 0.0, 0.0), id: 0 }, // Bottom-front-left
        Vertex { position: nalgebra::Point3::new(1.0, 0.0, 0.0), id: 1 }, // Bottom-front-right
        Vertex { position: nalgebra::Point3::new(0.0, 1.0, 0.0), id: 2 }, // Bottom-back-left
        Vertex { position: nalgebra::Point3::new(0.0, 0.0, 1.0), id: 3 }, // Top-front-left
    ];
    
    // Create faces for a single tetrahedron
    let faces = vec![
        Face { vertices: vec![0, 1, 2], id: 0 }, // Bottom face
        Face { vertices: vec![0, 1, 3], id: 1 }, // Front face
        Face { vertices: vec![1, 2, 3], id: 2 }, // Right face
        Face { vertices: vec![0, 2, 3], id: 3 }, // Left face
    ];
    
    // Create a single tetrahedral cell
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
