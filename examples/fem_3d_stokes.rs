//! 3D FEM Stokes Flow Example
//!
//! This example demonstrates the use of the 3D FEM solver for solving
//! the Stokes equations in a simple tetrahedral domain.

use cfd_3d::fem::{FemSolver, FemConfig, ElementType, StokesFlowProblem};
use cfd_core::fluid::Fluid;
use cfd_mesh::{Mesh, Vertex, Cell};
use nalgebra::{Point3, Vector3};
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("3D FEM Stokes Flow Example");
    println!("==========================");
    
    // Create fluid properties (water at 20°C)
    let fluid = Fluid::<f64>::water()?;
    println!("Fluid: Water at 20°C");
    println!("Density: {:.1} kg/m³", fluid.density());
    println!("Viscosity: {:.6} Pa·s", fluid.dynamic_viscosity());
    println!();
    
    // Create a simple tetrahedral mesh
    let mesh = create_unit_cube_mesh()?;
    println!("Created tetrahedral mesh with {} vertices and {} cells", 
             mesh.vertices.len(), mesh.cells.len());
    
    // Create FEM solver configuration
    let base = cfd_core::solver::SolverConfig::<f64>::builder()
        .tolerance(1e-6)
        .max_iterations(1000)
        .build();

    let config = FemConfig {
        base,
        use_stabilization: true,
        tau: 0.1,
        dt: Some(0.01),
        reynolds: Some(100.0),
        element_type: ElementType::Tetrahedral,
        quadrature_order: 2,
    };
    
    // Create FEM solver
    let mut solver = FemSolver::new(config);
    
    // Define boundary conditions for velocity at specific nodes
    let mut velocity_bcs = HashMap::new();
    
    // Top nodes: moving lid (velocity = (1, 0, 0))
    velocity_bcs.insert(4, Vector3::new(1.0, 0.0, 0.0));
    velocity_bcs.insert(5, Vector3::new(1.0, 0.0, 0.0));
    velocity_bcs.insert(6, Vector3::new(1.0, 0.0, 0.0));
    velocity_bcs.insert(7, Vector3::new(1.0, 0.0, 0.0));
    
    // Bottom nodes: no-slip (velocity = 0)
    velocity_bcs.insert(0, Vector3::new(0.0, 0.0, 0.0));
    velocity_bcs.insert(1, Vector3::new(0.0, 0.0, 0.0));
    velocity_bcs.insert(2, Vector3::new(0.0, 0.0, 0.0));
    velocity_bcs.insert(3, Vector3::new(0.0, 0.0, 0.0));
    
    println!("Boundary conditions set:");
    println!("  - Bottom nodes (0-3): No-slip walls");
    println!("  - Top nodes (4-7): Moving lid with velocity (1, 0, 0)");
    println!();
    
    // Create Stokes flow problem
    println!("Solving Stokes flow...");
    let problem = StokesFlowProblem::new(mesh, fluid, velocity_bcs);
    let solution = solver.solve(&problem)?;
    
    // Solution is now available in the returned StokesFlowSolution
    let velocity_solution = &solution.velocity;
    let pressure_solution = &solution.pressure;
    
    println!("\nSolution obtained!");
    
    // Calculate maximum velocity
    let max_velocity = velocity_solution.iter()
        .map(|v| v.abs())
        .fold(0.0, f64::max);
    println!("Maximum velocity component: {:.6} m/s", max_velocity);
    
    // Calculate maximum pressure
    let max_pressure = pressure_solution.iter()
        .map(|p| p.abs())
        .fold(0.0, f64::max);
    println!("Maximum pressure: {:.6} Pa", max_pressure);
    
    println!("\nFEM Stokes flow simulation completed successfully!");
    
    Ok(())
}

/// Create a simple unit cube mesh with tetrahedra
fn create_unit_cube_mesh() -> Result<Mesh<f64>, Box<dyn std::error::Error>> {
    // Create 8 vertices of a unit cube
    let vertices = vec![
        Vertex { position: Point3::new(0.0, 0.0, 0.0), id: 0 }, // 0
        Vertex { position: Point3::new(1.0, 0.0, 0.0), id: 1 }, // 1
        Vertex { position: Point3::new(1.0, 1.0, 0.0), id: 2 }, // 2
        Vertex { position: Point3::new(0.0, 1.0, 0.0), id: 3 }, // 3
        Vertex { position: Point3::new(0.0, 0.0, 1.0), id: 4 }, // 4
        Vertex { position: Point3::new(1.0, 0.0, 1.0), id: 5 }, // 5
        Vertex { position: Point3::new(1.0, 1.0, 1.0), id: 6 }, // 6
        Vertex { position: Point3::new(0.0, 1.0, 1.0), id: 7 }, // 7
    ];
    
    // Create tetrahedra by subdividing the cube
    // A cube can be divided into 5 or 6 tetrahedra
    let cells = vec![
        // Tetrahedron 1
        Cell { faces: vec![0, 1, 2, 4], id: 0, element_type: cfd_mesh::ElementType::Tetrahedron },
        // Tetrahedron 2
        Cell { faces: vec![1, 2, 4, 5], id: 1, element_type: cfd_mesh::ElementType::Tetrahedron },
        // Tetrahedron 3
        Cell { faces: vec![2, 4, 5, 6], id: 2, element_type: cfd_mesh::ElementType::Tetrahedron },
        // Tetrahedron 4
        Cell { faces: vec![0, 2, 3, 4], id: 3, element_type: cfd_mesh::ElementType::Tetrahedron },
        // Tetrahedron 5
        Cell { faces: vec![2, 3, 4, 7], id: 4, element_type: cfd_mesh::ElementType::Tetrahedron },
        // Tetrahedron 6
        Cell { faces: vec![2, 4, 6, 7], id: 5, element_type: cfd_mesh::ElementType::Tetrahedron },
    ];
    
    // Create faces (not strictly necessary for FEM, but good for completeness)
    let faces = vec![];
    
    let topology = cfd_mesh::MeshTopology {
        num_vertices: vertices.len(),
        num_edges: 0,
        num_faces: faces.len(),
        num_cells: cells.len(),
    };
    
    Ok(Mesh {
        vertices,
        edges: vec![],
        faces,
        cells,
        topology,
    })
}
