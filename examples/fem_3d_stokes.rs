//! 3D FEM Stokes Flow Example
//!
//! This example demonstrates the use of the 3D FEM solver for solving
//! the Stokes equations in a simple tetrahedral domain.

use cfd_3d::fem::{FemConfig, FemSolver};
use cfd_core::domain::mesh::ElementType;
use cfd_core::physics::fluid::{ConstantFluid, ConstantPropertyFluid};
use cfd_mesh::prelude::{Cell, Mesh, Vertex};
use nalgebra::Point3;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("3D FEM Stokes Flow Example");
    println!("==========================");

    // Create fluid properties (water at 20°C)
    let fluid = ConstantPropertyFluid::<f64>::water_20c()?;
    println!("Fluid: Water at 20°C");
    println!("Density: {:.1} kg/m³", fluid.density());
    println!("Viscosity: {:.6} Pa·s", fluid.dynamic_viscosity());
    println!();

    // Create a simple tetrahedral mesh
    let mesh = create_unit_cube_mesh()?;
    println!(
        "Created tetrahedral mesh with {} vertices and {} cells",
        mesh.vertices().len(),
        mesh.cells().len()
    );

    // Create FEM solver configuration
    let base = cfd_core::compute::solver::SolverConfig::builder()
        .tolerance(1e-6)
        .max_iterations(1000)
        .build();

    let config = FemConfig {
        base,
        use_stabilization: true,
        tau: 0.1,
        dt: Some(0.01),
        reynolds: Some(100.0),
        element_type: ElementType::Tetrahedron,
        quadrature_order: 2,
    };

    // Create FEM solver
    let _solver = FemSolver::new(config);

    println!("FEM solver configured with:");
    println!("  - Element type: Linear tetrahedron (Tet4)");
    println!("  - Stabilization: SUPG/PSPG");
    println!("  - Quadrature order: 2");
    println!("  - Reynolds number: 100");

    // Note: The actual problem setup and solving would require proper
    // boundary condition implementation that matches the API.
    // This example demonstrates the basic setup.

    println!("\nFEM solver setup completed successfully!");
    println!("Note: Full Stokes flow solving requires proper boundary condition setup.");

    Ok(())
}

/// Create a simple unit cube mesh with tetrahedra
fn create_unit_cube_mesh() -> Result<Mesh<f64>, Box<dyn std::error::Error>> {
    // Create 8 vertices of a unit cube
    let vertices = vec![
        Vertex::new(Point3::new(0.0, 0.0, 0.0)), // 0
        Vertex::new(Point3::new(1.0, 0.0, 0.0)), // 1
        Vertex::new(Point3::new(1.0, 1.0, 0.0)), // 2
        Vertex::new(Point3::new(0.0, 1.0, 0.0)), // 3
        Vertex::new(Point3::new(0.0, 0.0, 1.0)), // 4
        Vertex::new(Point3::new(1.0, 0.0, 1.0)), // 5
        Vertex::new(Point3::new(1.0, 1.0, 1.0)), // 6
        Vertex::new(Point3::new(0.0, 1.0, 1.0)), // 7
    ];

    // Create tetrahedra by subdividing the cube
    // A cube can be divided into 5 or 6 tetrahedra
    let cells = vec![
        // Tetrahedron 1
        Cell::tetrahedron(0, 1, 2, 3), // Tetrahedron with faces 0,1,2,3
        Cell::tetrahedron(4, 5, 6, 7), // Tetrahedron with faces 4,5,6,7
        Cell::tetrahedron(8, 9, 10, 11), // Tetrahedron with faces 8,9,10,11
        Cell::tetrahedron(12, 13, 14, 15), // Tetrahedron with faces 12,13,14,15
        Cell::tetrahedron(16, 17, 18, 19), // Tetrahedron with faces 16,17,18,19
        Cell::tetrahedron(20, 21, 22, 23), // Tetrahedron with faces 20,21,22,23
    ];

    // Create faces (not strictly necessary for FEM, but good for completeness)
    let _faces: Vec<u8> = vec![];

    let mut mesh = Mesh::new();

    // Add vertices to mesh
    for vertex in vertices {
        mesh.add_vertex(vertex);
    }

    // Add cells to mesh
    for cell in cells {
        mesh.add_cell(cell);
    }

    Ok(mesh)
}
