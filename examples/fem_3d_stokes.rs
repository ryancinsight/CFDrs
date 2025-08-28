//! 3D FEM Stokes Flow Example
//!
//! This example demonstrates the use of the 3D FEM solver for solving
//! the Stokes equations in a simple tetrahedral domain.

use cfd_3d::fem::{FemConfig, FemSolver};
use cfd_core::domains::mesh_operations::ElementType;
use cfd_core::fluid::Fluid;
use cfd_mesh::{Cell, Mesh, Vertex};
use nalgebra::Point3;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("3D FEM Stokes Flow Example");
    println!("==========================");

    // Create fluid properties (water at 20°C)
    let fluid = Fluid::<f64>::water_20c();
    println!("Fluid: Water at 20°C");
    println!("Density: {:.1} kg/m³", fluid.density);
    println!("Viscosity: {:.6} Pa·s", fluid.dynamic_viscosity());
    println!();

    // Create a simple tetrahedral mesh
    let mesh = create_unit_cube_mesh()?;
    println!(
        "Created tetrahedral mesh with {} vertices and {} cells",
        mesh.vertices.len(),
        mesh.cells.len()
    );

    // Create FEM solver configuration
    let base = cfd_core::solver::SolverConfig::builder()
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
    let solver = FemSolver::new(config);

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
        Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            id: 0,
        }, // 0
        Vertex {
            position: Point3::new(1.0, 0.0, 0.0),
            id: 1,
        }, // 1
        Vertex {
            position: Point3::new(1.0, 1.0, 0.0),
            id: 2,
        }, // 2
        Vertex {
            position: Point3::new(0.0, 1.0, 0.0),
            id: 3,
        }, // 3
        Vertex {
            position: Point3::new(0.0, 0.0, 1.0),
            id: 4,
        }, // 4
        Vertex {
            position: Point3::new(1.0, 0.0, 1.0),
            id: 5,
        }, // 5
        Vertex {
            position: Point3::new(1.0, 1.0, 1.0),
            id: 6,
        }, // 6
        Vertex {
            position: Point3::new(0.0, 1.0, 1.0),
            id: 7,
        }, // 7
    ];

    // Create tetrahedra by subdividing the cube
    // A cube can be divided into 5 or 6 tetrahedra
    let cells = vec![
        // Tetrahedron 1
        Cell {
            vertices: vec![0, 1, 2, 4],
            element_type: ElementType::Tetrahedron,
        },
        // Tetrahedron 2
        Cell {
            vertices: vec![1, 2, 4, 5],
            element_type: ElementType::Tetrahedron,
        },
        // Tetrahedron 3
        Cell {
            vertices: vec![2, 4, 5, 6],
            element_type: ElementType::Tetrahedron,
        },
        // Tetrahedron 4
        Cell {
            vertices: vec![0, 2, 3, 4],
            element_type: ElementType::Tetrahedron,
        },
        // Tetrahedron 5
        Cell {
            vertices: vec![2, 3, 4, 7],
            element_type: ElementType::Tetrahedron,
        },
        // Tetrahedron 6
        Cell {
            vertices: vec![2, 4, 6, 7],
            element_type: ElementType::Tetrahedron,
        },
    ];

    // Create faces (not strictly necessary for FEM, but good for completeness)
    let faces = vec![];

    let topology = cfd_mesh::MeshTopology::default();

    Ok(Mesh {
        vertices,
        edges: vec![],
        faces,
        cells,
        topology,
    })
}
