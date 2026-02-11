//! 3D FEM Stokes Flow Example
//!
//! This example demonstrates the use of the 3D FEM solver for solving
//! the Stokes equations in a simple tetrahedral domain.

use cfd_3d::fem::{FemConfig, FemSolver, StokesFlowProblem};
use cfd_core::geometry::ElementType;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::{ConstantFluid, ConstantPropertyFluid};
use cfd_mesh::prelude::{Cell, Face, Mesh, Vertex};
use nalgebra::{Point3, Vector3};
use std::collections::HashMap;

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

    // Identify boundary nodes and set boundary conditions
    // Lid-driven cavity:
    // Top face (y=1): Velocity inlet (u=1, v=0, w=0)
    // Other faces: No-slip wall
    let mut boundary_conditions = HashMap::new();
    let tolerance = 1e-6;

    for (i, vertex) in mesh.vertices().iter().enumerate() {
        let pos = vertex.position;
        let x = pos.x;
        let y = pos.y;
        let z = pos.z;

        // Check if on boundary
        let on_boundary = x.abs() < tolerance
            || (x - 1.0).abs() < tolerance
            || y.abs() < tolerance
            || (y - 1.0).abs() < tolerance
            || z.abs() < tolerance
            || (z - 1.0).abs() < tolerance;

        if on_boundary {
            if (y - 1.0).abs() < tolerance {
                // Top face: Lid velocity
                boundary_conditions.insert(
                    i,
                    BoundaryCondition::velocity_inlet(Vector3::new(1.0, 0.0, 0.0)),
                );
            } else {
                // Other faces: No-slip wall
                boundary_conditions.insert(i, BoundaryCondition::wall_no_slip());
            }
        }
    }

    println!(
        "Boundary conditions set for {} nodes",
        boundary_conditions.len()
    );

    // Create problem definition
    let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions);

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
    let mut solver = FemSolver::new(config);

    println!("FEM solver configured with:");
    println!("  - Element type: Linear tetrahedron (Tet4)");
    println!("  - Stabilization: SUPG/PSPG");
    println!("  - Quadrature order: 2");
    println!("  - Reynolds number: 100");

    println!("\nStarting solver...");
    let solution = solver.solve(&problem, None)?;

    println!("FEM solver setup and execution completed successfully!");
    println!("Solution computed with {} nodes", solution.n_nodes);

    Ok(())
}

/// Create a simple unit cube mesh with tetrahedra
fn create_unit_cube_mesh() -> Result<Mesh<f64>, Box<dyn std::error::Error>> {
    let mut mesh = Mesh::new();

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

    for vertex in vertices {
        mesh.add_vertex(vertex);
    }

    // 6-tetrahedron decomposition of a cube
    // Each tet is defined by 4 vertex indices
    let tet_indices = vec![
        vec![0, 1, 2, 6],
        vec![0, 2, 3, 6],
        vec![0, 5, 1, 6],
        vec![0, 4, 5, 6],
        vec![0, 3, 7, 6],
        vec![0, 7, 4, 6],
    ];

    // Deduplicate faces
    let mut face_map: HashMap<Vec<usize>, usize> = HashMap::new();

    for indices in tet_indices {
        let v0 = indices[0];
        let v1 = indices[1];
        let v2 = indices[2];
        let v3 = indices[3];

        // 4 faces per tet (vertex triplets)
        let faces_verts = vec![
            vec![v0, v1, v2],
            vec![v0, v1, v3],
            vec![v1, v2, v3],
            vec![v2, v0, v3],
        ];

        let mut face_indices = Vec::new();

        for f_verts in faces_verts {
            let mut sorted_verts = f_verts.clone();
            sorted_verts.sort_unstable();

            let face_idx = if let Some(&idx) = face_map.get(&sorted_verts) {
                idx
            } else {
                let idx = mesh.add_face(Face::triangle(f_verts[0], f_verts[1], f_verts[2]));
                face_map.insert(sorted_verts, idx);
                idx
            };
            face_indices.push(face_idx);
        }

        mesh.add_cell(Cell::tetrahedron(
            face_indices[0],
            face_indices[1],
            face_indices[2],
            face_indices[3],
        ));
    }

    Ok(mesh)
}
