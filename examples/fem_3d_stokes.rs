//! 3D FEM Stokes Flow Example
//!
//! This example demonstrates the use of the 3D FEM solver for solving
//! the Stokes equations in a simple tetrahedral domain (lid-driven cavity).
//!
//! # Geometry
//!
//! Unit cube [0,1]^3 discretized into tetrahedra via Kuhn triangulation.
//!
//! # Boundary Conditions
//!
//! - Top face (y=1): Lid velocity u = (1, 0, 0)
//! - All other faces: No-slip wall
//!
//! # Reference
//!
//! Ghia, Ghia & Shin (1982). *J. Comp. Physics*, 48, pp. 387-411.

use cfd_3d::fem::{FemConfig, FemSolver, StokesFlowProblem};
use cfd_core::geometry::ElementType;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::domain::core::index::VertexId;
use cfd_mesh::domain::topology::Cell;
use cfd_mesh::IndexedMesh;
use nalgebra::{Point3, Vector3};
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("3D FEM Stokes Flow Example");
    println!("==========================");

    // Create fluid properties (water at 20C)
    let fluid = ConstantPropertyFluid::<f64>::water_20c()?;
    println!("Fluid: Water at 20C");
    println!("Density: {:.1} kg/m^3", fluid.density);
    println!("Viscosity: {:.6} Pa.s", fluid.viscosity);
    println!();

    // Create a simple tetrahedral mesh
    let nx = 5; // Grid resolution (5x5x5 = 125 vertices)
    let mesh = create_refined_cube_mesh(nx)?;
    println!(
        "Created tetrahedral mesh with {} vertices and {} cells",
        mesh.vertices.len(),
        mesh.cells.len()
    );

    // Identify boundary nodes and set boundary conditions
    // Lid-driven cavity:
    // Top face (y=1): Velocity inlet (u=1, v=0, w=0)
    // All other faces: No-slip wall
    let mut boundary_conditions = HashMap::new();
    let tolerance = 1e-6;

    for (vid, vdata) in mesh.vertices.iter() {
        let i = vid.as_usize();
        let x = vdata.position.x;
        let y = vdata.position.y;
        let z = vdata.position.z;

        // Check if on any boundary face
        if x.abs() < tolerance
            || (x - 1.0).abs() < tolerance
            || y.abs() < tolerance
            || z.abs() < tolerance
            || (z - 1.0).abs() < tolerance
        {
            // Wall faces: No-slip
            boundary_conditions.insert(i, BoundaryCondition::wall_no_slip());
        } else if (y - 1.0).abs() < tolerance {
            // Top face (y=1): Lid velocity
            boundary_conditions.insert(
                i,
                BoundaryCondition::velocity_inlet(Vector3::new(1.0, 0.0, 0.0)),
            );
        }
    }

    let n_corner_nodes = mesh.vertices.len();
    let n_cells = mesh.cells.len();

    println!(
        "Boundary conditions set for {} nodes ({} total nodes, {} interior)",
        boundary_conditions.len(),
        n_corner_nodes,
        n_corner_nodes - boundary_conditions.len()
    );

    // Create problem definition
    let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, n_corner_nodes);

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
        grad_div_penalty: 0.0,
    };

    // Create FEM solver
    let mut solver = FemSolver::new(config);

    println!("FEM solver configured with:");
    println!("  - Element type: Linear tetrahedron (Tet4)");
    println!("  - Stabilization: SUPG/PSPG");
    println!("  - Quadrature order: 2");
    println!("  - Reynolds number: 100");
    println!("  - Total nodes: {}", n_corner_nodes);
    println!("  - Total cells: {}", n_cells);

    println!("\nStarting solver...");
    let solution = solver.solve(&problem, None)?;

    println!("FEM solver setup and execution completed successfully!");
    println!("Solution computed with {} nodes", solution.n_nodes);

    Ok(())
}

/// Create a refined unit cube mesh with interior nodes.
///
/// Creates an nx x nx x nx structured grid subdivided into tetrahedra
/// using Kuhn triangulation (5 tets per hex cell).
fn create_refined_cube_mesh(nx: usize) -> Result<IndexedMesh<f64>, Box<dyn std::error::Error>> {
    if nx < 2 {
        return Err("Need at least 2 nodes per direction".into());
    }

    let mut mesh = IndexedMesh::<f64>::new();

    // Create structured grid vertices
    let h = 1.0 / (nx as f64 - 1.0);
    for k in 0..nx {
        for j in 0..nx {
            for i in 0..nx {
                let x = i as f64 * h;
                let y = j as f64 * h;
                let z = k as f64 * h;
                mesh.vertices
                    .insert_unique(Point3::new(x, y, z), Vector3::z());
            }
        }
    }

    // Helper: convert (i,j,k) to vertex index
    let vertex_index = |i: usize, j: usize, k: usize| -> usize { k * nx * nx + j * nx + i };

    // Subdivide each hex into 5 tetrahedra (Kuhn triangulation)
    // Track faces for deduplication
    let mut face_map: HashMap<[u32; 3], usize> = HashMap::new();
    let mut face_counter: usize = 0;

    for k in 0..(nx - 1) {
        for j in 0..(nx - 1) {
            for i in 0..(nx - 1) {
                // 8 corners of hex cell
                let v000 = vertex_index(i, j, k);
                let v100 = vertex_index(i + 1, j, k);
                let v010 = vertex_index(i, j + 1, k);
                let v110 = vertex_index(i + 1, j + 1, k);
                let v001 = vertex_index(i, j, k + 1);
                let v101 = vertex_index(i + 1, j, k + 1);
                let v011 = vertex_index(i, j + 1, k + 1);
                let _v111 = vertex_index(i + 1, j + 1, k + 1);

                // Kuhn triangulation (5 tets)
                let tets = [
                    [v000, v100, v010, v001],
                    [v100, v110, v010, v101],
                    [v010, v110, v011, v101],
                    [v001, v101, v011, v010],
                    [v101, _v111, v110, v011],
                ];

                for tet_verts in &tets {
                    // Create 4 triangular faces for this tet
                    let face_triplets = [
                        [tet_verts[0], tet_verts[1], tet_verts[2]],
                        [tet_verts[0], tet_verts[1], tet_verts[3]],
                        [tet_verts[1], tet_verts[2], tet_verts[3]],
                        [tet_verts[2], tet_verts[0], tet_verts[3]],
                    ];

                    let mut cell_face_indices = Vec::with_capacity(4);
                    for f_verts in &face_triplets {
                        let mut sorted = [f_verts[0] as u32, f_verts[1] as u32, f_verts[2] as u32];
                        sorted.sort_unstable();

                        let face_idx = if let Some(&idx) = face_map.get(&sorted) {
                            idx
                        } else {
                            let v0 = VertexId(f_verts[0] as u32);
                            let v1 = VertexId(f_verts[1] as u32);
                            let v2 = VertexId(f_verts[2] as u32);
                            mesh.faces.add_triangle(v0, v1, v2);
                            let idx = face_counter;
                            face_map.insert(sorted, idx);
                            face_counter += 1;
                            idx
                        };
                        cell_face_indices.push(face_idx);
                    }

                    mesh.cells.push(Cell::tetrahedron(
                        cell_face_indices[0],
                        cell_face_indices[1],
                        cell_face_indices[2],
                        cell_face_indices[3],
                    ));
                }
            }
        }
    }

    Ok(mesh)
}
