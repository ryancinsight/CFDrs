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
    // NOTE: Need refined mesh with interior nodes for proper cavity flow simulation
    let nx = 5; // Grid resolution (5x5x5 = 125 vertices for validation)
    let mesh = create_refined_cube_mesh(nx)?;
    println!(
        "Created tetrahedral mesh with {} vertices and {} cells",
        mesh.vertices().len(),
        mesh.cells().len()
    );

    // Identify boundary nodes and set boundary conditions
    // Lid-driven cavity:
    // Top face (y=1): Velocity inlet (u=1, v=0, w=0)
    // All other faces: No-slip wall
    let mut boundary_conditions = HashMap::new();
    let tolerance = 1e-6;

    for (i, vertex) in mesh.vertices().iter().enumerate() {
        let pos = vertex.position;
        let x = pos.x;
        let y = pos.y;
        let z = pos.z;

        // Check if on any boundary face
        if x.abs() < tolerance {
            // Left face (x=0): No-slip wall
            boundary_conditions.insert(i, BoundaryCondition::wall_no_slip());
        } else if (x - 1.0).abs() < tolerance {
            // Right face (x=1): No-slip wall
            boundary_conditions.insert(i, BoundaryCondition::wall_no_slip());
        } else if y.abs() < tolerance {
            // Bottom face (y=0): No-slip wall
            boundary_conditions.insert(i, BoundaryCondition::wall_no_slip());
        } else if (y - 1.0).abs() < tolerance {
            // Top face (y=1): Lid velocity
            boundary_conditions.insert(
                i,
                BoundaryCondition::velocity_inlet(Vector3::new(1.0, 0.0, 0.0)),
            );
        } else if z.abs() < tolerance {
            // Front face (z=0): No-slip wall
            boundary_conditions.insert(i, BoundaryCondition::wall_no_slip());
        } else if (z - 1.0).abs() < tolerance {
            // Back face (z=1): No-slip wall
            boundary_conditions.insert(i, BoundaryCondition::wall_no_slip());
        }
    }

    println!(
        "Boundary conditions set for {} nodes ({} total nodes, {} interior)",
        boundary_conditions.len(),
        mesh.vertices().len(),
        mesh.vertices().len() - boundary_conditions.len()
    );

    // Get number of vertices before P2 conversion (for Taylor-Hood elements)
    let n_corner_nodes = mesh.vertices().len();
    let n_cells = mesh.cells().len();

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
    println!("  - Corner nodes (P1 pressure): {}", n_corner_nodes);
    println!("  - Total cells: {}", n_cells);

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

/// Create a refined unit cube mesh with interior nodes
/// 
/// This creates an nx × nx × nx structured grid subdivided into tetrahedra.
/// Each hex cell is split into 5 or 6 tetrahedra.
///
/// # Arguments
/// * `nx` - Number of nodes in each direction (e.g., 5 gives 125 nodes)
fn create_refined_cube_mesh(nx: usize) -> Result<Mesh<f64>, Box<dyn std::error::Error>> {
    if nx < 2 {
        return Err("Need at least 2 nodes per direction".into());
    }
    
    let mut mesh = Mesh::new();
    
    // Create structured grid vertices
    let h = 1.0 / (nx as f64 - 1.0); // Grid spacing
    for k in 0..nx {
        for j in 0..nx {
            for i in 0..nx {
                let x = i as f64 * h;
                let y = j as f64 * h;
                let z = k as f64 * h;
                mesh.add_vertex(Vertex::new(Point3::new(x, y, z)));
            }
        }
    }
    
    // Helper: convert (i,j,k) to vertex index
    let vertex_index = |i: usize, j: usize, k: usize| -> usize {
        k * nx * nx + j * nx + i
    };
    
    // Subdivide each hex into 5 tetrahedra (Kuhn triangulation)
    let mut face_map: HashMap<Vec<usize>, usize> = HashMap::new();
    
    for k in 0..(nx-1) {
        for j in 0..(nx-1) {
            for i in 0..(nx-1) {
                // 8 corners of hex cell
                let v000 = vertex_index(i,   j,   k);
                let v100 = vertex_index(i+1, j,   k);
                let v010 = vertex_index(i,   j+1, k);
                let v110 = vertex_index(i+1, j+1, k);
                let v001 = vertex_index(i,   j,   k+1);
                let v101 = vertex_index(i+1, j,   k+1);
                let v011 = vertex_index(i,   j+1, k+1);
                let v111 = vertex_index(i+1, j+1, k+1);
                
                // Kuhn triangulation (5 tets)
                let tets = vec![
                    vec![v000, v100, v010, v001],
                    vec![v100, v110, v010, v101],
                    vec![v010, v110, v011, v101],
                    vec![v001, v101, v011, v010],
                    vec![v101, v111, v110, v011],
                ];
                
                for tet_verts in tets {
                    // Create 4 triangular faces for this tet
                    let face_triplets = vec![
                        vec![tet_verts[0], tet_verts[1], tet_verts[2]],
                        vec![tet_verts[0], tet_verts[1], tet_verts[3]],
                        vec![tet_verts[1], tet_verts[2], tet_verts[3]],
                        vec![tet_verts[2], tet_verts[0], tet_verts[3]],
                    ];
                    
                    let mut face_indices = Vec::new();
                    for f_verts in face_triplets {
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
            }
        }
    }
    
    Ok(mesh)
}
