use cfd_3d::fem::{FemConfig, FemSolver, StokesFlowProblem};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::mesh::Mesh;
use cfd_mesh::topology::{Cell, Face, Vertex};
use nalgebra::Point3;
use std::collections::HashMap;

#[test]
fn test_robin_bc_assembly() {
    // Create a single tetrahedron mesh
    let mut mesh = Mesh::new();
    mesh.add_vertex(Vertex::new(Point3::new(0.0, 0.0, 0.0)));
    mesh.add_vertex(Vertex::new(Point3::new(1.0, 0.0, 0.0)));
    mesh.add_vertex(Vertex::new(Point3::new(0.0, 1.0, 0.0)));
    mesh.add_vertex(Vertex::new(Point3::new(0.0, 0.0, 1.0)));

    let f0 = mesh.add_face(Face::triangle(0, 1, 2));
    let f1 = mesh.add_face(Face::triangle(0, 1, 3));
    let f2 = mesh.add_face(Face::triangle(1, 2, 3));
    let f3 = mesh.add_face(Face::triangle(2, 0, 3));

    mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));

    // Setup Fluid
    let fluid = ConstantPropertyFluid::water_20c().unwrap();

    // Setup Boundary Conditions (Robin on all nodes)
    // Robin: 1.0 * u + 1.0 * du/dn = 2.0
    let bc = BoundaryCondition::Robin {
        alpha: 1.0,
        beta: 1.0,
        gamma: 2.0,
    };

    let mut boundary_conditions = HashMap::new();
    for i in 0..4 {
        boundary_conditions.insert(i, bc.clone());
    }

    let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions);

    // Setup Solver
    let config = FemConfig::default();
    let mut solver = FemSolver::new(config);

    // Solve (should not panic)
    let result = solver.solve(&problem);

    // We expect it to run. The physical result might be nonsense due to simple mesh/BC setup,
    // but the matrix assembly should work.
    assert!(result.is_ok(), "Solver failed with error: {:?}", result.err());
}
