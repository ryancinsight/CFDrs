use cfd_3d::fem::{FemConfig, FemSolver, StokesFlowProblem};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::IndexedMesh;
use cfd_mesh::domain::topology::{Cell, Face};
use nalgebra::Point3;
use std::collections::HashMap;

#[test]
fn test_robin_bc_assembly() {
    // Create a single tetrahedron mesh
    let mut mesh = IndexedMesh::new();
    let v0 = mesh.add_vertex_pos(Point3::new(0.0, 0.0, 0.0));
    let v1 = mesh.add_vertex_pos(Point3::new(1.0, 0.0, 0.0));
    let v2 = mesh.add_vertex_pos(Point3::new(0.0, 1.0, 0.0));
    let v3 = mesh.add_vertex_pos(Point3::new(0.0, 0.0, 1.0));

    let f0 = mesh.add_face(v0, v1, v2).0;
    let f1 = mesh.add_face(v0, v1, v3).0;
    let f2 = mesh.add_face(v1, v2, v3).0;
    let f3 = mesh.add_face(v2, v0, v3).0;

    mesh.add_cell(Cell::tetrahedron(f0 as usize, f1 as usize, f2 as usize, f3 as usize));

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

    let n_corner_nodes = mesh.vertex_count();
    let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, n_corner_nodes);

    // Setup Solver
    let config = FemConfig::default();
    let mut solver = FemSolver::new(config);

    // Solve (should not panic)
    let result = solver.solve(&problem, None);

    // We expect it to run. The physical result might be nonsense due to simple mesh/BC setup,
    // but the matrix assembly should work.
    assert!(
        result.is_ok(),
        "Solver failed with error: {:?}",
        result.err()
    );
}
