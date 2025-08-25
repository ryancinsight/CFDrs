//! Integration tests for CFD Suite
//!
//! These tests verify that different components work together correctly

use approx::assert_relative_eq;
use cfd_core::domains::fluid_dynamics::{
    FlowField, FlowOperations, KEpsilonModel, TurbulenceModel,
};
use cfd_core::values::{FlowGeometry, ReynoldsNumber};
use cfd_math::linear_solver::{ConjugateGradient, LinearSolver};
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::DVector;

#[test]
fn test_flow_field_complete_workflow() {
    // Create flow field
    let flow_field = FlowField::<f64>::new(32, 32, 32);

    // Initialize turbulence model
    let mut k_epsilon = KEpsilonModel::new();
    k_epsilon.initialize_state(&flow_field);

    // Calculate turbulent viscosity
    let nu_t = k_epsilon.turbulent_viscosity(&flow_field);
    assert_eq!(nu_t.len(), 32 * 32 * 32);

    // Calculate flow properties
    let divergence = FlowOperations::divergence(&flow_field.velocity);
    let vorticity = FlowOperations::vorticity(&flow_field.velocity);
    let kinetic_energy = FlowOperations::kinetic_energy(&flow_field.velocity);
    let enstrophy = FlowOperations::enstrophy(&flow_field.velocity);

    // Verify sizes
    assert_eq!(divergence.len(), 32 * 32 * 32);
    assert_eq!(vorticity.len(), 32 * 32 * 32);
    assert_eq!(kinetic_energy.len(), 32 * 32 * 32);
    assert_eq!(enstrophy.len(), 32 * 32 * 32);

    // For zero initial velocity, all should be zero
    for &div in &divergence {
        assert_relative_eq!(div, 0.0, epsilon = 1e-10);
    }
}

#[test]
fn test_reynolds_number_classification() {
    // Test laminar flow
    let re_laminar = ReynoldsNumber::new(1000.0, FlowGeometry::Pipe).unwrap();
    assert!(re_laminar.is_laminar());
    assert!(!re_laminar.is_turbulent());

    // Test transitional flow
    let re_transitional = ReynoldsNumber::new(2500.0, FlowGeometry::Pipe).unwrap();
    assert!(re_transitional.is_transitional());
    assert!(!re_transitional.is_laminar());
    assert!(!re_transitional.is_turbulent());

    // Test turbulent flow
    let re_turbulent = ReynoldsNumber::new(5000.0, FlowGeometry::Pipe).unwrap();
    assert!(re_turbulent.is_turbulent());
    assert!(!re_turbulent.is_laminar());
}

#[test]
fn test_linear_solver_integration() {
    // Create a simple Poisson problem
    let n = 50;
    let mut builder = SparseMatrixBuilder::new(n, n);

    // Build tridiagonal matrix (1D Poisson)
    for i in 0..n {
        builder.add_entry(i, i, 2.0).unwrap();
        if i > 0 {
            builder.add_entry(i, i - 1, -1.0).unwrap();
        }
        if i < n - 1 {
            builder.add_entry(i, i + 1, -1.0).unwrap();
        }
    }

    let matrix = builder.build().unwrap();
    let b = DVector::from_element(n, 1.0);

    // Solve with Conjugate Gradient
    let solver = ConjugateGradient::<f64>::default();
    let solution = solver.solve(&matrix, &b, None).unwrap();

    // Verify solution satisfies Ax = b
    let residual = &matrix * &solution - &b;
    let residual_norm = residual.norm();
    assert!(
        residual_norm < 1e-6,
        "Residual too large: {}",
        residual_norm
    );
}

#[test]
fn test_mesh_element_types() {
    use cfd_core::domains::mesh_operations::ElementType;
    use cfd_mesh::{Cell, Mesh, Vertex};
    use nalgebra::Point3;

    // Create vertices for a tetrahedron
    let vertices = vec![
        Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            id: 0,
        },
        Vertex {
            position: Point3::new(1.0, 0.0, 0.0),
            id: 1,
        },
        Vertex {
            position: Point3::new(0.0, 1.0, 0.0),
            id: 2,
        },
        Vertex {
            position: Point3::new(0.0, 0.0, 1.0),
            id: 3,
        },
    ];

    // Create a tetrahedral cell
    let cell = Cell {
        element_type: ElementType::Tetrahedron,
        vertices: vec![0, 1, 2, 3],
    };

    // Verify element properties
    assert_eq!(cell.element_type.num_nodes(), 4);
    assert_eq!(cell.element_type.dimension(), 3);
    assert!(cell.element_type.is_linear());
    assert!(!cell.element_type.is_quadratic());

    // Create mesh
    let mesh = Mesh {
        vertices,
        edges: vec![],
        faces: vec![],
        cells: vec![cell],
        topology: Default::default(),
    };

    assert_eq!(mesh.vertices.len(), 4);
    assert_eq!(mesh.cells.len(), 1);
}

#[test]
fn test_turbulence_model_workflow() {
    // Create flow field with non-zero velocity
    let mut flow_field = FlowField::<f64>::new(10, 10, 10);

    // Set some velocity values
    for i in 0..10 {
        for j in 0..10 {
            for k in 0..10 {
                let idx = i * 100 + j * 10 + k;
                flow_field.velocity.components[idx] =
                    nalgebra::Vector3::new(0.1 * i as f64, 0.05 * j as f64, 0.02 * k as f64);
            }
        }
    }

    // Test k-epsilon model
    let mut k_epsilon = KEpsilonModel::new();
    k_epsilon.initialize_state(&flow_field);

    let nu_t = k_epsilon.turbulent_viscosity(&flow_field);
    let tke = k_epsilon.turbulent_kinetic_energy(&flow_field);

    assert_eq!(nu_t.len(), 1000);
    assert_eq!(tke.len(), 1000);

    // Turbulent viscosity should be non-negative
    for &nu in &nu_t {
        assert!(nu >= 0.0, "Negative turbulent viscosity: {}", nu);
    }
}

#[test]
fn test_flow_operations_consistency() {
    // Create a simple flow field
    let flow_field = FlowField::<f64>::new(20, 20, 20);

    // Calculate properties multiple times - should be consistent
    let div1 = FlowOperations::divergence(&flow_field.velocity);
    let div2 = FlowOperations::divergence(&flow_field.velocity);

    for i in 0..div1.len() {
        assert_eq!(div1[i], div2[i], "Divergence calculation not consistent");
    }

    let ke1 = FlowOperations::kinetic_energy(&flow_field.velocity);
    let ke2 = FlowOperations::kinetic_energy(&flow_field.velocity);

    for i in 0..ke1.len() {
        assert_eq!(ke1[i], ke2[i], "Kinetic energy calculation not consistent");
    }
}

#[test]
fn test_end_to_end_cfd_workflow() {
    // This test simulates a complete CFD workflow

    // 1. Setup problem
    let grid_size = 16;
    let flow_field = FlowField::<f64>::new(grid_size, grid_size, grid_size);

    // 2. Calculate Reynolds number
    let velocity = 1.0;
    let length = 0.1;
    let viscosity = 1e-6;
    let reynolds = velocity * length / viscosity;
    let re = ReynoldsNumber::new(reynolds, FlowGeometry::Pipe).unwrap();

    assert!(re.is_turbulent());

    // 3. Setup turbulence model based on Reynolds number
    let mut turbulence_model = KEpsilonModel::new();
    turbulence_model.initialize_state(&flow_field);

    // 4. Calculate flow properties
    let divergence = FlowOperations::divergence(&flow_field.velocity);
    let vorticity = FlowOperations::vorticity(&flow_field.velocity);

    // 5. Setup and solve a linear system (pressure correction)
    let n = grid_size * grid_size;
    let mut builder = SparseMatrixBuilder::new(n, n);

    // Build a simple Laplacian
    for i in 0..n {
        let mut diagonal = 0.0f64;

        // Add off-diagonal entries
        if i >= grid_size {
            builder.add_entry(i, i - grid_size, -1.0f64).unwrap();
            diagonal += 1.0f64;
        }
        if i < n - grid_size {
            builder.add_entry(i, i + grid_size, -1.0f64).unwrap();
            diagonal += 1.0f64;
        }
        if i % grid_size != 0 {
            builder.add_entry(i, i - 1, -1.0f64).unwrap();
            diagonal += 1.0f64;
        }
        if i % grid_size != grid_size - 1 {
            builder.add_entry(i, i + 1, -1.0f64).unwrap();
            diagonal += 1.0f64;
        }

        // Add diagonal entry - ensure matrix is well-conditioned
        builder.add_entry(i, i, diagonal.max(4.0f64)).unwrap();
    }

    let matrix = builder.build().unwrap();
    let rhs = DVector::from_element(n, 0.1);

    let solver = ConjugateGradient::<f64>::default();
    let pressure_correction = solver.solve(&matrix, &rhs, None).unwrap();

    // 6. Verify results
    assert_eq!(pressure_correction.len(), n);
    assert_eq!(divergence.len(), grid_size * grid_size * grid_size);
    assert_eq!(vorticity.len(), grid_size * grid_size * grid_size);

    // Check that solver converged
    let residual = &matrix * &pressure_correction - &rhs;
    let residual_norm = residual.norm();
    let rhs_norm = rhs.norm();
    let relative_residual = if rhs_norm > 1e-10 {
        residual_norm / rhs_norm
    } else {
        residual_norm
    };

    // Use a more realistic tolerance for iterative solver
    assert!(
        relative_residual < 1e-3,
        "Solver did not converge: relative residual = {:.6e}",
        relative_residual
    );
}
