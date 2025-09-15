//! Integration test demonstrating working CFD functionality

use cfd_core::fluid::{ConstantFluid, ConstantPropertyFluid};

#[test]
fn test_fluid_properties() {
    // Create fluid
    let fluid = ConstantPropertyFluid {
        name: "water".to_string(),
        density: 1000.0,
        viscosity: 0.001,
        specific_heat: 4186.0,
        thermal_conductivity: 0.6,
    };

    assert_eq!(fluid.name, "water");
    let density_diff: f64 = fluid.density - 1000.0;
    assert!(density_diff.abs() < 1e-10);
    let viscosity_diff: f64 = fluid.kinematic_viscosity() - 1e-6;
    assert!(viscosity_diff.abs() < 1e-7);

    // Test Reynolds number calculation using the fluid properties
    let velocity = 1.0; // 1 m/s
    let length = 0.1; // 0.1 m
                      // Reynolds number = density * velocity * length / dynamic_viscosity
    let reynolds = fluid.density() * velocity * length / fluid.dynamic_viscosity();

    // Re = 1000 * 1 * 0.1 / 0.001 = 100,000
    assert!(
        (reynolds - 100_000.0).abs() < 1.0,
        "Reynolds number should be ~100,000"
    );
}

#[test]
fn test_2d_grid_creation() {
    use cfd_2d::grid::{Grid2D, StructuredGrid2D};

    let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)
        .expect("Grid creation should succeed");

    // Use Grid2D trait methods
    assert_eq!(grid.nx(), 10);
    assert_eq!(grid.ny(), 10);
    assert_eq!(grid.nx() * grid.ny(), 100);

    // Test grid spacing
    let (dx, dy) = grid.spacing();
    let dx_diff: f64 = dx - 0.1;
    assert!(dx_diff.abs() < 1e-10);
    let dy_diff: f64 = dy - 0.1;
    assert!(dy_diff.abs() < 1e-10);
}

#[test]
fn test_analytical_validation() {
    use cfd_validation::analytical::{AnalyticalSolution, PoiseuilleFlow};

    // Poiseuille flow parameters
    let channel_width = 0.1; // 10 cm channel width
    let viscosity = 0.001; // water
    let pressure_gradient = -100.0; // Pa/m
    let _length = 1.0; // 1 m length

    // Calculate max velocity for channel flow
    let u_max = -pressure_gradient * channel_width * channel_width / (8.0 * viscosity);

    use cfd_validation::analytical::PoiseuilleGeometry;

    let poiseuille = PoiseuilleFlow::create(
        u_max,
        channel_width,
        pressure_gradient,
        viscosity,
        PoiseuilleGeometry::Plates,
    );

    // Test velocity using the trait method
    let t = 0.0; // time
                 // At center of channel (y = 0 for symmetric channel)
    let v_center = poiseuille.velocity(0.0, 0.0, 0.0, t);

    // Center velocity should be close to u_max
    let center_velocity: f64 = v_center[0]; // x-component for channel flow
    assert!(
        (center_velocity - u_max).abs() / u_max < 0.01,
        "Center velocity should match analytical: got {}, expected {}",
        center_velocity,
        u_max
    );

    // Test velocity at wall (y = channel_width, the boundary)
    let v_wall = poiseuille.velocity(0.0, channel_width, 0.0, t);

    // At wall, velocity should be zero (no-slip)
    let wall_velocity: f64 = v_wall[0];
    assert!(wall_velocity.abs() < 1e-10, "Wall velocity should be zero");
}

#[test]
fn test_fem_matrix_assembly() {
    use cfd_3d::fem::{ElementMatrices, FluidElement};

    // Create a simple tetrahedral element
    let nodes = vec![0, 1, 2, 3];
    let element: FluidElement<f64> = FluidElement::new(nodes);

    // Element should have 4 nodes
    assert_eq!(element.nodes.len(), 4);

    // Create element matrices
    let n_dof = 16; // 4 nodes * (3 velocity + 1 pressure)
    let matrices: ElementMatrices<f64> = ElementMatrices::new(n_dof);

    // Matrices should be correct size
    assert_eq!(matrices.k_e.nrows(), 16);
    assert_eq!(matrices.k_e.ncols(), 16);
    assert_eq!(matrices.m_e.nrows(), 16);
}
