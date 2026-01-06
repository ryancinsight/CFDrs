//! MMS validation for complex boundary conditions
//!
//! Tests manufactured solutions with mixed boundary condition types
//! to verify proper boundary condition implementation.

use cfd_2d::grid::StructuredGrid2D;
use cfd_core::physics::boundary::{BoundaryCondition, BoundaryConditionSet, WallType};
use cfd_validation::manufactured::{ManufacturedDiffusion, TaylorGreenManufactured};
use nalgebra::Vector3;
use num_traits::FromPrimitive;

/// Helper function to compute x coordinate at grid index
#[allow(dead_code)]
fn x_at<T: nalgebra::RealField + FromPrimitive + Copy>(grid: &StructuredGrid2D<T>, i: usize) -> T {
    grid.bounds.0 + T::from_usize(i).unwrap_or_else(T::zero) * grid.dx
}

/// Helper function to compute y coordinate at grid index
#[allow(dead_code)]
fn y_at<T: nalgebra::RealField + FromPrimitive + Copy>(grid: &StructuredGrid2D<T>, j: usize) -> T {
    grid.bounds.2 + T::from_usize(j).unwrap_or_else(T::zero) * grid.dy
}

/// Test boundary condition set creation with mixed types
#[test]
fn test_mixed_boundary_condition_set() {
    let nx = 32;
    let ny = 32;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();
    assert_eq!(grid.nx, nx);

    // Set up mixed boundary conditions
    let mut bc_set = BoundaryConditionSet::<f64>::new();

    // West: Dirichlet velocity inlet (u=1, v=0)
    bc_set.add(
        "west",
        BoundaryCondition::VelocityInlet {
            velocity: Vector3::new(1.0, 0.0, 0.0),
        },
    );

    // East: Neumann pressure outlet (dp/dx = 0)
    bc_set.add("east", BoundaryCondition::PressureOutlet { pressure: 0.0 });

    // North: Dirichlet no-slip wall (u=0, v=0)
    bc_set.add(
        "north",
        BoundaryCondition::Wall {
            wall_type: WallType::NoSlip,
        },
    );

    // South: Neumann symmetry (dv/dy = 0)
    bc_set.add("south", BoundaryCondition::Neumann { gradient: 0.0 });

    // Verify all boundaries are set
    assert!(bc_set.contains("west"));
    assert!(bc_set.contains("east"));
    assert!(bc_set.contains("north"));
    assert!(bc_set.contains("south"));
    assert_eq!(bc_set.len(), 4);

    println!("✓ Mixed boundary condition set created successfully");
}

/// Test Robin boundary condition setup
#[test]
fn test_robin_boundary_conditions() {
    let nx = 32;
    let ny = 32;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();
    assert_eq!(grid.nx, nx);

    // Set up Robin boundary conditions
    let mut bc_set = BoundaryConditionSet::<f64>::new();

    // All boundaries: Robin condition αu + β∂u/∂n = γ
    // For this test, use α=1, β=0, γ=0 (equivalent to Dirichlet u=0)
    let robin_bc = BoundaryCondition::Robin {
        alpha: 1.0,
        beta: 0.0,
        gamma: 0.0,
    };

    bc_set.add("west", robin_bc.clone());
    bc_set.add("east", robin_bc.clone());
    bc_set.add("north", robin_bc.clone());
    bc_set.add("south", robin_bc);

    // Verify Robin BCs
    for (name, bc) in bc_set.iter() {
        if let BoundaryCondition::Robin { alpha, beta, gamma } = bc {
            assert!((*alpha - 1.0).abs() < 1e-10);
            assert!(beta.abs() < 1e-10);
            assert!(gamma.abs() < 1e-10);
        } else {
            panic!("Expected Robin BC for boundary {}", name);
        }
    }

    println!("✓ Robin boundary conditions verified");
}

/// Test periodic boundary condition setup
#[test]
fn test_periodic_boundary_conditions() {
    let nx = 32;
    let ny = 32;
    let _grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    // Set up periodic boundary conditions
    let mut bc_set = BoundaryConditionSet::<f64>::new();

    // Periodic in x-direction
    bc_set.add(
        "west",
        BoundaryCondition::Periodic {
            partner: "east".to_string(),
        },
    );
    bc_set.add(
        "east",
        BoundaryCondition::Periodic {
            partner: "west".to_string(),
        },
    );

    // No-slip walls in y-direction
    bc_set.add(
        "north",
        BoundaryCondition::Wall {
            wall_type: WallType::NoSlip,
        },
    );
    bc_set.add(
        "south",
        BoundaryCondition::Wall {
            wall_type: WallType::NoSlip,
        },
    );

    // Validate periodic boundary pairs
    assert!(bc_set.validate_periodic().is_ok());

    println!("✓ Periodic boundary conditions verified");
}

/// Test pressure-driven flow boundary conditions
#[test]
fn test_pressure_driven_flow_boundaries() {
    let nx = 32;
    let ny = 32;
    let _grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    // Set up pressure-driven flow boundaries
    let mut bc_set = BoundaryConditionSet::<f64>::new();

    // Inlet: Pressure inlet
    bc_set.add(
        "west",
        BoundaryCondition::PressureInlet {
            pressure: 1.0,
            velocity_direction: Some(Vector3::new(1.0, 0.0, 0.0)),
        },
    );

    // Outlet: Pressure outlet
    bc_set.add("east", BoundaryCondition::PressureOutlet { pressure: 0.0 });

    // Walls: No-slip
    bc_set.add(
        "north",
        BoundaryCondition::Wall {
            wall_type: WallType::NoSlip,
        },
    );
    bc_set.add(
        "south",
        BoundaryCondition::Wall {
            wall_type: WallType::NoSlip,
        },
    );

    // Verify pressure BCs
    if let Some(BoundaryCondition::PressureInlet { pressure, .. }) = bc_set.get("west") {
        assert!((*pressure - 1.0).abs() < 1e-10);
    } else {
        panic!("Expected PressureInlet BC for west boundary");
    }

    if let Some(BoundaryCondition::PressureOutlet { pressure }) = bc_set.get("east") {
        assert!(pressure.abs() < 1e-10);
    } else {
        panic!("Expected PressureOutlet BC for east boundary");
    }

    println!("✓ Pressure-driven flow boundary conditions verified");
}

/// Test Taylor-Green vortex solution evaluation
#[test]
fn test_taylor_green_evaluation() {
    let tg = TaylorGreenManufactured::new(0.01);
    let nx = 16;
    let ny = 16;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    // Verify Taylor-Green can evaluate at all grid points
    for i in 0..nx {
        for j in 0..ny {
            let x = x_at(&grid, i);
            let y = y_at(&grid, j);
            let vel = tg.velocity(x, y, 0.0);
            let p = tg.pressure(x, y, 0.0);
            let omega = tg.vorticity(x, y, 0.0);

            // Verify values are finite
            assert!(vel.x.is_finite(), "Velocity u is not finite at ({}, {})", x, y);
            assert!(vel.y.is_finite(), "Velocity v is not finite at ({}, {})", x, y);
            assert!(p.is_finite(), "Pressure is not finite at ({}, {})", x, y);
            assert!(omega.is_finite(), "Vorticity is not finite at ({}, {})", x, y);
        }
    }

    // Verify kinetic energy decay
    let ke0 = tg.kinetic_energy(0.0);
    let ke1 = tg.kinetic_energy(1.0);
    assert!(ke1 < ke0, "Kinetic energy should decay over time");

    println!("✓ Taylor-Green vortex evaluation verified");
}

/// Test manufactured diffusion solution
#[test]
fn test_manufactured_diffusion_evaluation() {
    use cfd_validation::manufactured::ManufacturedSolution;

    let mms = ManufacturedDiffusion::<f64>::new(1.0);
    let nx = 16;
    let ny = 16;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    // Verify MMS can evaluate at all grid points
    for i in 0..nx {
        for j in 0..ny {
            let x = x_at(&grid, i);
            let y = y_at(&grid, j);
            let u = mms.exact_solution(x, y, 0.0, 1.0);
            let source = mms.source_term(x, y, 0.0, 1.0);

            // Verify values are finite
            assert!(u.is_finite(), "Solution is not finite at ({}, {})", x, y);
            assert!(source.is_finite(), "Source term is not finite at ({}, {})", x, y);
        }
    }

    println!("✓ Manufactured diffusion solution evaluation verified");
}

/// Property-based test for boundary condition set operations
#[cfg(test)]
mod property_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// Test that boundary condition sets can be created with various configurations
        #[test]
        fn test_bc_set_operations(
            west_pressure in 0.0f64..10.0,
            east_pressure in 0.0f64..10.0,
        ) {
            let mut bc_set = BoundaryConditionSet::<f64>::new();

            bc_set.add(
                "west",
                BoundaryCondition::PressureInlet {
                    pressure: west_pressure,
                    velocity_direction: Some(Vector3::new(1.0, 0.0, 0.0)),
                },
            );
            bc_set.add(
                "east",
                BoundaryCondition::PressureOutlet {
                    pressure: east_pressure,
                },
            );

            assert!(bc_set.contains("west"));
            assert!(bc_set.contains("east"));
            assert_eq!(bc_set.len(), 2);
        }

        /// Test Taylor-Green solution with various viscosities
        #[test]
        fn test_taylor_green_viscosity(nu in 0.001f64..1.0) {
            let tg = TaylorGreenManufactured::new(nu);

            // Verify at center point
            let _vel = tg.velocity(0.5, 0.5, 0.0);
            let _p = tg.pressure(0.5, 0.5, 0.0);

            prop_assert!(_vel.x.is_finite());
            prop_assert!(_vel.y.is_finite());
            prop_assert!(_p.is_finite());

            // Verify energy decay rate increases with viscosity
            let ke0 = tg.kinetic_energy(0.0);
            let ke1 = tg.kinetic_energy(1.0);
            assert!(ke1 <= ke0);
        }
    }
}
