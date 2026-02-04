//! 2D Blood Flow in a Venturi (IBM + Non-Newtonian)
//!
//! This example simulates blood flow through a Venturi constriction using:
//! - 2D Incompressible Navier-Stokes (SIMPLEC)
//! - Immersed Boundary Method (IBM) for the geometry
//! - Carreau-Yasuda non-Newtonian viscosity model
//!
//! Domain: 4cm x 1cm
//! Throat: 0.5cm (50% stenosis)
//! Inlet Velocity: Parabolic profile, max 0.5 m/s

use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::physics::immersed_boundary::{BoundaryPoint, ImmersedBoundaryMethod, ImmersedBoundaryConfig};
use cfd_2d::simplec_pimple::{SimplecPimpleConfig, SimplecPimpleSolver};
use cfd_core::physics::boundary::{BoundaryCondition, WallType};
use cfd_core::physics::fluid::non_newtonian::CarreauYasuda;
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::Vector2;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Setup Domain
    let length = 0.04; // 4 cm
    let height = 0.01; // 1 cm
    let nx = 80;
    let ny = 20;
    let grid = StructuredGrid2D::new(nx, ny, 0.0, length, 0.0, height)?;

    // 2. Setup Fluid (Blood)
    let blood = CarreauYasuda::<f64>::blood();
    let rho = blood.density;

    // 3. Setup Solver
    let config = SimplecPimpleConfig {
        tolerance: 1e-3,
        max_inner_iterations: 5,
        n_outer_correctors: 2,
        n_inner_correctors: 1,
        alpha_u: 0.5, // Relax velocity
        alpha_p: 0.3, // Relax pressure
        ..Default::default()
    };
    let mut solver = SimplecPimpleSolver::new(grid.clone(), config)?;

    // Boundary Conditions
    // Inlet: Parabolic profile u(y) = 4 * U_max * y/H * (1 - y/H)
    // We implement this by setting fixed values at inlet nodes manually in the loop or using BC
    // Standard BC set for solver:
    solver.set_boundary("west".to_string(), BoundaryCondition::VelocityInlet { velocity: nalgebra::Vector3::new(0.5, 0.0, 0.0) }); // Average/Max
    solver.set_boundary("east".to_string(), BoundaryCondition::PressureOutlet { pressure: 0.0 });
    solver.set_boundary("north".to_string(), BoundaryCondition::Wall { wall_type: WallType::NoSlip });
    solver.set_boundary("south".to_string(), BoundaryCondition::Wall { wall_type: WallType::NoSlip });

    // 4. Setup Immersed Boundary (Venturi Throat)
    let mut ibm = ImmersedBoundaryMethod::new((nx, ny), (length, height));
    // Set IBM parameters
    let ibm_config = ImmersedBoundaryConfig {
        regularization: 1e-4, // Stiffer penalty
        ..Default::default()
    };
    ibm.set_config(ibm_config);

    // Define throat geometry: Cosine shape
    // y_wall(x) = H/2 - (H/2 - w/2) * (1 + cos(2*pi*(x - x_c)/L_c))/2
    // Throat at center, length 2cm
    let x_center = length / 2.0;
    let throat_length = 0.02;
    let throat_width = 0.005; // 0.5 cm
    let num_points = 100;

    for i in 0..num_points {
        let x = x_center - throat_length/2.0 + (throat_length * i as f64 / (num_points as f64 - 1.0));

        // Cosine constriction
        let phase = 2.0 * std::f64::consts::PI * (x - x_center) / throat_length;
        // Depth of constriction from wall
        let depth = (height - throat_width) / 4.0 * (1.0 + phase.cos());

        // Top wall point
        let p_top = Vector2::new(x, height - depth);
        ibm.add_boundary_point(BoundaryPoint {
            position: p_top,
            desired_velocity: Vector2::zeros(),
            force: Vector2::zeros(),
            segment_length: length / nx as f64, // Approx
        });

        // Bottom wall point
        let p_bot = Vector2::new(x, depth);
        ibm.add_boundary_point(BoundaryPoint {
            position: p_bot,
            desired_velocity: Vector2::zeros(),
            force: Vector2::zeros(),
            segment_length: length / nx as f64,
        });
    }

    // 5. Initialize Fields
    let mut fields = SimulationFields::new(nx, ny);
    // Initialize viscosity to constant high-shear limit
    let mu_inf = blood.viscosity_inf;
    for j in 0..ny {
        for i in 0..nx {
            fields.viscosity.set(i, j, mu_inf);
            fields.density.set(i, j, rho);
        }
    }

    // 6. Time Loop
    let dt = 0.0001;
    let end_time = 0.005; // Simulate 50 steps
    let mut time = 0.0;

    println!("Starting simulation...");
    println!("Grid: {}x{}, Length: {}m, Height: {}m", nx, ny, length, height);
    println!("Fluid: {}, Viscosity range: {:.4}-{:.4} Pa.s", blood.name(), blood.viscosity_inf, blood.viscosity_zero);

    while time < end_time {
        // A. Update Viscosity (Non-Newtonian)
        for j in 1..ny-1 {
            for i in 1..nx-1 {
                // Calculate strain rate magnitude: sqrt(2 * S_ij * S_ij)
                // S_11 = du/dx, S_22 = dv/dy, S_12 = 0.5(du/dy + dv/dx)
                let du_dx = (fields.u.at(i+1, j) - fields.u.at(i-1, j)) / (2.0 * grid.dx);
                let dv_dy = (fields.v.at(i, j+1) - fields.v.at(i, j-1)) / (2.0 * grid.dy);
                let du_dy = (fields.u.at(i, j+1) - fields.u.at(i, j-1)) / (2.0 * grid.dy);
                let dv_dx = (fields.v.at(i+1, j) - fields.v.at(i-1, j)) / (2.0 * grid.dx);

                let s11 = du_dx;
                let s22 = dv_dy;
                let s12 = 0.5 * (du_dy + dv_dx);

                let strain_rate = (2.0 * (s11*s11 + s22*s22 + 2.0*s12*s12)).sqrt();

                let mu = blood.viscosity_at_shear(strain_rate, 310.15, 101325.0).unwrap_or(mu_inf);
                fields.viscosity.set(i, j, mu);
            }
        }

        // B. Immersed Boundary Method
        // 1. Interpolate velocity to boundary points
        // Use separate vector for velocity field from fields
        let mut velocity_matrix = nalgebra::DMatrix::zeros(nx * ny * 2, 1);
        for j in 0..ny {
            for i in 0..nx {
                velocity_matrix[2 * (j * nx + i)] = fields.u.at(i, j);
                velocity_matrix[2 * (j * nx + i) + 1] = fields.v.at(i, j);
            }
        }

        let boundary_velocities = ibm.interpolate_velocities(&velocity_matrix)?;

        // 2. Compute forces
        ibm.update_forces(&boundary_velocities)?;

        // 3. Spread forces to grid
        let mut force_matrix = nalgebra::DMatrix::zeros(nx * ny * 2, 1);
        ibm.spread_forces(&mut force_matrix)?;

        // 4. Update fields.force
        for j in 0..ny {
            for i in 0..nx {
                let fx = force_matrix[2 * (j * nx + i)];
                let fy = force_matrix[2 * (j * nx + i) + 1];
                fields.force_u.set(i, j, fx);
                fields.force_v.set(i, j, fy);
            }
        }

        // C. Solve Flow
        let residual = solver.solve_time_step(&mut fields, dt, 0.0, rho)?;

        time += dt;
        if (time / dt) as usize % 10 == 0 {
            let u_max = fields.max_velocity_magnitude();
            println!("Time: {:.3}s, Residual: {:.2e}, Max U: {:.2} m/s", time, residual, u_max);
        }
    }

    println!("Simulation complete.");

    // Sample velocity at throat center
    let center_u = fields.u.at(nx/2, ny/2);
    println!("Velocity at throat center: {:.2} m/s", center_u);

    Ok(())
}
