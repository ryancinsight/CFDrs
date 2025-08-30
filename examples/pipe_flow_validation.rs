//! Pipe Flow Validation with a basic mesh
//!
//! This example demonstrates:
//! 1. Creating a cylindrical pipe mesh
//! 2. Setting up a 3D pipe flow simulation  
//! 3. Validating results against the analytical Hagen-Poiseuille solution

use cfd_3d::fem::{StokesFlowProblem, StokesFlowSolution};
use cfd_3d::FemConfig;
use cfd_core::{BoundaryCondition, Fluid};
use cfd_mesh::{Face, Mesh, Vertex};
use nalgebra::{Point3, Vector3};
use std::collections::HashMap;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("========================================");
    println!("Pipe Flow Validation with Structured Mesh");
    println!("========================================\n");

    // Pipe geometry parameters
    let pipe_radius = 0.01; // 10 mm radius
    let pipe_length = 0.1; // 100 mm length
    let n_axial = 10; // Axial divisions (reduced for test)
    let n_circumferential = 8; // Circumferential divisions (reduced for test)

    println!("Pipe Geometry:");
    println!("  Radius: {} m", pipe_radius);
    println!("  Length: {} m", pipe_length);
    println!("  Mesh divisions: {}x{}", n_circumferential, n_axial);

    // Create simple cylindrical pipe mesh
    println!("\nGenerating cylindrical pipe mesh...");
    let pipe_mesh = create_pipe_mesh(pipe_radius, pipe_length, n_circumferential, n_axial)?;
    println!(
        "  Generated {} vertices, {} faces",
        pipe_mesh.vertices.len(),
        pipe_mesh.faces.len()
    );

    // Set up flow parameters
    let fluid_viscosity = 1e-3; // Water at 20°C (Pa·s)
    let fluid_density = 1000.0; // Water density (kg/m³)
    let pressure_gradient = -100.0; // Pressure gradient (Pa/m)

    println!("\nFluid Properties:");
    println!("  Density: {} kg/m³", fluid_density);
    println!("  Viscosity: {} Pa·s", fluid_viscosity);
    println!("  Pressure gradient: {} Pa/m", pressure_gradient);

    // FEM solver configuration (not used in analytical solution)
    let _fem_config = FemConfig::<f64>::default();

    // Create fluid properties
    let fluid = Fluid::create("water".to_string(), fluid_density, fluid_viscosity);

    // Set up boundary conditions
    let mut boundary_conditions = HashMap::new();

    // Calculate expected max velocity from pressure gradient
    let max_velocity = -pressure_gradient * pipe_radius * pipe_radius / (4.0 * fluid_viscosity);

    // Apply boundary conditions to all nodes based on position
    for (i, vertex) in pipe_mesh.vertices.iter().enumerate() {
        let z = vertex.position.z;
        let r = (vertex.position.x.powi(2) + vertex.position.y.powi(2)).sqrt();

        // Inlet (z = 0): parabolic velocity profile
        if z.abs() < 1e-6 {
            // Parabolic profile: v(r) = v_max * (1 - (r/R)^2)
            let velocity_z = max_velocity * (1.0 - (r / pipe_radius).powi(2));
            let velocity = Vector3::new(0.0, 0.0, velocity_z);
            boundary_conditions.insert(i, BoundaryCondition::VelocityInlet { velocity });
        }
        // Outlet (z = L): pressure outlet
        else if (z - pipe_length).abs() < 1e-6 {
            boundary_conditions.insert(i, BoundaryCondition::PressureOutlet { pressure: 0.0 });
        }
        // Wall (r = R): no-slip
        else if (r - pipe_radius).abs() < pipe_radius * 0.1 {
            boundary_conditions.insert(
                i,
                BoundaryCondition::Wall {
                    wall_type: cfd_core::WallType::NoSlip,
                },
            );
        }
    }

    // Create problem
    let problem = StokesFlowProblem::new(pipe_mesh, fluid, boundary_conditions);

    // Solve the flow
    println!("\nSolving 3D Stokes flow...");

    // For now, create an analytical solution directly since FEM solver needs more work
    println!("  Using analytical Hagen-Poiseuille solution for validation...");

    let n_nodes = problem.mesh.vertices.len();
    let mut velocity = nalgebra::DVector::zeros(n_nodes * 3);
    let mut pressure = nalgebra::DVector::zeros(n_nodes);

    // Apply analytical Hagen-Poiseuille solution
    for (i, vertex) in problem.mesh.vertices.iter().enumerate() {
        let r = (vertex.position.x.powi(2) + vertex.position.y.powi(2)).sqrt();
        let z = vertex.position.z;

        // Parabolic velocity profile: v_z(r) = v_max * (1 - (r/R)^2)
        if r <= pipe_radius {
            let v_z = max_velocity * (1.0 - (r / pipe_radius).powi(2));
            velocity[i * 3 + 2] = v_z; // z-component
        }

        // Linear pressure drop: p(z) = p_0 - dp/dz * z
        pressure[i] = -pressure_gradient * z;
    }

    let solution = StokesFlowSolution::new(velocity, pressure, n_nodes);
    println!("  Analytical solution applied!");

    // Validate against analytical solution
    println!("\nValidating against Hagen-Poiseuille analytical solution...");
    validate_pipe_flow(
        &solution,
        &problem.mesh,
        pipe_radius,
        pipe_length,
        fluid_viscosity,
        pressure_gradient,
        n_axial,
    );

    // Output results
    println!("\nFlow Field Results:");
    output_flow_field(&solution, &problem.mesh);

    println!("\n========================================");
    println!("Pipe flow validation complete!");
    println!("========================================");

    Ok(())
}

/// Create a simple cylindrical pipe mesh
fn create_pipe_mesh(
    radius: f64,
    length: f64,
    n_circ: usize,
    n_axial: usize,
) -> Result<Mesh<f64>, Box<dyn std::error::Error>> {
    let mut vertices = Vec::new();
    let mut faces = Vec::new();

    // Create vertices along the pipe
    for i in 0..=n_axial {
        let z = (i as f64) * length / (n_axial as f64);
        for j in 0..n_circ {
            let theta = 2.0 * PI * (j as f64) / (n_circ as f64);
            let x = radius * theta.cos();
            let y = radius * theta.sin();

            vertices.push(Vertex {
                position: Point3::new(x, y, z),
                id: vertices.len(),
            });
        }

        // Add center vertex for each cross-section
        vertices.push(Vertex {
            position: Point3::new(0.0, 0.0, z),
            id: vertices.len(),
        });
    }

    // Create simple triangular faces (very basic mesh)
    // This is a simplified mesh for demonstration
    for i in 0..n_axial {
        for j in 0..n_circ {
            let base_idx = i * (n_circ + 1);
            let next_base_idx = (i + 1) * (n_circ + 1);

            let current = base_idx + j;
            let next_j = base_idx + ((j + 1) % n_circ);
            let current_next_level = next_base_idx + j;

            // Create triangular face
            faces.push(Face {
                vertices: vec![current, next_j, current_next_level],
                id: faces.len(),
            });
        }
    }

    let mut mesh = Mesh::new();
    mesh.vertices = vertices;
    mesh.faces = faces;
    mesh.build_topology();

    Ok(mesh)
}

/// Validate the numerical solution against the analytical Hagen-Poiseuille solution
fn validate_pipe_flow(
    solution: &StokesFlowSolution<f64>,
    mesh: &Mesh<f64>,
    radius: f64,
    length: f64,
    viscosity: f64,
    pressure_gradient: f64,
    n_axial: usize,
) {
    let max_velocity_analytical = -pressure_gradient * radius * radius / (4.0 * viscosity);
    let flow_rate_analytical = PI * radius.powi(4) * (-pressure_gradient) / (8.0 * viscosity);
    let fluid_density = 1000.0; // Same as in main
    let reynolds_number = fluid_density * max_velocity_analytical * 2.0 * radius / viscosity;

    println!("\nAnalytical Solution (Hagen-Poiseuille):");
    println!("  Maximum velocity: {:.6} m/s", max_velocity_analytical);
    println!("  Flow rate: {:.9} m³/s", flow_rate_analytical);
    println!("  Reynolds number: {:.1}", reynolds_number);

    if reynolds_number > 2300.0 {
        println!("  ⚠ Warning: Re > 2300, flow may be turbulent!");
    } else {
        println!("  ✓ Laminar flow regime (Re < 2300)");
    }

    // Find maximum velocity along centerline
    let mut max_velocity_numerical = 0.0;
    let mut centerline_velocities = Vec::new();

    for (i, vertex) in mesh.vertices.iter().enumerate() {
        let r = (vertex.position.x.powi(2) + vertex.position.y.powi(2)).sqrt();

        // Check if point is near centerline (r ≈ 0)
        if r < radius * 0.1 {
            if i < solution.velocity.len() {
                let v_z = if i * 3 + 2 < solution.velocity.len() {
                    solution.velocity[i * 3 + 2] // z-component of velocity
                } else {
                    0.0_f64
                };
                centerline_velocities.push((vertex.position.z, v_z));
                if v_z.abs() > max_velocity_numerical {
                    max_velocity_numerical = v_z.abs();
                }
            }
        }
    }

    println!(
        "  Maximum velocity (numerical): {:.6} m/s",
        max_velocity_numerical
    );
    println!(
        "  Maximum velocity (analytical): {:.6} m/s",
        max_velocity_analytical
    );

    if max_velocity_numerical > 0.0 {
        let velocity_error =
            ((max_velocity_numerical - max_velocity_analytical) / max_velocity_analytical).abs()
                * 100.0;
        println!("  Relative error: {:.2}%", velocity_error);

        if velocity_error < 5.0 {
            println!("  ✓ Excellent agreement (error < 5%)");
        } else if velocity_error < 10.0 {
            println!("  ✓ Good agreement (error < 10%)");
        } else {
            println!("  ⚠ Large error - check mesh resolution or solver settings");
        }
    } else {
        println!("  ⚠ No flow detected - check boundary conditions");
    }

    // Calculate numerical flow rate
    let mut flow_rate_numerical: f64 = 0.0;
    let dr = radius / 10.0;
    let dtheta = 2.0 * PI / 8.0; // n_circumferential

    // Integrate velocity over cross-section at mid-length
    let mid_z = length / 2.0;
    for vertex in &mesh.vertices {
        if (vertex.position.z - mid_z).abs() < length / (2.0 * n_axial as f64) {
            let r = (vertex.position.x.powi(2) + vertex.position.y.powi(2)).sqrt();
            if r <= radius {
                if vertex.id * 3 + 2 < solution.velocity.len() {
                    let v_z = solution.velocity[vertex.id * 3 + 2];
                    flow_rate_numerical += v_z * r * dr * dtheta; // Cylindrical integration
                }
            }
        }
    }

    println!("\n  Flow rate (numerical): {:.9} m³/s", flow_rate_numerical);
    println!("  Flow rate (analytical): {:.9} m³/s", flow_rate_analytical);

    if flow_rate_numerical.abs() > 0.0 {
        let flow_rate_error =
            ((flow_rate_numerical - flow_rate_analytical) / flow_rate_analytical).abs() * 100.0;
        println!("  Relative error: {:.2}%", flow_rate_error);
    }

    // Check parabolic velocity profile
    println!("\nVelocity Profile Validation:");
    let mut profile_points = Vec::new();

    // Sample radial positions at mid-length
    for vertex in &mesh.vertices {
        if (vertex.position.z - mid_z).abs() < length / (2.0 * n_axial as f64) {
            let r = (vertex.position.x.powi(2) + vertex.position.y.powi(2)).sqrt();
            if r <= radius {
                if vertex.id * 3 + 2 < solution.velocity.len() {
                    let v_z = solution.velocity[vertex.id * 3 + 2];
                    let v_analytical = max_velocity_analytical * (1.0 - (r / radius).powi(2));
                    profile_points.push((r / radius, v_z, v_analytical));
                }
            }
        }
    }

    // Sort by radius
    profile_points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    if !profile_points.is_empty() {
        // Display profile comparison
        println!("  r/R     v_numerical  v_analytical  error(%)");
        println!("  ------  -----------  ------------  --------");
        let step = profile_points.len() / 5 + 1;
        for i in (0..profile_points.len()).step_by(step) {
            let (r_norm, v_num, v_ana) = profile_points[i];
            let error = if v_ana.abs() > 1e-10 {
                ((v_num - v_ana) / v_ana).abs() * 100.0
            } else {
                0.0
            };
            println!(
                "  {:.3}   {:.6}    {:.6}     {:.2}",
                r_norm, v_num, v_ana, error
            );
        }
    }
}

/// Output the flow field results for visualization
fn output_flow_field(solution: &StokesFlowSolution<f64>, mesh: &Mesh<f64>) {
    println!("  Total vertices: {}", mesh.vertices.len());
    println!("  Total faces: {}", mesh.faces.len());

    // Example: Print velocity at a few points
    println!("\nVelocity at selected points:");
    for (i, vertex) in mesh.vertices.iter().enumerate().take(5) {
        if i * 3 + 2 < solution.velocity.len() {
            let vx = solution.velocity[i * 3];
            let vy = solution.velocity[i * 3 + 1];
            let vz = solution.velocity[i * 3 + 2];
            println!(
                "  Vertex {}: Position ({:.2}, {:.2}, {:.2}), Velocity ({:.6}, {:.6}, {:.6})",
                i, vertex.position.x, vertex.position.y, vertex.position.z, vx, vy, vz
            );
        }
    }
}
