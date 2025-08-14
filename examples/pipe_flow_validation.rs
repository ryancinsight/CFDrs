//! Pipe Flow Validation with CSGrs-generated Mesh
//!
//! This example demonstrates:
//! 1. Creating a cylindrical pipe mesh using CSGrs
//! 2. Setting up a 3D pipe flow simulation
//! 3. Validating results against the analytical Hagen-Poiseuille solution

use cfd_mesh::{Mesh, Vertex, Face, MeshTopology, csg::CsgMeshAdapter};
use cfd_3d::{FemSolver, FemConfig, FluidProperties};
use nalgebra::{Point3, Vector3};
use std::f64::consts::PI;
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("========================================");
    println!("Pipe Flow Validation with CSGrs");
    println!("========================================\n");
    
    // Pipe geometry parameters
    let pipe_radius = 0.01; // 10 mm radius
    let pipe_length = 0.1;  // 100 mm length
    let n_axial = 20;       // Axial divisions
    let n_circumferential = 16; // Circumferential divisions
    
    println!("Pipe Geometry:");
    println!("  Radius: {} m", pipe_radius);
    println!("  Length: {} m", pipe_length);
    println!("  Mesh divisions: {}x{}", n_circumferential, n_axial);
    
    // Create cylindrical pipe mesh
    println!("\nGenerating cylindrical pipe mesh...");
    let pipe_mesh = create_pipe_mesh(pipe_radius, pipe_length, n_circumferential, n_axial)?;
    println!("  Generated {} vertices, {} faces", pipe_mesh.vertices.len(), pipe_mesh.faces.len());
    
    // Use CSGrs to create end caps
    let csg_adapter = CsgMeshAdapter::<f64>::new();
    
        // Create a simple closed pipe mesh using cylinder primitive from CSG
    println!("\nCreating complete pipe mesh using CSG cylinder primitive...");
    let complete_pipe = csg_adapter.create_cylinder(pipe_radius, pipe_length, n_circumferential)?;
    println!("  Complete pipe mesh: {} vertices, {} faces", 
             complete_pipe.vertices.len(), complete_pipe.faces.len());
    
    // Note: CSG boolean operations are not implemented, so we use basic primitives
    
    // Set up flow parameters
    let fluid_viscosity = 1e-3;  // Water at 20°C (Pa·s)
    let fluid_density = 1000.0;   // Water density (kg/m³)
    let pressure_gradient = -100.0; // Pressure gradient (Pa/m)
    
    println!("\nFlow Parameters:");
    println!("  Fluid viscosity: {} Pa·s", fluid_viscosity);
    println!("  Fluid density: {} kg/m³", fluid_density);
    println!("  Pressure gradient: {} Pa/m", pressure_gradient);
    
    // Calculate analytical solution (Hagen-Poiseuille)
    let max_velocity_analytical = -pressure_gradient * pipe_radius * pipe_radius / (4.0 * fluid_viscosity);
    let flow_rate_analytical = PI * pipe_radius.powi(4) * (-pressure_gradient) / (8.0 * fluid_viscosity);
    let reynolds_number = fluid_density * max_velocity_analytical * 2.0 * pipe_radius / fluid_viscosity;
    
    println!("\nAnalytical Solution (Hagen-Poiseuille):");
    println!("  Maximum velocity: {:.6} m/s", max_velocity_analytical);
    println!("  Flow rate: {:.9} m³/s", flow_rate_analytical);
    println!("  Reynolds number: {:.1}", reynolds_number);
    
    if reynolds_number > 2300.0 {
        println!("  ⚠ Warning: Re > 2300, flow may be turbulent!");
    } else {
        println!("  ✓ Laminar flow regime (Re < 2300)");
    }
    
    // Create FEM solver
    println!("\nSetting up FEM solver...");
    let base_config = cfd_core::SolverConfig::<f64>::builder()
        .tolerance(1e-6)
        .max_iterations(1000)
        .verbosity(1)
        .build();
    
    let fem_config = FemConfig {
        base: base_config,
        use_stabilization: true,
        tau: 0.1,
        dt: None, // Steady state
        reynolds: Some(reynolds_number),
    };
    
    // Create fluid properties
    let fluid_props = FluidProperties {
        density: fluid_density,
        viscosity: fluid_viscosity,
        body_force: None, // No body forces for horizontal pipe
    };
    
    let mut fem_solver = FemSolver::new(fem_config, complete_pipe.clone(), fluid_props);
    
    // Apply boundary conditions
    println!("\nApplying boundary conditions:");
    let mut velocity_bcs = HashMap::new();
    
    // Inlet: parabolic velocity profile (Poiseuille flow)
    // u_z(r) = u_max * (1 - (r/R)^2)
    let u_max = 2.0 * max_velocity_analytical; // Maximum velocity at centerline
    println!("  Inlet: Parabolic profile with u_max = {:.4} m/s", u_max);
    
    // Apply parabolic profile at inlet nodes
    for (i, vertex) in complete_pipe.vertices.iter().enumerate() {
        if vertex.position.z.abs() < 1e-6 { // Inlet at z=0
            let r = (vertex.position.x * vertex.position.x + 
                    vertex.position.y * vertex.position.y).sqrt();
            let u_z = u_max * (1.0 - (r/pipe_radius).powi(2));
            velocity_bcs.insert(i, Vector3::new(0.0, 0.0, u_z));
        }
    }
    
    // Walls: no-slip condition
    println!("  Walls: No-slip (u = v = w = 0)");
    for (i, vertex) in complete_pipe.vertices.iter().enumerate() {
        let r = (vertex.position.x * vertex.position.x + 
                vertex.position.y * vertex.position.y).sqrt();
        if (r - pipe_radius).abs() < 1e-6 { // On pipe wall
            velocity_bcs.insert(i, Vector3::new(0.0, 0.0, 0.0));
        }
    }
    
    // Outlet: stress-free (natural BC, no explicit constraint needed)
    println!("  Outlet: Stress-free (natural BC)");
    
    // Solve for steady-state flow
    println!("\nSolving steady-state Stokes flow...");
    fem_solver.solve_stokes(&velocity_bcs)?;
    
    // Get velocity solution
    let velocity = fem_solver.get_velocity_field();
    
    // Validate results
    println!("\nValidation Results:");
    
    // Find maximum velocity along centerline
    let mut max_velocity_numerical = 0.0;
    let mut centerline_velocities = Vec::new();
    
    for (i, vertex) in complete_pipe.vertices.iter().enumerate() {
        let r = (vertex.position.x.powi(2) + vertex.position.y.powi(2)).sqrt();
        
        // Check if point is near centerline (r ≈ 0)
        if r < pipe_radius * 0.1 {
            // Velocity is stored as [u0, v0, w0, u1, v1, w1, ...]
            let vel_idx = i * 3;
            if vel_idx + 2 < velocity.len() {
                let v_z = velocity[vel_idx + 2]; // z-component of velocity
                centerline_velocities.push((vertex.position.z, v_z));
                if v_z.abs() > max_velocity_numerical {
                    max_velocity_numerical = v_z.abs();
                }
            }
        }
    }
    
    println!("  Maximum velocity (numerical): {:.6} m/s", max_velocity_numerical);
    println!("  Maximum velocity (analytical): {:.6} m/s", max_velocity_analytical);
    
    if max_velocity_numerical > 0.0 {
        let velocity_error = ((max_velocity_numerical - max_velocity_analytical) / max_velocity_analytical).abs() * 100.0;
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
    let dr = pipe_radius / 10.0;
    let dtheta = 2.0 * PI / n_circumferential as f64;
    
    // Integrate velocity over cross-section at mid-length
    let mid_z = pipe_length / 2.0;
    for vertex in &complete_pipe.vertices {
        if (vertex.position.z - mid_z).abs() < pipe_length / (2.0 * n_axial as f64) {
            let r = (vertex.position.x.powi(2) + vertex.position.y.powi(2)).sqrt();
            if r <= pipe_radius {
                            let vel_idx = vertex.id * 3;
            if vel_idx + 2 < velocity.len() {
                let v_z = velocity[vel_idx + 2];
                flow_rate_numerical += v_z * r * dr * dtheta; // Cylindrical integration
            }
            }
        }
    }
    
    println!("\n  Flow rate (numerical): {:.9} m³/s", flow_rate_numerical);
    println!("  Flow rate (analytical): {:.9} m³/s", flow_rate_analytical);
    
    if flow_rate_numerical.abs() > 0.0 {
        let flow_rate_error = ((flow_rate_numerical - flow_rate_analytical) / flow_rate_analytical).abs() * 100.0;
        println!("  Relative error: {:.2}%", flow_rate_error);
    }
    
    // Check parabolic velocity profile
    println!("\nVelocity Profile Validation:");
    let mut profile_points = Vec::new();
    
    // Sample radial positions at mid-length
    for vertex in &complete_pipe.vertices {
        if (vertex.position.z - mid_z).abs() < pipe_length / (2.0 * n_axial as f64) {
            let r = (vertex.position.x.powi(2) + vertex.position.y.powi(2)).sqrt();
            if r <= pipe_radius {
                let vel_idx = vertex.id * 3;
                if vel_idx + 2 < velocity.len() {
                    let v_z = velocity[vel_idx + 2];
                    let v_analytical = max_velocity_analytical * (1.0 - (r / pipe_radius).powi(2));
                    profile_points.push((r / pipe_radius, v_z, v_analytical));
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
            println!("  {:.3}   {:.6}    {:.6}     {:.2}", r_norm, v_num, v_ana, error);
        }
    }
    
    println!("\n========================================");
    println!("Pipe Flow Validation Complete!");
    println!("========================================");
    
    Ok(())
}

/// Create a cylindrical pipe mesh
fn create_pipe_mesh(
    radius: f64,
    length: f64,
    n_circumferential: usize,
    n_axial: usize,
) -> Result<Mesh<f64>, Box<dyn std::error::Error>> {
    let mut mesh = Mesh::new();
    let mut vertex_id = 0;
    
    // Generate vertices in cylindrical coordinates
    for k in 0..=n_axial {
        let z = k as f64 * length / n_axial as f64;
        
        for j in 0..n_circumferential {
            let theta = j as f64 * 2.0 * PI / n_circumferential as f64;
            
            // Only outer surface for hollow pipe
            let r = radius;
            let x = r * theta.cos();
            let y = r * theta.sin();
            
            mesh.vertices.push(Vertex {
                id: vertex_id,
                position: Point3::new(x, y, z),
            });
            vertex_id += 1;
        }
    }
    
    // Generate faces (quadrilaterals split into triangles)
    let mut face_id = 0;
    for k in 0..n_axial {
        for j in 0..n_circumferential {
            let j_next = (j + 1) % n_circumferential;
            
            // Current ring
            let v0 = k * n_circumferential + j;
            let v1 = k * n_circumferential + j_next;
            
            // Next ring
            let v2 = (k + 1) * n_circumferential + j;
            let v3 = (k + 1) * n_circumferential + j_next;
            
            // First triangle
            mesh.faces.push(Face {
                id: face_id,
                vertices: vec![v0, v1, v2],
            });
            face_id += 1;
            
            // Second triangle
            mesh.faces.push(Face {
                id: face_id,
                vertices: vec![v1, v3, v2],
            });
            face_id += 1;
        }
    }
    
    // Update topology
    mesh.topology = MeshTopology {
        num_vertices: mesh.vertices.len(),
        num_edges: 0, // Not computed
        num_faces: mesh.faces.len(),
        num_cells: 0, // Surface mesh only
    };
    
    Ok(mesh)
}

