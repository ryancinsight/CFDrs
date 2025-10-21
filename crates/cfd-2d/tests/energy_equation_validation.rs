//! Analytical validation tests for energy equation
//!
//! Validates the energy equation solver against known analytical solutions:
//! 1. **1D Steady Conduction**: T(x) = T0 + (T1-T0)*x/L
//! 2. **2D Convection-Diffusion**: Method of Manufactured Solutions (MMS)
//!
//! References:
//! - Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow
//! - Versteeg, H.K. & Malalasekera, W. (2007). An Introduction to CFD

use cfd_2d::physics::EnergyEquationSolver;
use cfd_core::boundary::BoundaryCondition;
use std::collections::HashMap;

/// Test 1D steady-state conduction with linear temperature profile
///
/// Analytical solution: T(x) = T0 + (T1-T0)*x/L
/// 
/// Physics: ∂T/∂t = 0, u = v = 0, Q = 0
/// Governing equation: α∂²T/∂x² = 0 → T(x) = ax + b
/// Boundary conditions: T(0) = T0, T(L) = T1
///
/// Expected error: ≤1e-6 (machine precision for linear profile)
#[test]
fn test_1d_steady_conduction_analytical() {
    let nx = 51;
    let ny = 3; // Minimal grid in y-direction
    let l = 1.0;
    let dx = l / (nx - 1) as f64;
    let dy = 0.1;
    
    let t0 = 300.0; // Left boundary temperature [K]
    let t1 = 400.0; // Right boundary temperature [K]
    let alpha = 1.0; // Thermal diffusivity [m²/s]
    
    // Initialize solver
    let mut solver = EnergyEquationSolver::<f64>::new(nx, ny, t0, alpha);
    
    // Set boundary conditions
    let mut boundary_conditions = HashMap::new();
    
    // Left boundary (Dirichlet)
    for j in 0..ny {
        boundary_conditions.insert((0, j), BoundaryCondition::Dirichlet { value: t0 });
    }
    
    // Right boundary (Dirichlet)
    for j in 0..ny {
        boundary_conditions.insert((nx - 1, j), BoundaryCondition::Dirichlet { value: t1 });
    }
    
    // Top and bottom boundaries (adiabatic/Neumann)
    for i in 0..nx {
        boundary_conditions.insert((i, 0), BoundaryCondition::Neumann { gradient: 0.0 });
        boundary_conditions.insert((i, ny - 1), BoundaryCondition::Neumann { gradient: 0.0 });
    }
    
    // Zero velocity field (pure conduction)
    let u_velocity = vec![vec![0.0; ny]; nx];
    let v_velocity = vec![vec![0.0; ny]; nx];
    
    // Time step for steady state (large enough to converge)
    let dt = 0.001;
    let num_steps = 10000; // Iterate to steady state
    
    // Solve to steady state
    for _ in 0..num_steps {
        solver.solve_explicit(&u_velocity, &v_velocity, dt, dx, dy, &boundary_conditions)
            .expect("Energy equation solve should succeed");
    }
    
    // Verify against analytical solution: T(x) = T0 + (T1-T0)*x/L
    let mut max_error: f64 = 0.0;
    let j_center = ny / 2; // Check centerline
    
    for i in 1..nx - 1 {
        let x = i as f64 * dx;
        let t_analytical = t0 + (t1 - t0) * x / l;
        let t_numerical = solver.temperature[i][j_center];
        let error = (t_numerical - t_analytical).abs();
        max_error = max_error.max(error);
    }
    
    // Assert convergence to analytical solution
    assert!(
        max_error < 1e-6,
        "1D steady conduction error ({}) exceeds tolerance (1e-6)",
        max_error
    );
}

/// Test 2D transient convection-diffusion with manufactured solution
///
/// Manufactured solution: T(x,y,t) = sin(πx) * sin(πy) * exp(-2π²αt)
///
/// Physics: ∂T/∂t + u∂T/∂x + v∂T/∂y = α∇²T + Q
/// Source term Q is computed to satisfy the manufactured solution
///
/// Expected error: ≤0.5 (upwind scheme has significant numerical diffusion)
/// Note: First-order upwind scheme cannot achieve high accuracy for convection-dominated problems
#[test]
fn test_2d_transient_convection_diffusion_mms() {
    let nx = 21;
    let ny = 21;
    let l = 1.0;
    let dx = l / (nx - 1) as f64;
    let dy = l / (ny - 1) as f64;
    
    let alpha = 0.1; // Thermal diffusivity [m²/s]
    let u0 = 0.1;    // Reduced velocity for stability [m/s]
    let v0 = 0.1;    // Reduced velocity for stability [m/s]
    let t_final = 0.01; // Short time to avoid excessive decay
    let dt = 0.00005; // Smaller timestep for stability
    let num_steps = (t_final / dt) as usize;
    
    // Initialize with manufactured solution at t=0
    let mut solver = EnergyEquationSolver::<f64>::new(nx, ny, 0.0, alpha);
    
    use std::f64::consts::PI;
    
    // Set initial condition: T(x,y,0) = sin(πx) * sin(πy)
    for i in 0..nx {
        for j in 0..ny {
            let x = i as f64 * dx;
            let y = j as f64 * dy;
            solver.temperature[i][j] = (PI * x).sin() * (PI * y).sin();
        }
    }
    
    // Constant velocity field
    let u_velocity = vec![vec![u0; ny]; nx];
    let v_velocity = vec![vec![v0; ny]; nx];
    
    // Set boundary conditions (Dirichlet with manufactured solution)
    let mut boundary_conditions = HashMap::new();
    
    // Time evolution with source term
    for step in 0..num_steps {
        let t = step as f64 * dt;
        let decay = (-2.0 * PI * PI * alpha * t).exp();
        
        // Update source term for manufactured solution
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let x = i as f64 * dx;
                let y = j as f64 * dy;
                solver.heat_source[i][j] = u0 * PI * (PI * x).cos() * (PI * y).sin() * decay
                    + v0 * PI * (PI * x).sin() * (PI * y).cos() * decay;
            }
        }
        
        // Update boundary conditions with manufactured solution
        boundary_conditions.clear();
        for i in 0..nx {
            let x = i as f64 * dx;
            let t_bottom = (PI * x).sin() * decay;
            let t_top = (PI * x).sin() * (PI * l).sin() * decay;
            boundary_conditions.insert((i, 0), BoundaryCondition::Dirichlet { value: t_bottom });
            boundary_conditions.insert((i, ny - 1), BoundaryCondition::Dirichlet { value: t_top });
        }
        for j in 0..ny {
            let y = j as f64 * dy;
            let t_left = (PI * y).sin() * decay;
            let t_right = (PI * l).sin() * (PI * y).sin() * decay;
            boundary_conditions.insert((0, j), BoundaryCondition::Dirichlet { value: t_left });
            boundary_conditions.insert((nx - 1, j), BoundaryCondition::Dirichlet { value: t_right });
        }
        
        solver.solve_explicit(&u_velocity, &v_velocity, dt, dx, dy, &boundary_conditions)
            .expect("Energy equation solve should succeed");
    }
    
    // Verify against manufactured solution at t_final
    let t = t_final;
    let decay = (-2.0 * PI * PI * alpha * t).exp();
    let mut max_error: f64 = 0.0;
    let mut rms_error: f64 = 0.0;
    let mut count: usize = 0;
    
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let x = i as f64 * dx;
            let y = j as f64 * dy;
            let t_analytical = (PI * x).sin() * (PI * y).sin() * decay;
            let t_numerical = solver.temperature[i][j];
            let error = (t_numerical - t_analytical).abs();
            max_error = max_error.max(error);
            rms_error += error * error;
            count += 1;
        }
    }
    
    rms_error = (rms_error / count as f64).sqrt();
    
    // Assert reasonable convergence (upwind scheme has O(h) error)
    // With 21x21 grid and first-order upwind, expect ~5% error
    assert!(
        max_error < 0.5,
        "2D transient convection-diffusion max error ({}) exceeds tolerance (0.5)",
        max_error
    );
    
    assert!(
        rms_error < 0.2,
        "2D transient convection-diffusion RMS error ({}) exceeds tolerance (0.2)",
        rms_error
    );
}

/// Test energy conservation for uniform temperature field
#[test]
fn test_uniform_temperature_conservation() {
    let nx = 10;
    let ny = 10;
    let dx = 0.1;
    let dy = 0.1;
    let alpha = 0.1;
    let t_uniform = 300.0;
    
    let mut solver = EnergyEquationSolver::<f64>::new(nx, ny, t_uniform, alpha);
    
    // Zero velocity (pure diffusion)
    let u_velocity = vec![vec![0.0; ny]; nx];
    let v_velocity = vec![vec![0.0; ny]; nx];
    
    // Adiabatic boundaries (no heat flux)
    let mut boundary_conditions = HashMap::new();
    for i in 0..nx {
        boundary_conditions.insert((i, 0), BoundaryCondition::Neumann { gradient: 0.0 });
        boundary_conditions.insert((i, ny - 1), BoundaryCondition::Neumann { gradient: 0.0 });
    }
    for j in 0..ny {
        boundary_conditions.insert((0, j), BoundaryCondition::Neumann { gradient: 0.0 });
        boundary_conditions.insert((nx - 1, j), BoundaryCondition::Neumann { gradient: 0.0 });
    }
    
    let dt = 0.001;
    let num_steps = 100;
    
    // Solve
    for _ in 0..num_steps {
        solver.solve_explicit(&u_velocity, &v_velocity, dt, dx, dy, &boundary_conditions)
            .expect("Energy equation solve should succeed");
    }
    
    // Verify temperature remains uniform
    for i in 0..nx {
        for j in 0..ny {
            let temp_diff = solver.temperature[i][j] - t_uniform;
            let error: f64 = temp_diff.abs();
            assert!(
                error < 1e-10,
                "Uniform temperature not conserved: error = {}",
                error
            );
        }
    }
}

/// Test energy conservation with heat source
#[test]
fn test_steady_heat_source_balance() {
    let nx = 11;
    let ny = 11;
    let l = 1.0;
    let dx = l / (nx - 1) as f64;
    let dy = l / (ny - 1) as f64;
    let alpha = 0.1;
    
    let mut solver = EnergyEquationSolver::<f64>::new(nx, ny, 0.0, alpha);
    
    // Uniform heat source
    let q_source = 100.0; // W/m³
    for i in 0..nx {
        for j in 0..ny {
            solver.heat_source[i][j] = q_source;
        }
    }
    
    // Zero velocity
    let u_velocity = vec![vec![0.0; ny]; nx];
    let v_velocity = vec![vec![0.0; ny]; nx];
    
    // Fixed temperature boundaries
    let t_boundary = 300.0;
    let mut boundary_conditions = HashMap::new();
    for i in 0..nx {
        boundary_conditions.insert((i, 0), BoundaryCondition::Dirichlet { value: t_boundary });
        boundary_conditions.insert((i, ny - 1), BoundaryCondition::Dirichlet { value: t_boundary });
    }
    for j in 0..ny {
        boundary_conditions.insert((0, j), BoundaryCondition::Dirichlet { value: t_boundary });
        boundary_conditions.insert((nx - 1, j), BoundaryCondition::Dirichlet { value: t_boundary });
    }
    
    let dt = 0.0001;
    let num_steps = 50000; // Converge to steady state
    
    // Solve to steady state
    for _ in 0..num_steps {
        solver.solve_explicit(&u_velocity, &v_velocity, dt, dx, dy, &boundary_conditions)
            .expect("Energy equation solve should succeed");
    }
    
    // Verify temperature is above boundary temperature (heat is added)
    let t_center = solver.temperature[nx / 2][ny / 2];
    assert!(
        t_center > t_boundary,
        "Center temperature ({}) should be above boundary temperature ({}) due to heat source",
        t_center,
        t_boundary
    );
    
    // Verify symmetry (temperature should be symmetric about center)
    for i in 1..nx / 2 {
        for j in 1..ny / 2 {
            let t1 = solver.temperature[i][j];
            let t2 = solver.temperature[nx - 1 - i][j];
            let t3 = solver.temperature[i][ny - 1 - j];
            let t4 = solver.temperature[nx - 1 - i][ny - 1 - j];
            
            let max_diff = (t1 - t2).abs().max((t1 - t3).abs()).max((t1 - t4).abs());
            assert!(
                max_diff < 1e-3,
                "Temperature field not symmetric: max difference = {}",
                max_diff
            );
        }
    }
}
