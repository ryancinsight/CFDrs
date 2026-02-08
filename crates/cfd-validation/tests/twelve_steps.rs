use cfd_2d::grid::{StructuredGrid2D, traits::Grid2D};
use cfd_2d::fields::{Field2D, SimulationFields};
use approx::assert_relative_eq;
use cfd_math::sparse::SparseMatrixBuilder;
use cfd_math::linear_solver::{ConjugateGradient, IterativeSolverConfig, preconditioners::IncompleteLU};
use nalgebra::{DVector, Vector3};
use cfd_2d::simplec_pimple::{SimplecPimpleSolver, config::{SimplecPimpleConfig, AlgorithmType}};
use cfd_core::physics::boundary::BoundaryCondition;

/// Step 1: 1D Linear Convection
/// Equation: du/dt + c*du/dx = 0
/// Method: 1st-order Upwind FVM
#[test]
fn test_step_1_linear_convection() {
    type T = f64; // Use f64 for validation
    
    // 1. Setup Grid (1D strip)
    let nx = 50;
    let ny = 1; // 1D equivalent
    let length = 2.0;
    let dx = length / nx as f64;
    let c = 1.0; // Wave speed
    let dt = 0.025; // CFL = 0.625 (stable for upwind)
    let t_end = 0.5;
    
    let grid = StructuredGrid2D::new(
        nx, ny,
        0.0, length,
        0.0, 0.1, // Small y-height
    ).expect("Grid creation failed");

    // 2. Initialize Fields
    let mut u = Field2D::new(nx, ny, 0.0);
    let mut u_new = Field2D::new(nx, ny, 0.0);

    // Initial Condition: Square wave from 0.5 to 1.0
    for i in 0..nx {
        let x = grid.cell_center(i, 0).unwrap().x;
        if x >= 0.5 && x <= 1.0 {
            u.set(i, 0, 2.0);
        } else {
            u.set(i, 0, 1.0);
        }
    }

    // 3. Time Stepping Loop (Explicit First-Order Upwind)
    let mut t = 0.0;
    while t < t_end {
        // Apply Boundary Conditions (Periodic)
        // For i=0, 'west' is i=nx-1
        // For i=nx-1, 'east' is i=0
        
        for i in 0..nx {
            let u_curr = u.at(i, 0);
            
            // Upwind neighbor (west)
            let u_west = if i == 0 {
                u.at(nx - 1, 0) // Periodic
            } else {
                u.at(i - 1, 0)
            };

            // Discretization: u_new = u - c * dt/dx * (u - u_west)
            // Flux out (east face) - Flux in (west face)
            // F_e = c * u_curr
            // F_w = c * u_west (upwind)
            // du/dt = - (F_e - F_w) / dx
             
            let val = u_curr - c * (dt / dx) * (u_curr - u_west);
            u_new.set(i, 0, val);
        }

        // Update time and field
        u = u_new.clone();
        t += dt;
    }

    // 4. Verification against Analytical Solution
    // Analytical: Wave moves distance c*t = 1.0 * 0.5 = 0.5
    // Original pulse: [0.5, 1.0] -> New pulse: [1.0, 1.5]
    
    // Check center of pulse
    let mut max_val = 0.0;
    let mut max_x = 0.0;
    for i in 0..nx {
        let val = u.at(i, 0);
        if val > max_val {
            max_val = val;
            max_x = grid.cell_center(i, 0).unwrap().x;
        }
    }

    // Numerical diffusion will smear the square wave, but the "center of mass" 
    // or the location of the peak (if it were a gaussian) should translate.
    // For a square wave with upwind, the leading edge smears forward and trailing edge counts.
    // Let's check mass conservation as a strict metric first.
    
    let initial_mass: f64 = (0.5 * 2.0) + (1.5 * 1.0); // (width 0.5 * height 2) + (rest 1.5 * bg 1) = 1.0 + 1.5 = 2.5
    // Actually sum of discrete cells
    let mut final_mass = 0.0;
    for i in 0..nx {
        final_mass += u.at(i, 0) * dx;
    }
    
    // Discrete Initial Mass explicitly
    let mut init_mass_discrete = 0.0;
    for i in 0..nx {
        let x = grid.cell_center(i, 0).unwrap().x;
        let val = if x >= 0.5 && x <= 1.0 { 2.0 } else { 1.0 };
        init_mass_discrete += val * dx;
    }

    assert_relative_eq!(final_mass, init_mass_discrete, epsilon = 1e-6);
    println!("Step 1 (Linear Convection): Mass Conserved. Init: {}, Final: {}, Diff: {}", init_mass_discrete, final_mass, (final_mass - init_mass_discrete).abs());
}

/// Step 2: 1D Nonlinear Convection (Burgers' Inviscid)
/// Equation: du/dt + u*du/dx = 0
/// Method: Conservative Upwind FVM
#[test]
fn test_step_2_nonlinear_convection() {
    type T = f64;
    
    // 1. Setup Grid (1D strip)
    let nx = 50;
    let ny = 1;
    let length = 2.0;
    let dx = length / nx as f64;
    let dt = 0.01; // Smaller dt for stability with u=2 (CFL ~ 0.5 using max(u))
    let t_end = 0.5;
    
    let grid = StructuredGrid2D::new(nx, ny, 0.0, length, 0.0, 0.1).expect("Grid failed");

    // 2. Initialize Fields
    let mut u = Field2D::new(nx, ny, 0.0);
    let mut u_new = Field2D::new(nx, ny, 0.0);

    // Initial Condition: Square wave u=2 inside [0.5, 1.0], u=1 elsewhere
    for i in 0..nx {
        let x = grid.cell_center(i, 0).unwrap().x;
        if x >= 0.5 && x <= 1.0 {
            u.set(i, 0, 2.0);
        } else {
            u.set(i, 0, 1.0);
        }
    }

    // 3. Time Stepping Loop
    let mut t = 0.0;
    while t < t_end {
        for i in 0..nx {
            // Upwind neighbor (Periodic)
            let u_curr = u.at(i, 0);
            let u_west = if i == 0 { u.at(nx - 1, 0) } else { u.at(i - 1, 0) };

            // Conservative Flux: F = u^2 / 2
            // Upwind Flux (assuming u > 0 always): F_face = F(u_upwind)
            // East Face of cell i: F_{i+1/2}. Upwind is cell i.
            // West Face of cell i: F_{i-1/2}. Upwind is cell i-1 (west).
            
            let f_e = 0.5 * u_curr * u_curr;
            let f_w = 0.5 * u_west * u_west;

            // Update: du/dt = - (F_e - F_w) / dx
            let val = u_curr - (dt / dx) * (f_e - f_w);
            u_new.set(i, 0, val);
        }
        u = u_new.clone();
        t += dt;
    }

    // 4. Verification
    // Check Mass Conservation
    let mut final_mass = 0.0;
    for i in 0..nx {
        final_mass += u.at(i, 0) * dx;
    }
    
    // Initial mass calculation
    let mut init_mass_discrete = 0.0;
    for i in 0..nx {
        let x = grid.cell_center(i, 0).unwrap().x;
        let val = if x >= 0.5 && x <= 1.0 { 2.0 } else { 1.0 };
        init_mass_discrete += val * dx;
    }

    assert_relative_eq!(final_mass, init_mass_discrete, epsilon = 1e-6);
    println!("Step 2 (Nonlinear Convection): Mass Conserved. Init: {}, Final: {}, Diff: {}", init_mass_discrete, final_mass, (final_mass - init_mass_discrete).abs());
}

/// Step 3: 1D Diffusion
/// Equation: du/dt = nu * d^2u/dx^2
/// Method: 1st-order explicit time, Central Difference space
#[test]
fn test_step_3_diffusion() {
    type T = f64;
    
    // 1. Setup Grid
    let nx = 100;
    let length = 4.0;
    let dx = length / nx as f64;
    let nu = 0.1;
    let sigma = 0.25;
    let dt = 0.001; // Stability: nu*dt/dx^2 <= 0.5. 0.1*0.001/(0.04^2) = 0.0001/0.0016 = 0.0625 < 0.5 OK
    let t_end = 0.5;
    
    let grid = StructuredGrid2D::new(nx, 1, -2.0, 2.0, 0.0, 0.1).expect("Grid failed");

    // Analytical Solution function
    let t0 = (sigma * sigma) / (2.0 * nu);
    let analytical = |x: f64, t: f64| -> f64 {
        let effective_t = t + t0;
        (1.0 / (4.0 * std::f64::consts::PI * nu * effective_t).sqrt()) * (-x*x / (4.0 * nu * effective_t)).exp()
    };

    // 2. Initialize Fields (at t=0)
    let mut u = Field2D::new(nx, 1, 0.0);
    let mut u_new = Field2D::new(nx, 1, 0.0);
    
    for i in 0..nx {
        let x = grid.cell_center(i, 0).unwrap().x;
        u.set(i, 0, analytical(x, 0.0));
    }

    // 3. Time Stepping Loop
    let mut t = 0.0;
    while t < t_end {
        for i in 0..nx {
            let u_curr = u.at(i, 0);
            
            // Neighbors with Fixed BC u=0 at walls (Dirichlet for infinite domain approx)
            // Or use analytical BC
            let x_w = grid.cell_center(0, 0).unwrap().x - dx;
            let val_w_bc = analytical(x_w, t); // Analytical BC at west ghost
            let x_e = grid.cell_center(nx-1, 0).unwrap().x + dx;
            let val_e_bc = analytical(x_e, t); // Analytical BC at east ghost

            let u_west = if i == 0 { val_w_bc } else { u.at(i - 1, 0) };
            let u_east = if i == nx - 1 { val_e_bc } else { u.at(i + 1, 0) };

            // Diffusion: nu * (u_E - 2u_P + u_W) / dx^2
            let diffusion = nu * (u_east - 2.0 * u_curr + u_west) / (dx * dx);
            
            // Update
            let val = u_curr + dt * diffusion;
            u_new.set(i, 0, val);
        }
        u = u_new.clone();
        t += dt;
    }

    // 4. Verification
    let mut l2_error = 0.0;
    for i in 0..nx {
        let x = grid.cell_center(i, 0).unwrap().x;
        let u_num = u.at(i, 0);
        let u_ana = analytical(x, t);
        l2_error += (u_num - u_ana).powi(2);
    }
    l2_error = (l2_error / nx as f64).sqrt();

    println!("Step 3 (Diffusion): L2 Error = {:.6e}", l2_error);
    assert!(l2_error < 1e-3, "Diffusion accuracy failed");
}

/// Step 4: 1D Burgers' Equation
/// Equation: du/dt + u*du/dx = nu * d^2u/dx^2
/// Method: Conservative Upwind (Convection) + Central Diff (Diffusion)
#[test]
fn test_step_4_burgers() {
    type T = f64;
    
    // 1. Setup Grid (Periodic)
    let nx = 100;
    let length = 2.0 * std::f64::consts::PI;
    let dx = length / nx as f64;
    let nu = 0.07;
    let dt = 0.005; // Stability: Combined condition
    let t_end = 2.0;
    
    let grid = StructuredGrid2D::new(nx, 1, 0.0, length, 0.0, 0.1).expect("Grid failed");

    // 2. Initialize Fields (Sawtooth)
    let mut u = Field2D::new(nx, 1, 0.0);
    let mut u_new = Field2D::new(nx, 1, 0.0);
    
    // Analytical validation for Burgers is complex, but we can check mass conservation with Periodic BC.
    // Init: u = - sin(x) for shock formation verification? 
    // Or u = 4.0 - 2.0*nu * dphi/dx / phi, where phi = exp ...
    // Let's use Sawtooth wave: u = 4 - x inside [0, 2pi].
    // Better: Sinusoidal. u(x,0) = sin(x).
    // It will steepen into a N-wave shock and then decay.
    
    for i in 0..nx {
        let x = grid.cell_center(i, 0).unwrap().x;
        u.set(i, 0, (x).sin());
    }

    // 3. Time Stepping Loop
    let mut t = 0.0;
    while t < t_end {
        for i in 0..nx {
            let u_curr = u.at(i, 0);
            
            // Periodic BC
            let u_west = if i == 0 { u.at(nx-1, 0) } else { u.at(i-1, 0) };
            let u_east = if i == nx-1 { u.at(0, 0) } else { u.at(i+1, 0) };

            // Fluxes
            // Convection (Upwind for u>0, Downwind for u<0? No, just Upwind or average).
            // Conservative Flux F_c = u^2/2.
            // Face East (i+1/2): F_c_e = 0.5 * u_curr^2 (if u_curr > 0). If u_east < 0?
            // Simple upwind for Burgers u u_x usually handles signs.
            // Let's use simple non-conservative for internal check, OR strict conservative flux.
            // F_{i+1/2} = 0.5 * u_{upwind}^2.
            // But 'upwind' depends on velocity direction at face.
            // Approx face velocity u_{i+1/2} = 0.5*(u_i + u_{i+1}).
            // If u_{i+1/2} > 0, upwind is i. Else i+1.
            
            // Let's implement: F_{i+1/2} = 0.5 * u_L^2 if u_avg > 0 else 0.5 * u_R^2.
            
            let u_face_e = 0.5 * (u_curr + u_east);
            let f_c_e = if u_face_e >= 0.0 {
                0.5 * u_curr * u_curr
            } else {
                0.5 * u_east * u_east
            };
            
            let u_west_cell = if i == 0 { u.at(nx-1, 0) } else { u.at(i-1, 0) }; // Re-get for consistency
            // Wait, I have u_west. Need u_west_west for face west?
            // Face West (i-1/2): Same logic.
            // u_{i-1/2} = 0.5*(u_west + u_curr).
            let u_face_w = 0.5 * (u_west + u_curr);
            let f_c_w = if u_face_w >= 0.0 {
                0.5 * u_west * u_west
            } else {
                0.5 * u_curr * u_curr
            };
            
            // Diffusion Flux (Central)
            // F_d_e = nu * (u_east - u_curr) / dx
            // F_d_w = nu * (u_curr - u_west) / dx
            let f_d_e = nu * (u_east - u_curr) / dx;
            let f_d_w = nu * (u_curr - u_west) / dx;
            
            // Total Flux F = F_c - F_d
            let f_total_e = f_c_e - f_d_e;
            let f_total_w = f_c_w - f_d_w;
            
            // Update: du/dt = - (F_e - F_w) / dx
            let val = u_curr - (dt / dx) * (f_total_e - f_total_w);
            u_new.set(i, 0, val);
        }
        u = u_new.clone();
        t += dt;
    }

    // 4. Verification
    // Check Mass Conservation (periodic)
    let mut final_mass = 0.0;
    for i in 0..nx {
        final_mass += u.at(i, 0) * dx;
    }
    
    // Initial mass
    let mut init_mass = 0.0;
    for i in 0..nx {
        let x = grid.cell_center(i, 0).unwrap().x;
        init_mass += (x).sin() * dx;
    }
    
    // Note: Integral of sin(x) over 0..2pi is 0.
    // Discrete sum might be close to 0 but float error.
    assert!(init_mass.abs() < 1e-10); // Check valid init
    
    assert_relative_eq!(final_mass, init_mass, epsilon = 1e-6);
    println!("Step 4 (Burgers): Mass Conserved. Init: {}, Final: {}, Diff: {}", init_mass, final_mass, (final_mass - init_mass).abs());
}

/// Step 5: 2D Linear Convection
/// Equation: du/dt + cx*du/dx + cy*du/dy = 0
/// Method: 1st-order Upwind FVM on 2D grid
#[test]
fn test_step_5_2d_linear_convection() {
    type T = f64;
    
    // 1. Setup Grid (Square Domain)
    let nx = 50;
    let ny = 50;
    let length = 2.0;
    let dx = length / nx as f64;
    let dy = length / ny as f64;
    let cx = 1.0;
    let cy = 1.0;
    let dt = 0.01; // CFL = 0.5 < 1.0
    let t_end = 0.5;
    
    let grid = StructuredGrid2D::new(nx, ny, 0.0, length, 0.0, length).expect("Grid failed");

    // 2. Initialize Fields (Square Wave in center)
    let mut u = Field2D::new(nx, ny, 0.0);
    let mut u_new = Field2D::new(nx, ny, 0.0);

    // Initial Condition: Square block u=2 in [0.5, 1.0] x [0.5, 1.0]
    for i in 0..nx {
        for j in 0..ny {
            let center = grid.cell_center(i, j).unwrap();
            if center.x >= 0.5 && center.x <= 1.0 && center.y >= 0.5 && center.y <= 1.0 {
                u.set(i, j, 2.0);
            } else {
                u.set(i, j, 1.0);
            }
        }
    }

    // 3. Time Stepping Loop
    let mut t = 0.0;
    while t < t_end {
        for i in 0..nx {
            for j in 0..ny {
                let u_curr = u.at(i, j);
                
                // Upwind Neighbors (Periodic BC for simplicity)
                let u_west = if i == 0 { u.at(nx-1, j) } else { u.at(i-1, j) };
                let u_south = if j == 0 { u.at(i, ny-1) } else { u.at(i, j-1) };

                // Discretization: du/dt = - [ cx*(u - u_w)/dx + cy*(u - u_s)/dy ]
                let adv_x = cx * (u_curr - u_west) / dx;
                let adv_y = cy * (u_curr - u_south) / dy;
                
                let val = u_curr - dt * (adv_x + adv_y);
                u_new.set(i, j, val);
            }
        }
        u = u_new.clone();
        t += dt;
    }

    // 4. Verification
    // Check Mass Conservation
    let mut final_mass = 0.0;
    for i in 0..nx {
        for j in 0..ny {
            final_mass += u.at(i, j) * dx * dy;
        }
    }
    
    // Initial Mass
    let mut init_mass_discrete = 0.0;
    for i in 0..nx {
        for j in 0..ny {
            let center = grid.cell_center(i, j).unwrap();
            let val = if center.x >= 0.5 && center.x <= 1.0 && center.y >= 0.5 && center.y <= 1.0 { 2.0 } else { 1.0 };
            init_mass_discrete += val * dx * dy;
        }
    }
    
    assert_relative_eq!(final_mass, init_mass_discrete, epsilon = 1e-6);
    println!("Step 5 (2D Linear Convection): Mass Conserved. Init: {}, Final: {}, Diff: {}", init_mass_discrete, final_mass, (final_mass - init_mass_discrete).abs());
}

/// Step 6: 2D Nonlinear Convection
/// Equation: du/dt + u*du/dx + u*du/dy = 0 (Scalar Burgers-like)
/// Method: Conservative Upwind FVM on 2D grid
#[test]
fn test_step_6_nonlinear_convection() {
    type T = f64;
    
    // 1. Setup Grid (Square Domain)
    let nx = 50;
    let ny = 50;
    let length = 2.0;
    let dx = length / nx as f64;
    let dy = length / ny as f64;
    let dt = 0.01;
    let t_end = 0.5;
    
    let grid = StructuredGrid2D::new(nx, ny, 0.0, length, 0.0, length).expect("Grid failed");

    // 2. Initialize Fields (Square Wave)
    let mut u = Field2D::new(nx, ny, 0.0);
    let mut u_new = Field2D::new(nx, ny, 0.0);

    // Initial Condition
    for i in 0..nx {
        for j in 0..ny {
            let center = grid.cell_center(i, j).unwrap();
            if center.x >= 0.5 && center.x <= 1.0 && center.y >= 0.5 && center.y <= 1.0 {
                u.set(i, j, 2.0);
            } else {
                u.set(i, j, 1.0);
            }
        }
    }

    // 3. Time Stepping Loop
    let mut t = 0.0;
    while t < t_end {
        for i in 0..nx {
            for j in 0..ny {
                let u_curr = u.at(i, j);
                
                // Periodic BC
                let u_west = if i == 0 { u.at(nx-1, j) } else { u.at(i-1, j) };
                let u_south = if j == 0 { u.at(i, ny-1) } else { u.at(i, j-1) };

                // Conservative Fluxes F = u^2/2 (in x and y)
                let f_e_x = 0.5 * u_curr * u_curr;
                let f_w_x = 0.5 * u_west * u_west;
                
                let f_n_y = 0.5 * u_curr * u_curr;
                let f_s_y = 0.5 * u_south * u_south;

                // Update: du/dt = - [ (F_e - F_w)/dx + (F_n - F_s)/dy ]
                let val = u_curr - dt * ((f_e_x - f_w_x)/dx + (f_n_y - f_s_y)/dy);
                u_new.set(i, j, val);
            }
        }
        u = u_new.clone();
        t += dt;
    }

    // 4. Verification
    let mut final_mass = 0.0;
    for i in 0..nx {
        for j in 0..ny {
            final_mass += u.at(i, j) * dx * dy;
        }
    }
    
    let mut init_mass_discrete = 0.0;
    for i in 0..nx {
        for j in 0..ny {
            let center = grid.cell_center(i, j).unwrap();
            let val = if center.x >= 0.5 && center.x <= 1.0 && center.y >= 0.5 && center.y <= 1.0 { 2.0 } else { 1.0 };
            init_mass_discrete += val * dx * dy;
        }
    }
    
    assert_relative_eq!(final_mass, init_mass_discrete, epsilon = 1e-6);
    println!("Step 6 (2D Nonlinear Convection): Mass Conserved. Init: {}, Final: {}", init_mass_discrete, final_mass);
}

/// Step 7: 2D Diffusion
/// Equation: du/dt = nu * (d^2u/dx^2 + d^2u/dy^2)
/// Method: 1st-order explicit time, Central Difference space
#[test]
fn test_step_7_diffusion() {
    type T = f64;
    
    // 1. Setup Grid
    let nx = 31;
    let ny = 31;
    let length = 2.0;
    let dx = length / (nx as f64 - 1.0);
    let dy = length / (ny as f64 - 1.0);
    // Grid from 0 to 2
    let nu = 0.05;
    let sigma = 0.25;
    let dt = 0.001;
    let t_end = 0.25;
    
    let grid = StructuredGrid2D::new(nx, ny, 0.0, length, 0.0, length).expect("Grid failed");

    // Analytical Solution (2D Gaussian)
    // u(x,y,t) = exp( - ((x-xc)^2 + (y-yc)^2) / (4 nu (t+t0)) ) / (4 pi nu (t+t0))
    // Note: To match peak height logic better:
    let t0 = (sigma * sigma) / (4.0 * nu); // t0 such that at t=0, width is sigma
    let xc = 1.0;
    let yc = 1.0;
    
    let analytical = |x: f64, y: f64, t: f64| -> f64 {
        let effective_t = t + t0;
        // 2D Heat kernel
        let r2 = (x - xc).powi(2) + (y - yc).powi(2);
        (1.0 / (4.0 * std::f64::consts::PI * nu * effective_t)) * (-r2 / (4.0 * nu * effective_t)).exp()
    };

    // 2. Initialize Fields
    let mut u = Field2D::new(nx, ny, 0.0);
    let mut u_new = Field2D::new(nx, ny, 0.0);
    
    for i in 0..nx {
        for j in 0..ny {
            let center = grid.cell_center(i, j).unwrap();
            u.set(i, j, analytical(center.x, center.y, 0.0));
        }
    }

    // 3. Time Stepping Loop
    let mut t = 0.0;
    while t < t_end {
        for i in 0..nx {
            for j in 0..ny {
                let u_curr = u.at(i, j);
                let center = grid.cell_center(i, j).unwrap();
                
                // Boundaries fixed to analytical solution
                if i == 0 || i == nx-1 || j == 0 || j == ny-1 {
                    u_new.set(i, j, analytical(center.x, center.y, t + dt));
                    continue;
                }
                
                let u_west = u.at(i-1, j);
                let u_east = u.at(i+1, j);
                let u_south = u.at(i, j-1);
                let u_north = u.at(i, j+1);

                let d2u_dx2 = (u_east - 2.0*u_curr + u_west) / (dx*dx);
                let d2u_dy2 = (u_north - 2.0*u_curr + u_south) / (dy*dy);
                
                let val = u_curr + dt * nu * (d2u_dx2 + d2u_dy2);
                u_new.set(i, j, val);
            }
        }
        u = u_new.clone();
        t += dt;
    }

    // 4. Verification
    let mut l2_error = 0.0;
    for i in 1..nx-1 {
        for j in 1..ny-1 {
            let center = grid.cell_center(i, j).unwrap();
            let u_num = u.at(i, j);
            let u_ana = analytical(center.x, center.y, t);
            l2_error += (u_num - u_ana).powi(2);
        }
    }
    l2_error = (l2_error / ((nx-2)*(ny-2)) as f64).sqrt();

    println!("Step 7 (2D Diffusion): L2 Error = {:.6e}", l2_error);
    assert!(l2_error < 5e-3, "2D Diffusion accuracy failed");
}

/// Step 8: 2D Burgers' Equation
/// Equation: du/dt + u*du/dx + v*du/dy = nu * (d^2u/dx^2 + d^2u/dy^2)
/// Simplified: u_t + u u_x + u u_y = nu Lap(u) (Scalar)
#[test]
fn test_step_8_burgers_2d() {
    type T = f64;
    
    // 1. Setup
    let nx = 50;
    let ny = 50;
    let length = 2.0;
    let dx = length / nx as f64;
    let dy = length / ny as f64;
    let nu = 0.01;
    let dt = 0.002;
    let t_end = 0.5;
    
    let grid = StructuredGrid2D::new(nx, ny, 0.0, length, 0.0, length).expect("Grid failed");

    // 2. Init (Gaussian bump to diffuse and advect)
    let mut u = Field2D::new(nx, ny, 0.0);
    let mut u_new = Field2D::new(nx, ny, 0.0);
    
    for i in 0..nx {
        for j in 0..ny {
            let center = grid.cell_center(i, j).unwrap();
            let r2 = (center.x - 1.0).powi(2) + (center.y - 1.0).powi(2);
            u.set(i, j, (-r2 / 0.1).exp());
        }
    }
    
    let init_mass = u.data().iter().sum::<f64>() * dx * dy;

    // 3. Time Stepping
    let mut t = 0.0;
    while t < t_end {
        for i in 0..nx {
             for j in 0..ny {
                let u_curr = u.at(i, j);
                
                // Periodic BC
                let u_w = if i == 0 { u.at(nx-1, j) } else { u.at(i-1, j) };
                let u_e = if i == nx-1 { u.at(0, j) } else { u.at(i+1, j) };
                let u_s = if j == 0 { u.at(i, ny-1) } else { u.at(i, j-1) };
                let u_n = if j == ny-1 { u.at(i, 0) } else { u.at(i, j+1) }; // Periodic y
                
                // Advection (Scalar Burgers: u u_x + u u_y)
                // Conservative fluxes approx
                let f_e_x = 0.5 * u_curr * u_curr;
                let f_w_x = 0.5 * u_w * u_w;
                let f_n_y = 0.5 * u_curr * u_curr;
                let f_s_y = 0.5 * u_s * u_s;
                
                let adv = (f_e_x - f_w_x)/dx + (f_n_y - f_s_y)/dy;
                
                // Diffusion
                let diff = nu * ((u_e - 2.0*u_curr + u_w)/(dx*dx) + (u_n - 2.0*u_curr + u_s)/(dy*dy));
                
                let val = u_curr - dt * adv + dt * diff;
                u_new.set(i, j, val);
            }
        }
        u = u_new.clone();
        t += dt;
    }
    
    // 4. Verification
    let final_mass = u.data().iter().sum::<f64>() * dx * dy;
    assert_relative_eq!(final_mass, init_mass, epsilon = 2e-4);
    println!("Step 8 (2D Burgers): Mass Conserved. Init: {:.6}, Final: {:.6}, Diff: {:.6e}", init_mass, final_mass, (final_mass - init_mass).abs());
}


#[test]
fn test_step_11_lid_driven_cavity() {
    // Lid Driven Cavity (Re=100)
    let nx = 32;
    let ny = 32;
    let width = 1.0;
    let height = 1.0;
    
    let grid = StructuredGrid2D::new(nx, ny, 0.0, width, 0.0, height).expect("Grid creation failed");
    
    // Config for Re = 100
    // L = 1, U = 1, nu = 0.01 -> Re = 100
    let rho = 1.0;
    let nu = 0.01;
    let dt = 0.005; 
    
    let config = SimplecPimpleConfig {
        algorithm: AlgorithmType::Simplec,
        tolerance: 1e-4, 
        max_inner_iterations: 5,
        alpha_u: 0.7,
        alpha_p: 0.3,
        n_outer_correctors: 1,
        n_inner_correctors: 2,
        ..Default::default()
    };
    
    let mut solver = SimplecPimpleSolver::new(grid.clone(), config).expect("Solver creation failed");
    
    // Set Boundary Conditions
    solver.set_boundary("north".to_string(), BoundaryCondition::Wall { 
        wall_type: cfd_core::physics::boundary::WallType::Moving { velocity: Vector3::new(1.0, 0.0, 0.0) } 
    });
    solver.set_boundary("south".to_string(), BoundaryCondition::Wall { 
        wall_type: cfd_core::physics::boundary::WallType::NoSlip 
    });
    solver.set_boundary("east".to_string(), BoundaryCondition::Wall { 
        wall_type: cfd_core::physics::boundary::WallType::NoSlip 
    });
    solver.set_boundary("west".to_string(), BoundaryCondition::Wall { 
        wall_type: cfd_core::physics::boundary::WallType::NoSlip 
    });
    
    let mut fields = SimulationFields::new(nx, ny);
    
    // Run until convergence or max steps
    let max_steps = 1000;
    let mut converged = false;
    
    for _step in 0..max_steps {
        let (_effective_dt, residual) = solver.solve_adaptive(&mut fields, dt, nu, rho, 1, 1e-5).expect("Solve failed");
        if residual < 1e-4 {
            converged = true;
            break;
        }
    }
    
    assert!(converged, "Solver did not converge for Lid-Driven Cavity");
    
    // Verification: Center of vortex check and Top Lid Velocity
    let u_center = fields.u.at(nx/2, ny/2); 
    assert!(u_center < 0.0, "Center u-velocity should be negative (vortex recirculation). Got {}", u_center);
    
    let u_top = fields.u.at(nx/2, ny-1);
    assert!(u_top > 0.5, "Top cell velocity should be positive. Got {}", u_top);
}

#[test]
fn test_step_12_channel_flow() {
    // Channel Flow (Poiseuille)
    let nx = 50;
    let ny = 10;
    let width = 5.0;  // Length
    let height = 1.0; // Height
    
    let grid = StructuredGrid2D::new(nx, ny, 0.0, width, 0.0, height).expect("Grid creation failed");
    
    let rho = 1.0;
    let nu = 0.1; // Re = 10 (Laminar)
    let dt = 0.01;
    
    let config = SimplecPimpleConfig {
        algorithm: AlgorithmType::Simplec,
        tolerance: 1e-4,
        ..Default::default()
    };
    
    let mut solver = SimplecPimpleSolver::new(grid.clone(), config).unwrap();
    
    // Inlet (West): Moving Wall trick for Plug Flow (u=1, v=0)
    solver.set_boundary("west".to_string(), BoundaryCondition::Wall { 
        wall_type: cfd_core::physics::boundary::WallType::Moving { velocity: Vector3::new(1.0, 0.0, 0.0) } 
    });
    
    // Outlet (East): Neumann (Gradient = 0)
    solver.set_boundary("east".to_string(), BoundaryCondition::Neumann { gradient: 0.0 });
    
    // Walls (North/South): No Slip
    solver.set_boundary("north".to_string(), BoundaryCondition::Wall { 
        wall_type: cfd_core::physics::boundary::WallType::NoSlip 
    });
    solver.set_boundary("south".to_string(), BoundaryCondition::Wall { 
        wall_type: cfd_core::physics::boundary::WallType::NoSlip 
    });
    
    let mut fields = SimulationFields::new(nx, ny);
    
    let mut converged = false;
    for _ in 0..1000 {
        let (_, res) = solver.solve_adaptive(&mut fields, dt, nu, rho, 1, 1e-5).unwrap();
        if res < 1e-4 { converged = true; break; }
    }
    
    assert!(converged, "Channel flow solver did not converge");
    
    // Check parabolic profile at outlet
    let u_outlet_center = fields.u.at(nx-1, ny/2);
    // 1.5 is theoretical max for Poiseuille
    assert!(u_outlet_center > 1.2, "Outlet centerline velocity should be > 1.2. Got {}", u_outlet_center);
    
    let u_outlet_wall = fields.u.at(nx-1, 0);
    assert!(u_outlet_wall < 0.1, "Outlet wall velocity should be near 0. Got {}", u_outlet_wall);
}

/// Helper: Solve Poisson Equation on 2D Grid with Dirichlet BC = 0 (or custom via rhs adjustment)
/// Solves: nabla^2 p = source
fn solve_poisson_2d(
    grid: &StructuredGrid2D<f64>,
    source: &Field2D<f64>, // RHS field
    boundary_val_func: impl Fn(f64, f64) -> f64, // Analytical Dirichlet BC
) -> Field2D<f64> {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx;
    let dy = grid.dy;
    
    // Build Matrix
    // Stencil: (u_{i+1} + u_{i-1} - 2u_i)/dx^2 + ...
    let mut matrix_builder = SparseMatrixBuilder::new(nx * ny, nx * ny);
    let mut rhs_vec = DVector::zeros(nx * ny);
    
    for j in 0..ny {
        for i in 0..nx {
            let idx = j * nx + i;
            let center = grid.cell_center(i, j).unwrap();
            
            // Check boundary
            if i == 0 || i == nx-1 || j == 0 || j == ny-1 {
                // Dirichlet BC
                matrix_builder.add_entry(idx, idx, 1.0).unwrap();
                rhs_vec[idx] = boundary_val_func(center.x, center.y);
            } else {
                // Interior
                let coeff_x = 1.0/(dx*dx);
                let coeff_y = 1.0/(dy*dy);
                let center_coeff = -2.0*coeff_x - 2.0*coeff_y;
                
                matrix_builder.add_entry(idx, idx, center_coeff).unwrap();
                matrix_builder.add_entry(idx, j*nx + (i-1), coeff_x).unwrap(); // West
                matrix_builder.add_entry(idx, j*nx + (i+1), coeff_x).unwrap(); // East
                matrix_builder.add_entry(idx, (j-1)*nx + i, coeff_y).unwrap(); // South
                matrix_builder.add_entry(idx, (j+1)*nx + i, coeff_y).unwrap(); // North
                
                // RHS
                rhs_vec[idx] = source.at(i, j);
            }
        }
    }
    
    // Solve
    let matrix = matrix_builder.build().unwrap();
    let solver = ConjugateGradient::new(IterativeSolverConfig {
        tolerance: 1e-10,
        max_iterations: 2000,
        ..Default::default()
    });
    
    // Use IncompleteLU preconditioner
    let preconditioner = IncompleteLU::new(&matrix).unwrap();
    
    let mut x_vec = DVector::zeros(nx * ny);
    solver.solve_preconditioned(&matrix, &rhs_vec, &preconditioner, &mut x_vec).unwrap();
    
    let mut p = Field2D::new(nx, ny, 0.0);
    for j in 0..ny {
        for i in 0..nx {
            p.set(i, j, x_vec[j*nx + i]);
        }
    }
    p
}

/// Step 9: Laplace Equation
/// nabla^2 p = 0
/// BCs: p = y on top, p = 0 on bottom, periodic or linear on sides?
/// Let's use analytical solution p(x,y) = x + y. nabla^2 = 0.
#[test]
fn test_step_9_laplace() {
    let nx = 30;
    let ny = 30;
    let grid = StructuredGrid2D::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();
    
    let analytical = |x: f64, y: f64| x * x - y * y; // p = x^2 - y^2 => p_xx = 2, p_yy = -2 => sum = 0.
    
    let mut source = Field2D::new(nx, ny, 0.0); // Source = 0
    let p_num = solve_poisson_2d(&grid, &source, analytical);
    
    // Check Error
    let mut l2_error = 0.0;
    for i in 0..nx {
        for j in 0..ny {
            let center = grid.cell_center(i, j).unwrap();
            let p_ana = analytical(center.x, center.y);
            l2_error += (p_num.at(i, j) - p_ana).powi(2);
        }
    }
    l2_error = (l2_error / ((nx*ny) as f64)).sqrt();
    
    println!("Step 9 (Laplace): L2 Error = {:.6e}", l2_error);
    assert!(l2_error < 1e-4);
}

/// Step 10: Poisson Equation
/// nabla^2 p = b
/// Analytical: p = sin(pi*x)*sin(pi*y).
/// nabla^2 p = -pi^2 sin sin - pi^2 sin sin = -2pi^2 sin(pi*x)sin(pi*y).
/// BCs: p=0 on boundaries.
#[test]
fn test_step_10_poisson() {
    let nx = 40;
    let ny = 40;
    let grid = StructuredGrid2D::new(nx, ny, 0.0, 2.0, 0.0, 1.0).unwrap();
    // Use domain [0,1]x[0,1] actually to match sin(pi*x) logic easily
    // Or just map analytical
    
    let pi = std::f64::consts::PI;
    let analytical_func = |x: f64, y: f64| (pi*x).sin() * (pi*y).sin();
    
    let mut source = Field2D::new(nx, ny, 0.0);
    // Source b = -2*pi^2 * sin(pi*x)*sin(pi*y)
    for i in 0..nx {
        for j in 0..ny {
            let center = grid.cell_center(i, j).unwrap();
            source.set(i, j, -2.0 * pi * pi * (pi*center.x).sin() * (pi*center.y).sin());
        }
    }
    
    // BCs are 0 for this function at x=0,1 y=0,1 (if domain 0..1)
    // Actually domain width 2.0 -> sin(2pi) = 0. Correct.
    // Wait, grid above is 0..2.0 x.
    // If x=2.0 -> sin(2pi) = 0.
    // If y=1.0 -> sin(pi) = 0.
    // So BCs are naturally 0.
    
    let p_num = solve_poisson_2d(&grid, &source, analytical_func);
    
    let mut l2_error = 0.0;
    for i in 0..nx {
        for j in 0..ny {
            let center = grid.cell_center(i, j).unwrap();
            let p_ana = analytical_func(center.x, center.y);
            l2_error += (p_num.at(i, j) - p_ana).powi(2);
        }
    }
    l2_error = (l2_error / ((nx*ny) as f64)).sqrt();
    
    println!("Step 10 (Poisson): L2 Error = {:.6e}", l2_error);
    assert!(l2_error < 5e-4);
}
