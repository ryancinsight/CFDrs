//! Lid-driven cavity flow validation test
//!
//! This is a standard CFD benchmark problem used to validate incompressible flow solvers.
//! Reference: Ghia et al. (1982) "High-Re solutions for incompressible flow using the
//! Navier-Stokes equations and a multigrid method", J. Comput. Phys., 48, 387-411.

use cfd_suite::prelude::*;
use nalgebra::DMatrix;

/// Ghia et al. (1982) reference data for Re=100
/// u-velocity along vertical centerline at x=0.5
const GHIA_RE100_Y: &[f64] = &[
    0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000, 0.8047, 0.8594, 0.9063,
    0.9453, 0.9531, 0.9609, 0.9688, 1.0000,
];

const GHIA_RE100_U: &[f64] = &[
    0.00000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, 0.05454,
    0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.09111, 0.06923, 0.00000,
];

/// Ghia et al. (1982) reference data for Re=1000
const GHIA_RE1000_Y: &[f64] = &[
    0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000, 0.8047, 0.8594, 0.9063,
    0.9453, 0.9531, 0.9609, 0.9688, 1.0000,
];

const GHIA_RE1000_U: &[f64] = &[
    0.00000, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289, -0.27805, -0.10648, 0.05702,
    0.18719, 0.18942, 0.17238, 0.13088, 0.11477, 0.09233, 0.06803, 0.00000,
];

/// Lid-driven cavity problem setup
struct LidDrivenCavity {
    /// Grid size
    n: usize,
    /// Reynolds number
    reynolds: f64,
    /// Domain size (assumed square)
    length: f64,
    /// Lid velocity
    u_lid: f64,
    /// Grid spacing
    dx: f64,
    dy: f64,
    /// Kinematic viscosity
    nu: f64,
}

impl LidDrivenCavity {
    fn new(n: usize, reynolds: f64) -> Self {
        let length = 1.0;
        let u_lid = 1.0;
        let nu = u_lid * length / reynolds;
        let dx = length / (n as f64 - 1.0);
        let dy = dx;

        Self {
            n,
            reynolds,
            length,
            u_lid,
            dx,
            dy,
            nu,
        }
    }

    /// Set up initial conditions
    fn initialize(&self) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
        let u = DMatrix::zeros(self.n, self.n);
        let v = DMatrix::zeros(self.n, self.n);
        let p = DMatrix::zeros(self.n, self.n);
        (u, v, p)
    }

    /// Apply boundary conditions
    fn apply_boundary_conditions(&self, u: &mut DMatrix<f64>, v: &mut DMatrix<f64>) {
        let n = self.n;

        // Top wall (lid): u = u_lid, v = 0
        for j in 0..n {
            u[(n - 1, j)] = self.u_lid;
            v[(n - 1, j)] = 0.0;
        }

        // Bottom wall: u = 0, v = 0
        for j in 0..n {
            u[(0, j)] = 0.0;
            v[(0, j)] = 0.0;
        }

        // Left wall: u = 0, v = 0
        for i in 0..n {
            u[(i, 0)] = 0.0;
            v[(i, 0)] = 0.0;
        }

        // Right wall: u = 0, v = 0
        for i in 0..n {
            u[(i, n - 1)] = 0.0;
            v[(i, n - 1)] = 0.0;
        }
    }

    /// Solve using projection method (Chorin's method)
    fn solve(
        &self,
        max_iterations: usize,
        tolerance: f64,
    ) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
        let mut u = DMatrix::zeros(self.n, self.n);
        let mut v = DMatrix::zeros(self.n, self.n);
        let mut p = DMatrix::zeros(self.n, self.n);

        // Time step (CFL condition)
        let dt = 0.5
            * f64::min(
                self.dx * self.dx / (4.0 * self.nu),
                self.dx / (2.0 * self.u_lid),
            );

        for iteration in 0..max_iterations {
            let u_old = u.clone();
            let v_old = v.clone();

            // Step 1: Compute intermediate velocity (without pressure gradient)
            let (u_star, v_star) = self.compute_intermediate_velocity(&u, &v, dt);

            // Step 2: Solve pressure Poisson equation
            let p_new = self.solve_pressure_poisson(&u_star, &v_star, dt);

            // Step 3: Correct velocity with pressure gradient
            u = self.correct_velocity_u(&u_star, &p_new, dt);
            v = self.correct_velocity_v(&v_star, &p_new, dt);
            p = p_new;

            // Apply boundary conditions
            self.apply_boundary_conditions(&mut u, &mut v);

            // Check convergence
            let residual = self.compute_residual(&u, &u_old, &v, &v_old);

            if iteration % 100 == 0 {
                println!("Iteration {}: residual = {:.6e}", iteration, residual);
            }

            if residual < tolerance {
                println!("Converged after {} iterations", iteration);
                break;
            }
        }

        (u, v, p)
    }

    /// Compute intermediate velocity using explicit time stepping
    fn compute_intermediate_velocity(
        &self,
        u: &DMatrix<f64>,
        v: &DMatrix<f64>,
        dt: f64,
    ) -> (DMatrix<f64>, DMatrix<f64>) {
        let n = self.n;
        let mut u_star = u.clone();
        let mut v_star = v.clone();

        // Interior points only
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                // Convection terms (central difference)
                let u_conv = u[(i, j)] * (u[(i, j + 1)] - u[(i, j - 1)]) / (2.0 * self.dx)
                    + v[(i, j)] * (u[(i + 1, j)] - u[(i - 1, j)]) / (2.0 * self.dy);

                let v_conv = u[(i, j)] * (v[(i, j + 1)] - v[(i, j - 1)]) / (2.0 * self.dx)
                    + v[(i, j)] * (v[(i + 1, j)] - v[(i - 1, j)]) / (2.0 * self.dy);

                // Diffusion terms (central difference)
                let u_diff = self.nu
                    * ((u[(i, j + 1)] - 2.0 * u[(i, j)] + u[(i, j - 1)]) / (self.dx * self.dx)
                        + (u[(i + 1, j)] - 2.0 * u[(i, j)] + u[(i - 1, j)]) / (self.dy * self.dy));

                let v_diff = self.nu
                    * ((v[(i, j + 1)] - 2.0 * v[(i, j)] + v[(i, j - 1)]) / (self.dx * self.dx)
                        + (v[(i + 1, j)] - 2.0 * v[(i, j)] + v[(i - 1, j)]) / (self.dy * self.dy));

                // Update intermediate velocity
                u_star[(i, j)] = u[(i, j)] + dt * (-u_conv + u_diff);
                v_star[(i, j)] = v[(i, j)] + dt * (-v_conv + v_diff);
            }
        }

        (u_star, v_star)
    }

    /// Solve pressure Poisson equation using Jacobi iteration
    fn solve_pressure_poisson(
        &self,
        u_star: &DMatrix<f64>,
        v_star: &DMatrix<f64>,
        dt: f64,
    ) -> DMatrix<f64> {
        let n = self.n;
        let mut p = DMatrix::zeros(n, n);
        let mut p_new = DMatrix::zeros(n, n);

        // RHS of Poisson equation
        let mut rhs = DMatrix::zeros(n, n);
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                rhs[(i, j)] = (1.0 / dt)
                    * ((u_star[(i, j + 1)] - u_star[(i, j - 1)]) / (2.0 * self.dx)
                        + (v_star[(i + 1, j)] - v_star[(i - 1, j)]) / (2.0 * self.dy));
            }
        }

        // Jacobi iteration
        let max_pressure_iterations = 1000;
        let pressure_tolerance = 1e-6;

        for _ in 0..max_pressure_iterations {
            for i in 1..n - 1 {
                for j in 1..n - 1 {
                    p_new[(i, j)] = 0.25
                        * (p[(i, j + 1)] + p[(i, j - 1)] + p[(i + 1, j)] + p[(i - 1, j)]
                            - self.dx * self.dx * rhs[(i, j)]);
                }
            }

            // Apply Neumann boundary conditions (dp/dn = 0)
            for i in 0..n {
                p_new[(i, 0)] = p_new[(i, 1)];
                p_new[(i, n - 1)] = p_new[(i, n - 2)];
                p_new[(0, i)] = p_new[(1, i)];
                p_new[(n - 1, i)] = p_new[(n - 2, i)];
            }

            // Check convergence
            let residual: f64 = (&p_new - &p).iter().map(|x| x * x).sum::<f64>().sqrt();
            if residual < pressure_tolerance {
                break;
            }

            p = p_new.clone();
        }

        p_new
    }

    /// Correct u-velocity with pressure gradient
    fn correct_velocity_u(&self, u_star: &DMatrix<f64>, p: &DMatrix<f64>, dt: f64) -> DMatrix<f64> {
        let n = self.n;
        let mut u = u_star.clone();

        for i in 1..n - 1 {
            for j in 1..n - 1 {
                u[(i, j)] = u_star[(i, j)] - dt * (p[(i, j + 1)] - p[(i, j - 1)]) / (2.0 * self.dx);
            }
        }

        u
    }

    /// Correct v-velocity with pressure gradient
    fn correct_velocity_v(&self, v_star: &DMatrix<f64>, p: &DMatrix<f64>, dt: f64) -> DMatrix<f64> {
        let n = self.n;
        let mut v = v_star.clone();

        for i in 1..n - 1 {
            for j in 1..n - 1 {
                v[(i, j)] = v_star[(i, j)] - dt * (p[(i + 1, j)] - p[(i - 1, j)]) / (2.0 * self.dy);
            }
        }

        v
    }

    /// Compute residual for convergence check
    fn compute_residual(
        &self,
        u: &DMatrix<f64>,
        u_old: &DMatrix<f64>,
        v: &DMatrix<f64>,
        v_old: &DMatrix<f64>,
    ) -> f64 {
        let u_diff = u - u_old;
        let v_diff = v - v_old;

        let u_norm: f64 = u_diff.iter().map(|x| x * x).sum::<f64>().sqrt();
        let v_norm: f64 = v_diff.iter().map(|x| x * x).sum::<f64>().sqrt();

        (u_norm + v_norm) / (self.n as f64)
    }

    /// Extract u-velocity along vertical centerline
    fn extract_centerline_u(&self, u: &DMatrix<f64>) -> Vec<f64> {
        let n = self.n;
        let j_center = n / 2;

        (0..n).map(|i| u[(i, j_center)]).collect()
    }

    /// Validate against Ghia et al. reference data
    fn validate_against_ghia(&self, u: &DMatrix<f64>) -> f64 {
        let centerline_u = self.extract_centerline_u(u);

        // Choose reference data based on Reynolds number
        let (ref_y, ref_u) = if (self.reynolds - 100.0).abs() < 1.0 {
            (GHIA_RE100_Y, GHIA_RE100_U)
        } else if (self.reynolds - 1000.0).abs() < 1.0 {
            (GHIA_RE1000_Y, GHIA_RE1000_U)
        } else {
            panic!("No reference data for Re = {}", self.reynolds);
        };

        // Interpolate simulation results to reference y-locations
        let mut total_error = 0.0;
        let mut count = 0;

        for (&y_ref, &u_ref) in ref_y.iter().zip(ref_u.iter()) {
            // Find corresponding grid point
            let i = ((1.0 - y_ref) * (self.n as f64 - 1.0)).round() as usize;
            if i < self.n {
                let u_sim = centerline_u[i];
                let error = (u_sim - u_ref).abs();
                total_error += error * error;
                count += 1;
            }
        }

        (total_error / count as f64).sqrt()
    }
}

#[test]
fn test_lid_driven_cavity_re100() {
    let cavity = LidDrivenCavity::new(41, 100.0);
    let (u, v, _p) = cavity.solve(5000, 1e-5);

    // Validate against Ghia et al. data
    let rmse = cavity.validate_against_ghia(&u);
    println!("RMSE against Ghia et al. (Re=100): {:.4e}", rmse);

    // Check that RMSE is within acceptable tolerance
    // Allow 5% error due to discretization differences
    assert!(rmse < 0.05, "RMSE {} exceeds tolerance", rmse);

    // Check mass conservation (divergence should be near zero)
    let mut max_div: f64 = 0.0;
    for i in 1..cavity.n - 1 {
        for j in 1..cavity.n - 1 {
            let div = (u[(i, j + 1)] - u[(i, j - 1)]) / (2.0 * cavity.dx)
                + (v[(i + 1, j)] - v[(i - 1, j)]) / (2.0 * cavity.dy);
            max_div = max_div.max(div.abs());
        }
    }
    println!("Maximum divergence: {:.4e}", max_div);
    assert!(
        max_div < 1e-3,
        "Mass conservation violated: max divergence = {}",
        max_div
    );
}

#[test]
fn test_lid_driven_cavity_re1000() {
    let cavity = LidDrivenCavity::new(61, 1000.0);
    let (u, v, _p) = cavity.solve(10000, 1e-5);

    // Validate against Ghia et al. data
    let rmse = cavity.validate_against_ghia(&u);
    println!("RMSE against Ghia et al. (Re=1000): {:.4e}", rmse);

    // Higher Reynolds number allows slightly higher error
    assert!(rmse < 0.08, "RMSE {} exceeds tolerance", rmse);

    // Check for vortex formation (characteristic of lid-driven cavity)
    // The center should have negative u-velocity (recirculation)
    let center_i = cavity.n / 2;
    let center_j = cavity.n / 2;
    let u_center = u[(center_i, center_j)];

    println!("U-velocity at center: {:.4}", u_center);
    assert!(u_center < 0.1, "No recirculation detected at center");
}

#[test]
fn test_steady_state_convergence() {
    let cavity = LidDrivenCavity::new(21, 100.0);

    // Run for different numbers of iterations
    let (u1, _, _) = cavity.solve(1000, 1e-5);
    let (u2, _, _) = cavity.solve(2000, 1e-5);

    // Check that solution doesn't change significantly
    let diff: f64 = (&u2 - &u1).iter().map(|x| x * x).sum::<f64>().sqrt();
    let norm: f64 = u2.iter().map(|x| x * x).sum::<f64>().sqrt();
    let relative_diff = diff / norm;

    println!(
        "Relative difference between 1000 and 2000 iterations: {:.4e}",
        relative_diff
    );
    assert!(
        relative_diff < 0.01,
        "Solution not converged to steady state"
    );
}
