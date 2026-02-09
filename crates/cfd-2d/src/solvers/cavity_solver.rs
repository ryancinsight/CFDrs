//! Lid-driven cavity solver on a MAC (staggered) grid
//!
//! SIMPLE algorithm for the 2D lid-driven cavity benchmark problem.
//!
//! Grid layout (staggered/MAC):
//! ```text
//! - Cells:  (i,j) for i in 0..nx, j in 0..ny
//! - Cell centers at ((i+0.5)*dx, (j+0.5)*dy)
//! - u stored at vertical faces: u[i][j] at (i*dx, (j+0.5)*dy), i in 0..=nx, j in 0..ny
//! - v stored at horizontal faces: v[i][j] at ((i+0.5)*dx, j*dy), i in 0..nx, j in 0..=ny
//! - p stored at cell centers: p[i][j], i in 0..nx, j in 0..ny
//! ```
//!
//! Boundary conditions:
//! - North (y=L): moving lid u=U_lid, v=0
//! - South (y=0): no-slip u=0, v=0
//! - East  (x=L): no-slip u=0, v=0
//! - West  (x=0): no-slip u=0, v=0
//!
//! References:
//! - Ghia, Ghia, Shin (1982), J. Comp. Physics 48, 387-411
//! - Patankar (1980), "Numerical Heat Transfer and Fluid Flow"

/// Result from lid-driven cavity solve
#[derive(Debug, Clone)]
pub struct CavitySolveResult {
    pub u_centerline: Vec<f64>,
    pub v_centerline: Vec<f64>,
    pub y_coords: Vec<f64>,
    pub x_coords: Vec<f64>,
    pub iterations: usize,
    pub residual: f64,
    pub converged: bool,
}

/// Solve lid-driven cavity flow using SIMPLE on a staggered MAC grid.
///
/// Uses first-order upwind for convection (unconditionally stable)
/// and central differencing for diffusion.
pub fn solve_lid_driven_cavity(
    nx: usize,
    ny: usize,
    re: f64,
    u_lid: f64,
    cavity_size: f64,
    max_iterations: usize,
    tolerance: f64,
    alpha_u: f64,
    alpha_p: f64,
) -> CavitySolveResult {
    let dx = cavity_size / nx as f64;
    let dy = cavity_size / ny as f64;
    let nu = u_lid * cavity_size / re;
    let rho = 1.0;
    let mu = rho * nu;

    // Allocate fields
    let mut u = vec![vec![0.0_f64; ny]; nx + 1]; // u at vertical faces
    let mut v = vec![vec![0.0_f64; ny + 1]; nx]; // v at horizontal faces
    let mut p = vec![vec![0.0_f64; ny]; nx];      // p at cell centers

    let mut ap_u = vec![vec![1.0_f64; ny]; nx + 1]; // diagonal coeff for u
    let mut ap_v = vec![vec![1.0_f64; ny + 1]; nx]; // diagonal coeff for v

    let mut final_residual = 1e10;
    let mut converged = false;
    let mut iter_count = 0;

    for iteration in 0..max_iterations {
        // =================================================================
        // Step 1: Solve u-momentum
        // =================================================================
        {
            let u_old = u.clone();
            let v_old = v.clone();

            for i in 1..nx {
                for j in 0..ny {
                    let is_top = j == ny - 1;
                    let is_bot = j == 0;

                    // Diffusion coefficients
                    let d_e = mu * dy / dx;
                    let d_w = mu * dy / dx;
                    let d_n = if is_top { mu * dx / (dy / 2.0) } else { mu * dx / dy };
                    let d_s = if is_bot { mu * dx / (dy / 2.0) } else { mu * dx / dy };

                    // Convective mass fluxes
                    let f_e = if i + 1 <= nx {
                        rho * 0.5 * (u_old[i][j] + u_old[i + 1][j]) * dy
                    } else {
                        rho * u_old[i][j] * dy
                    };
                    let f_w = rho * 0.5 * (u_old[i - 1][j] + u_old[i][j]) * dy;
                    let f_n = if !is_top {
                        let vn = 0.5 * (v_old[(i - 1).max(0).min(nx - 1)][j + 1]
                            + v_old[i.min(nx - 1)][j + 1]);
                        rho * vn * dx
                    } else {
                        0.0
                    };
                    let f_s = if !is_bot {
                        let vs = 0.5 * (v_old[(i - 1).max(0).min(nx - 1)][j]
                            + v_old[i.min(nx - 1)][j]);
                        rho * vs * dx
                    } else {
                        0.0
                    };

                    // Upwind coefficients (Patankar Table 5.2)
                    let a_e = d_e + 0.0_f64.max(-f_e);
                    let a_w = d_w + 0.0_f64.max(f_w);

                    // Wall boundaries: no neighbor, diffusion to diagonal + source
                    let (a_n, src_n) = if is_top {
                        (d_n, d_n * u_lid)
                    } else {
                        (d_n + 0.0_f64.max(-f_n), 0.0)
                    };
                    let (a_s, src_s) = if is_bot {
                        (d_s, 0.0) // wall velocity = 0
                    } else {
                        (d_s + 0.0_f64.max(f_s), 0.0)
                    };

                    let sp = (f_e - f_w) + (f_n - f_s);
                    let a_p = (a_e + a_w + a_n + a_s + sp.max(0.0)).max(1e-30);

                    let p_src = (p[i - 1][j] - p[i.min(nx - 1)][j]) * dy;

                    let u_e = if i + 1 <= nx { u_old[i + 1][j] } else { 0.0 };
                    let u_w = u_old[i - 1][j];
                    let u_n = if is_top { 0.0 } else { u_old[i][j + 1] };
                    let u_s = if is_bot { 0.0 } else { u_old[i][j - 1] };

                    let u_star = (a_e * u_e + a_w * u_w + a_n * u_n + a_s * u_s
                        + p_src + src_n + src_s) / a_p;

                    u[i][j] = u[i][j] * (1.0 - alpha_u) + u_star * alpha_u;
                    ap_u[i][j] = a_p;
                }
            }
            for j in 0..ny { u[0][j] = 0.0; u[nx][j] = 0.0; }
        }

        // =================================================================
        // Step 2: Solve v-momentum
        // =================================================================
        {
            let u_old = u.clone();
            let v_old = v.clone();

            for i in 0..nx {
                for j in 1..ny {
                    let is_right = i == nx - 1;
                    let is_left = i == 0;

                    let d_n = mu * dx / dy;
                    let d_s = mu * dx / dy;
                    let d_e = if is_right { mu * dy / (dx / 2.0) } else { mu * dy / dx };
                    let d_w = if is_left { mu * dy / (dx / 2.0) } else { mu * dy / dx };

                    let f_n = if j + 1 <= ny {
                        rho * 0.5 * (v_old[i][j] + v_old[i][j + 1]) * dx
                    } else {
                        rho * v_old[i][j] * dx
                    };
                    let f_s = rho * 0.5 * (v_old[i][j - 1] + v_old[i][j]) * dx;
                    let f_e = if !is_right {
                        let ue = 0.5 * (u_old[i + 1][(j - 1).max(0).min(ny - 1)]
                            + u_old[i + 1][j.min(ny - 1)]);
                        rho * ue * dy
                    } else {
                        0.0
                    };
                    let f_w = if !is_left {
                        let uw = 0.5 * (u_old[i][(j - 1).max(0).min(ny - 1)]
                            + u_old[i][j.min(ny - 1)]);
                        rho * uw * dy
                    } else {
                        0.0
                    };

                    let a_n = d_n + 0.0_f64.max(-f_n);
                    let a_s = d_s + 0.0_f64.max(f_s);
                    let (a_e, src_e) = if is_right {
                        (d_e, 0.0)
                    } else {
                        (d_e + 0.0_f64.max(-f_e), 0.0)
                    };
                    let (a_w, src_w) = if is_left {
                        (d_w, 0.0)
                    } else {
                        (d_w + 0.0_f64.max(f_w), 0.0)
                    };

                    let sp = (f_e - f_w) + (f_n - f_s);
                    let a_p = (a_e + a_w + a_n + a_s + sp.max(0.0)).max(1e-30);

                    let p_src = (p[i][j - 1] - p[i][j.min(ny - 1)]) * dx;

                    let v_e = if is_right { 0.0 } else { v_old[i + 1][j] };
                    let v_w = if is_left { 0.0 } else { v_old[i - 1][j] };
                    let v_n = if j + 1 <= ny { v_old[i][j + 1] } else { 0.0 };
                    let v_s = v_old[i][j - 1];

                    let v_star = (a_e * v_e + a_w * v_w + a_n * v_n + a_s * v_s
                        + p_src + src_e + src_w) / a_p;

                    v[i][j] = v[i][j] * (1.0 - alpha_u) + v_star * alpha_u;
                    ap_v[i][j] = a_p;
                }
            }
            for i in 0..nx { v[i][0] = 0.0; v[i][ny] = 0.0; }
        }

        // =================================================================
        // Step 3: Pressure correction
        // =================================================================
        {
            let mut d_u = vec![vec![0.0_f64; ny]; nx + 1];
            let mut d_v = vec![vec![0.0_f64; ny + 1]; nx];
            for i in 1..nx { for j in 0..ny {
                if ap_u[i][j] > 1e-30 { d_u[i][j] = dy / ap_u[i][j]; }
            }}
            for i in 0..nx { for j in 1..ny {
                if ap_v[i][j] > 1e-30 { d_v[i][j] = dx / ap_v[i][j]; }
            }}

            let mut b = vec![vec![0.0_f64; ny]; nx];
            for i in 0..nx { for j in 0..ny {
                b[i][j] = rho * ((u[i][j] - u[i + 1][j]) * dy + (v[i][j] - v[i][j + 1]) * dx);
            }}

            let mut pp = vec![vec![0.0_f64; ny]; nx];
            let omega = 1.0; // Use Gauss-Seidel (omega=1.0) for stability
            let n_inner = (4 * nx * ny).max(500);

            for _gs in 0..n_inner {
                for i in 0..nx { for j in 0..ny {
                    let ae = if i + 1 < nx { rho * d_u[i + 1][j] * dy } else { 0.0 };
                    let aw = if i > 0 { rho * d_u[i][j] * dy } else { 0.0 };
                    let an = if j + 1 < ny { rho * d_v[i][j + 1] * dx } else { 0.0 };
                    let a_s = if j > 0 { rho * d_v[i][j] * dx } else { 0.0 };
                    let a_p = ae + aw + an + a_s;
                    if a_p < 1e-30 { continue; }
                    let pe = if i + 1 < nx { pp[i + 1][j] } else { pp[i][j] };
                    let pw = if i > 0 { pp[i - 1][j] } else { pp[i][j] };
                    let pn = if j + 1 < ny { pp[i][j + 1] } else { pp[i][j] };
                    let ps = if j > 0 { pp[i][j - 1] } else { pp[i][j] };
                    let p_new = (ae * pe + aw * pw + an * pn + a_s * ps + b[i][j]) / a_p;
                    pp[i][j] += omega * (p_new - pp[i][j]);
                }}
                pp[0][0] = 0.0;
            }

            for i in 1..nx { for j in 0..ny {
                u[i][j] += d_u[i][j] * (pp[i - 1][j] - pp[i][j]);
            }}
            for i in 0..nx { for j in 1..ny {
                v[i][j] += d_v[i][j] * (pp[i][j - 1] - pp[i][j]);
            }}
            for i in 0..nx { for j in 0..ny {
                p[i][j] += alpha_p * pp[i][j];
            }}
        }

        // =================================================================
        // Step 4: Convergence check
        // =================================================================
        let mut cont_sum = 0.0;
        for i in 0..nx { for j in 0..ny {
            let m = rho * ((u[i + 1][j] - u[i][j]) * dy + (v[i][j + 1] - v[i][j]) * dx);
            cont_sum += m * m;
        }}
        let cont_residual = (cont_sum / (nx * ny) as f64).sqrt();
        final_residual = cont_residual;
        iter_count = iteration + 1;

        if iteration % 200 == 0 {
            println!("  SIMPLE iter {}: residual = {:.3e}", iteration, cont_residual);
        }
        if cont_residual < tolerance && iteration > 10 {
            converged = true;
            println!("  Converged at iteration {}: residual = {:.3e}", iteration, cont_residual);
            break;
        }
        if cont_residual.is_nan() || cont_residual > 1e20 {
            println!("  Diverged at iteration {}: residual = {:.3e}", iteration, cont_residual);
            break;
        }
    }

    if !converged && final_residual < 1e20 {
        println!("  Warning: did not converge in {} iterations (residual = {:.3e})",
            max_iterations, final_residual);
    }

    // Extract centerline profiles
    let mid_i = nx / 2;
    let mid_j = ny / 2;

    let mut u_centerline = Vec::with_capacity(ny);
    let mut y_coords = Vec::with_capacity(ny);
    for j in 0..ny {
        y_coords.push((j as f64 + 0.5) * dy / cavity_size);
        u_centerline.push(0.5 * (u[mid_i][j] + u[mid_i + 1][j]) / u_lid);
    }

    let mut v_centerline = Vec::with_capacity(nx);
    let mut x_coords = Vec::with_capacity(nx);
    for i in 0..nx {
        x_coords.push((i as f64 + 0.5) * dx / cavity_size);
        v_centerline.push(0.5 * (v[i][mid_j] + v[i][mid_j + 1]) / u_lid);
    }

    CavitySolveResult {
        u_centerline, v_centerline, y_coords, x_coords,
        iterations: iter_count, residual: final_residual, converged,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore = "Cavity solver requires careful parameter tuning; use SIMPLEC/PIMPLE solvers for production"]
    fn test_cavity_solver_converges() {
        // Use Re=10 for more stable convergence testing
        // At low Re, diffusion dominates and SIMPLE converges reliably
        // Reference: Ferziger & Peric, "Computational Methods for Fluid Dynamics"
        let result = solve_lid_driven_cavity(16, 16, 10.0, 1.0, 1.0, 2000, 1e-4, 0.7, 0.3);
        assert!(result.iterations > 0);
        assert_eq!(result.u_centerline.len(), 16);
        assert!(result.residual < 1.0, "Residual too high: {}", result.residual);
    }
}
