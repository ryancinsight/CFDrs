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
/// and central differencing for diffusion.  Boundary diffusion uses
/// the half-cell distance where the first interior node is offset
/// from the wall (y-direction for u, x-direction for v), following
/// Patankar (1980) §5.3.
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
    let mut p = vec![vec![0.0_f64; ny]; nx]; // p at cell centres

    let mut ap_u = vec![vec![1.0_f64; ny]; nx + 1];
    let mut ap_v = vec![vec![1.0_f64; ny + 1]; nx];

    let mut final_residual = 1e10;
    let mut converged = false;
    let mut iter_count = 0;

    for iteration in 0..max_iterations {
        // =================================================================
        // Step 1: Solve u-momentum  (faces i=1..nx-1, j=0..ny-1)
        // =================================================================
        //
        //  u[i][j] lives at (i*dx, (j+0.5)*dy).
        //  Boundary distances:
        //    East/West neighbours are other u nodes at full dx apart.
        //    North:  When j==ny-1, the lid (u=u_lid) is at y=L, and
        //            the u-node is at (ny-0.5)*dy → distance dy/2.
        //    South:  When j==0, the wall (u=0) is at y=0, and the
        //            u-node is at 0.5*dy → distance dy/2.
        {
            let u_old = u.clone();
            let v_old = v.clone();

            for i in 1..nx {
                for j in 0..ny {
                    let is_top = j == ny - 1;
                    let is_bot = j == 0;

                    // --- diffusion coefficients ---
                    // x-direction: node-to-node distance is always dx
                    let d_e = mu * dy / dx;
                    let d_w = mu * dy / dx;
                    // y-direction: half-cell at top/bottom walls
                    let d_n = if is_top { mu * dx / (dy / 2.0) } else { mu * dx / dy };
                    let d_s = if is_bot { mu * dx / (dy / 2.0) } else { mu * dx / dy };

                    // --- convective mass fluxes through CV faces ---
                    // East face at x=(i+0.5)*dx
                    let f_e = if i + 1 <= nx {
                        rho * 0.5 * (u_old[i][j] + u_old[i + 1][j]) * dy
                    } else {
                        rho * u_old[i][j] * dy
                    };
                    // West face at x=(i-0.5)*dx
                    let f_w = rho * 0.5 * (u_old[i - 1][j] + u_old[i][j]) * dy;
                    // North face at y=(j+1)*dy – v interpolated at x=i*dx
                    let f_n = if !is_top {
                        let il = (i - 1).min(nx - 1);
                        let ir = i.min(nx - 1);
                        rho * 0.5 * (v_old[il][j + 1] + v_old[ir][j + 1]) * dx
                    } else {
                        0.0 // lid is solid → v=0
                    };
                    // South face at y=j*dy
                    let f_s = if !is_bot {
                        let il = (i - 1).min(nx - 1);
                        let ir = i.min(nx - 1);
                        rho * 0.5 * (v_old[il][j] + v_old[ir][j]) * dx
                    } else {
                        0.0 // wall → v=0
                    };

                    // --- upwind coefficients (Patankar Table 5.2) ---
                    let a_e = d_e + 0.0_f64.max(-f_e);
                    let a_w = d_w + 0.0_f64.max(f_w);
                    let (a_n, src_n) = if is_top {
                        // Wall treatment: known u_lid goes into source
                        (d_n, d_n * u_lid)
                    } else {
                        (d_n + 0.0_f64.max(-f_n), 0.0)
                    };
                    let (a_s, src_s) = if is_bot {
                        (d_s, 0.0) // wall u=0
                    } else {
                        (d_s + 0.0_f64.max(f_s), 0.0)
                    };

                    // Mass source imbalance (goes → 0 as continuity is satisfied)
                    let delta_f = (f_e - f_w) + (f_n - f_s);
                    let a_p = (a_e + a_w + a_n + a_s + delta_f.max(0.0)).max(1e-30);

                    // Pressure gradient force
                    let p_src = (p[i - 1][j] - p[i.min(nx - 1)][j]) * dy;

                    // Neighbour velocities (wall → 0 except lid → handled via src_n)
                    let u_e = if i + 1 <= nx { u_old[i + 1][j] } else { 0.0 };
                    let u_w = u_old[i - 1][j];
                    let u_n = if is_top { 0.0 } else { u_old[i][j + 1] };
                    let u_s = if is_bot { 0.0 } else { u_old[i][j - 1] };

                    let u_star =
                        (a_e * u_e + a_w * u_w + a_n * u_n + a_s * u_s + p_src + src_n + src_s)
                            / a_p;

                    u[i][j] = u[i][j] * (1.0 - alpha_u) + u_star * alpha_u;
                    ap_u[i][j] = a_p;
                }
            }
            // Enforce wall BCs on u
            for j in 0..ny {
                u[0][j] = 0.0;
                u[nx][j] = 0.0;
            }
        }

        // =================================================================
        // Step 2: Solve v-momentum  (faces i=0..nx-1, j=1..ny-1)
        // =================================================================
        //
        //  v[i][j] lives at ((i+0.5)*dx, j*dy).
        //  Boundary distances:
        //    North/South neighbours are other v nodes at full dy apart.
        //    East:  When i==nx-1 the wall (v=0) is at x=L, and the
        //           v-node is at (nx-0.5)*dx → distance dx/2.
        //    West:  When i==0 the wall (v=0) is at x=0, and the
        //           v-node is at 0.5*dx → distance dx/2.
        {
            let u_old = u.clone();
            let v_old = v.clone();

            for i in 0..nx {
                for j in 1..ny {
                    let is_right = i == nx - 1;
                    let is_left = i == 0;

                    // --- diffusion coefficients ---
                    // y-direction: node-to-node distance is always dy
                    let d_n = mu * dx / dy;
                    let d_s = mu * dx / dy;
                    // x-direction: half-cell at left/right walls
                    let d_e = if is_right { mu * dy / (dx / 2.0) } else { mu * dy / dx };
                    let d_w = if is_left { mu * dy / (dx / 2.0) } else { mu * dy / dx };

                    // --- convective mass fluxes ---
                    let f_n = if j + 1 <= ny {
                        rho * 0.5 * (v_old[i][j] + v_old[i][j + 1]) * dx
                    } else {
                        rho * v_old[i][j] * dx
                    };
                    let f_s = rho * 0.5 * (v_old[i][j - 1] + v_old[i][j]) * dx;
                    let f_e = if !is_right {
                        let jb = (j - 1).min(ny - 1);
                        let jt = j.min(ny - 1);
                        rho * 0.5 * (u_old[i + 1][jb] + u_old[i + 1][jt]) * dy
                    } else {
                        0.0
                    };
                    let f_w = if !is_left {
                        let jb = (j - 1).min(ny - 1);
                        let jt = j.min(ny - 1);
                        rho * 0.5 * (u_old[i][jb] + u_old[i][jt]) * dy
                    } else {
                        0.0
                    };

                    let a_n = d_n + 0.0_f64.max(-f_n);
                    let a_s = d_s + 0.0_f64.max(f_s);
                    let (a_e, _src_e) = if is_right {
                        (d_e, 0.0)
                    } else {
                        (d_e + 0.0_f64.max(-f_e), 0.0)
                    };
                    let (a_w, _src_w) = if is_left {
                        (d_w, 0.0)
                    } else {
                        (d_w + 0.0_f64.max(f_w), 0.0)
                    };

                    let delta_f = (f_e - f_w) + (f_n - f_s);
                    let a_p = (a_e + a_w + a_n + a_s + delta_f.max(0.0)).max(1e-30);

                    let p_src = (p[i][j - 1] - p[i][j.min(ny - 1)]) * dx;

                    let v_e = if is_right { 0.0 } else { v_old[i + 1][j] };
                    let v_w = if is_left { 0.0 } else { v_old[i - 1][j] };
                    let v_n = if j + 1 <= ny { v_old[i][j + 1] } else { 0.0 };
                    let v_s = v_old[i][j - 1];

                    let v_star =
                        (a_e * v_e + a_w * v_w + a_n * v_n + a_s * v_s + p_src) / a_p;

                    v[i][j] = v[i][j] * (1.0 - alpha_u) + v_star * alpha_u;
                    ap_v[i][j] = a_p;
                }
            }
            // Enforce wall BCs on v
            for i in 0..nx {
                v[i][0] = 0.0;
                v[i][ny] = 0.0;
            }
        }

        // =================================================================
        // Step 3: Pressure correction (Poisson via SOR)
        // =================================================================
        {
            // d-coefficients: relate velocity correction to pressure correction
            let mut d_u = vec![vec![0.0_f64; ny]; nx + 1];
            let mut d_v = vec![vec![0.0_f64; ny + 1]; nx];
            for i in 1..nx {
                for j in 0..ny {
                    if ap_u[i][j] > 1e-30 {
                        d_u[i][j] = dy / ap_u[i][j];
                    }
                }
            }
            for i in 0..nx {
                for j in 1..ny {
                    if ap_v[i][j] > 1e-30 {
                        d_v[i][j] = dx / ap_v[i][j];
                    }
                }
            }

            // Mass source (negative of divergence of starred velocity)
            let mut b = vec![vec![0.0_f64; ny]; nx];
            for i in 0..nx {
                for j in 0..ny {
                    b[i][j] = rho
                        * ((u[i][j] - u[i + 1][j]) * dy + (v[i][j] - v[i][j + 1]) * dx);
                }
            }

            // SOR for pressure correction pp
            let mut pp = vec![vec![0.0_f64; ny]; nx];
            let omega = 1.2_f64; // conservative SOR (< optimal for stability)
            let n_inner = (8 * nx * ny).max(1000);

            for _gs in 0..n_inner {
                for i in 0..nx {
                    for j in 0..ny {
                        // Skip the pinned reference cell
                        if i == 0 && j == 0 {
                            continue;
                        }
                        let ae = if i + 1 < nx { rho * d_u[i + 1][j] * dy } else { 0.0 };
                        let aw = if i > 0 { rho * d_u[i][j] * dy } else { 0.0 };
                        let an = if j + 1 < ny { rho * d_v[i][j + 1] * dx } else { 0.0 };
                        let a_s = if j > 0 { rho * d_v[i][j] * dx } else { 0.0 };
                        let a_p = ae + aw + an + a_s;
                        if a_p < 1e-30 {
                            continue;
                        }
                        let pe = if i + 1 < nx { pp[i + 1][j] } else { pp[i][j] };
                        let pw = if i > 0 { pp[i - 1][j] } else { pp[i][j] };
                        let pn = if j + 1 < ny { pp[i][j + 1] } else { pp[i][j] };
                        let ps = if j > 0 { pp[i][j - 1] } else { pp[i][j] };
                        let p_new = (ae * pe + aw * pw + an * pn + a_s * ps + b[i][j]) / a_p;
                        pp[i][j] += omega * (p_new - pp[i][j]);
                    }
                }
                // Reference pressure pinned at (0,0)
                pp[0][0] = 0.0;
            }

            // Correct velocities
            for i in 1..nx {
                for j in 0..ny {
                    u[i][j] += d_u[i][j] * (pp[i - 1][j] - pp[i][j]);
                }
            }
            for i in 0..nx {
                for j in 1..ny {
                    v[i][j] += d_v[i][j] * (pp[i][j - 1] - pp[i][j]);
                }
            }

            // Update pressure
            for i in 0..nx {
                for j in 0..ny {
                    p[i][j] += alpha_p * pp[i][j];
                }
            }
        }

        // =================================================================
        // Step 4: Convergence check (RMS of mass imbalance)
        // =================================================================
        let mut cont_sum = 0.0;
        for i in 0..nx {
            for j in 0..ny {
                let m = rho
                    * ((u[i + 1][j] - u[i][j]) * dy + (v[i][j + 1] - v[i][j]) * dx);
                cont_sum += m * m;
            }
        }
        let cont_residual = (cont_sum / (nx * ny) as f64).sqrt();
        final_residual = cont_residual;
        iter_count = iteration + 1;

        if cont_residual < tolerance && iteration > 10 {
            converged = true;
            break;
        }
        if cont_residual.is_nan() || cont_residual > 1e20 {
            break;
        }
    }

    // Extract centerline profiles (cell-centred averages)
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
        u_centerline,
        v_centerline,
        y_coords,
        x_coords,
        iterations: iter_count,
        residual: final_residual,
        converged,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cavity_solver_converges() {
        // Re=1, 16×16 grid – diffusion-dominated, converges rapidly.
        // Validated against Ghia et al. (1982) in Python scripts for
        // higher Re values.
        let result = solve_lid_driven_cavity(16, 16, 1.0, 1.0, 1.0, 5000, 1e-4, 0.7, 0.3);
        assert!(result.iterations > 0);
        assert_eq!(result.u_centerline.len(), 16);
        assert!(
            result.residual < 1e-3,
            "Residual too high: {}",
            result.residual
        );
    }

    #[test]
    fn test_cavity_solver_re100_32x32() {
        // Re=100, 32×32 grid with conservative under-relaxation.
        // Cell Peclet number Pe ≈ 3.1 → first-order upwind is stable.
        let result = solve_lid_driven_cavity(32, 32, 100.0, 1.0, 1.0, 10_000, 1e-4, 0.5, 0.1);
        assert!(result.iterations > 0);
        assert_eq!(result.u_centerline.len(), 32);
        assert!(
            !result.residual.is_nan() && result.residual < 1.0,
            "Solver diverged: residual = {}",
            result.residual
        );
    }
}
