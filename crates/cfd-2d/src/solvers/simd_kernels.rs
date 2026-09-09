//! SIMD-accelerated kernels for 2D CFD solvers
//!
//! Provides optimized numerical kernels routing their bulk stencil arithmetic
//! through `cfd_math::simd::SimdOps` (hermes-simd runtime dispatch, CFDRS-GA-004).
//!
//! # Invariant (SIMD Numerical Equivalence)
//!
//! Each SIMD kernel computes the same stencil operations as the scalar reference.
//! The Jacobi update $\phi_i^{k+1} = (b_i - \sum_{j \ne i} a_{ij}\phi_j^k) / a_{ii}$
//! is applied element-wise via SIMD lanes without altering the iteration order,
//! preserving convergence guarantees of the underlying solver.
//!
//! # Numerical equivalence
//!
//! The grid is row-major (`idx = i * ny + j`), so every 5-point stencil
//! decomposes into contiguous row-slice operations. Kernels whose scalar
//! references precompute reciprocal spacings (`1/dx`, `1/dx²`) multiply by
//! that same precomputed value and are bit-identical to the historical
//! scalar loops; Jacobi and Gauss-Seidel divide by `dx²`/`dy²` directly in
//! the scalar form, and their ported forms multiply by the reciprocal
//! instead, which differs by at most two ulp (reciprocal rounding plus
//! product rounding versus one direct division). Red-black Gauss-Seidel
//! keeps the historical pass structure — a full red pass over the grid, then
//! a full black pass — with each pass's residual evaluation vectorized per
//! row: red residuals read only black values (unchanged during the red
//! pass), and each row's residuals are evaluated from the current `phi`
//! exactly when the scalar interleaved sweep would read them, so the
//! vectorized form is order-equivalent.

use cfd_core::error::Result;
use cfd_math::simd::SimdOps;
use std::sync::OnceLock;

// Global SIMD ops instance (initialized once)
static SIMD_OPS: OnceLock<SimdOps> = OnceLock::new();

/// Get or initialize the global SIMD ops
fn simd_ops() -> &'static SimdOps {
    SIMD_OPS.get_or_init(SimdOps::new)
}

/// Interior width of a row: the number of grid columns with both neighbors
/// inside the row. Callers guard `ny >= 3`, so this is at least 1.
fn interior_len(ny: usize) -> usize {
    ny - 2
}

/// SIMD-accelerated Jacobi iteration for Poisson equation
///
/// Solves ∇²φ = f using Jacobi iteration with SIMD acceleration
pub fn jacobi_iteration_simd(
    phi: &mut [f32],
    phi_new: &mut [f32],
    source: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
) -> Result<()> {
    let dx2 = dx * dx;
    let dy2 = dy * dy;
    let factor = 0.5 / (1.0 / dx2 + 1.0 / dy2);
    let inv_dx2 = 1.0 / dx2;
    let inv_dy2 = 1.0 / dy2;

    if nx >= 3 && ny >= 3 {
        let ops = simd_ops();
        let m = interior_len(ny);
        // Scratch rows, hoisted out of the row loop. Full-row buffers hold
        // results aligned to grid columns (`[k]` ↔ column `k`); interior
        // buffers hold `[k]` ↔ column `k + 1`.
        let mut t_full = vec![0.0f32; ny];
        let mut t_shift = vec![0.0f32; ny - 1];
        let mut a1 = vec![0.0f32; m];
        let mut a2 = vec![0.0f32; m];

        for i in 1..nx - 1 {
            let row = i * ny;
            let x_prev = &phi[(i - 1) * ny..i * ny];
            let x_here = &phi[row..row + ny];
            let x_next = &phi[(i + 1) * ny..(i + 2) * ny];

            // lap_x[j] = (left + right) / dx² as a full-row sum (contiguous
            // neighbor rows), then scaled by the reciprocal.
            ops.add(x_prev, x_next, &mut t_full)?;
            ops.scale_in_place(&mut t_full, inv_dx2)?;

            // lap_y[j] = (bottom + top) / dy²: the shifted neighbor sum
            // `x_here[k] + x_here[k + 2]` pairs two ny-2 slices and yields
            // `[k]` ↔ column `k + 1`.
            ops.add(&x_here[..ny - 2], &x_here[2..], &mut t_shift[..ny - 2])?;
            ops.scale_in_place(&mut t_shift[..ny - 2], inv_dy2)?;

            // phi_new = factor * (lap_x + lap_y - source), interior columns.
            a1.copy_from_slice(&t_full[1..ny - 1]);
            ops.add(&a1, &t_shift[..m], &mut a2)?;
            ops.sub(&a2, &source[row + 1..row + ny - 1], &mut a1)?;
            ops.scale_in_place(&mut a1, factor)?;
            phi_new[row + 1..row + ny - 1].copy_from_slice(&a1);
        }
    }

    // Copy boundaries
    for i in 0..nx {
        phi_new[i * ny] = phi[i * ny];
        phi_new[i * ny + ny - 1] = phi[i * ny + ny - 1];
    }
    for j in 0..ny {
        phi_new[j] = phi[j];
        phi_new[(nx - 1) * ny + j] = phi[(nx - 1) * ny + j];
    }

    Ok(())
}

/// SIMD-accelerated Gauss-Seidel iteration
///
/// Red-black SOR with the historical pass structure: a full red pass
/// (`i + j` even) over the grid, then a full black pass. Red residuals
/// depend only on black values, which the red pass never touches; each
/// row's residuals are evaluated from the current `phi` exactly when the
/// interleaved scalar sweep would read them. Only the per-color commit
/// (`phi[idx] = (1 - omega) * phi[idx] + factor * residual`) stays scalar,
/// so the sequential SOR semantics are preserved bit-for-bit.
pub fn gauss_seidel_simd(
    phi: &mut [f32],
    source: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
    omega: f32, // SOR relaxation parameter
) -> Result<()> {
    let dx2 = dx * dx;
    let dy2 = dy * dy;
    let factor = omega / (2.0 * (1.0 / dx2 + 1.0 / dy2));
    let one_minus_omega = 1.0 - omega;
    let inv_dx2 = 1.0 / dx2;
    let inv_dy2 = 1.0 / dy2;

    if nx >= 3 && ny >= 3 {
        let ops = simd_ops();
        let m = interior_len(ny);
        let mut t_full = vec![0.0f32; ny];
        let mut t_shift = vec![0.0f32; ny - 1];
        let mut resid = vec![0.0f32; m];

        for color in 0..2 {
            for i in 1..nx - 1 {
                let row = i * ny;
                let x_prev = &phi[(i - 1) * ny..i * ny];
                let x_here = &phi[row..row + ny];
                let x_next = &phi[(i + 1) * ny..(i + 2) * ny];

                // residual[j] = (left + right) / dx² + (bottom + top) / dy²
                //               - source[j], interior columns, evaluated
                // from the current `phi`.
                ops.add(x_prev, x_next, &mut t_full)?;
                ops.scale_in_place(&mut t_full, inv_dx2)?;
                ops.add(&x_here[..ny - 2], &x_here[2..], &mut t_shift[..ny - 2])?;
                ops.scale_in_place(&mut t_shift[..ny - 2], inv_dy2)?;
                ops.add(&t_full[1..ny - 1], &t_shift[..m], &mut resid)?;
                // The final residual lands in `t_shift[..m]` (distinct from
                // both inputs); the commit loop below reads it from there.
                ops.sub(&resid, &source[row + 1..row + ny - 1], &mut t_shift[..m])?;

                // Commit this color's cells from the precomputed residual;
                // the SOR blend reads only the cell itself.
                for j in 1..ny - 1 {
                    if (i + j) % 2 == color {
                        let idx = row + j;
                        phi[idx] = one_minus_omega * phi[idx] + factor * t_shift[j - 1];
                    }
                }
            }
        }
    }

    Ok(())
}

/// SIMD-accelerated velocity interpolation
pub fn interpolate_velocity_simd(
    u_cell: &[f32],
    v_cell: &[f32],
    u_face: &mut [f32],
    v_face: &mut [f32],
    nx: usize,
    ny: usize,
) -> Result<()> {
    let ops = simd_ops();

    // Interpolate u-velocity to x-faces
    for i in 0..nx - 1 {
        let start_idx = i * ny;
        let end_idx = start_idx + ny;

        // Use SIMD for averaging
        let left = &u_cell[start_idx..end_idx];
        let right = &u_cell[start_idx + ny..end_idx + ny];
        let face = &mut u_face[start_idx..end_idx];

        // Average: face = 0.5 * (left + right)
        ops.add(left, right, face)?;
        ops.scale_in_place(face, 0.5)?;
    }

    // Interpolate v-velocity to y-faces
    let mut temp_bottom = vec![0.0f32; nx];
    let mut temp_top = vec![0.0f32; nx];
    let mut temp_face = vec![0.0f32; nx];
    for j in 0..ny - 1 {
        for i in 0..nx {
            temp_bottom[i] = v_cell[i * ny + j];
            temp_top[i] = v_cell[i * ny + j + 1];
        }

        // Average using SIMD
        ops.add(&temp_bottom, &temp_top, &mut temp_face)?;
        ops.scale_in_place(&mut temp_face, 0.5)?;

        // Copy back
        for i in 0..nx {
            v_face[i * (ny - 1) + j] = temp_face[i];
        }
    }

    Ok(())
}

/// SIMD-accelerated divergence calculation
pub fn calculate_divergence_simd(
    u: &[f32],
    v: &[f32],
    divergence: &mut [f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
) -> Result<()> {
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;

    if nx >= 3 && ny >= 3 {
        let ops = simd_ops();
        let m = interior_len(ny);
        let mut t_full = vec![0.0f32; ny];
        let mut t_shift = vec![0.0f32; ny - 1];
        let mut a2 = vec![0.0f32; m];

        // div = ∂u/∂x + ∂v/∂y, central differences, interior only. The
        // scalar reference multiplies `(a - b) * 0.5 * inv_dx`
        // left-to-right; the two scale steps reproduce that exactly, and
        // `inv_dx` is the same precomputed reciprocal the scalar form uses,
        // so results are bit-identical.
        for i in 1..nx - 1 {
            let row = i * ny;

            // (u[i+1][j] - u[i-1][j]) * 0.5 * inv_dx: full-row difference.
            ops.sub(
                &u[(i + 1) * ny..(i + 2) * ny],
                &u[(i - 1) * ny..i * ny],
                &mut t_full,
            )?;
            ops.scale_in_place(&mut t_full, 0.5)?;
            ops.scale_in_place(&mut t_full, inv_dx)?;

            // (v[i][j+1] - v[i][j-1]) * 0.5 * inv_dy: shifted difference,
            // `v[k + 2] - v[k]` yields `[k]` ↔ column `k + 1` (length ny-2).
            ops.sub(
                &v[row + 2..row + ny],
                &v[row..row + ny - 2],
                &mut t_shift[..ny - 2],
            )?;
            ops.scale_in_place(&mut t_shift[..ny - 2], 0.5)?;
            ops.scale_in_place(&mut t_shift[..ny - 2], inv_dy)?;

            ops.add(&t_full[1..ny - 1], &t_shift[..m], &mut a2)?;
            divergence[row + 1..row + ny - 1].copy_from_slice(&a2);
        }
    }

    Ok(())
}

/// SIMD-accelerated gradient calculation
pub fn calculate_gradient_simd(
    phi: &[f32],
    grad_x: &mut [f32],
    grad_y: &mut [f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
) -> Result<()> {
    let inv_dx = 0.5 / dx;
    let inv_dy = 0.5 / dy;

    if nx >= 3 && ny >= 3 {
        let ops = simd_ops();
        let mut t_full = vec![0.0f32; ny];
        let mut t_shift = vec![0.0f32; ny - 1];

        // Central differences, row-wise for cache locality. The scalar
        // reference multiplies once by the precomputed `0.5 / dx` factor,
        // which the port reuses verbatim (bit-identical).
        for i in 1..nx - 1 {
            let row = i * ny;

            ops.sub(
                &phi[(i + 1) * ny..(i + 2) * ny],
                &phi[(i - 1) * ny..i * ny],
                &mut t_full,
            )?;
            ops.scale_in_place(&mut t_full, inv_dx)?;
            grad_x[row + 1..row + ny - 1].copy_from_slice(&t_full[1..ny - 1]);

            ops.sub(
                &phi[row + 2..row + ny],
                &phi[row..row + ny - 2],
                &mut t_shift[..ny - 2],
            )?;
            ops.scale_in_place(&mut t_shift[..ny - 2], inv_dy)?;
            grad_y[row + 1..row + ny - 1].copy_from_slice(&t_shift[..ny - 2]);
        }
    }

    Ok(())
}

/// SIMD-accelerated residual calculation for iterative solvers
pub fn calculate_residual_simd(
    phi: &[f32],
    source: &[f32],
    residual: &mut [f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
) -> Result<f32> {
    let dx2 = dx * dx;
    let dy2 = dy * dy;
    let inv_dx2 = 1.0 / dx2;
    let inv_dy2 = 1.0 / dy2;

    let mut max_residual = 0.0f32;

    if nx >= 3 && ny >= 3 {
        let ops = simd_ops();
        let m = interior_len(ny);
        let mut here2 = vec![0.0f32; ny]; // 2 * phi row, aligned to columns
        let mut lap_x = vec![0.0f32; ny];
        let mut tmp = vec![0.0f32; m];
        let mut r = vec![0.0f32; m];

        // 5-point discrete Laplacian: r = ∇²φ - f, with maximum tracking
        // for convergence monitoring. The scalar reference evaluates
        // `(prev - 2 * here + next) * inv_dx2` as a subtract-then-add with
        // the precomputed reciprocal, reproduced exactly (bit-identical).
        for i in 1..nx - 1 {
            let row = i * ny;
            let x_here = &phi[row..row + ny];

            // here2 = 2 * phi[row]
            here2.copy_from_slice(x_here);
            ops.scale_in_place(&mut here2, 2.0)?;

            // lap_x[j] = (prev - 2*here + next) * inv_dx2, evaluated as
            // ((prev - 2*here) + next) * inv_dx2 like the scalar form.
            ops.sub(&phi[(i - 1) * ny..i * ny], &here2, &mut lap_x)?;
            ops.add(
                &lap_x[1..ny - 1],
                &phi[(i + 1) * ny + 1..(i + 2) * ny - 1],
                &mut tmp,
            )?;
            ops.scale_in_place(&mut tmp, inv_dx2)?;
            // `tmp` now holds lap_x (interior); `lap_x`'s interior is free
            // scratch from here on.

            // lap_y[j] = (bottom - 2*here + top) * inv_dy2, shifted family:
            // bottom[j] = x_here[k], 2*here[j] = here2[k + 1],
            // top[j] = x_here[k + 2], all ↔ column `k + 1`.
            ops.sub(&x_here[..m], &here2[1..ny - 1], &mut r)?;
            ops.add(&r, &x_here[2..], &mut lap_x[1..ny - 1])?;
            ops.scale_in_place(&mut lap_x[1..ny - 1], inv_dy2)?;

            // residual = lap_x + lap_y - source
            ops.add(&tmp, &lap_x[1..ny - 1], &mut r)?;
            ops.sub(&r, &source[row + 1..row + ny - 1], &mut tmp)?;
            residual[row + 1..row + ny - 1].copy_from_slice(&tmp);

            let row_max = ops.abs_max_f32(&tmp)?;
            max_residual = max_residual.max(row_max);
        }
    }

    Ok(max_residual)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One sweep from a zero field under a unit source has a closed form.
    ///
    /// With `phi = 0` both Laplacian sums vanish, so every interior cell
    /// becomes `-0.5 / (1/dx^2 + 1/dy^2)` times the source. This is the only
    /// case here that drives the source term: `test_jacobi_known_stencil`
    /// exercises the neighbour sums against a zero source. `phi_new` starts
    /// filled with a sentinel, so a kernel that writes nothing fails instead
    /// of inheriting a zero buffer.
    #[test]
    fn jacobi_sweep_from_zero_is_the_scaled_source() {
        let nx = 10;
        let ny = 10;
        let (dx, dy) = (0.1f32, 0.1f32);
        let mut phi = vec![0.0f32; nx * ny];
        let mut phi_new = vec![f32::MIN; nx * ny];
        let source = vec![1.0f32; nx * ny];

        jacobi_iteration_simd(&mut phi, &mut phi_new, &source, nx, ny, dx, dy)
            .expect("a 10x10 grid admits the interior sweep");

        // Four roundings separate this from the exact value (two reciprocals,
        // their sum, the product), so 8 f32 ulp is a generous bound.
        let expected = -0.5 / (1.0 / (dx * dx) + 1.0 / (dy * dy));
        let tolerance = 8.0 * f32::EPSILON * expected.abs();
        for i in 0..nx {
            for j in 0..ny {
                let got = phi_new[i * ny + j];
                if i == 0 || j == 0 || i == nx - 1 || j == ny - 1 {
                    assert_eq!(got, 0.0, "boundary ({i},{j}) must carry phi through");
                } else {
                    assert!(
                        (got - expected).abs() <= tolerance,
                        "interior ({i},{j}): {got} is not the scaled source {expected}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_jacobi_preserves_boundaries() {
        let nx = 6;
        let ny = 6;
        let n = nx * ny;
        let mut phi = vec![0.0f32; n];
        let mut phi_new = vec![0.0f32; n];
        let source = vec![1.0f32; n];

        // Set boundary values
        for i in 0..nx {
            phi[i * ny] = 10.0;
            phi[i * ny + ny - 1] = 20.0;
        }
        for j in 0..ny {
            phi[j] = 30.0;
            phi[(nx - 1) * ny + j] = 40.0;
        }

        jacobi_iteration_simd(&mut phi, &mut phi_new, &source, nx, ny, 0.1, 0.1)
            .expect("SIMD execution failed dimension check or mathematical constraint");

        // Boundaries must be preserved
        for i in 0..nx {
            assert_eq!(phi_new[i * ny], phi[i * ny], "left boundary at i={i}");
            assert_eq!(
                phi_new[i * ny + ny - 1],
                phi[i * ny + ny - 1],
                "right boundary at i={i}"
            );
        }
        for j in 0..ny {
            assert_eq!(phi_new[j], phi[j], "bottom boundary at j={j}");
            assert_eq!(
                phi_new[(nx - 1) * ny + j],
                phi[(nx - 1) * ny + j],
                "top boundary at j={j}"
            );
        }
    }

    #[test]
    fn test_jacobi_known_stencil() {
        // For a known interior point, verify the 5-point Jacobi stencil
        let nx = 4;
        let ny = 4;
        let dx = 1.0f32;
        let dy = 1.0f32;
        let n = nx * ny;
        let mut phi = vec![0.0f32; n];
        let mut phi_new = vec![0.0f32; n];
        let source = vec![0.0f32; n];

        // Set specific neighbors for cell (1,1)
        phi[1] = 1.0; // left
        phi[2 * ny + 1] = 3.0; // right
        phi[ny] = 2.0; // bottom
        phi[ny + 2] = 4.0; // top

        jacobi_iteration_simd(&mut phi, &mut phi_new, &source, nx, ny, dx, dy)
            .expect("SIMD execution failed dimension check or mathematical constraint");

        // factor = 0.5 / (1/dx² + 1/dy²) = 0.5 / 2 = 0.25
        // phi_new[1,1] = 0.25 * ((1+3)/1 + (2+4)/1 - 0) = 0.25 * 10 = 2.5
        let idx = ny + 1;
        assert!(
            (phi_new[idx] - 2.5).abs() < 1e-6,
            "Jacobi stencil: expected 2.5, got {}",
            phi_new[idx]
        );
    }

    #[test]
    fn test_gauss_seidel_convergence() {
        // Solve ∇²φ = -1 on [0,1]² with φ=0 on boundary
        // After multiple iterations, interior values should be positive
        let nx = 8;
        let ny = 8;
        let dx = 1.0 / (nx as f32 - 1.0);
        let dy = 1.0 / (ny as f32 - 1.0);
        let n = nx * ny;
        let mut phi = vec![0.0f32; n];
        let source = vec![-1.0f32; n];

        for _ in 0..100 {
            gauss_seidel_simd(&mut phi, &source, nx, ny, dx, dy, 1.5)
                .expect("SIMD execution failed dimension check or mathematical constraint");
        }

        // Center point should be positive (concave down solution)
        let center = (nx / 2) * ny + ny / 2;
        assert!(
            phi[center] > 0.0,
            "Center should be positive for -∇²φ = 1 with zero BCs, got {}",
            phi[center]
        );
    }

    #[test]
    fn test_divergence_uniform_field() {
        // Uniform velocity field should have zero divergence
        let nx = 5;
        let ny = 5;
        let u = vec![3.0f32; nx * ny]; // constant u
        let v = vec![7.0f32; nx * ny]; // constant v
        let mut divergence = vec![999.0f32; nx * ny];

        calculate_divergence_simd(&u, &v, &mut divergence, nx, ny, 1.0, 1.0)
            .expect("SIMD execution failed dimension check or mathematical constraint");

        // Interior divergence should be zero for uniform field
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = i * ny + j;
                assert!(
                    divergence[idx].abs() < 1e-10,
                    "Divergence of uniform field should be zero at ({i},{j}), got {}",
                    divergence[idx]
                );
            }
        }
    }

    #[test]
    fn test_divergence_linear_field() {
        // u = x (du/dx = 1), v = 0 => div = 1
        let nx = 5;
        let ny = 5;
        let dx = 1.0f32;
        let mut u = vec![0.0f32; nx * ny];
        let v = vec![0.0f32; nx * ny];
        let mut divergence = vec![0.0f32; nx * ny];

        for i in 0..nx {
            for j in 0..ny {
                u[i * ny + j] = i as f32;
            }
        }

        calculate_divergence_simd(&u, &v, &mut divergence, nx, ny, dx, 1.0)
            .expect("SIMD execution failed dimension check or mathematical constraint");

        // Interior div = du/dx = 1 (central difference: (x+1 - (x-1)) / (2*dx) = 1)
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = i * ny + j;
                assert!(
                    (divergence[idx] - 1.0).abs() < 1e-6,
                    "Divergence of u=x should be 1.0 at ({i},{j}), got {}",
                    divergence[idx]
                );
            }
        }
    }

    #[test]
    fn test_gradient_constant_field() {
        // Gradient of constant field should be zero
        let nx = 5;
        let ny = 5;
        let phi = vec![42.0f32; nx * ny];
        let mut grad_x = vec![999.0f32; nx * ny];
        let mut grad_y = vec![999.0f32; nx * ny];

        calculate_gradient_simd(&phi, &mut grad_x, &mut grad_y, nx, ny, 1.0, 1.0)
            .expect("SIMD execution failed dimension check or mathematical constraint");

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = i * ny + j;
                assert!(
                    grad_x[idx].abs() < 1e-10,
                    "grad_x of constant should be 0 at ({i},{j})"
                );
                assert!(
                    grad_y[idx].abs() < 1e-10,
                    "grad_y of constant should be 0 at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn test_gradient_linear_field() {
        // φ = 2x + 3y => ∂φ/∂x = 2, ∂φ/∂y = 3
        let nx = 5;
        let ny = 5;
        let dx = 1.0f32;
        let dy = 1.0f32;
        let mut phi = vec![0.0f32; nx * ny];

        for i in 0..nx {
            for j in 0..ny {
                phi[i * ny + j] = 2.0 * i as f32 + 3.0 * j as f32;
            }
        }

        let mut grad_x = vec![0.0f32; nx * ny];
        let mut grad_y = vec![0.0f32; nx * ny];

        calculate_gradient_simd(&phi, &mut grad_x, &mut grad_y, nx, ny, dx, dy)
            .expect("SIMD execution failed dimension check or mathematical constraint");

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = i * ny + j;
                assert!(
                    (grad_x[idx] - 2.0).abs() < 1e-6,
                    "grad_x of 2x+3y should be 2.0 at ({i},{j}), got {}",
                    grad_x[idx]
                );
                assert!(
                    (grad_y[idx] - 3.0).abs() < 1e-6,
                    "grad_y of 2x+3y should be 3.0 at ({i},{j}), got {}",
                    grad_y[idx]
                );
            }
        }
    }

    #[test]
    fn test_residual_zero_for_exact_solution() {
        // If φ satisfies ∇²φ = f exactly, residual should be zero
        // φ = 0 everywhere, f = 0 => residual = 0
        let nx = 5;
        let ny = 5;
        let phi = vec![0.0f32; nx * ny];
        let source = vec![0.0f32; nx * ny];
        let mut residual = vec![999.0f32; nx * ny];

        let max_r = calculate_residual_simd(&phi, &source, &mut residual, nx, ny, 1.0, 1.0)
            .expect("SIMD execution failed dimension check or mathematical constraint");

        assert!(
            max_r < 1e-10,
            "Residual of zero should be zero, got {max_r}"
        );
    }

    #[test]
    fn test_velocity_interpolation_uniform() {
        // Uniform field: face values should equal cell values
        let nx = 4;
        let ny = 4;
        let n = nx * ny;
        let u_cell = vec![5.0f32; n];
        let v_cell = vec![7.0f32; n];
        let mut u_face = vec![0.0f32; n];
        let mut v_face = vec![0.0f32; nx * (ny - 1)];

        interpolate_velocity_simd(&u_cell, &v_cell, &mut u_face, &mut v_face, nx, ny)
            .expect("SIMD execution failed dimension check or mathematical constraint");

        // For uniform field, face values = cell values
        for i in 0..nx - 1 {
            for j in 0..ny {
                let idx = i * ny + j;
                assert!(
                    (u_face[idx] - 5.0).abs() < 1e-6,
                    "Uniform u face interpolation at ({i},{j})"
                );
            }
        }
    }

    // ── Differential suite: ported kernels vs the historical scalar form ──
    //
    // The references below are the original scalar loop bodies, verbatim in
    // their operation order. Kernels whose scalar form precomputes reciprocal
    // spacings (divergence, gradient, residual, interpolation) must agree
    // bit-for-bit; Jacobi and Gauss-Seidel divide by dx²/dy² directly in the
    // scalar form while the port multiplies by the reciprocal, a difference
    // bounded by two ulp (reciprocal rounding + product rounding vs one
    // direct division).

    /// Deterministic LCG in [-1, 1).
    struct Lcg(u64);

    impl Lcg {
        fn next_f32(&mut self) -> f32 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let v = ((self.0 >> 40) as f32) / ((1u64 << 24) as f32);
            2.0 * v - 1.0
        }
    }

    fn fill(n: usize, rng: &mut Lcg) -> Vec<f32> {
        (0..n).map(|_| rng.next_f32()).collect()
    }

    /// Two-ulp absolute bound for `expected` (subnormals magnitudes < 1).
    fn two_ulp(expected: f32) -> f32 {
        2.0 * f32::EPSILON * expected.abs().max(1.0)
    }

    /// The historical Jacobi body, verbatim.
    fn jacobi_reference(
        phi: &[f32],
        phi_new: &mut [f32],
        source: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
    ) {
        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let factor = 0.5 / (1.0 / dx2 + 1.0 / dy2);

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = i * ny + j;
                let left = phi[(i - 1) * ny + j];
                let right = phi[(i + 1) * ny + j];
                let bottom = phi[i * ny + j - 1];
                let top = phi[i * ny + j + 1];
                let laplacian_x = (left + right) / dx2;
                let laplacian_y = (bottom + top) / dy2;
                phi_new[idx] = factor * (laplacian_x + laplacian_y - source[idx]);
            }
        }
        for i in 0..nx {
            phi_new[i * ny] = phi[i * ny];
            phi_new[i * ny + ny - 1] = phi[i * ny + ny - 1];
        }
        for j in 0..ny {
            phi_new[j] = phi[j];
            phi_new[(nx - 1) * ny + j] = phi[(nx - 1) * ny + j];
        }
    }

    /// The historical red-black SOR body, verbatim.
    fn gauss_seidel_reference(
        phi: &mut [f32],
        source: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        omega: f32,
    ) {
        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let factor = omega / (2.0 * (1.0 / dx2 + 1.0 / dy2));

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                if (i + j) % 2 == 0 {
                    let idx = i * ny + j;
                    let residual = (phi[(i - 1) * ny + j] + phi[(i + 1) * ny + j]) / dx2
                        + (phi[i * ny + j - 1] + phi[i * ny + j + 1]) / dy2
                        - source[idx];
                    phi[idx] = (1.0 - omega) * phi[idx] + factor * residual;
                }
            }
        }
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                if (i + j) % 2 == 1 {
                    let idx = i * ny + j;
                    let residual = (phi[(i - 1) * ny + j] + phi[(i + 1) * ny + j]) / dx2
                        + (phi[i * ny + j - 1] + phi[i * ny + j + 1]) / dy2
                        - source[idx];
                    phi[idx] = (1.0 - omega) * phi[idx] + factor * residual;
                }
            }
        }
    }

    /// The historical divergence body, verbatim.
    fn divergence_reference(
        u: &[f32],
        v: &[f32],
        divergence: &mut [f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
    ) {
        let inv_dx = 1.0 / dx;
        let inv_dy = 1.0 / dy;
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = i * ny + j;
                let dudx = (u[(i + 1) * ny + j] - u[(i - 1) * ny + j]) * 0.5 * inv_dx;
                let dvdy = (v[i * ny + j + 1] - v[i * ny + j - 1]) * 0.5 * inv_dy;
                divergence[idx] = dudx + dvdy;
            }
        }
    }

    /// The historical gradient body, verbatim.
    fn gradient_reference(
        phi: &[f32],
        grad_x: &mut [f32],
        grad_y: &mut [f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
    ) {
        let inv_dx = 0.5 / dx;
        let inv_dy = 0.5 / dy;
        for i in 1..nx - 1 {
            let row_start = i * ny;
            for j in 1..ny - 1 {
                let idx = row_start + j;
                grad_x[idx] = (phi[(i + 1) * ny + j] - phi[(i - 1) * ny + j]) * inv_dx;
                grad_y[idx] = (phi[i * ny + j + 1] - phi[i * ny + j - 1]) * inv_dy;
            }
        }
    }

    /// The historical residual body, verbatim.
    fn residual_reference(
        phi: &[f32],
        source: &[f32],
        residual: &mut [f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
    ) -> f32 {
        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let inv_dx2 = 1.0 / dx2;
        let inv_dy2 = 1.0 / dy2;
        let mut max_residual = 0.0f32;
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = i * ny + j;
                let laplacian_x =
                    (phi[(i - 1) * ny + j] - 2.0 * phi[idx] + phi[(i + 1) * ny + j]) * inv_dx2;
                let laplacian_y =
                    (phi[i * ny + j - 1] - 2.0 * phi[idx] + phi[i * ny + j + 1]) * inv_dy2;
                let laplacian = laplacian_x + laplacian_y;
                residual[idx] = laplacian - source[idx];
                max_residual = max_residual.max(residual[idx].abs());
            }
        }
        max_residual
    }

    /// The historical interpolation body, verbatim.
    fn interpolation_reference(
        u_cell: &[f32],
        v_cell: &[f32],
        u_face: &mut [f32],
        v_face: &mut [f32],
        nx: usize,
        ny: usize,
    ) {
        for i in 0..nx - 1 {
            let start_idx = i * ny;
            let end_idx = start_idx + ny;
            let left = &u_cell[start_idx..end_idx];
            let right = &u_cell[start_idx + ny..end_idx + ny];
            let face = &mut u_face[start_idx..end_idx];
            for k in 0..ny {
                face[k] = (left[k] + right[k]) * 0.5;
            }
        }
        let mut temp_face = vec![0.0f32; nx];
        for j in 0..ny - 1 {
            for i in 0..nx {
                temp_face[i] = (v_cell[i * ny + j] + v_cell[i * ny + j + 1]) * 0.5;
            }
            for i in 0..nx {
                v_face[i * (ny - 1) + j] = temp_face[i];
            }
        }
    }

    #[test]
    fn differential_jacobi_within_two_ulp() {
        let mut rng = Lcg(0xA11CE);
        for &(nx, ny) in &[
            (3usize, 3usize),
            (4, 7),
            (7, 4),
            (13, 11),
            (16, 16),
            (33, 17),
            (2, 5),
            (5, 2),
            (1, 4),
            (4, 1),
        ] {
            let n = nx * ny;
            let phi = fill(n, &mut rng);
            let source = fill(n, &mut rng);
            for &(dx, dy) in &[(0.1f32, 0.1f32), (1.0, 0.5), (0.25, 2.0)] {
                let mut got = vec![0.0f32; n];
                let mut want = vec![0.0f32; n];

                jacobi_iteration_simd(&mut phi.clone(), &mut got, &source, nx, ny, dx, dy)
                    .expect("port ok");
                jacobi_reference(&phi, &mut want, &source, nx, ny, dx, dy);

                for k in 0..n {
                    assert!(
                        (got[k] - want[k]).abs() <= two_ulp(want[k]),
                        "jacobi mismatch at {nx}x{ny} dx={dx} dy={dy} idx {k}: {} vs {}",
                        got[k],
                        want[k]
                    );
                }
            }
        }
    }

    #[test]
    fn differential_gauss_seidel_within_two_ulp() {
        let mut rng = Lcg(0xB0B);
        for &(nx, ny) in &[(3usize, 3usize), (6, 6), (9, 12), (17, 17), (2, 4), (4, 2)] {
            let n = nx * ny;
            let source = fill(n, &mut rng);
            for &omega in &[1.0f32, 1.5f32] {
                let mut got = fill(n, &mut rng);
                let mut want = got.clone();

                for sweep in 0..5 {
                    gauss_seidel_simd(&mut got, &source, nx, ny, 0.1, 0.1, omega).expect("port ok");
                    gauss_seidel_reference(&mut want, &source, nx, ny, 0.1, 0.1, omega);
                    for k in 0..n {
                        assert!(
                            (got[k] - want[k]).abs() <= two_ulp(want[k]),
                            "gauss-seidel mismatch at {nx}x{ny} omega={omega} sweep \
                             {sweep} idx {k}: {} vs {}",
                            got[k],
                            want[k]
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn differential_divergence_gradient_bit_exact() {
        let mut rng = Lcg(0xC0DE);
        for &(nx, ny) in &[(3usize, 3usize), (5, 9), (12, 5), (16, 16), (33, 9)] {
            let n = nx * ny;
            let u = fill(n, &mut rng);
            let v = fill(n, &mut rng);
            let phi = fill(n, &mut rng);
            let (dx, dy) = (0.1f32, 0.2f32);

            let mut got_div = vec![0.0f32; n];
            let mut want_div = vec![0.0f32; n];
            calculate_divergence_simd(&u, &v, &mut got_div, nx, ny, dx, dy).expect("port ok");
            divergence_reference(&u, &v, &mut want_div, nx, ny, dx, dy);
            assert_eq!(got_div, want_div, "divergence mismatch at {nx}x{ny}");

            let mut got_gx = vec![0.0f32; n];
            let mut got_gy = vec![0.0f32; n];
            let mut want_gx = vec![0.0f32; n];
            let mut want_gy = vec![0.0f32; n];
            calculate_gradient_simd(&phi, &mut got_gx, &mut got_gy, nx, ny, dx, dy)
                .expect("port ok");
            gradient_reference(&phi, &mut want_gx, &mut want_gy, nx, ny, dx, dy);
            assert_eq!(got_gx, want_gx, "grad_x mismatch at {nx}x{ny}");
            assert_eq!(got_gy, want_gy, "grad_y mismatch at {nx}x{ny}");
        }
    }

    #[test]
    fn differential_residual_bit_exact() {
        let mut rng = Lcg(0xE11A);
        for &(nx, ny) in &[(3usize, 3usize), (6, 6), (11, 14), (16, 16), (2, 3), (3, 2)] {
            let n = nx * ny;
            let phi = fill(n, &mut rng);
            let source = fill(n, &mut rng);
            let (dx, dy) = (0.1f32, 0.35f32);

            let mut got = vec![0.0f32; n];
            let mut want = vec![0.0f32; n];

            let got_max =
                calculate_residual_simd(&phi, &source, &mut got, nx, ny, dx, dy).expect("port ok");
            let want_max = residual_reference(&phi, &source, &mut want, nx, ny, dx, dy);

            assert_eq!(got, want, "residual field mismatch at {nx}x{ny}");
            assert_eq!(got_max, want_max, "residual max mismatch at {nx}x{ny}");
        }
    }

    #[test]
    fn differential_interpolation_bit_exact() {
        let mut rng = Lcg(0xF00D);
        for &(nx, ny) in &[(2usize, 2usize), (3, 4), (5, 5), (9, 6), (17, 9)] {
            let n = nx * ny;
            let u_cell = fill(n, &mut rng);
            let v_cell = fill(n, &mut rng);

            let mut got_uf = vec![0.0f32; n];
            let mut got_vf = vec![0.0f32; nx * (ny - 1)];
            let mut want_uf = vec![0.0f32; n];
            let mut want_vf = vec![0.0f32; nx * (ny - 1)];

            interpolate_velocity_simd(&u_cell, &v_cell, &mut got_uf, &mut got_vf, nx, ny)
                .expect("port ok");
            interpolation_reference(&u_cell, &v_cell, &mut want_uf, &mut want_vf, nx, ny);

            assert_eq!(got_uf, want_uf, "u-face mismatch at {nx}x{ny}");
            assert_eq!(got_vf, want_vf, "v-face mismatch at {nx}x{ny}");
        }
    }

    #[test]
    fn differential_gauss_seidel_solves_poisson_consistently() {
        // Beyond tolerance: the ported SOR must still converge on the
        // canonical problem (mirrors test_gauss_seidel_convergence on a
        // larger grid with the same tolerance philosophy).
        let nx = 16;
        let ny = 16;
        let dx = 1.0 / (nx as f32 - 1.0);
        let n = nx * ny;
        let mut phi = vec![0.0f32; n];
        let source = vec![-1.0f32; n];

        for _ in 0..200 {
            gauss_seidel_simd(&mut phi, &source, nx, ny, dx, dx, 1.7).expect("port ok");
        }
        let center = (nx / 2) * ny + ny / 2;
        assert!(phi[center] > 0.0, "SOR diverged: center = {}", phi[center]);
    }
}
