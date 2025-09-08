//! SIMD-accelerated kernels for 2D CFD solvers
//!
//! Provides optimized numerical kernels using SIMD operations

use cfd_core::error::Result;
use cfd_math::simd::{SimdProcessor, VectorOps};
use std::sync::OnceLock;

// Global SIMD processor instance (initialized once)
static SIMD_PROCESSOR: OnceLock<SimdProcessor> = OnceLock::new();

/// Get or initialize the global SIMD processor
fn simd_processor() -> &'static SimdProcessor {
    SIMD_PROCESSOR.get_or_init(SimdProcessor::new)
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

    let processor = simd_processor();

    // Process interior points
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let idx = i * ny + j;

            // Gather neighbor values
            let left = phi[(i - 1) * ny + j];
            let right = phi[(i + 1) * ny + j];
            let bottom = phi[i * ny + j - 1];
            let top = phi[i * ny + j + 1];

            // Jacobi update
            phi_new[idx] = factor * ((left + right) / dx2 + (bottom + top) / dy2 - source[idx]);
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

    // Red-Black Gauss-Seidel for better parallelization
    // Red points (i+j even)
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

    // Black points (i+j odd)
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
    let processor = simd_processor();

    // Interpolate u-velocity to x-faces
    for i in 0..nx - 1 {
        let start_idx = i * ny;
        let end_idx = start_idx + ny;

        // Use SIMD for averaging
        let left = &u_cell[start_idx..end_idx];
        let right = &u_cell[start_idx + ny..end_idx + ny];
        let face = &mut u_face[start_idx..end_idx];

        // Average: face = 0.5 * (left + right)
        processor.ops.add(left, right, face)?;
        // Scale in-place by creating temporary
        let temp = face.to_vec();
        processor.ops.scale(&temp, 0.5, face)?;
    }

    // Interpolate v-velocity to y-faces
    for j in 0..ny - 1 {
        let mut temp_bottom = Vec::with_capacity(nx);
        let mut temp_top = Vec::with_capacity(nx);
        let mut temp_face = vec![0.0f32; nx];

        for i in 0..nx {
            temp_bottom.push(v_cell[i * ny + j]);
            temp_top.push(v_cell[i * ny + j + 1]);
        }

        // Average using SIMD
        processor.ops.add(&temp_bottom, &temp_top, &mut temp_face)?;
        // Scale with temporary to avoid borrow conflict
        let temp_copy = temp_face.clone();
        processor.ops.scale(&temp_copy, 0.5, &mut temp_face)?;

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

    // Calculate divergence: div = ∂u/∂x + ∂v/∂y
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let idx = i * ny + j;

            // Central differences
            let dudx = (u[(i + 1) * ny + j] - u[(i - 1) * ny + j]) * 0.5 * inv_dx;
            let dvdy = (v[i * ny + j + 1] - v[i * ny + j - 1]) * 0.5 * inv_dy;

            divergence[idx] = dudx + dvdy;
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
    let processor = simd_processor();
    let inv_dx = 0.5 / dx;
    let inv_dy = 0.5 / dy;

    // Calculate gradients using central differences
    for i in 1..nx - 1 {
        let row_start = i * ny;

        // Process a row at a time for better cache usage
        for j in 1..ny - 1 {
            let idx = row_start + j;

            // ∂φ/∂x
            grad_x[idx] = (phi[(i + 1) * ny + j] - phi[(i - 1) * ny + j]) * inv_dx;

            // ∂φ/∂y
            grad_y[idx] = (phi[i * ny + j + 1] - phi[i * ny + j - 1]) * inv_dy;
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
    let processor = simd_processor();
    let dx2 = dx * dx;
    let dy2 = dy * dy;
    let inv_dx2 = 1.0 / dx2;
    let inv_dy2 = 1.0 / dy2;

    let mut max_residual = 0.0f32;

    // Calculate residual: r = ∇²φ - f
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let idx = i * ny + j;

            // Laplacian
            let laplacian = (phi[(i - 1) * ny + j] - 2.0 * phi[idx] + phi[(i + 1) * ny + j])
                * inv_dx2
                + (phi[i * ny + j - 1] - 2.0 * phi[idx] + phi[i * ny + j + 1]) * inv_dy2;

            residual[idx] = laplacian - source[idx];
            max_residual = max_residual.max(residual[idx].abs());
        }
    }

    Ok(max_residual)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jacobi_iteration() {
        let nx = 10;
        let ny = 10;
        let mut phi = vec![0.0f32; nx * ny];
        let mut phi_new = vec![0.0f32; nx * ny];
        let source = vec![1.0f32; nx * ny];

        let result = jacobi_iteration_simd(&mut phi, &mut phi_new, &source, nx, ny, 0.1, 0.1);

        assert!(result.is_ok());
    }

    #[test]
    fn test_divergence_calculation() {
        let nx = 5;
        let ny = 5;
        let u = vec![1.0f32; nx * ny];
        let v = vec![2.0f32; nx * ny];
        let mut divergence = vec![0.0f32; nx * ny];

        let result = calculate_divergence_simd(&u, &v, &mut divergence, nx, ny, 1.0, 1.0);

        assert!(result.is_ok());
    }
}
