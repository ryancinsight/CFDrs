//! `x86/x86_64` SIMD implementations using AVX2 and SSE4.1

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    _mm256_add_ps, _mm256_blendv_ps, _mm256_cmp_ps, _mm256_div_ps, _mm256_loadu_ps, _mm256_mul_ps,
    _mm256_set1_ps, _mm256_setzero_ps, _mm256_storeu_ps, _mm256_sub_ps, _mm_add_ps, _mm_blendv_ps,
    _mm_cmpgt_ps, _mm_div_ps, _mm_loadu_ps, _mm_mul_ps, _mm_set1_ps, _mm_setzero_ps, _mm_storeu_ps,
    _mm_sub_ps, _CMP_GT_OQ,
};

/// AVX2 implementation for advection kernel (256-bit vectors, 8 floats)
///
/// # Safety
///
/// This function requires AVX2 instruction set support. Caller must ensure
/// the target hardware supports AVX2 and all slice inputs have compatible lengths
/// for the given grid dimensions (nx, ny).
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
pub unsafe fn advection_avx2(
    input: &[f32],
    output: &mut [f32],
    velocity_x: &[f32],
    velocity_y: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
    dt: f32,
) {
    // Process 8 elements at a time with AVX2
    let x_spacing_vec = _mm256_set1_ps(dx);
    let y_spacing_vec = _mm256_set1_ps(dy);
    let time_step_vec = _mm256_set1_ps(dt);
    let zero = _mm256_setzero_ps();

    for j in 1..ny - 1 {
        let mut i = 1;

        // Process 8 elements at a time
        while i + 7 < nx - 1 {
            let idx = j * nx + i;

            // Load current values
            let u = _mm256_loadu_ps(&raw const input[idx]);
            let vx = _mm256_loadu_ps(&raw const velocity_x[idx]);
            let vy = _mm256_loadu_ps(&raw const velocity_y[idx]);

            // Compute upwind differences for x-direction
            let u_left = _mm256_loadu_ps(&raw const input[idx - 1]);
            let u_right = _mm256_loadu_ps(&raw const input[idx + 1]);

            // vx > 0 ? (u - u_left) : (u_right - u)
            let mask_x = _mm256_cmp_ps(vx, zero, _CMP_GT_OQ);
            let x_gradient_forward = _mm256_sub_ps(u, u_left);
            let x_gradient_backward = _mm256_sub_ps(u_right, u);
            let gradient_x = _mm256_blendv_ps(x_gradient_backward, x_gradient_forward, mask_x);
            let gradient_x = _mm256_div_ps(gradient_x, x_spacing_vec);

            // Compute upwind differences for y-direction
            let u_bottom = _mm256_loadu_ps(&raw const input[idx - nx]);
            let u_top = _mm256_loadu_ps(&raw const input[idx + nx]);

            // vy > 0 ? (u - u_bottom) : (u_top - u)
            let mask_y = _mm256_cmp_ps(vy, zero, _CMP_GT_OQ);
            let y_gradient_forward = _mm256_sub_ps(u, u_bottom);
            let y_gradient_backward = _mm256_sub_ps(u_top, u);
            let gradient_y = _mm256_blendv_ps(y_gradient_backward, y_gradient_forward, mask_y);
            let gradient_y = _mm256_div_ps(gradient_y, y_spacing_vec);

            // Compute advection term: u - dt * (vx * gradient_x + vy * gradient_y)
            let advection_x = _mm256_mul_ps(vx, gradient_x);
            let advection_y = _mm256_mul_ps(vy, gradient_y);
            let advection = _mm256_add_ps(advection_x, advection_y);
            let dt_advection = _mm256_mul_ps(time_step_vec, advection);
            let result = _mm256_sub_ps(u, dt_advection);

            // Store result
            _mm256_storeu_ps(&raw mut output[idx], result);

            i += 8;
        }

        // Handle remaining elements with scalar code
        while i < nx - 1 {
            let idx = j * nx + i;
            let vx = velocity_x[idx];
            let vy = velocity_y[idx];

            let grad_x = if vx > 0.0 {
                (input[idx] - input[idx - 1]) / dx
            } else {
                (input[idx + 1] - input[idx]) / dx
            };

            let grad_y = if vy > 0.0 {
                (input[idx] - input[idx - nx]) / dy
            } else {
                (input[idx + nx] - input[idx]) / dy
            };

            output[idx] = input[idx] - dt * (vx * grad_x + vy * grad_y);
            i += 1;
        }
    }
}

/// SSE4.1 implementation for advection kernel (128-bit vectors, 4 floats)
///
/// # Safety
///
/// This function requires SSE4.1 instruction set support. Caller must ensure
/// the target hardware supports SSE4.1 and all slice inputs have compatible lengths
/// for the given grid dimensions (nx, ny).
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
#[allow(clippy::too_many_arguments)]
pub unsafe fn advection_sse41(
    input: &[f32],
    output: &mut [f32],
    velocity_x: &[f32],
    velocity_y: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
    dt: f32,
) {
    // Process 4 elements at a time with SSE4.1
    let x_spacing_vec = _mm_set1_ps(dx);
    let y_spacing_vec = _mm_set1_ps(dy);
    let time_step_vec = _mm_set1_ps(dt);
    let zero = _mm_setzero_ps();

    for j in 1..ny - 1 {
        let mut i = 1;

        // Process 4 elements at a time
        while i + 3 < nx - 1 {
            let idx = j * nx + i;

            // Load current values
            let u = _mm_loadu_ps(&raw const input[idx]);
            let vx = _mm_loadu_ps(&raw const velocity_x[idx]);
            let vy = _mm_loadu_ps(&raw const velocity_y[idx]);

            // Compute upwind differences for x-direction
            let u_left = _mm_loadu_ps(&raw const input[idx - 1]);
            let u_right = _mm_loadu_ps(&raw const input[idx + 1]);

            // vx > 0 ? (u - u_left) : (u_right - u)
            let mask_x = _mm_cmpgt_ps(vx, zero);
            let x_gradient_forward = _mm_sub_ps(u, u_left);
            let x_gradient_backward = _mm_sub_ps(u_right, u);
            let gradient_x = _mm_blendv_ps(x_gradient_backward, x_gradient_forward, mask_x);
            let gradient_x = _mm_div_ps(gradient_x, x_spacing_vec);

            // Compute upwind differences for y-direction
            let u_bottom = _mm_loadu_ps(&raw const input[idx - nx]);
            let u_top = _mm_loadu_ps(&raw const input[idx + nx]);

            // vy > 0 ? (u - u_bottom) : (u_top - u)
            let mask_y = _mm_cmpgt_ps(vy, zero);
            let y_gradient_forward = _mm_sub_ps(u, u_bottom);
            let y_gradient_backward = _mm_sub_ps(u_top, u);
            let gradient_y = _mm_blendv_ps(y_gradient_backward, y_gradient_forward, mask_y);
            let gradient_y = _mm_div_ps(gradient_y, y_spacing_vec);

            // Compute advection term
            let advection_x = _mm_mul_ps(vx, gradient_x);
            let advection_y = _mm_mul_ps(vy, gradient_y);
            let advection = _mm_add_ps(advection_x, advection_y);
            let dt_advection = _mm_mul_ps(time_step_vec, advection);
            let result = _mm_sub_ps(u, dt_advection);

            // Store result
            _mm_storeu_ps(&raw mut output[idx], result);

            i += 4;
        }

        // Handle remaining elements with scalar code
        while i < nx - 1 {
            let idx = j * nx + i;
            let vx = velocity_x[idx];
            let vy = velocity_y[idx];

            let grad_x = if vx > 0.0 {
                (input[idx] - input[idx - 1]) / dx
            } else {
                (input[idx + 1] - input[idx]) / dx
            };

            let grad_y = if vy > 0.0 {
                (input[idx] - input[idx - nx]) / dy
            } else {
                (input[idx + nx] - input[idx]) / dy
            };

            output[idx] = input[idx] - dt * (vx * grad_x + vy * grad_y);
            i += 1;
        }
    }
}

/// AVX2 implementation for diffusion kernel
///
/// # Safety
///
/// This function requires AVX2 instruction set support. Caller must ensure
/// the target hardware supports AVX2 and all slice inputs have compatible lengths
/// for the given grid dimensions (nx, ny).
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
pub unsafe fn diffusion_avx2(
    input: &[f32],
    output: &mut [f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
    dt: f32,
    nu: f32, // Kinematic viscosity
) {
    let x_inv_squared = 1.0 / (dx * dx);
    let y_inv_squared = 1.0 / (dy * dy);
    let time_viscosity = dt * nu;

    let x_inv_sq_vec = _mm256_set1_ps(x_inv_squared);
    let y_inv_sq_vec = _mm256_set1_ps(y_inv_squared);
    let time_visc_vec = _mm256_set1_ps(time_viscosity);
    let two = _mm256_set1_ps(2.0);

    for j in 1..ny - 1 {
        let mut i = 1;

        // Process 8 elements at a time
        while i + 7 < nx - 1 {
            let idx = j * nx + i;

            // Load stencil values
            let u_center = _mm256_loadu_ps(&raw const input[idx]);
            let u_left = _mm256_loadu_ps(&raw const input[idx - 1]);
            let u_right = _mm256_loadu_ps(&raw const input[idx + 1]);
            let u_bottom = _mm256_loadu_ps(&raw const input[idx - nx]);
            let u_top = _mm256_loadu_ps(&raw const input[idx + nx]);

            // Compute Laplacian: (u_left - 2*u + u_right)/dx^2 + (u_bottom - 2*u + u_top)/dy^2
            let laplacian_x = _mm256_sub_ps(u_left, _mm256_mul_ps(two, u_center));
            let laplacian_x = _mm256_add_ps(laplacian_x, u_right);
            let laplacian_x = _mm256_mul_ps(laplacian_x, x_inv_sq_vec);

            let laplacian_y = _mm256_sub_ps(u_bottom, _mm256_mul_ps(two, u_center));
            let laplacian_y = _mm256_add_ps(laplacian_y, u_top);
            let laplacian_y = _mm256_mul_ps(laplacian_y, y_inv_sq_vec);

            let laplacian = _mm256_add_ps(laplacian_x, laplacian_y);

            // Update: u + dt * nu * laplacian
            let update = _mm256_mul_ps(laplacian, time_visc_vec);
            let result = _mm256_add_ps(u_center, update);

            // Store result
            _mm256_storeu_ps(&raw mut output[idx], result);

            i += 8;
        }

        // Handle remaining elements
        while i < nx - 1 {
            let idx = j * nx + i;

            let laplacian_x = (input[idx - 1] - 2.0 * input[idx] + input[idx + 1]) * x_inv_squared;
            let laplacian_y =
                (input[idx - nx] - 2.0 * input[idx] + input[idx + nx]) * y_inv_squared;
            let laplacian = laplacian_x + laplacian_y;

            output[idx] = input[idx] + time_viscosity * laplacian;
            i += 1;
        }
    }
}
