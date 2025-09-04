//! x86/x86_64 SIMD implementations using AVX2 and SSE4.1

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// AVX2 implementation for advection kernel (256-bit vectors, 8 floats)
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
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
    let spacing_x_vec = _mm256_set1_ps(dx);
    let spacing_y_vec = _mm256_set1_ps(dy);
    let time_step_vec = _mm256_set1_ps(dt);
    let zero = _mm256_setzero_ps();

    for j in 1..ny - 1 {
        let mut i = 1;

        // Process 8 elements at a time
        while i + 7 < nx - 1 {
            let idx = j * nx + i;

            // Load current values
            let u = _mm256_loadu_ps(&input[idx]);
            let vx = _mm256_loadu_ps(&velocity_x[idx]);
            let vy = _mm256_loadu_ps(&velocity_y[idx]);

            // Compute upwind differences for x-direction
            let u_left = _mm256_loadu_ps(&input[idx - 1]);
            let u_right = _mm256_loadu_ps(&input[idx + 1]);

            // vx > 0 ? (u - u_left) : (u_right - u)
            let mask_x = _mm256_cmp_ps(vx, zero, _CMP_GT_OQ);
            let gradient_x_pos = _mm256_sub_ps(u, u_left);
            let gradient_x_neg = _mm256_sub_ps(u_right, u);
            let gradient_x = _mm256_blendv_ps(gradient_x_neg, gradient_x_pos, mask_x);
            let gradient_x = _mm256_div_ps(gradient_x, spacing_x_vec);

            // Compute upwind differences for y-direction
            let u_bottom = _mm256_loadu_ps(&input[idx - nx]);
            let u_top = _mm256_loadu_ps(&input[idx + nx]);

            // vy > 0 ? (u - u_bottom) : (u_top - u)
            let mask_y = _mm256_cmp_ps(vy, zero, _CMP_GT_OQ);
            let gradient_y_pos = _mm256_sub_ps(u, u_bottom);
            let gradient_y_neg = _mm256_sub_ps(u_top, u);
            let gradient_y = _mm256_blendv_ps(gradient_y_neg, gradient_y_pos, mask_y);
            let gradient_y = _mm256_div_ps(gradient_y, spacing_y_vec);

            // Compute advection term: u - dt * (vx * gradient_x + vy * gradient_y)
            let advection_x = _mm256_mul_ps(vx, gradient_x);
            let advection_y = _mm256_mul_ps(vy, gradient_y);
            let advection = _mm256_add_ps(advection_x, advection_y);
            let advection = _mm256_mul_ps(advection, time_step_vec);
            let result = _mm256_sub_ps(u, advection);

            // Store result
            _mm256_storeu_ps(&mut output[idx], result);

            i += 8;
        }

        // Handle remaining elements with scalar code
        while i < nx - 1 {
            let idx = j * nx + i;
            let vx = velocity_x[idx];
            let vy = velocity_y[idx];

            let du_dx = if vx > 0.0 {
                (input[idx] - input[idx - 1]) / dx
            } else {
                (input[idx + 1] - input[idx]) / dx
            };

            let du_dy = if vy > 0.0 {
                (input[idx] - input[idx - nx]) / dy
            } else {
                (input[idx + nx] - input[idx]) / dy
            };

            output[idx] = input[idx] - dt * (vx * du_dx + vy * du_dy);
            i += 1;
        }
    }
}

/// SSE4.1 implementation for advection kernel (128-bit vectors, 4 floats)
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
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
    let spacing_x_vec = _mm_set1_ps(dx);
    let spacing_y_vec = _mm_set1_ps(dy);
    let time_step_vec = _mm_set1_ps(dt);
    let zero = _mm_setzero_ps();

    for j in 1..ny - 1 {
        let mut i = 1;

        // Process 4 elements at a time
        while i + 3 < nx - 1 {
            let idx = j * nx + i;

            // Load current values
            let u = _mm_loadu_ps(&input[idx]);
            let vx = _mm_loadu_ps(&velocity_x[idx]);
            let vy = _mm_loadu_ps(&velocity_y[idx]);

            // Compute upwind differences for x-direction
            let u_left = _mm_loadu_ps(&input[idx - 1]);
            let u_right = _mm_loadu_ps(&input[idx + 1]);

            // vx > 0 ? (u - u_left) : (u_right - u)
            let mask_x = _mm_cmpgt_ps(vx, zero);
            let gradient_x_pos = _mm_sub_ps(u, u_left);
            let gradient_x_neg = _mm_sub_ps(u_right, u);
            let gradient_x = _mm_blendv_ps(gradient_x_neg, gradient_x_pos, mask_x);
            let gradient_x = _mm_div_ps(gradient_x, spacing_x_vec);

            // Compute upwind differences for y-direction
            let u_bottom = _mm_loadu_ps(&input[idx - nx]);
            let u_top = _mm_loadu_ps(&input[idx + nx]);

            // vy > 0 ? (u - u_bottom) : (u_top - u)
            let mask_y = _mm_cmpgt_ps(vy, zero);
            let gradient_y_pos = _mm_sub_ps(u, u_bottom);
            let gradient_y_neg = _mm_sub_ps(u_top, u);
            let gradient_y = _mm_blendv_ps(gradient_y_neg, gradient_y_pos, mask_y);
            let gradient_y = _mm_div_ps(gradient_y, spacing_y_vec);

            // Compute advection term
            let advection_x = _mm_mul_ps(vx, gradient_x);
            let advection_y = _mm_mul_ps(vy, gradient_y);
            let advection = _mm_add_ps(advection_x, advection_y);
            let advection = _mm_mul_ps(advection, time_step_vec);
            let result = _mm_sub_ps(u, advection);

            // Store result
            _mm_storeu_ps(&mut output[idx], result);

            i += 4;
        }

        // Handle remaining elements with scalar code
        while i < nx - 1 {
            let idx = j * nx + i;
            let vx = velocity_x[idx];
            let vy = velocity_y[idx];

            let du_dx = if vx > 0.0 {
                (input[idx] - input[idx - 1]) / dx
            } else {
                (input[idx + 1] - input[idx]) / dx
            };

            let du_dy = if vy > 0.0 {
                (input[idx] - input[idx - nx]) / dy
            } else {
                (input[idx + nx] - input[idx]) / dy
            };

            output[idx] = input[idx] - dt * (vx * du_dx + vy * du_dy);
            i += 1;
        }
    }
}

/// AVX2 implementation for diffusion kernel
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
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
    let dx2_inv = 1.0 / (dx * dx);
    let dy2_inv = 1.0 / (dy * dy);
    let dt_nu = dt * nu;

    let dx2_inv_vec = _mm256_set1_ps(dx2_inv);
    let dy2_inv_vec = _mm256_set1_ps(dy2_inv);
    let dt_nu_vec = _mm256_set1_ps(dt_nu);
    let two = _mm256_set1_ps(2.0);

    for j in 1..ny - 1 {
        let mut i = 1;

        // Process 8 elements at a time
        while i + 7 < nx - 1 {
            let idx = j * nx + i;

            // Load stencil values
            let u_center = _mm256_loadu_ps(&input[idx]);
            let u_left = _mm256_loadu_ps(&input[idx - 1]);
            let u_right = _mm256_loadu_ps(&input[idx + 1]);
            let u_bottom = _mm256_loadu_ps(&input[idx - nx]);
            let u_top = _mm256_loadu_ps(&input[idx + nx]);

            // Compute Laplacian: (u_left - 2*u + u_right)/dx^2 + (u_bottom - 2*u + u_top)/dy^2
            let d2u_dx2 = _mm256_sub_ps(u_left, _mm256_mul_ps(two, u_center));
            let d2u_dx2 = _mm256_add_ps(d2u_dx2, u_right);
            let d2u_dx2 = _mm256_mul_ps(d2u_dx2, dx2_inv_vec);

            let d2u_dy2 = _mm256_sub_ps(u_bottom, _mm256_mul_ps(two, u_center));
            let d2u_dy2 = _mm256_add_ps(d2u_dy2, u_top);
            let d2u_dy2 = _mm256_mul_ps(d2u_dy2, dy2_inv_vec);

            let laplacian = _mm256_add_ps(d2u_dx2, d2u_dy2);

            // Update: u + dt * nu * laplacian
            let update = _mm256_mul_ps(laplacian, dt_nu_vec);
            let result = _mm256_add_ps(u_center, update);

            // Store result
            _mm256_storeu_ps(&mut output[idx], result);

            i += 8;
        }

        // Handle remaining elements
        while i < nx - 1 {
            let idx = j * nx + i;

            let d2u_dx2 = (input[idx - 1] - 2.0 * input[idx] + input[idx + 1]) * dx2_inv;
            let d2u_dy2 = (input[idx - nx] - 2.0 * input[idx] + input[idx + nx]) * dy2_inv;
            let laplacian = d2u_dx2 + d2u_dy2;

            output[idx] = input[idx] + dt_nu * laplacian;
            i += 1;
        }
    }
}
