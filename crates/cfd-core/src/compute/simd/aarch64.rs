//! AArch64 SIMD implementations using NEON

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

/// NEON implementation for advection kernel (128-bit vectors, 4 floats)
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub unsafe fn advection_neon(
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
    // Process 4 elements at a time with NEON
    let dx_vec = vdupq_n_f32(dx);
    let dy_vec = vdupq_n_f32(dy);
    let dt_vec = vdupq_n_f32(dt);
    let zero = vdupq_n_f32(0.0);

    for j in 1..ny - 1 {
        let mut i = 1;

        // Process 4 elements at a time
        while i + 3 < nx - 1 {
            let idx = j * nx + i;

            // Load current values
            let u = vld1q_f32(&input[idx]);
            let vx = vld1q_f32(&velocity_x[idx]);
            let vy = vld1q_f32(&velocity_y[idx]);

            // Compute upwind differences for x-direction
            let u_left = vld1q_f32(&input[idx - 1]);
            let u_right = vld1q_f32(&input[idx + 1]);

            // vx > 0 ? (u - u_left) : (u_right - u)
            let mask_x = vcgtq_f32(vx, zero);
            let du_dx_pos = vsubq_f32(u, u_left);
            let du_dx_neg = vsubq_f32(u_right, u);
            let du_dx = vbslq_f32(mask_x, du_dx_pos, du_dx_neg);
            let du_dx = vdivq_f32(du_dx, dx_vec);

            // Compute upwind differences for y-direction
            let u_bottom = vld1q_f32(&input[idx - nx]);
            let u_top = vld1q_f32(&input[idx + nx]);

            // vy > 0 ? (u - u_bottom) : (u_top - u)
            let mask_y = vcgtq_f32(vy, zero);
            let du_dy_pos = vsubq_f32(u, u_bottom);
            let du_dy_neg = vsubq_f32(u_top, u);
            let du_dy = vbslq_f32(mask_y, du_dy_pos, du_dy_neg);
            let du_dy = vdivq_f32(du_dy, dy_vec);

            // Compute advection term: u - dt * (vx * du_dx + vy * du_dy)
            let advection_x = vmulq_f32(vx, du_dx);
            let advection_y = vmulq_f32(vy, du_dy);
            let advection = vaddq_f32(advection_x, advection_y);
            let advection = vmulq_f32(advection, dt_vec);
            let result = vsubq_f32(u, advection);

            // Store result
            vst1q_f32(&mut output[idx], result);

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

/// NEON implementation for diffusion kernel
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub unsafe fn diffusion_neon(
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

    let dx2_inv_vec = vdupq_n_f32(dx2_inv);
    let dy2_inv_vec = vdupq_n_f32(dy2_inv);
    let dt_nu_vec = vdupq_n_f32(dt_nu);
    let two = vdupq_n_f32(2.0);

    for j in 1..ny - 1 {
        let mut i = 1;

        // Process 4 elements at a time
        while i + 3 < nx - 1 {
            let idx = j * nx + i;

            // Load stencil values
            let u_center = vld1q_f32(&input[idx]);
            let u_left = vld1q_f32(&input[idx - 1]);
            let u_right = vld1q_f32(&input[idx + 1]);
            let u_bottom = vld1q_f32(&input[idx - nx]);
            let u_top = vld1q_f32(&input[idx + nx]);

            // Compute Laplacian
            let d2u_dx2 = vsubq_f32(u_left, vmulq_f32(two, u_center));
            let d2u_dx2 = vaddq_f32(d2u_dx2, u_right);
            let d2u_dx2 = vmulq_f32(d2u_dx2, dx2_inv_vec);

            let d2u_dy2 = vsubq_f32(u_bottom, vmulq_f32(two, u_center));
            let d2u_dy2 = vaddq_f32(d2u_dy2, u_top);
            let d2u_dy2 = vmulq_f32(d2u_dy2, dy2_inv_vec);

            let laplacian = vaddq_f32(d2u_dx2, d2u_dy2);

            // Update: u + dt * nu * laplacian
            let update = vmulq_f32(laplacian, dt_nu_vec);
            let result = vaddq_f32(u_center, update);

            // Store result
            vst1q_f32(&mut output[idx], result);

            i += 4;
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

/// NEON implementation for vector dot product
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub unsafe fn dot_product_neon(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    let n = a.len();

    let mut sum = vdupq_n_f32(0.0);
    let mut i = 0;

    // Process 4 elements at a time
    while i + 3 < n {
        let a_vec = vld1q_f32(&a[i]);
        let b_vec = vld1q_f32(&b[i]);
        let prod = vmulq_f32(a_vec, b_vec);
        sum = vaddq_f32(sum, prod);
        i += 4;
    }

    // Sum the vector elements
    let mut result = vaddvq_f32(sum);

    // Handle remaining elements
    while i < n {
        result += a[i] * b[i];
        i += 1;
    }

    result
}
