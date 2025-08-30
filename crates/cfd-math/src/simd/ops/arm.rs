//! ARM NEON SIMD implementations

use crate::error::Result;

// NEON implementations for f32
#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
#[inline]
pub unsafe fn add_neon_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::aarch64::*;

    let len = a.len();
    let simd_len = len & !3; // Process 4 elements at a time

    for i in (0..simd_len).step_by(4) {
        let va = vld1q_f32(a.as_ptr().add(i));
        let vb = vld1q_f32(b.as_ptr().add(i));
        let vr = vaddq_f32(va, vb);
        vst1q_f32(result.as_mut_ptr().add(i), vr);
    }

    // Handle remainder
    for i in simd_len..len {
        result[i] = a[i] + b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
#[inline]
pub unsafe fn sub_neon_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::aarch64::*;

    let len = a.len();
    let simd_len = len & !3;

    for i in (0..simd_len).step_by(4) {
        let va = vld1q_f32(a.as_ptr().add(i));
        let vb = vld1q_f32(b.as_ptr().add(i));
        let vr = vsubq_f32(va, vb);
        vst1q_f32(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] - b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
#[inline]
pub unsafe fn mul_neon_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::aarch64::*;

    let len = a.len();
    let simd_len = len & !3;

    for i in (0..simd_len).step_by(4) {
        let va = vld1q_f32(a.as_ptr().add(i));
        let vb = vld1q_f32(b.as_ptr().add(i));
        let vr = vmulq_f32(va, vb);
        vst1q_f32(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] * b[i];
    }

    Ok(())
}
