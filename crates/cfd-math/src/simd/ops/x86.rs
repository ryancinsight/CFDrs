//! x86/x86_64 SIMD implementations (SSE, AVX2)


// AVX2 implementations for f32
#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn add_avx2_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !7; // Process 8 elements at a time

    for i in (0..simd_len).step_by(8) {
        let va = _mm256_loadu_ps(a.as_ptr().add(i));
        let vb = _mm256_loadu_ps(b.as_ptr().add(i));
        let vr = _mm256_add_ps(va, vb);
        _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
    }

    // Handle remainder
    for i in simd_len..len {
        result[i] = a[i] + b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn sub_avx2_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !7;

    for i in (0..simd_len).step_by(8) {
        let va = _mm256_loadu_ps(a.as_ptr().add(i));
        let vb = _mm256_loadu_ps(b.as_ptr().add(i));
        let vr = _mm256_sub_ps(va, vb);
        _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] - b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn mul_avx2_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !7;

    for i in (0..simd_len).step_by(8) {
        let va = _mm256_loadu_ps(a.as_ptr().add(i));
        let vb = _mm256_loadu_ps(b.as_ptr().add(i));
        let vr = _mm256_mul_ps(va, vb);
        _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] * b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "fma"))]
#[inline]
pub unsafe fn fma_avx2_f32(a: &[f32], b: &[f32], c: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !7;

    for i in (0..simd_len).step_by(8) {
        let va = _mm256_loadu_ps(a.as_ptr().add(i));
        let vb = _mm256_loadu_ps(b.as_ptr().add(i));
        let vc = _mm256_loadu_ps(c.as_ptr().add(i));
        let vr = _mm256_fmadd_ps(va, vb, vc);
        _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i].mul_add(b[i], c[i]);
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn scale_avx2_f32(input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = input.len();
    let simd_len = len & !7;
    let vs = _mm256_set1_ps(scalar);

    for i in (0..simd_len).step_by(8) {
        let vi = _mm256_loadu_ps(input.as_ptr().add(i));
        let vr = _mm256_mul_ps(vi, vs);
        _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = input[i] * scalar;
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn dot_avx2_f32(a: &[f32], b: &[f32]) -> f32 {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !7;
    let mut sum = _mm256_setzero_ps();

    for i in (0..simd_len).step_by(8) {
        let va = _mm256_loadu_ps(a.as_ptr().add(i));
        let vb = _mm256_loadu_ps(b.as_ptr().add(i));
        sum = _mm256_fmadd_ps(va, vb, sum);
    }

    // Horizontal sum
    let sum_array = std::mem::transmute::<__m256, [f32; 8]>(sum);
    let mut result = sum_array.iter().sum::<f32>();

    // Handle remainder
    for i in simd_len..len {
        result += a[i] * b[i];
    }

    result
}

// SSE4.2 implementations for f32
#[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
#[inline]
pub unsafe fn add_sse42_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !3; // Process 4 elements at a time

    for i in (0..simd_len).step_by(4) {
        let va = _mm_loadu_ps(a.as_ptr().add(i));
        let vb = _mm_loadu_ps(b.as_ptr().add(i));
        let vr = _mm_add_ps(va, vb);
        _mm_storeu_ps(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] + b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
#[inline]
pub unsafe fn sub_sse42_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !3;

    for i in (0..simd_len).step_by(4) {
        let va = _mm_loadu_ps(a.as_ptr().add(i));
        let vb = _mm_loadu_ps(b.as_ptr().add(i));
        let vr = _mm_sub_ps(va, vb);
        _mm_storeu_ps(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] - b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
#[inline]
pub unsafe fn mul_sse42_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !3;

    for i in (0..simd_len).step_by(4) {
        let va = _mm_loadu_ps(a.as_ptr().add(i));
        let vb = _mm_loadu_ps(b.as_ptr().add(i));
        let vr = _mm_mul_ps(va, vb);
        _mm_storeu_ps(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] * b[i];
    }

    Ok(())
}

// AVX2 implementations for f64
#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn add_avx2_f64(a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !3; // Process 4 elements at a time

    for i in (0..simd_len).step_by(4) {
        let va = _mm256_loadu_pd(a.as_ptr().add(i));
        let vb = _mm256_loadu_pd(b.as_ptr().add(i));
        let vr = _mm256_add_pd(va, vb);
        _mm256_storeu_pd(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] + b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn sub_avx2_f64(a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !3;

    for i in (0..simd_len).step_by(4) {
        let va = _mm256_loadu_pd(a.as_ptr().add(i));
        let vb = _mm256_loadu_pd(b.as_ptr().add(i));
        let vr = _mm256_sub_pd(va, vb);
        _mm256_storeu_pd(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] - b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn mul_avx2_f64(a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !3;

    for i in (0..simd_len).step_by(4) {
        let va = _mm256_loadu_pd(a.as_ptr().add(i));
        let vb = _mm256_loadu_pd(b.as_ptr().add(i));
        let vr = _mm256_mul_pd(va, vb);
        _mm256_storeu_pd(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = a[i] * b[i];
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn scale_avx2_f64(input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
    use std::arch::x86_64::*;

    let len = input.len();
    let simd_len = len & !3;
    let vs = _mm256_set1_pd(scalar);

    for i in (0..simd_len).step_by(4) {
        let vi = _mm256_loadu_pd(input.as_ptr().add(i));
        let vr = _mm256_mul_pd(vi, vs);
        _mm256_storeu_pd(result.as_mut_ptr().add(i), vr);
    }

    for i in simd_len..len {
        result[i] = input[i] * scalar;
    }

    Ok(())
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline]
pub unsafe fn dot_avx2_f64(a: &[f64], b: &[f64]) -> f64 {
    use std::arch::x86_64::*;

    let len = a.len();
    let simd_len = len & !3;
    let mut sum = _mm256_setzero_pd();

    for i in (0..simd_len).step_by(4) {
        let va = _mm256_loadu_pd(a.as_ptr().add(i));
        let vb = _mm256_loadu_pd(b.as_ptr().add(i));
        sum = _mm256_fmadd_pd(va, vb, sum);
    }

    // Horizontal sum
    let sum_array = std::mem::transmute::<__m256d, [f64; 4]>(sum);
    let mut result = sum_array.iter().sum::<f64>();

    // Handle remainder
    for i in simd_len..len {
        result += a[i] * b[i];
    }

    result
}
