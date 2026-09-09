//! Tests for SIMD operations backed by hermes-simd

#[cfg(test)]
use super::*;
use eunomia::assert_relative_eq;
use hermes_simd::dispatch::SimdOps as HermesOps;

#[test]
fn test_simd_add_f32() {
    let simd = SimdOps::new();
    let a = vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let b = vec![8.0f32, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
    let mut result = vec![0.0f32; 8];

    simd.add(&a, &b, &mut result).expect("expected value");

    for i in 0..8 {
        assert_relative_eq!(result[i], 9.0, epsilon = 1e-6);
    }
}

#[test]
fn test_simd_mul_f32() {
    let simd = SimdOps::new();
    let a = vec![1.0f32, 2.0, 3.0, 4.0];
    let b = vec![2.0f32, 3.0, 4.0, 5.0];
    let mut result = vec![0.0f32; 4];

    simd.mul(&a, &b, &mut result).expect("expected value");

    assert_relative_eq!(result[0], 2.0, epsilon = 1e-6);
    assert_relative_eq!(result[1], 6.0, epsilon = 1e-6);
    assert_relative_eq!(result[2], 12.0, epsilon = 1e-6);
    assert_relative_eq!(result[3], 20.0, epsilon = 1e-6);
}

#[test]
fn test_simd_scale_f32() {
    let simd = SimdOps::new();
    let input = vec![1.0f32, 2.0, 3.0, 4.0];
    let scalar = 2.5;
    let mut result = vec![0.0f32; 4];

    simd.scale(&input, scalar, &mut result)
        .expect("expected value");

    assert_relative_eq!(result[0], 2.5, epsilon = 1e-6);
    assert_relative_eq!(result[1], 5.0, epsilon = 1e-6);
    assert_relative_eq!(result[2], 7.5, epsilon = 1e-6);
    assert_relative_eq!(result[3], 10.0, epsilon = 1e-6);
}

#[test]
fn test_simd_dot_f32() {
    let simd = SimdOps::new();
    let a = vec![1.0f32, 2.0, 3.0, 4.0];
    let b = vec![4.0f32, 3.0, 2.0, 1.0];

    let dot = simd.dot(&a, &b).expect("expected value");

    // 1*4 + 2*3 + 3*2 + 4*1 = 4 + 6 + 6 + 4 = 20
    assert_relative_eq!(dot, 20.0, epsilon = 1e-6);
}

#[test]
fn test_simd_unaligned_lengths() {
    let simd = SimdOps::new();

    // Test with length not divisible by vector width
    let a = vec![1.0f32; 13];
    let b = vec![2.0f32; 13];
    let mut result = vec![0.0f32; 13];

    simd.add(&a, &b, &mut result).expect("expected value");

    for i in 0..13 {
        assert_relative_eq!(result[i], 3.0, epsilon = 1e-6);
    }
}

#[test]
fn test_simd_dot_f64() {
    let simd = SimdOps::new();
    let a = vec![1.0f64, 2.0, 3.0, 4.0];
    let b = vec![4.0f64, 3.0, 2.0, 1.0];

    let dot = simd.dot_f64(&a, &b).expect("expected value");

    assert_relative_eq!(dot, 20.0, epsilon = 1e-12);
}

#[test]
fn test_simd_sum_f32() {
    let simd = SimdOps::new();
    let input = vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];

    let sum = simd.sum_f32(&input).expect("expected value");

    assert_relative_eq!(sum, 45.0, epsilon = 1e-6);
}

#[test]
fn test_simd_max_f32() {
    let simd = SimdOps::new();
    let input = vec![3.0f32, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];

    let max = simd.max_f32(&input).expect("expected value");

    assert_relative_eq!(max, 9.0, epsilon = 1e-6);
}

#[test]
fn test_simd_add_u32() {
    let simd = SimdOps::new();
    let a = vec![1u32, 2, 3, 4, 5];
    let b = vec![5u32, 4, 3, 2, 1];
    let mut result = vec![0u32; 5];

    simd.add_u32(&a, &b, &mut result).expect("expected value");

    for i in 0..5 {
        assert_eq!(result[i], 6);
    }
}

#[test]
fn test_dimension_mismatch() {
    let simd = SimdOps::new();
    let a = vec![1.0f32; 5];
    let b = vec![1.0f32; 4];
    let mut result = vec![0.0f32; 5];

    assert!(simd.add(&a, &b, &mut result).is_err());
}

#[test]
fn test_empty_vectors() {
    let simd = SimdOps::new();
    let a: Vec<f32> = vec![];
    let b: Vec<f32> = vec![];
    let mut result: Vec<f32> = vec![];

    assert!(simd.add(&a, &b, &mut result).is_ok());
}

#[test]
fn test_simd_performance_characteristics() {
    let simd = SimdOps::new();

    // Create large vectors
    let size = 1024;
    let a = vec![1.0f32; size];
    let b = vec![2.0f32; size];
    let mut result = vec![0.0f32; size];

    // Should complete without error
    simd.add(&a, &b, &mut result).expect("expected value");

    // Verify correctness
    for i in 0..size {
        assert_relative_eq!(result[i], 3.0, epsilon = 1e-6);
    }
}

#[test]
fn test_convolution() {
    use crate::simd::vectorization::VectorizedOps;
    let ops = VectorizedOps::new();
    let signal = vec![1.0, 2.0, 3.0];
    let kernel = vec![0.5, 1.0];
    // Expected result len: 3 + 2 - 1 = 4
    let mut result = vec![0.0; 4];

    ops.convolution(&signal, &kernel, &mut result)
        .expect("expected value");

    // n=0: signal[0]*kernel[0] = 1*0.5 = 0.5
    // n=1: signal[1]*kernel[0] + signal[0]*kernel[1] = 2*0.5 + 1*1.0 = 2.0
    // n=2: signal[2]*kernel[0] + signal[1]*kernel[1] = 3*0.5 + 2*1.0 = 3.5
    // n=3: signal[2]*kernel[1] = 3*1.0 = 3.0

    assert_relative_eq!(result[0], 0.5, epsilon = 1e-6);
    assert_relative_eq!(result[1], 2.0, epsilon = 1e-6);
    assert_relative_eq!(result[2], 3.5, epsilon = 1e-6);
    assert_relative_eq!(result[3], 3.0, epsilon = 1e-6);
}

// ── Differential suite: hermes-dispatched kernels vs scalar references ───
//
// Every test is deterministic (fixed LCG seed) and sweeps lengths that hit
// each dispatch regime — empty, sub-lane tails, lane boundaries, unrolled
// chunk boundaries, and large arrays — so masked-tail handling is exercised
// on every path regardless of the host's widest architecture.

/// Deterministic LCG producing values in [-1, 1).
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let v = ((self.0 >> 40) as f64) / ((1u64 << 24) as f64);
        2.0 * v - 1.0
    }

    fn next_f32(&mut self) -> f32 {
        self.next_f64() as f32
    }

    fn fill_f32(&mut self, len: usize) -> Vec<f32> {
        (0..len).map(|_| self.next_f32()).collect()
    }

    fn fill_f64(&mut self, len: usize) -> Vec<f64> {
        (0..len).map(|_| self.next_f64()).collect()
    }
}

/// Lengths covering every dispatch regime: empty, 1..=3 (sub-lane tails),
/// lane widths (4/8/16/32), unrolled chunk boundaries, and large arrays.
const SWEEP: &[usize] = &[
    0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 255, 256, 257, 511,
    512, 1000, 1024,
];

fn ref_add<T: Copy + std::ops::Add<Output = T>>(a: &[T], b: &[T]) -> Vec<T> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect()
}

fn ref_sub<T: Copy + std::ops::Sub<Output = T>>(a: &[T], b: &[T]) -> Vec<T> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x - y).collect()
}

fn ref_mul<T: Copy + std::ops::Mul<Output = T>>(a: &[T], b: &[T]) -> Vec<T> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).collect()
}

fn ref_scale<T: Copy + std::ops::Mul<Output = T>>(a: &[T], s: T) -> Vec<T> {
    a.iter().map(|&x| x * s).collect()
}

#[test]
fn differential_elementwise_f32_bit_exact() {
    let simd = SimdOps::new();
    let mut rng = Lcg::new(0xC0FFEE);

    for &len in SWEEP {
        let a = rng.fill_f32(len);
        let b = rng.fill_f32(len);
        let scalar = rng.next_f32();

        let mut out = vec![0.0f32; len];
        simd.add(&a, &b, &mut out).expect("add ok");
        assert_eq!(out, ref_add(&a, &b), "add mismatch at len {len}");

        simd.sub(&a, &b, &mut out).expect("sub ok");
        assert_eq!(out, ref_sub(&a, &b), "sub mismatch at len {len}");

        simd.mul(&a, &b, &mut out).expect("mul ok");
        assert_eq!(out, ref_mul(&a, &b), "mul mismatch at len {len}");

        // Divisor shifted into [1, 2) so division is total (no NaN/inf).
        let b_safe: Vec<f32> = b.iter().map(|&v| v + 2.0).collect();
        simd.div(&a, &b_safe, &mut out).expect("div ok");
        let expected: Vec<f32> = a.iter().zip(b_safe.iter()).map(|(&x, &y)| x / y).collect();
        assert_eq!(out, expected, "div mismatch at len {len}");

        simd.scale(&a, scalar, &mut out).expect("scale ok");
        assert_eq!(out, ref_scale(&a, scalar), "scale mismatch at len {len}");
    }
}

#[test]
fn differential_elementwise_f64_bit_exact() {
    let simd = SimdOps::new();
    let mut rng = Lcg::new(0xBADF00D);

    for &len in SWEEP {
        let a = rng.fill_f64(len);
        let b = rng.fill_f64(len);
        let scalar = rng.next_f64();

        let mut out = vec![0.0f64; len];
        simd.add_f64(&a, &b, &mut out).expect("add ok");
        assert_eq!(out, ref_add(&a, &b), "add mismatch at len {len}");

        simd.sub_f64(&a, &b, &mut out).expect("sub ok");
        assert_eq!(out, ref_sub(&a, &b), "sub mismatch at len {len}");

        simd.mul_f64(&a, &b, &mut out).expect("mul ok");
        assert_eq!(out, ref_mul(&a, &b), "mul mismatch at len {len}");

        let b_safe: Vec<f64> = b.iter().map(|&v| v + 2.0).collect();
        simd.div_f64(&a, &b_safe, &mut out).expect("div ok");
        let expected: Vec<f64> = a.iter().zip(b_safe.iter()).map(|(&x, &y)| x / y).collect();
        assert_eq!(out, expected, "div mismatch at len {len}");

        simd.scale_f64(&a, scalar, &mut out).expect("scale ok");
        assert_eq!(out, ref_scale(&a, scalar), "scale mismatch at len {len}");
    }
}

#[test]
fn differential_fma_f32_matches_scalar_reference() {
    let simd = SimdOps::new();
    let mut rng = Lcg::new(0xFACEB00C);

    for &len in SWEEP {
        let a = rng.fill_f32(len);
        let b = rng.fill_f32(len);
        let c = rng.fill_f32(len);
        let mut out = vec![0.0f32; len];

        simd.fma(&a, &b, &c, &mut out).expect("fma ok");

        for i in 0..len {
            let expected = a[i].mul_add(b[i], c[i]);
            // Unit-alpha FMADD is single-rounding on FMA hardware; a
            // scalar fallback may double-round (separate mul then add),
            // so allow one ulp of disagreement.
            let ulp = f32::EPSILON * expected.abs().max(1.0);
            assert!(
                (out[i] - expected).abs() <= ulp,
                "fma mismatch at len {len}, idx {i}: {} vs {expected}",
                out[i]
            );
        }
    }
}

#[test]
fn differential_fma_f64_matches_scalar_reference() {
    let simd = SimdOps::new();
    let mut rng = Lcg::new(0xDEADBEEF);

    for &len in SWEEP {
        let a = rng.fill_f64(len);
        let b = rng.fill_f64(len);
        let c = rng.fill_f64(len);
        let mut out = vec![0.0f64; len];

        simd.fma_f64(&a, &b, &c, &mut out).expect("fma ok");

        for i in 0..len {
            let expected = a[i].mul_add(b[i], c[i]);
            let ulp = f64::EPSILON * expected.abs().max(1.0);
            assert!(
                (out[i] - expected).abs() <= ulp,
                "fma mismatch at len {len}, idx {i}: {} vs {expected}",
                out[i]
            );
        }
    }
}

#[test]
fn differential_dot_f32_against_f64_oracle() {
    let simd = SimdOps::new();
    let mut rng = Lcg::new(0x5EED5EED);

    for &len in SWEEP {
        // Strictly positive inputs keep the sum bounded away from zero, so
        // the relative tolerance is meaningful (no cancellation blowup).
        let a: Vec<f32> = (0..len).map(|_| rng.next_f32() * 0.5 + 1.0).collect();
        let b: Vec<f32> = (0..len).map(|_| rng.next_f32() * 0.5 + 1.0).collect();

        let got = simd.dot(&a, &b).expect("dot ok");

        // f64 accumulation of f32 products is essentially exact here:
        // f32*f32 products are exactly representable in f64.
        let oracle: f64 = a
            .iter()
            .zip(b.iter())
            .map(|(&x, &y)| f64::from(x) * f64::from(y))
            .sum();

        // Multi-accumulator partial sums reorder additions; with n <= 1024
        // and magnitudes <= 2 the accumulated relative drift stays orders
        // of magnitude below 1e-3.
        let tol = 1e-3 * oracle.abs().max(1.0);
        assert!(
            (f64::from(got) - oracle).abs() <= f64::from(tol),
            "dot drift at len {len}: {got} vs oracle {oracle}"
        );
    }
}

#[test]
fn differential_dot_f64_accumulation_order() {
    let simd = SimdOps::new();
    let mut rng = Lcg::new(0x1234ABCD);

    for &len in SWEEP {
        let a: Vec<f64> = (0..len).map(|_| rng.next_f64() * 0.5 + 1.0).collect();
        let b: Vec<f64> = (0..len).map(|_| rng.next_f64() * 0.5 + 1.0).collect();

        let got = simd.dot_f64(&a, &b).expect("dot ok");
        let sequential: f64 = a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum();

        let tol = 1e-12 * sequential.abs().max(1.0);
        assert!(
            (got - sequential).abs() <= tol,
            "dot drift at len {len}: {got} vs sequential {sequential}"
        );
    }
}

#[test]
fn differential_sum_and_max() {
    let simd = SimdOps::new();
    let mut rng = Lcg::new(0x0DDBA11);

    for &len in SWEEP {
        // f32 sum against an exact f64 oracle (positive inputs).
        let input: Vec<f32> = (0..len).map(|_| rng.next_f32() * 0.5 + 1.0).collect();
        let sum = simd.sum_f32(&input).expect("sum ok");
        let oracle: f64 = input.iter().map(|&v| f64::from(v)).sum();
        assert!(
            (f64::from(sum) - oracle).abs() <= 1e-3 * oracle.abs().max(1.0),
            "sum drift at len {len}: {sum} vs oracle {oracle}"
        );

        if len > 0 {
            let max = simd.max_f32(&input).expect("max ok");
            assert_eq!(max, input.iter().copied().fold(f32::MIN, f32::max));
        }

        // f64: reductions are order-sensitive only for sum; max is exact.
        let input64 = rng.fill_f64(len);
        let sum64 = <f64 as HermesOps>::sum(&input64);
        let sequential64: f64 = input64.iter().sum();
        let tol64 = 1e-12 * sequential64.abs().max(1.0);
        assert!(
            (sum64 - sequential64).abs() <= tol64,
            "f64 sum drift at len {len}: {sum64} vs {sequential64}"
        );
        if len > 0 {
            let max64 = <f64 as HermesOps>::max(&input64);
            assert_eq!(max64, input64.iter().copied().fold(f64::MIN, f64::max));
        }
    }
}

#[test]
fn fma_and_scale_empty_are_noops() {
    let simd = SimdOps::new();
    let mut out: Vec<f32> = vec![];
    simd.fma(&[], &[], &[], &mut out).expect("empty fma ok");
    simd.scale(&[], 2.0, &mut out).expect("empty scale ok");
}

#[test]
fn fma_length_mismatch_is_error() {
    let simd = SimdOps::new();
    let a = vec![1.0f32; 4];
    let b = vec![1.0f32; 4];
    let c = vec![1.0f32; 3];
    let mut out = vec![0.0f32; 4];
    // Input-side mismatch is caught by cfd-math validation...
    assert!(simd.fma(&a, &b, &c, &mut out).is_err());
    // ...and the output side is caught by the same validation before
    // any kernel runs.
    let mut short_out = vec![0.0f32; 3];
    assert!(simd.fma(&a, &b, &a, &mut short_out).is_err());
}

#[test]
fn scale_length_mismatch_is_error() {
    let simd = SimdOps::new();
    let input = vec![1.0f32; 4];
    let mut short_out = vec![0.0f32; 3];
    assert!(simd.scale(&input, 2.0, &mut short_out).is_err());
}
