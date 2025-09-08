//! Tests for SIMD and SWAR implementations

#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::simd::arch_detect::ArchDetect;
    use approx::assert_relative_eq;

    #[test]
    fn test_arch_detection() {
        let arch = ArchDetect::new();
        let capability = arch.capability();

        // Should detect something
        assert!(matches!(
            capability,
            SimdCapability::Avx2
                | SimdCapability::Sse42
                | SimdCapability::Neon
                | SimdCapability::Swar
        ));

        // Vector widths should be sensible
        assert!(arch.vector_width_f32() >= 1);
        assert!(arch.vector_width_f64() >= 1);
    }

    #[test]
    fn test_simd_add_f32() {
        let simd = SimdOps::new();
        let a = vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let b = vec![8.0f32, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        let mut result = vec![0.0f32; 8];

        simd.add(&a, &b, &mut result).unwrap();

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

        simd.mul(&a, &b, &mut result).unwrap();

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

        simd.scale(&input, scalar, &mut result).unwrap();

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

        let dot = simd.dot(&a, &b).unwrap();

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

        simd.add(&a, &b, &mut result).unwrap();

        for i in 0..13 {
            assert_relative_eq!(result[i], 3.0, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_swar_add_f32() {
        let swar = SwarOps::new();
        let a = vec![1.0f32, 2.0, 3.0, 4.0, 5.0];
        let b = vec![5.0f32, 4.0, 3.0, 2.0, 1.0];
        let mut result = vec![0.0f32; 5];

        swar.process_binary_f32(&a, &b, &mut result, |a, b, r| {
            for i in 0..a.len() {
                r[i] = a[i] + b[i];
            }
        })
        .unwrap();

        for i in 0..5 {
            assert_relative_eq!(result[i], 6.0, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_swar_mul_f32() {
        let swar = SwarOps::new();
        let a = vec![1.0f32, 2.0, 3.0, 4.0];
        let b = vec![2.0f32, 2.0, 2.0, 2.0];
        let mut result = vec![0.0f32; 4];

        swar.mul_f32(&a, &b, &mut result).unwrap();

        assert_relative_eq!(result[0], 2.0, epsilon = 1e-6);
        assert_relative_eq!(result[1], 4.0, epsilon = 1e-6);
        assert_relative_eq!(result[2], 6.0, epsilon = 1e-6);
        assert_relative_eq!(result[3], 8.0, epsilon = 1e-6);
    }

    #[test]
    fn test_swar_dot_f32() {
        let swar = SwarOps::new();
        let a = vec![1.0f32, 2.0, 3.0, 4.0];
        let b = vec![4.0f32, 3.0, 2.0, 1.0];

        let dot = swar.dot_f32(&a, &b).unwrap();

        assert_relative_eq!(dot, 20.0, epsilon = 1e-6);
    }

    #[test]
    fn test_swar_sum_f32() {
        let swar = SwarOps::new();
        let input = vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];

        let sum = swar.sum_f32(&input).unwrap();

        assert_relative_eq!(sum, 45.0, epsilon = 1e-6);
    }

    #[test]
    fn test_swar_max_f32() {
        let swar = SwarOps::new();
        let input = vec![3.0f32, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];

        let max = swar.max_f32(&input).unwrap();

        assert_relative_eq!(max, 9.0, epsilon = 1e-6);
    }

    #[test]
    fn test_swar_add_u32() {
        let swar = SwarOps::new();
        let a = vec![1u32, 2, 3, 4, 5];
        let b = vec![5u32, 4, 3, 2, 1];
        let mut result = vec![0u32; 5];

        swar.add_u32(&a, &b, &mut result).unwrap();

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

    // Benchmark-style test for performance validation
    #[test]
    fn test_simd_performance_characteristics() {
        let simd = SimdOps::new();
        let arch = ArchDetect::new();

        // Create large vectors
        let size = 1024;
        let a = vec![1.0f32; size];
        let b = vec![2.0f32; size];
        let mut result = vec![0.0f32; size];

        // Should complete without error
        simd.add(&a, &b, &mut result).unwrap();

        // Verify correctness
        for i in 0..size {
            assert_relative_eq!(result[i], 3.0, epsilon = 1e-6);
        }

        println!("SIMD capability detected: {:?}", arch.capability());
        println!("Vector width (f32): {}", arch.vector_width_f32());
        println!("Vector width (f64): {}", arch.vector_width_f64());
    }
}
