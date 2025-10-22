//! Comprehensive edge case tests for TVD limiters
//!
//! This module provides rigorous testing of TVD limiter functions with
//! property-based testing and edge case validation per CFD standards.
//!
//! References:
//! - Sweby, P.K. (1984). "High Resolution Schemes Using Flux Limiters"
//! - Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"
//! - Hirsch, C. (2007). "Numerical Computation of Internal and External Flows"

#[cfg(test)]
mod tvd_limiter_edge_cases {
    use super::super::*;
    use crate::physics::momentum::tvd_limiters::UpwindLimiter;
    use proptest::prelude::*;
    
    /// Edge case: Zero gradient ratio (extrema point)
    /// All limiters should return 0 at extrema (Sweby 1984)
    #[test]
    fn test_zero_gradient_ratio() {
        let limiters: Vec<(&str, Box<dyn TvdLimiter<f64>>)> = vec![
            ("Superbee", Box::new(Superbee)),
            ("VanLeer", Box::new(VanLeer)),
            ("Minmod", Box::new(Minmod)),
            ("MC", Box::new(MonotonizedCentral)),
            ("Upwind", Box::new(UpwindLimiter)),
        ];
        
        for (name, limiter) in limiters {
            let psi = limiter.limit(0.0);
            assert_eq!(psi, 0.0, "{}: ψ(0) must be 0 (upwind at extrema)", name);
        }
    }
    
    /// Edge case: Negative gradient ratios
    /// Negative r indicates extremum, should return 0
    #[test]
    fn test_negative_gradient_ratios() {
        let limiters: Vec<Box<dyn TvdLimiter<f64>>> = vec![
            Box::new(Superbee),
            Box::new(VanLeer),
            Box::new(Minmod),
            Box::new(MonotonizedCentral),
        ];
        
        for limiter in limiters {
            for r in [-10.0, -1.0, -0.5, -0.1, -1e-10] {
                let psi = limiter.limit(r);
                assert_eq!(
                    psi, 0.0,
                    "{}: ψ({}) should be 0 for negative r (downwind of extremum)",
                    limiter.name(), r
                );
            }
        }
    }
    
    /// Edge case: Very large gradient ratios
    /// For r >> 1, ψ should not exceed 2 (TVD constraint)
    #[test]
    fn test_large_gradient_ratios() {
        let limiters: Vec<Box<dyn TvdLimiter<f64>>> = vec![
            Box::new(Superbee),
            Box::new(VanLeer),
            Box::new(Minmod),
            Box::new(MonotonizedCentral),
        ];
        
        for limiter in limiters {
            for r in [10.0, 100.0, 1000.0, 1e6] {
                let psi = limiter.limit(r);
                assert!(
                    psi <= 2.0 + 1e-10,
                    "{}: ψ({}) = {} exceeds TVD limit of 2",
                    limiter.name(), r, psi
                );
            }
        }
    }
    
    /// Edge case: r = 1 (smooth gradient)
    /// At r=1, gradient is uniform across three points
    #[test]
    fn test_uniform_gradient() {
        let superbee = Superbee;
        let van_leer = VanLeer;
        let minmod = Minmod;
        let mc = MonotonizedCentral;
        
        // Superbee: ψ(1) = max(min(2*1, 1), min(1, 2*1)) = 1
        assert!((superbee.limit(1.0_f64) - 1.0).abs() < 1e-10);
        
        // Van Leer: ψ(1) = 2*1/(1+1) = 1
        assert!((van_leer.limit(1.0_f64) - 1.0).abs() < 1e-10);
        
        // Minmod: ψ(1) = min(1, 2*1) = 1
        assert!((minmod.limit(1.0_f64) - 1.0).abs() < 1e-10);
        
        // MC: ψ(1) = min(2*1, (1+1)/2, 2) = 1
        assert!((mc.limit(1.0_f64) - 1.0).abs() < 1e-10);
    }
    
    /// Edge case: Very small positive gradient ratios
    /// For r → 0+, ψ should approach 0
    #[test]
    fn test_small_positive_gradient_ratios() {
        let limiters: Vec<Box<dyn TvdLimiter<f64>>> = vec![
            Box::new(Superbee),
            Box::new(VanLeer),
            Box::new(Minmod),
            Box::new(MonotonizedCentral),
        ];
        
        for limiter in limiters {
            for r in [1e-10, 1e-8, 1e-6, 1e-4] {
                let psi = limiter.limit(r);
                assert!(
                    psi >= 0.0,
                    "{}: ψ({}) should be non-negative",
                    limiter.name(), r
                );
                // For small r, ψ should be small (but can be O(r) depending on limiter)
                assert!(
                    psi < 1.0,
                    "{}: ψ({}) = {} should be small for small r",
                    limiter.name(), r, psi
                );
            }
        }
    }
    
    /// Edge case: Uniform field (zero denominator in interpolation)
    /// Should return central value without error
    #[test]
    fn test_uniform_field_interpolation() {
        let limiters: Vec<Box<dyn TvdLimiter<f64>>> = vec![
            Box::new(Superbee),
            Box::new(VanLeer),
            Box::new(Minmod),
            Box::new(MonotonizedCentral),
        ];
        
        for limiter in limiters {
            let face = limiter.interpolate_face(5.0, 5.0, 5.0);
            assert_eq!(
                face, 5.0,
                "{}: Uniform field should return central value",
                limiter.name()
            );
        }
    }
    
    /// Edge case: Nearly uniform field (very small denominator)
    /// Should handle gracefully without division by near-zero
    #[test]
    fn test_nearly_uniform_field() {
        let limiters: Vec<Box<dyn TvdLimiter<f64>>> = vec![
            Box::new(Superbee),
            Box::new(VanLeer),
            Box::new(Minmod),
            Box::new(MonotonizedCentral),
        ];
        
        let phi_c = 1.0;
        let phi_d = phi_c + 1e-12; // Very small difference
        
        for limiter in limiters {
            let face = limiter.interpolate_face(phi_c - 1e-12, phi_c, phi_d);
            assert!(
                face.is_finite(),
                "{}: Should handle near-uniform field",
                limiter.name()
            );
            assert!(
                (face - phi_c).abs() < 1e-6,
                "{}: Face value should be close to central for near-uniform field",
                limiter.name()
            );
        }
    }
    
    /// Edge case: Shock-like discontinuity (r >> 1)
    /// Limiters should compress to prevent oscillations
    #[test]
    fn test_shock_discontinuity() {
        let superbee = Superbee;
        
        // Sharp gradient: r = (10-0)/(100-10) = 10/90 ≈ 0.11
        let face = superbee.interpolate_face(0.0, 10.0, 100.0);
        
        // Should be between upwind and downwind
        assert!(face >= 10.0, "Face should be ≥ upwind value");
        assert!(face <= 100.0, "Face should be ≤ downwind value");
    }
    
    /// Edge case: Reverse gradient (local maximum)
    /// r < 0 should result in first-order upwind
    #[test]
    fn test_local_maximum() {
        let limiters: Vec<Box<dyn TvdLimiter<f64>>> = vec![
            Box::new(Superbee),
            Box::new(VanLeer),
            Box::new(Minmod),
        ];
        
        // Local max: values 10, 20, 10 => r = (20-10)/(10-20) = -1
        for limiter in limiters {
            let face = limiter.interpolate_face(10.0, 20.0, 10.0);
            // With ψ(-1) = 0, face = 20 + 0.5*0*(10-20) = 20
            assert_eq!(
                face, 20.0,
                "{}: Local maximum should use upwind value",
                limiter.name()
            );
        }
    }
    
    /// Edge case: Monotonic increasing sequence
    /// Should apply high-order correction
    #[test]
    fn test_monotonic_increase() {
        let superbee = Superbee;
        
        // Monotonic: 1, 2, 3 => r = (2-1)/(3-2) = 1
        // ψ(1) = 1, face = 2 + 0.5*1*(3-2) = 2.5
        let face = superbee.interpolate_face(1.0_f64, 2.0, 3.0);
        assert!((face - 2.5).abs() < 1e-10, "Monotonic sequence should use high-order");
    }
    
    /// Edge case: Upwind limiter should always be first-order
    /// ψ(r) = 0 for all r
    #[test]
    fn test_upwind_always_first_order() {
        let upwind = UpwindLimiter;
        
        for r in [-10.0, -1.0, 0.0, 0.5, 1.0, 2.0, 10.0] {
            let psi = upwind.limit(r);
            assert_eq!(psi, 0.0, "Upwind limiter should always return 0");
        }
        
        // Interpolation should always return central value
        let face = upwind.interpolate_face(1.0_f64, 2.0, 3.0);
        assert_eq!(face, 2.0, "Upwind interpolation should return central value");
    }
    
    /// Property test: TVD region constraint for r ∈ [0, 1]
    /// ψ(r) ≤ 2r for 0 ≤ r ≤ 1
    proptest! {
        #[test]
        fn prop_tvd_region_0_to_1(r in 0.0..1.0f64) {
            let limiters: Vec<(&str, Box<dyn TvdLimiter<f64>>)> = vec![
                ("Superbee", Box::new(Superbee)),
                ("VanLeer", Box::new(VanLeer)),
                ("Minmod", Box::new(Minmod)),
                ("MC", Box::new(MonotonizedCentral)),
            ];
            
            for (name, limiter) in limiters {
                let psi = limiter.limit(r);
                
                // Property: 0 ≤ ψ(r) ≤ 2r for r ∈ [0, 1]
                assert!(psi >= 0.0, "{}: ψ({}) = {} < 0", name, r, psi);
                assert!(
                    psi <= 2.0 * r + 1e-10,
                    "{}: ψ({}) = {} > 2r = {}",
                    name, r, psi, 2.0 * r
                );
            }
        }
    }
    
    /// Property test: TVD region constraint for r > 1
    /// 0 ≤ ψ(r) ≤ 2 for r > 1
    proptest! {
        #[test]
        fn prop_tvd_region_above_1(r in 1.0..10.0f64) {
            let limiters: Vec<(&str, Box<dyn TvdLimiter<f64>>)> = vec![
                ("Superbee", Box::new(Superbee)),
                ("VanLeer", Box::new(VanLeer)),
                ("Minmod", Box::new(Minmod)),
                ("MC", Box::new(MonotonizedCentral)),
            ];
            
            for (name, limiter) in limiters {
                let psi = limiter.limit(r);
                
                // Property: 0 ≤ ψ(r) ≤ 2 for r > 1
                assert!(psi >= 0.0, "{}: ψ({}) = {} < 0", name, r, psi);
                assert!(
                    psi <= 2.0 + 1e-10,
                    "{}: ψ({}) = {} > 2",
                    name, r, psi
                );
            }
        }
    }
    
    /// Property test: Negative r always returns 0
    /// ψ(r) = 0 for r < 0 (downwind of extremum)
    proptest! {
        #[test]
        fn prop_negative_r_zero(r in -10.0..0.0f64) {
            let limiters: Vec<(&str, Box<dyn TvdLimiter<f64>>)> = vec![
                ("Superbee", Box::new(Superbee)),
                ("VanLeer", Box::new(VanLeer)),
                ("Minmod", Box::new(Minmod)),
                ("MC", Box::new(MonotonizedCentral)),
            ];
            
            for (name, limiter) in limiters {
                let psi = limiter.limit(r);
                assert_eq!(psi, 0.0, "{}: ψ({}) must be 0 for negative r", name, r);
            }
        }
    }
    
    /// Property test: Interpolated face value is bounded
    /// φ_min ≤ φ_face ≤ φ_max
    proptest! {
        #[test]
        fn prop_face_value_bounded(
            phi_u in -100.0..100.0f64,
            phi_c in -100.0..100.0f64,
            phi_d in -100.0..100.0f64
        ) {
            let limiters: Vec<Box<dyn TvdLimiter<f64>>> = vec![
                Box::new(Superbee),
                Box::new(VanLeer),
                Box::new(Minmod),
                Box::new(MonotonizedCentral),
            ];
            
            let phi_min = phi_u.min(phi_c).min(phi_d);
            let phi_max = phi_u.max(phi_c).max(phi_d);
            
            for limiter in limiters {
                let face = limiter.interpolate_face(phi_u, phi_c, phi_d);
                
                // Property: Face value should be within bounds (with tolerance)
                assert!(
                    face >= phi_min - 1e-6,
                    "{}: face = {} < min = {}",
                    limiter.name(), face, phi_min
                );
                assert!(
                    face <= phi_max + 1e-6,
                    "{}: face = {} > max = {}",
                    limiter.name(), face, phi_max
                );
            }
        }
    }
    
    /// Property test: Monotonicity preservation
    /// For monotonic data, face value should maintain monotonicity
    proptest! {
        #[test]
        fn prop_monotonicity_preservation(
            phi_u in 0.0..100.0f64,
            delta1 in 0.0..10.0f64,
            delta2 in 0.0..10.0f64
        ) {
            let phi_c = phi_u + delta1;
            let phi_d = phi_c + delta2;
            
            let limiters: Vec<Box<dyn TvdLimiter<f64>>> = vec![
                Box::new(Superbee),
                Box::new(VanLeer),
                Box::new(Minmod),
                Box::new(MonotonizedCentral),
            ];
            
            for limiter in limiters {
                let face = limiter.interpolate_face(phi_u, phi_c, phi_d);
                
                // Property: For monotonic increasing data, φ_u ≤ face ≤ φ_d
                assert!(
                    face >= phi_u - 1e-6,
                    "{}: face = {} < phi_u = {}",
                    limiter.name(), face, phi_u
                );
                assert!(
                    face <= phi_d + 1e-6,
                    "{}: face = {} > phi_d = {}",
                    limiter.name(), face, phi_d
                );
            }
        }
    }
    
    /// Edge case: Limiter compressiveness ordering
    /// Superbee > MC > Van Leer > Minmod (for most r values)
    #[test]
    fn test_compressiveness_ordering() {
        let superbee = Superbee;
        let mc = MonotonizedCentral;
        let van_leer = VanLeer;
        let minmod = Minmod;
        
        // Test at r = 0.5 (typical smooth gradient)
        let r = 0.5;
        let psi_superbee = superbee.limit(r);
        let psi_mc = mc.limit(r);
        let psi_van_leer = van_leer.limit(r);
        let psi_minmod = minmod.limit(r);
        
        // Superbee should be most compressive (largest ψ)
        assert!(psi_superbee >= psi_mc, "Superbee should be ≥ MC");
        assert!(psi_mc >= psi_van_leer, "MC should be ≥ Van Leer");
        assert!(psi_minmod <= psi_van_leer, "Minmod should be ≤ Van Leer");
    }
    
    /// Edge case: Limiter names are correct
    #[test]
    fn test_limiter_names() {
        assert_eq!(<Superbee as TvdLimiter<f64>>::name(&Superbee), "Superbee");
        assert_eq!(<VanLeer as TvdLimiter<f64>>::name(&VanLeer), "VanLeer");
        assert_eq!(<Minmod as TvdLimiter<f64>>::name(&Minmod), "Minmod");
        assert_eq!(<MonotonizedCentral as TvdLimiter<f64>>::name(&MonotonizedCentral), "MonotonizedCentral");
        assert_eq!(<UpwindLimiter as TvdLimiter<f64>>::name(&UpwindLimiter), "Upwind");
    }
}
