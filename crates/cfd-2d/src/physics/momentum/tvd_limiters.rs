//! Total Variation Diminishing (TVD) flux limiters for high-resolution schemes
//!
//! TVD limiters prevent spurious oscillations in convection-dominated flows
//! while maintaining high-order accuracy away from discontinuities.
//!
//! # Theory
//!
//! TVD limiters control the face value interpolation based on the local gradient ratio:
//! ```text
//! r = (φ_C - φ_U) / (φ_D - φ_C)
//! ```
//! where:
//! * φ_U = upwind value
//! * φ_C = central value  
//! * φ_D = downwind value
//!
//! The limiter function ψ(r) determines the face value:
//! ```text
//! φ_face = φ_C + 0.5 * ψ(r) * (φ_D - φ_C)
//! ```
//!
//! # TVD Regions
//!
//! For TVD property (Sweby 1984), ψ(r) must satisfy:
//! * 0 ≤ ψ(r) ≤ 2r for 0 ≤ r ≤ 1
//! * 0 ≤ ψ(r) ≤ 2 for r > 1
//! * ψ(0) = 0 (upwind at extrema)
//!
//! # Implemented Limiters
//!
//! * **Superbee** (Roe 1986): Most compressive, excellent shock capturing
//! * **Van Leer** (Van Leer 1974): Smooth, less compressive
//! * **Minmod** (Roe 1986): Most diffusive, very stable
//! * **MC** (Monotonized Central): Balance between Superbee and Van Leer
//!
//! # References
//!
//! * Sweby, P.K. (1984). "High Resolution Schemes Using Flux Limiters"
//! * Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"
//! * Van Leer, B. (1974). "Towards the Ultimate Conservative Difference Scheme"
//! * Hirsch, C. (2007). "Numerical Computation of Internal and External Flows"
//!
//! # Example
//!
//! ```ignore
//! use cfd_2d::physics::momentum::tvd_limiters::{TvdLimiter, Superbee};
//!
//! let limiter = Superbee;
//! let r = 2.0; // Gradient ratio
//! let psi = limiter.limit(r); // Returns limited value
//! ```

use nalgebra::RealField;

/// Trait for TVD flux limiters
///
/// Implements the limiter function ψ(r) that controls high-order
/// interpolation to maintain TVD property.
pub trait TvdLimiter<T: RealField + Copy> {
    /// Compute limiter function ψ(r)
    ///
    /// # Arguments
    /// * `r` - Gradient ratio: (φ_C - φ_U) / (φ_D - φ_C)
    ///
    /// # Returns
    /// Limiter value in [0, 2] satisfying TVD constraints
    fn limit(&self, r: T) -> T;

    /// Get limiter name for diagnostics
    fn name(&self) -> &'static str;

    /// Compute face value using limited interpolation
    ///
    /// # Arguments
    /// * `phi_u` - Upwind value
    /// * `phi_c` - Central value
    /// * `phi_d` - Downwind value
    ///
    /// # Returns
    /// Limited face value: φ_face = φ_C + 0.5 * ψ(r) * (φ_D - φ_C)
    fn interpolate_face(&self, phi_u: T, phi_c: T, phi_d: T) -> T {
        let half = T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()));
        
        // Compute gradient ratio r = (φ_C - φ_U) / (φ_D - φ_C)
        let denominator = phi_d - phi_c;
        
        // Handle zero or near-zero denominator (uniform region)
        if denominator.abs() < T::default_epsilon() * T::from_i32(10).unwrap_or(T::one()) {
            return phi_c;
        }
        
        let r = (phi_c - phi_u) / denominator;
        let psi = self.limit(r);
        
        // φ_face = φ_C + 0.5 * ψ(r) * (φ_D - φ_C)
        phi_c + half * psi * denominator
    }
}

/// Superbee limiter (Roe 1986)
///
/// Most compressive limiter, excellent for shock capturing.
/// ψ(r) = max(0, min(2r, 1), min(r, 2))
///
/// Properties:
/// * Sharpest discontinuity resolution
/// * May introduce slight overshoots
/// * Recommended for shock-like flows
///
/// # Reference
/// Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"
#[derive(Debug, Clone, Copy)]
pub struct Superbee;

impl<T: RealField + Copy> TvdLimiter<T> for Superbee {
    fn limit(&self, r: T) -> T {
        let zero = T::zero();
        let one = T::one();
        let two = one + one;
        
        if r <= zero {
            zero
        } else {
            // max(0, min(2r, 1), min(r, 2))
            let term1 = (two * r).min(one);
            let term2 = r.min(two);
            term1.max(term2)
        }
    }

    fn name(&self) -> &'static str {
        "Superbee"
    }
}

/// Van Leer limiter (Van Leer 1974)
///
/// Smooth limiter with good balance of accuracy and stability.
/// ψ(r) = (r + |r|) / (1 + |r|)
///
/// Properties:
/// * Smooth transition between schemes
/// * No overshoots or undershoots
/// * Good for general flows
///
/// # Reference
/// Van Leer, B. (1974). "Towards the Ultimate Conservative Difference Scheme"
#[derive(Debug, Clone, Copy)]
pub struct VanLeer;

impl<T: RealField + Copy> TvdLimiter<T> for VanLeer {
    fn limit(&self, r: T) -> T {
        let zero = T::zero();
        let one = T::one();
        
        if r <= zero {
            zero
        } else {
            // (r + |r|) / (1 + |r|) = 2r / (1 + r) for r > 0
            let two = one + one;
            two * r / (one + r)
        }
    }

    fn name(&self) -> &'static str {
        "VanLeer"
    }
}

/// Minmod limiter (Roe 1986)
///
/// Most diffusive limiter, extremely stable.
/// ψ(r) = max(0, min(r, 1))
///
/// Properties:
/// * Most dissipative TVD limiter
/// * No oscillations guaranteed
/// * Recommended for very difficult flows
///
/// # Reference
/// Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"
#[derive(Debug, Clone, Copy)]
pub struct Minmod;

impl<T: RealField + Copy> TvdLimiter<T> for Minmod {
    fn limit(&self, r: T) -> T {
        let zero = T::zero();
        let one = T::one();
        
        if r <= zero {
            zero
        } else {
            r.min(one)
        }
    }

    fn name(&self) -> &'static str {
        "Minmod"
    }
}

/// Monotonized Central (MC) limiter
///
/// Balance between Superbee and Van Leer.
/// ψ(r) = max(0, min(2r, (1+r)/2, 2))
///
/// Properties:
/// * More compressive than Van Leer
/// * More stable than Superbee
/// * Good general-purpose limiter
///
/// # Reference
/// Van Leer, B. (1977). "Towards the Ultimate Conservative Difference Scheme III"
#[derive(Debug, Clone, Copy)]
pub struct MonotonizedCentral;

impl<T: RealField + Copy> TvdLimiter<T> for MonotonizedCentral {
    fn limit(&self, r: T) -> T {
        let zero = T::zero();
        let one = T::one();
        let two = one + one;
        let half = T::from_f64(0.5).unwrap_or(one / two);
        
        if r <= zero {
            zero
        } else {
            // max(0, min(2r, (1+r)/2, 2))
            let term1 = two * r;
            let term2 = (one + r) * half;
            let term3 = two;
            term1.min(term2).min(term3)
        }
    }

    fn name(&self) -> &'static str {
        "MonotonizedCentral"
    }
}

/// Upwind limiter (first-order)
///
/// No high-order correction, pure upwind.
/// ψ(r) = 0 for all r
///
/// Properties:
/// * Most stable
/// * First-order accurate
/// * Use for comparison or very difficult flows
#[derive(Debug, Clone, Copy)]
pub struct UpwindLimiter;

impl<T: RealField + Copy> TvdLimiter<T> for UpwindLimiter {
    fn limit(&self, _r: T) -> T {
        T::zero()
    }

    fn name(&self) -> &'static str {
        "Upwind"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_superbee_limiter() {
        let limiter = Superbee;
        
        // Test TVD region for r in [0, 1]
        assert_eq!(limiter.limit(0.0), 0.0); // Upwind at extrema
        assert_eq!(limiter.limit(0.5), 1.0); // min(2*0.5, 1) = 1.0
        assert_eq!(limiter.limit(1.0), 1.0); // min(2, 1) max min(1, 2) = 1
        
        // Test TVD region for r > 1
        assert_eq!(limiter.limit(2.0), 2.0); // min(4, 1) max min(2, 2) = 2
        
        // Test negative r (downwind of extremum)
        assert_eq!(limiter.limit(-1.0), 0.0);
    }

    #[test]
    fn test_van_leer_limiter() {
        let limiter = VanLeer;
        
        assert_eq!(limiter.limit(0.0_f64), 0.0);
        assert!((limiter.limit(1.0_f64) - 1.0).abs() < 1e-10); // 2*1/(1+1) = 1
        assert!((limiter.limit(2.0_f64) - 4.0/3.0).abs() < 1e-10); // 2*2/(1+2) = 4/3
        assert_eq!(limiter.limit(-1.0_f64), 0.0);
    }

    #[test]
    fn test_minmod_limiter() {
        let limiter = Minmod;
        
        assert_eq!(limiter.limit(0.0), 0.0);
        assert_eq!(limiter.limit(0.5), 0.5);
        assert_eq!(limiter.limit(1.0), 1.0);
        assert_eq!(limiter.limit(2.0), 1.0); // Most diffusive
        assert_eq!(limiter.limit(-1.0), 0.0);
    }

    #[test]
    fn test_mc_limiter() {
        let limiter = MonotonizedCentral;
        
        assert_eq!(limiter.limit(0.0_f64), 0.0);
        assert!((limiter.limit(1.0_f64) - 1.0).abs() < 1e-10); // min(2, 1, 2) = 1
        assert!((limiter.limit(2.0_f64) - 1.5).abs() < 1e-10); // min(4, 1.5, 2) = 1.5
        assert_eq!(limiter.limit(-1.0_f64), 0.0);
    }

    #[test]
    fn test_face_interpolation() {
        let limiter = Superbee;
        
        // Uniform region: should return central value
        let face = limiter.interpolate_face(1.0_f64, 1.0, 1.0);
        assert!((face - 1.0).abs() < 1e-10);
        
        // Smooth gradient: r = (2-1)/(3-2) = 1, ψ(1) = 1
        // φ_face = 2 + 0.5 * 1 * (3-2) = 2.5
        let face = limiter.interpolate_face(1.0_f64, 2.0, 3.0);
        assert!((face - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_tvd_property() {
        // Verify all limiters satisfy TVD constraints
        let limiters: Vec<Box<dyn TvdLimiter<f64>>> = vec![
            Box::new(Superbee),
            Box::new(VanLeer),
            Box::new(Minmod),
            Box::new(MonotonizedCentral),
        ];
        
        for limiter in limiters {
            for r_int in 0..100 {
                let r = f64::from(r_int) / 10.0; // Test r from 0 to 10
                let psi = limiter.limit(r);
                
                // TVD constraint: 0 ≤ ψ(r) ≤ 2
                assert!(psi >= 0.0, "{}: ψ({}) = {} < 0", limiter.name(), r, psi);
                assert!(psi <= 2.0, "{}: ψ({}) = {} > 2", limiter.name(), r, psi);
                
                // TVD constraint for r ∈ [0,1]: ψ(r) ≤ 2r
                if r <= 1.0 {
                    assert!(
                        psi <= 2.0 * r + 1e-10,
                        "{}: ψ({}) = {} > 2r = {}",
                        limiter.name(),
                        r,
                        psi,
                        2.0 * r
                    );
                }
            }
        }
    }
}
