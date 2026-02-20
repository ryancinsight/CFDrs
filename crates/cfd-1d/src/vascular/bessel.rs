//! Exact evaluation of complex Bessel functions
//!
//! # Mathematical Foundation
//! The Bessel function of the first kind of order $\nu$, $J_\nu(z)$, is defined by the power series:
//! ```text
//! J_ν(z) = \sum_{m=0}^∞ \frac{(-1)^m}{m! Γ(m+ν+1)} (\frac{z}{2})^{2m+ν}
//! ```
//!
//! For Womersley flow, we require $J_0(z)$ and $J_1(z)$ for complex arguments $z \in \mathbb{C}$.
//!
//! ## Convergence and Precision Proof
//! The power series converges absolutely for all $z \in \mathbb{C}$. For $|z| \le 25$ (which covers typical physiological
//! Womersley numbers $\alpha \le 25$), maximal intermediate terms are on the order of $|z/2|^{20} / 10! \approx 10^7$.
//! Using double precision (`f64`), which carries ~15-17 significant decimal digits, catastrophic cancellation will
//! strip at most ~7-8 digits. This guarantees at least 7-9 significant digits of physical precision in the worst
//! physiological case ($\alpha = 25$), and near full precision for typical arteries ($\alpha \approx 3-10$).
//!
//! # Implementations
//!
//! - `bessel_j0(z)` computes $J_0(z)$
//! - `bessel_j1(z)` computes $J_1(z)$

use nalgebra::{Complex, ComplexField, RealField};
use num_traits::FromPrimitive;

/// Bessel function of the first kind, order zero: $J_0(z)$
///
/// Evaluates the exact infinite series:
/// $J_0(z) = \sum_{m=0}^\infty \frac{(-1)^m}{(m!)^2} \left(\frac{z}{2}\right)^{2m}$
pub fn bessel_j0<T: RealField + FromPrimitive + Copy>(z: Complex<T>) -> Complex<T> {
    let mut sum = Complex::new(T::one(), T::zero());
    let mut term = Complex::new(T::one(), T::zero());
    let z_half = z / T::from_f64(2.0).unwrap();
    let z_half_sq = z_half * z_half;

    let mut m = 1_u32;
    let max_iter = 150_u32; // Guaranteed convergence for |z| < 30
    let tolerance = T::from_f64(1e-15).unwrap();

    while m < max_iter {
        let m_t = T::from_u32(m).unwrap();
        // term(m) = term(m-1) * (-z^2 / 4) / m^2
        term = term * (-z_half_sq) / (m_t * m_t);
        sum += term;

        if term.modulus() < tolerance {
            break;
        }
        m += 1;
    }
    sum
}

/// Bessel function of the first kind, order one: $J_1(z)$
///
/// Evaluates the exact infinite series:
/// $J_1(z) = \sum_{m=0}^\infty \frac{(-1)^m}{m!(m+1)!} \left(\frac{z}{2}\right)^{2m+1}$
pub fn bessel_j1<T: RealField + FromPrimitive + Copy>(z: Complex<T>) -> Complex<T> {
    let z_half = z / T::from_f64(2.0).unwrap();
    let mut sum = z_half;
    let mut term = z_half;
    let z_half_sq = z_half * z_half;

    let mut m = 1_u32;
    let max_iter = 150_u32;
    let tolerance = T::from_f64(1e-15).unwrap();

    while m < max_iter {
        let m_t = T::from_u32(m).unwrap();
        let m_plus_1 = T::from_u32(m + 1).unwrap();
        // term(m) = term(m-1) * (-z^2 / 4) / (m * (m+1))
        term = term * (-z_half_sq) / (m_t * m_plus_1);
        sum += term;

        if term.modulus() < tolerance {
            break;
        }
        m += 1;
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bessel_j0_real() {
        // J_0(0) = 1
        let z0 = Complex::new(0.0, 0.0);
        assert_relative_eq!(bessel_j0(z0).re, 1.0, epsilon = 1e-12);
        assert_relative_eq!(bessel_j0(z0).im, 0.0, epsilon = 1e-12);

        // J_0(2.4048255576) approx 0 (first root)
        let z_root = Complex::new(2.404825557695773, 0.0);
        assert_relative_eq!(bessel_j0(z_root).re, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_bessel_j1_real() {
        // J_1(0) = 0
        let z0 = Complex::new(0.0, 0.0);
        assert_relative_eq!(bessel_j1(z0).re, 0.0, epsilon = 1e-12);
        
        // J_1(3.8317059702) approx 0 (first non-zero root)
        let z_root = Complex::new(3.831705970207512, 0.0);
        assert_relative_eq!(bessel_j1(z_root).re, 0.0, epsilon = 1e-5);
    }
}
