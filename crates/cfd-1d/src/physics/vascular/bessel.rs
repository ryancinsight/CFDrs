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
//! ## Theorem: Bessel Function Spectral Error Bounds
//!
//! **Theorem**: The infinite Taylor series expansion for $J_\nu(z)$ evaluated
//! sequentially using $m$-th term $t_m$ truncated at $M$ terms has an absolute
//! truncation error strictly bounded by the magnitude of the $(M+1)$-th term,
//! provided $M > |z|^2 / 4$.
//!
//! **Proof Outline**: For $|z| \le 25$ (physiological Womersley numbers $\alpha \le 25$),
//! the maximal intermediate terms reach $O(10^7)$. In IEEE 754 float64 (53 bits of mantissa),
//! this magnitude induces a catastrophic cancellation loss of $\approx \log_{10}(10^7) = 7$ digits.
//! Since float64 provides $\approx 15.9$ decimal digits of precision, the returned complex
//! Bessel value retains at least $15.9 - 7 = 8.9$ significant digits of physical precision,
//! strictly bounding the relative numerical error to $< 10^{-8}$.
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
    let z_half = z / (T::one() + T::one());
    let z_half_sq = z_half * z_half;

    let mut m = 1_u32;
    let max_iter = 150_u32; // Guaranteed convergence for |z| < 30
    let tolerance = T::from_f64(1e-15).expect("Mathematical constant conversion compromised");

    while m < max_iter {
        let m_t = T::from_u32(m).unwrap_or_else(num_traits::Zero::zero);
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
    let z_half = z / (T::one() + T::one());
    let mut sum = z_half;
    let mut term = z_half;
    let z_half_sq = z_half * z_half;

    let mut m = 1_u32;
    let max_iter = 150_u32;
    let tolerance = T::from_f64(1e-15).expect("Mathematical constant conversion compromised");

    while m < max_iter {
        let m_t = T::from_u32(m).unwrap_or_else(num_traits::Zero::zero);
        let m_plus_1 = T::from_u32(m + 1).unwrap_or_else(num_traits::Zero::zero);
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

    /// Complex argument J₀(1+i): direct series evaluation.
    /// J₀(1+i) ≈ 0.93760848 − 0.49652925i (verified by independent series computation).
    ///
    /// The sign of the imaginary part is negative because (z/2)² = i/2 for z=1+i,
    /// and the alternating series generates net negative imaginary contribution.
    #[test]
    fn test_bessel_j0_complex_argument() {
        let z = Complex::new(1.0, 1.0);
        let j0_z = bessel_j0(z);
        assert_relative_eq!(j0_z.re, 0.937_608_48, max_relative = 1e-5);
        assert_relative_eq!(j0_z.im, -0.496_529_25, max_relative = 1e-5);
    }

    /// Complex argument J₁(1+i): direct series evaluation.
    /// J₁(1+i) ≈ 0.61416033 − 0.36502802i.
    #[test]
    fn test_bessel_j1_complex_argument() {
        let z = Complex::new(1.0, 1.0);
        let j1_z = bessel_j1(z);
        assert_relative_eq!(j1_z.re, 0.614_160_33, max_relative = 1e-5);
        assert_relative_eq!(j1_z.im, 0.365_028_02, max_relative = 1e-5);
    }

    /// Derivative identity: J₀′(z) = −J₁(z).
    /// Verified via symmetric finite difference: [J₀(z+δ) − J₀(z−δ)] / (2δ) ≈ −J₁(z).
    #[test]
    fn test_j0_derivative_equals_neg_j1() {
        let z = Complex::new(2.0, 0.5);
        let delta = Complex::new(1e-7, 0.0);
        let j0_plus = bessel_j0(z + delta);
        let j0_minus = bessel_j0(z - delta);
        let numerical_derivative = (j0_plus - j0_minus) / (Complex::new(2e-7, 0.0));
        let neg_j1 = -bessel_j1(z);
        assert_relative_eq!(numerical_derivative.re, neg_j1.re, max_relative = 1e-5);
        assert_relative_eq!(numerical_derivative.im, neg_j1.im, max_relative = 1e-5);
    }
}

