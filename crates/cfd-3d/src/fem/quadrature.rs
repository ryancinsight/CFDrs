//! Quadrature rules for finite elements
//!
//! # Theorem — Keast Quadrature Exactness (Keast 1986)
//!
//! A degree-$d$ quadrature rule with $n$ points $\{(\mathbf{x}_k, w_k)\}$
//! integrates all polynomials of total degree $\leq d$ exactly over the
//! reference tetrahedron:
//!
//! ```text
//! ∫_T p(x) dx = Σ_k w_k p(x_k)    ∀ p ∈ P_d
//! ```
//!
//! The 5-point Keast rule used here has degree 3 (integrates cubics exactly),
//! which is sufficient for the bilinear forms arising from P2 Lagrange elements
//! (the integrand $\nabla N_i \cdot \nabla N_j$ is at most degree 2 for P2).
//!
//! **Reference:** Keast, P., "Moderate-degree tetrahedral quadrature formulas",
//! CMAME 55(3), 1986, pp. 339–348.

use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Numerical integration for tetrahedra
pub struct TetrahedronQuadrature<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    points: Vec<Vector3<T>>,
    weights: Vec<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive>
    TetrahedronQuadrature<T>
{
    /// Keast degree 3 quadrature rule (5 points)
    /// Precision O(h^4), enough for quadratic elements
    pub fn keast_degree_3() -> Self {
        let a = <T as FromPrimitive>::from_f64(0.25)
            .expect("0.25 is exactly representable in IEEE 754");
        let b =
            <T as FromPrimitive>::from_f64(0.5).expect("0.5 is exactly representable in IEEE 754");
        let c = <T as FromPrimitive>::from_f64(1.0 / 6.0)
            .expect("1/6 is an IEEE 754 representable f64 constant");

        let p1 = Vector3::new(a, a, a);
        let w1 = <T as FromPrimitive>::from_f64(-0.8)
            .expect("-0.8 is an IEEE 754 representable f64 constant")
            / <T as FromPrimitive>::from_f64(6.0)
                .expect("6.0 is representable in all IEEE 754 types"); // Normalized volume = 1/6

        // Other 4 points are permutations of (1/2, 1/6, 1/6)
        let points = vec![
            p1,
            Vector3::new(b, c, c),
            Vector3::new(c, b, c),
            Vector3::new(c, c, b),
            Vector3::new(c, c, c),
        ];

        let w2 = <T as FromPrimitive>::from_f64(0.45)
            .expect("0.45 is an IEEE 754 representable f64 constant")
            / <T as FromPrimitive>::from_f64(6.0)
                .expect("6.0 is representable in all IEEE 754 types");
        let weights = vec![w1, w2, w2, w2, w2];

        Self { points, weights }
    }

    /// Return the quadrature point coordinates
    pub fn points(&self) -> &[Vector3<T>] {
        &self.points
    }
    /// Return the quadrature weights
    pub fn weights(&self) -> &[T] {
        &self.weights
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Weight sum: Σw_k = 1/6 (volume of the reference tetrahedron).
    ///
    /// Theorem (Keast 1986): A valid quadrature rule on the reference tet
    /// must integrate the constant function 1 exactly, yielding
    /// Σw_k = Vol(T_ref) = 1/6.
    #[test]
    fn weight_sum_equals_reference_tet_volume() {
        let quad = TetrahedronQuadrature::<f64>::keast_degree_3();
        let sum: f64 = quad.weights().iter().sum();
        assert_relative_eq!(sum, 1.0 / 6.0, epsilon = 1e-14);
    }

    /// Linear polynomial exactness: the Keast degree-3 rule integrates
    /// linear functions exactly.
    ///
    /// For f(x,y,z) = x on the reference tet {x+y+z ≤ 1, x,y,z ≥ 0}:
    /// ∫_T x dV = 1/24 (analytically derived from triple integral).
    #[test]
    fn linear_polynomial_exactness() {
        let quad = TetrahedronQuadrature::<f64>::keast_degree_3();
        let mut integral = 0.0_f64;
        for (pt, w) in quad.points().iter().zip(quad.weights().iter()) {
            integral += w * pt.x; // f(x,y,z) = x
        }
        assert_relative_eq!(integral, 1.0 / 24.0, epsilon = 1e-14);
    }
}
