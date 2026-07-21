//! Tetrahedral quadrature rules (Keast rules)
//!
//! References:
//! Keast, P. (1986). "Moderate-degree quadrature rules for simplices".
//! Computer Methods in Applied Mechanics and Engineering.

use crate::integration::traits::Quadrature3D;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, RealField};

fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

/// Tetrahedral quadrature rule
pub struct TetrahedronQuadrature<T: RealField + Copy> {
    points: Vec<[T; 4]>,
    weights: Vec<T>,
    degree: usize,
}

impl<T: RealField + Copy + FloatElement> TetrahedronQuadrature<T> {
    /// Create a new tetrahedral quadrature rule of specified degree
    /// Supported degrees: 1, 2, 5
    pub fn new(degree: usize) -> Result<Self> {
        match degree {
            1 => Ok(Self::degree_1()),
            2 => Ok(Self::degree_2()),
            3 => Ok(Self::degree_3()),
            5 => Ok(Self::degree_5()),
            _ => Err(Error::InvalidInput(format!(
                "Unsupported tetrahedral quadrature degree: {degree}. Use 1, 2, 3, or 5."
            ))),
        }
    }

    /// Degree 1 (1 point) - exact for linear polynomials
    fn degree_1() -> Self {
        let q = from_f64(0.25);
        Self {
            points: vec![[q, q, q, q]],
            weights: vec![T::ONE],
            degree: 1,
        }
    }

    /// Degree 2 (4 points) - exact for quadratic polynomials
    fn degree_2() -> Self {
        let alpha = from_f64(0.5854101966249685);
        let beta = from_f64(0.1381966011250105);
        let w = from_f64(0.25);

        Self {
            points: vec![
                [alpha, beta, beta, beta],
                [beta, alpha, beta, beta],
                [beta, beta, alpha, beta],
                [beta, beta, beta, alpha],
            ],
            weights: vec![w, w, w, w],
            degree: 2,
        }
    }

    /// Degree 3 (5 points)
    fn degree_3() -> Self {
        let w1 = from_f64(-0.8);
        let w2 = from_f64(0.45);
        let q = from_f64(0.25);
        let a = from_f64(1.0 / 6.0);
        let b = from_f64(0.5);

        Self {
            points: vec![
                [q, q, q, q],
                [b, a, a, a],
                [a, b, a, a],
                [a, a, b, a],
                [a, a, a, b],
            ],
            weights: vec![w1, w2, w2, w2, w2],
            degree: 3,
        }
    }

    /// Degree 5 (15 points) - Keast Rule 5
    fn degree_5() -> Self {
        // Point 1 (1)
        let w1 = from_f64(0.030283678097089 * 6.0);
        let p1 = from_f64(1.0 / 4.0);

        // Point 2 (4)
        let w2 = from_f64(0.010959110309130 * 6.0);
        let x2_a = from_f64(0.0543467785602507);
        let x2_b = from_f64(0.3152177404799164);

        // Point 3 (4)
        let w3 = from_f64(0.013436499589410 * 6.0);
        let x3_a = from_f64(0.702284731154388);
        let x3_b = from_f64(0.0992384229485373);

        // Point 4 (6)
        let w4 = from_f64(0.007022847311543 * 6.0);
        let x4_a = from_f64(0.454763304449814);
        let x4_b = from_f64(0.0452366955501859);

        let mut points = Vec::with_capacity(15);
        let mut weights = Vec::with_capacity(15);

        // 1. One point (1/4, 1/4, 1/4, 1/4)
        points.push([p1, p1, p1, p1]);
        weights.push(w1);

        // 2. Four points (a, b, b, b) and perms
        points.push([x2_a, x2_b, x2_b, x2_b]);
        points.push([x2_b, x2_a, x2_b, x2_b]);
        points.push([x2_b, x2_b, x2_a, x2_b]);
        points.push([x2_b, x2_b, x2_b, x2_a]);
        for _ in 0..4 {
            weights.push(w2);
        }

        // 3. Four points (c, d, d, d) and perms
        points.push([x3_a, x3_b, x3_b, x3_b]);
        points.push([x3_b, x3_a, x3_b, x3_b]);
        points.push([x3_b, x3_b, x3_a, x3_b]);
        points.push([x3_b, x3_b, x3_b, x3_a]);
        for _ in 0..4 {
            weights.push(w3);
        }

        // 4. Six points (e, e, f, f) and perms
        points.push([x4_a, x4_a, x4_b, x4_b]);
        points.push([x4_a, x4_b, x4_a, x4_b]);
        points.push([x4_a, x4_b, x4_b, x4_a]);
        points.push([x4_b, x4_a, x4_a, x4_b]);
        points.push([x4_b, x4_a, x4_b, x4_a]);
        points.push([x4_b, x4_b, x4_a, x4_a]);
        for _ in 0..6 {
            weights.push(w4);
        }

        Self {
            points,
            weights,
            degree: 5,
        }
    }
}

impl<T: RealField + Copy> Quadrature3D<T> for TetrahedronQuadrature<T> {
    fn points(&self) -> &[[T; 4]] {
        &self.points
    }

    fn weights(&self) -> &[T] {
        &self.weights
    }

    fn order(&self) -> usize {
        self.degree
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn degree_one_rule_has_centroid_and_unit_weight() {
        let rule = TetrahedronQuadrature::<f64>::new(1).unwrap();

        assert_eq!(rule.num_points(), 1);
        assert_eq!(rule.order(), 1);
        assert_eq!(rule.points()[0], [0.25, 0.25, 0.25, 0.25]);
        assert_eq!(rule.weights(), &[1.0]);
    }

    #[test]
    fn degree_two_rule_weights_sum_to_one() {
        let rule = TetrahedronQuadrature::<f64>::new(2).unwrap();
        let weight_sum: f64 = rule.weights().iter().copied().sum();

        assert_relative_eq!(weight_sum, 1.0, epsilon = 1e-15);
        assert_eq!(rule.num_points(), 4);
    }

    #[test]
    fn tetrahedron_quadrature_rejects_unsupported_degree() {
        let Err(error) = TetrahedronQuadrature::<f64>::new(4) else {
            panic!("expected degree 4 to be rejected");
        };

        assert!(
            error.to_string().contains("Use 1, 2, 3, or 5"),
            "unexpected error: {error}"
        );
    }
}
