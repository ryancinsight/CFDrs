//! Quadrature rules for finite elements

use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Numerical integration for tetrahedra
pub struct TetrahedronQuadrature<T: RealField + Copy> {
    points: Vec<Vector3<T>>,
    weights: Vec<T>,
}

impl<T: RealField + Copy + FromPrimitive> TetrahedronQuadrature<T> {
    /// Keast degree 3 quadrature rule (5 points)
    /// Precision O(h^4), enough for quadratic elements
    pub fn keast_degree_3() -> Self {
        let a = T::from_f64(0.25).unwrap();
        let b = T::from_f64(0.5).unwrap();
        let c = T::from_f64(1.0/6.0).unwrap();
        
        let p1 = Vector3::new(a, a, a);
        let w1 = T::from_f64(-0.8).unwrap() / T::from_f64(6.0).unwrap(); // Normalized volume = 1/6
        
        // Other 4 points are permutations of (1/2, 1/6, 1/6)
        let points = vec![
            p1,
            Vector3::new(b, c, c),
            Vector3::new(c, b, c),
            Vector3::new(c, c, b),
            Vector3::new(c, c, c),
        ];
        
        let w2 = T::from_f64(0.45).unwrap() / T::from_f64(6.0).unwrap();
        let weights = vec![w1, w2, w2, w2, w2];
        
        Self { points, weights }
    }

    pub fn points(&self) -> &[Vector3<T>] { &self.points }
    pub fn weights(&self) -> &[T] { &self.weights }
}
