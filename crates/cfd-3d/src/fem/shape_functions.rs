//! Shape functions for high-order finite elements

use nalgebra::{Matrix3x4, RealField, DMatrix};
use num_traits::{Float, FromPrimitive};

/// Quadratic Lagrange shape functions for 10-node tetrahedra (P2)
pub struct LagrangeTet10<T: RealField + Copy> {
    /// Corner node gradients (∇L_i) - constant over the tet
    p1_gradients: Matrix3x4<T>,
}

impl<T: RealField + FromPrimitive + Copy + Float> LagrangeTet10<T> {
    /// Create new P2 shape functions from elemental P1 gradients
    pub fn new(p1_gradients: Matrix3x4<T>) -> Self {
        Self { p1_gradients }
    }

    /// Calculate shape function values N_i at barycentric coordinates L
    /// L = [L1, L2, L3, L4] where sum(L) = 1.0
    pub fn values(&self, l: &[T; 4]) -> [T; 10] {
        let mut n = [T::zero(); 10];
        let two = T::from_f64(2.0).unwrap();
        let four = T::from_f64(4.0).unwrap();

        // Corner nodes (0-3)
        for i in 0..4 {
            n[i] = (two * l[i] - T::one()) * l[i];
        }

        // Mid-edge nodes (4-9)
        // 4: Mid(0, 1)
        n[4] = four * l[0] * l[1];
        // 5: Mid(1, 2)
        n[5] = four * l[1] * l[2];
        // 6: Mid(2, 0)
        n[6] = four * l[2] * l[0];
        // 7: Mid(0, 3)
        n[7] = four * l[0] * l[3];
        // 8: Mid(1, 3)
        n[8] = four * l[1] * l[3];
        // 9: Mid(2, 3)
        n[9] = four * l[2] * l[3];

        n
    }

    /// Calculate shape function gradients ∇N_i at barycentric coordinates L
    /// Returns 3x10 matrix (rows for x,y,z; columns for nodes 0-9)
    pub fn gradients(&self, l: &[T; 4]) -> DMatrix<T> {
        let mut grad = DMatrix::zeros(3, 10);
        let four = T::from_f64(4.0).unwrap();

        // Corner nodes (0-3): ∇Ni = (4Li - 1) ∇Li
        for i in 0..4 {
            let factor = four * l[i] - T::one();
            let g_i = self.p1_gradients.column(i) * factor;
            grad.set_column(i, &g_i);
        }

        // Mid-edge nodes (4-9): ∇Nij = 4(Li ∇Lj + Lj ∇Li)
        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        for (idx, &(i, j)) in edges.iter().enumerate() {
            let g_ij = (self.p1_gradients.column(j) * l[i] + self.p1_gradients.column(i) * l[j]) * four;
            grad.set_column(4 + idx, &g_ij);
        }

        grad
    }
}
