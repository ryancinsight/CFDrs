//! Shape functions for high-order finite elements
//!
//! # Theorem — P2 Lagrange Basis Properties (Ciarlet 1978)
//!
//! The 10-node quadratic Lagrange basis on a tetrahedron satisfies:
//!
//! 1. **Partition of unity:** $\sum_{i=0}^{9} N_i(\mathbf{x}) = 1$ for all $\mathbf{x}$ in the element.
//! 2. **Kronecker delta:** $N_i(\mathbf{x}_j) = \delta_{ij}$ at each node.
//! 3. **Completeness:** The basis spans all polynomials of degree $\leq 2$ on the tetrahedron.
//!
//! **Proof sketch.** Corner functions $N_i = L_i(2L_i - 1)$ and mid-edge functions
//! $N_{ij} = 4 L_i L_j$ form a unisolvent set on the 10-node principal lattice of
//! $\mathbb{P}_2$. Partition of unity follows from the identity
//! $(\sum L_i)^2 = \sum L_i(2L_i - 1) + 2\sum_{i<j} 2L_i L_j = 1$.
//!
//! # Theorem — Quadratic Element Convergence Rate
//!
//! For a P2 finite element solution $u_h$ of a smooth solution $u \in H^3(\Omega)$:
//!
//! ```text
//! ‖u − u_h‖_{H^1} ≤ C h² |u|_{H^3}
//! ```
//!
//! where $h$ is the maximum element diameter and $C$ depends only on the mesh shape regularity.
//!
//! **Reference:** Brenner & Scott, "Math. Theory of FEM", 3rd Ed., Thm. 4.4.20.

use nalgebra::{DMatrix, Matrix3x4, RealField};
use num_traits::{Float, FromPrimitive};

/// Quadratic Lagrange shape functions for 10-node tetrahedra (P2)
pub struct LagrangeTet10<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Corner node gradients (∇L_i) - constant over the tet
    p1_gradients: Matrix3x4<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy + Float>
    LagrangeTet10<T>
{
    /// Create new P2 shape functions from elemental P1 gradients
    pub fn new(p1_gradients: Matrix3x4<T>) -> Self {
        Self { p1_gradients }
    }

    /// Calculate shape function values N_i at barycentric coordinates L
    /// L = [L1, L2, L3, L4] where sum(L) = 1.0
    pub fn values(&self, l: &[T; 4]) -> [T; 10] {
        let mut n = [T::zero(); 10];
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");
        let four = <T as FromPrimitive>::from_f64(4.0)
            .expect("4.0 is representable in all IEEE 754 types");

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
        let four = <T as FromPrimitive>::from_f64(4.0)
            .expect("4.0 is representable in all IEEE 754 types");

        // Corner nodes (0-3): ∇Ni = (4Li - 1) ∇Li
        for (i, &li) in l.iter().enumerate() {
            let factor = four * li - T::one();
            let g_i = self.p1_gradients.column(i) * factor;
            grad.set_column(i, &g_i);
        }

        // Mid-edge nodes (4-9): ∇Nij = 4(Li ∇Lj + Lj ∇Li)
        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        for (idx, &(i, j)) in edges.iter().enumerate() {
            let g_ij =
                (self.p1_gradients.column(j) * l[i] + self.p1_gradients.column(i) * l[j]) * four;
            grad.set_column(4 + idx, &g_ij);
        }

        grad
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Reference tetrahedron P1 gradients: vertices at
    /// v0=(0,0,0), v1=(1,0,0), v2=(0,1,0), v3=(0,0,1).
    /// Barycentric coords: L1=1-x-y-z, L2=x, L3=y, L4=z.
    /// ∇L1=(-1,-1,-1), ∇L2=(1,0,0), ∇L3=(0,1,0), ∇L4=(0,0,1).
    fn reference_tet_shape() -> LagrangeTet10<f64> {
        let p1_gradients = Matrix3x4::new(
            -1.0, 1.0, 0.0, 0.0, // row x
            -1.0, 0.0, 1.0, 0.0, // row y
            -1.0, 0.0, 0.0, 1.0, // row z
        );
        LagrangeTet10::new(p1_gradients)
    }

    /// Partition of unity: ΣN_i(L) = 1 for any L with Σ L_k = 1.
    ///
    /// Theorem (Ciarlet 1978): The P2 Lagrange basis satisfies
    /// Σ_{i=0}^{9} N_i(x) = 1 for all x in the element, which follows
    /// from (ΣL_i)² = 1.
    #[test]
    fn partition_of_unity() {
        let shape = reference_tet_shape();

        // Interior point
        let l = [0.1, 0.2, 0.3, 0.4];
        let n = shape.values(&l);
        let sum: f64 = n.iter().sum();
        assert_relative_eq!(sum, 1.0, epsilon = 1e-14);

        // Centroid
        let l_centroid = [0.25, 0.25, 0.25, 0.25];
        let n_c = shape.values(&l_centroid);
        let sum_c: f64 = n_c.iter().sum();
        assert_relative_eq!(sum_c, 1.0, epsilon = 1e-14);

        // Near-vertex
        let l_near = [0.9, 0.05, 0.03, 0.02];
        let n_nv = shape.values(&l_near);
        let sum_nv: f64 = n_nv.iter().sum();
        assert_relative_eq!(sum_nv, 1.0, epsilon = 1e-14);
    }

    /// Kronecker delta: N_i(x_j) = δ_{ij} at the 10 nodal points.
    ///
    /// The 10 nodes of the P2 tet are: 4 corners (where exactly one L_k=1)
    /// and 6 mid-edge points (where exactly two L_k = 0.5).
    #[test]
    fn kronecker_delta_at_nodes() {
        let shape = reference_tet_shape();

        // Corner nodes: L = (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)
        let corners: [[f64; 4]; 4] = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];

        for (node_idx, l) in corners.iter().enumerate() {
            let n = shape.values(l);
            for i in 0..10 {
                let expected = if i == node_idx { 1.0 } else { 0.0 };
                assert_relative_eq!(n[i], expected, epsilon = 1e-14,);
            }
        }

        // Mid-edge nodes: edges (0,1),(1,2),(2,0),(0,3),(1,3),(2,3)
        // with L at midpoint = 0.5 for each endpoint
        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        for (edge_idx, &(a, b)) in edges.iter().enumerate() {
            let mut l = [0.0_f64; 4];
            l[a] = 0.5;
            l[b] = 0.5;
            let n = shape.values(&l);
            let mid_node = 4 + edge_idx;
            for i in 0..10 {
                let expected = if i == mid_node { 1.0 } else { 0.0 };
                assert_relative_eq!(n[i], expected, epsilon = 1e-14,);
            }
        }
    }

    /// Gradient consistency: Σ∇N_i = 0 (constant completeness).
    ///
    /// Since Σ N_i = 1 everywhere, differentiating gives Σ ∇N_i = 0.
    #[test]
    fn gradient_sum_is_zero() {
        let shape = reference_tet_shape();
        let l = [0.15, 0.25, 0.35, 0.25];
        let grad = shape.gradients(&l);

        // Sum over all 10 shape function gradients (columns)
        for row in 0..3 {
            let sum: f64 = (0..10).map(|col| grad[(row, col)]).sum();
            assert_relative_eq!(sum, 0.0, epsilon = 1e-13);
        }
    }
}
