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
