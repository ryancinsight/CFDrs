//! Element-level operations for FEM
//!
//! # Theorem — Tetrahedron Volume (Euler 1752)
//!
//! Given a tetrahedron with vertices $\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3$,
//! its signed volume is
//!
//! ```text
//! 6V = det(v₁ − v₀, v₂ − v₀, v₃ − v₀)
//! ```
//!
//! The absolute value gives the positive volume $|V|$. A positive determinant
//! indicates right-hand orientation of the vertex ordering.
//!
//! **Proof sketch.** The volume of a parallelepiped spanned by three edge
//! vectors is the scalar triple product. A tetrahedron occupies exactly
//! $1/6$ of the parallelepiped, giving the factor of 6.
//!
//! # Theorem — Shape Function Gradients via Cofactor Matrix
//!
//! For linear (P1) tetrahedral elements, the shape function gradient $\nabla N_i$
//! is constant over the element and equals the cofactor of the $i$-th row of the
//! 4×4 coordinate matrix divided by $6V$:
//!
//! ```text
//! ∂N_i/∂x_j = (-1)^{i+j+1} M_{i,j} / (6V)
//! ```
//!
//! where $M_{i,j}$ is the $(i,j)$ minor of the coordinate matrix.
//!
//! **Reference:** Zienkiewicz & Taylor, "The Finite Element Method", Vol. 1, 6th Ed., §7.3.

use crate::fem::constants;
use nalgebra::{DMatrix, DVector, Matrix3, RealField, Vector3};
use num_traits::{Float, FromPrimitive};

/// Element matrices for FEM assembly
#[derive(Debug, Clone)]
pub struct ElementMatrices<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Element stiffness matrix
    pub k_e: DMatrix<T>,
    /// Element mass matrix
    pub m_e: DMatrix<T>,
    /// Element force vector
    pub f_e: DMatrix<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> ElementMatrices<T> {
    /// Create new element matrices
    #[must_use]
    pub fn new(n_dof: usize) -> Self {
        Self {
            k_e: DMatrix::zeros(n_dof, n_dof),
            m_e: DMatrix::zeros(n_dof, n_dof),
            f_e: DMatrix::zeros(n_dof, 1),
        }
    }
}

/// Fluid element for FEM calculations
#[derive(Debug, Clone)]
pub struct FluidElement<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Node indices
    pub nodes: Vec<usize>,
    /// Element volume
    pub volume: T,
    /// Shape function derivatives (3 rows for x/y/z, 4 columns for nodes 1-4 for P1, or 10 columns for P2)
    pub shape_derivatives: DMatrix<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy + Float> FluidElement<T> {
    /// Create new fluid element
    #[must_use]
    pub fn new(nodes: Vec<usize>) -> Self {
        let n = nodes.len();
        Self {
            nodes,
            volume: T::zero(),
            shape_derivatives: DMatrix::zeros(3, n),
        }
    }

    /// Calculate element volume (for tetrahedron) using nodal coordinates.
    ///
    /// # Theorem (Tetrahedron Volume from Jacobian)
    ///
    /// For a tetrahedron with vertices $v_0, v_1, v_2, v_3$, the signed volume is
    ///
    /// $$V = \frac{1}{6} \det\bigl[\,v_1 - v_0 \;\big|\; v_2 - v_0 \;\big|\; v_3 - v_0\,\bigr]$$
    ///
    /// **Proof sketch**: The three edge vectors emanating from $v_0$ form a
    /// parallelepiped whose volume is $|\det J|$. A tetrahedron occupies exactly
    /// $1/6$ of this parallelepiped (combinatorial identity for simplices in
    /// $\mathbb{R}^3$). The sign of $\det J$ encodes orientation.
    pub fn calculate_volume(&mut self, vertices: &[Vector3<T>]) -> T {
        if vertices.len() < 4 {
            return T::zero();
        }

        let v0 = &vertices[0];
        let v1 = &vertices[1];
        let v2 = &vertices[2];
        let v3 = &vertices[3];

        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let e3 = v3 - v0;

        let det_j = e1.dot(&e2.cross(&e3)); // Signed 6V
        self.volume = num_traits::Float::abs(det_j)
            / <T as FromPrimitive>::from_f64(6.0)
                .expect("6.0 is representable in all IEEE 754 types");
        det_j // Return signed 6V for orientation-correct assembly
    }

    /// Calculate shape function derivatives for tetrahedral element.
    ///
    /// # Theorem (Isoparametric Shape Function Gradients)
    ///
    /// For linear tetrahedral (P1) elements, the shape function gradients are
    /// constant within each element and given by $\nabla N_i = J^{-T} \hat{\nabla} \hat{N}_i$
    /// where $J$ is the Jacobian of the reference-to-physical mapping. Since
    /// $\hat{N}_i$ are barycentric coordinates with constant gradients on the
    /// reference tetrahedron, $\nabla N_i$ is also constant per element.
    ///
    /// **Proof sketch**: The affine map $x = v_0 + J\hat{x}$ with
    /// $J = [v_1 - v_0 | v_2 - v_0 | v_3 - v_0]$ transforms the reference
    /// simplex to the physical element. Applying the chain rule:
    /// $\partial N_i / \partial x_j = (J^{-T})_{jk}\,\partial\hat{N}_i / \partial\hat{x}_k$.
    /// The cofactor formula gives $J^{-T}$ from cross products of edge vectors
    /// divided by $6V$.
    ///
    /// **Reference**: Zienkiewicz & Taylor (2005), *The Finite Element Method*, Vol. 1, 6th Ed., §8.5.
    pub fn calculate_shape_derivatives(&mut self, vertices: &[Vector3<T>]) {
        if vertices.len() < 4 {
            return;
        }

        let v1 = vertices[0];
        let v2 = vertices[1];
        let v3 = vertices[2];
        let v4 = vertices[3];

        // Determinant of Jacobian matrix (signed 6V)
        let d1_vec = v2 - v1;
        let d2_vec = v3 - v1;
        let d3_vec = v4 - v1;
        let six_v = d1_vec.dot(&d2_vec.cross(&d3_vec));

        if num_traits::Float::abs(six_v)
            < <T as FromPrimitive>::from_f64(1e-24)
                .expect("1e-24 is an IEEE 754 representable f64 constant")
        {
            self.shape_derivatives = DMatrix::zeros(3, 4);
            return;
        }

        // Calculate derivatives using signed cofactor matrix of the 4×4 coordinate
        // determinant (Zienkiewicz & Taylor, Vol 1, 6th Ed, §7.3).
        //
        // For node i of a tet with 0-based indexing in the 4×4 coordinate matrix:
        //   ∂N_i/∂x = B_i / (6V),  where B_i = (-1)^{i+1} * M_{i,1}
        //   ∂N_i/∂y = C_i / (6V),  where C_i = (-1)^{i+2} * M_{i,2}
        //   ∂N_i/∂z = D_i / (6V),  where D_i = (-1)^{i+3} * M_{i,3}
        //
        // The raw minor expressions below compute M_{i,j}. The alternating cofactor
        // sign (-1)^{i+j} is applied via the leading negation on each vector.

        // Face opposite to node 0 (v1) -> grad_N0: cofactor signs are (-,+,-)
        let grad_n0 = -Vector3::new(
            (v2.y * (v3.z - v4.z) + v3.y * (v4.z - v2.z) + v4.y * (v2.z - v3.z)) / six_v,
            (v2.x * (v4.z - v3.z) + v3.x * (v2.z - v4.z) + v4.x * (v3.z - v2.z)) / six_v,
            (v2.x * (v3.y - v4.y) + v3.x * (v4.y - v2.y) + v4.x * (v2.y - v3.y)) / six_v,
        );
        // Face opposite to node 1 (v2) -> grad_N1: cofactor signs are (+,-,+)
        let grad_n1 = -Vector3::new(
            (v1.y * (v4.z - v3.z) + v3.y * (v1.z - v4.z) + v4.y * (v3.z - v1.z)) / six_v,
            (v1.x * (v3.z - v4.z) + v3.x * (v4.z - v1.z) + v4.x * (v1.z - v3.z)) / six_v,
            (v1.x * (v4.y - v3.y) + v3.x * (v1.y - v4.y) + v4.x * (v3.y - v1.y)) / six_v,
        );
        // Face opposite to node 2 (v3) -> grad_N2: cofactor signs are (-,+,-)
        let grad_n2 = -Vector3::new(
            (v1.y * (v2.z - v4.z) + v2.y * (v4.z - v1.z) + v4.y * (v1.z - v2.z)) / six_v,
            (v1.x * (v4.z - v2.z) + v2.x * (v1.z - v4.z) + v4.x * (v2.z - v1.z)) / six_v,
            (v1.x * (v2.y - v4.y) + v2.x * (v4.y - v1.y) + v4.x * (v1.y - v2.y)) / six_v,
        );

        // Sum of gradients is zero: grad_N3 = -(grad_N0 + grad_N1 + grad_N2)
        let grad_n3 = -(grad_n0 + grad_n1 + grad_n2);

        // CRITICAL FIX: Columns must correspond to nodes [0, 1, 2, 3]
        let c0 = DVector::from_iterator(3, grad_n0.iter().copied());
        let c1 = DVector::from_iterator(3, grad_n1.iter().copied());
        let c2 = DVector::from_iterator(3, grad_n2.iter().copied());
        let c3 = DVector::from_iterator(3, grad_n3.iter().copied());
        self.shape_derivatives = DMatrix::from_columns(&[c0, c1, c2, c3]);
    }

    /// Calculate element stiffness matrix contribution for Stokes flow
    /// Based on Hughes, "The Finite Element Method", Dover Publications, 2000
    pub fn stiffness_contribution(&self, viscosity: T) -> DMatrix<T> {
        let n_nodes = self.nodes.len();
        let n_dof = n_nodes * constants::VELOCITY_COMPONENTS;
        let mut k_e = DMatrix::zeros(n_dof, n_dof);

        if n_nodes == 4 {
            // Tetrahedral element: Single-point integration is exact for linear elements
            let factor = viscosity * self.volume;

            for i in 0..4 {
                for j in 0..4 {
                    // Viscous stiffness block (velocity-velocity coupling)
                    for d in 0..3 {
                        let row = i * 3 + d;
                        let col = j * 3 + d;

                        // Viscous contribution: μ ∫ (∇Ni : ∇Nj) dΩ
                        let mut visc_term = T::zero();
                        for k in 0..3 {
                            visc_term +=
                                self.shape_derivatives[(k, i)] * self.shape_derivatives[(k, j)];
                        }
                        k_e[(row, col)] += factor * visc_term;

                        // Additional viscous terms from strain rate tensor: μ ∫ (∂Ni/∂xk)(∂Nj/∂xk)
                        for k in 0..3 {
                            if k != d {
                                let cross_term = self.shape_derivatives[(d, i)]
                                    * self.shape_derivatives[(k, j)]
                                    + self.shape_derivatives[(k, i)]
                                        * self.shape_derivatives[(d, j)];
                                k_e[(row, col)] += factor
                                    * cross_term
                                    * <T as FromPrimitive>::from_f64(0.5)
                                        .expect("0.5 is exactly representable in IEEE 754");
                            }
                        }
                    }
                }
            }
        }
        k_e
    }

    /// Calculate strain rate from velocity gradient
    pub fn strain_rate(&self, velocity_gradient: &Matrix3<T>) -> Matrix3<T> {
        let half =
            <T as FromPrimitive>::from_f64(0.5).expect("0.5 is exactly representable in IEEE 754");
        (velocity_gradient + velocity_gradient.transpose()) * half
    }
}
