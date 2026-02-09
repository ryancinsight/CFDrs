//! Element-level operations for FEM

use crate::fem::constants;
use nalgebra::{DMatrix, Matrix3, Matrix3x4, RealField, Vector3};
use num_traits::{Float, FromPrimitive};

/// Element matrices for FEM assembly
#[derive(Debug, Clone)]
pub struct ElementMatrices<T: RealField + Copy> {
    /// Element stiffness matrix
    pub k_e: DMatrix<T>,
    /// Element mass matrix
    pub m_e: DMatrix<T>,
    /// Element force vector
    pub f_e: DMatrix<T>,
}

impl<T: RealField + FromPrimitive + Copy> ElementMatrices<T> {
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
pub struct FluidElement<T: RealField + Copy> {
    /// Node indices
    pub nodes: Vec<usize>,
    /// Element volume
    pub volume: T,
    /// Shape function derivatives (3 rows for x/y/z, 4 columns for nodes 1-4)
    pub shape_derivatives: Matrix3x4<T>,
}

impl<T: RealField + FromPrimitive + Copy + Float> FluidElement<T> {
    /// Create new fluid element
    #[must_use]
    pub fn new(nodes: Vec<usize>) -> Self {
        Self {
            nodes,
            volume: T::zero(),
            shape_derivatives: Matrix3x4::zeros(),
        }
    }

    /// Calculate element volume (for tetrahedron) using nodal coordinates
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
        self.volume = num_traits::Float::abs(det_j) / T::from_f64(6.0).unwrap_or_else(T::one);
        det_j // Return signed 6V for orientation-correct assembly
    }

    /// Calculate shape function derivatives for tetrahedral element
    /// Based on Zienkiewicz & Taylor, "The Finite Element Method", Vol. 1, 6th Ed.
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

        if num_traits::Float::abs(six_v) < T::from_f64(1e-24).unwrap_or_else(T::zero) {
            self.shape_derivatives = Matrix3x4::zeros();
            return;
        }

        // Calculate derivatives using signed cofactor matrix
        // d1 corresponding to face opposite to node 0 (v1) -> grad_N0
        let grad_n0 = Vector3::new(
            (v2.y * (v3.z - v4.z) + v3.y * (v4.z - v2.z) + v4.y * (v2.z - v3.z)) / six_v,
            (v2.x * (v4.z - v3.z) + v3.x * (v2.z - v4.z) + v4.x * (v3.z - v2.z)) / six_v,
            (v2.x * (v3.y - v4.y) + v3.x * (v4.y - v2.y) + v4.x * (v2.y - v3.y)) / six_v,
        );
        // d2 corresponding to face opposite to node 1 (v2) -> grad_N1
        let grad_n1 = Vector3::new(
            (v1.y * (v4.z - v3.z) + v3.y * (v1.z - v4.z) + v4.y * (v3.z - v1.z)) / six_v,
            (v1.x * (v3.z - v4.z) + v3.x * (v4.z - v1.z) + v4.x * (v1.z - v3.z)) / six_v,
            (v1.x * (v4.y - v3.y) + v3.x * (v1.y - v4.y) + v4.x * (v3.y - v1.y)) / six_v,
        );
        // d3 corresponding to face opposite to node 2 (v3) -> grad_N2
        let grad_n2 = Vector3::new(
            (v1.y * (v2.z - v4.z) + v2.y * (v4.z - v1.z) + v4.y * (v1.z - v2.z)) / six_v,
            (v1.x * (v4.z - v2.z) + v2.x * (v1.z - v4.z) + v4.x * (v2.z - v1.z)) / six_v,
            (v1.x * (v2.y - v4.y) + v2.x * (v4.y - v1.y) + v4.x * (v1.y - v2.y)) / six_v,
        );
        
        // Sum of gradients is zero: grad_N3 = -(grad_N0 + grad_N1 + grad_N2)
        let grad_n3 = -(grad_n0 + grad_n1 + grad_n2);
        
        // CRITICAL FIX: Columns must correspond to nodes [0, 1, 2, 3]
        self.shape_derivatives = Matrix3x4::from_columns(&[grad_n0, grad_n1, grad_n2, grad_n3]);
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
                                    * T::from_f64(0.5).unwrap_or_else(T::one);
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
        let half = T::from_f64(0.5).unwrap_or_else(T::zero);
        (velocity_gradient + velocity_gradient.transpose()) * half
    }
}
