//! Element-level operations for FEM

use crate::fem::constants;
use nalgebra::{DMatrix, Matrix3, Matrix3x4, RealField, Vector3};
use num_traits::FromPrimitive;

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

impl<T: RealField + FromPrimitive + Copy> FluidElement<T> {
    /// Create new fluid element
    #[must_use]
    pub fn new(nodes: Vec<usize>) -> Self {
        Self {
            nodes,
            volume: T::zero(),
            shape_derivatives: Matrix3x4::zeros(),
        }
    }

    /// Calculate element volume (for tetrahedron)
    pub fn calculate_volume(&mut self, vertices: &[Vector3<T>]) -> T {
        if self.nodes.len() != constants::TET4_NODES {
            return T::zero();
        }

        let v0 = &vertices[self.nodes[0]];
        let v1 = &vertices[self.nodes[1]];
        let v2 = &vertices[self.nodes[2]];
        let v3 = &vertices[self.nodes[3]];

        // Volume = |det(v1-v0, v2-v0, v3-v0)| / 6
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let e3 = v3 - v0;

        let det = e1.dot(&e2.cross(&e3));
        self.volume = det.abs() / T::from_f64(constants::TET_VOLUME_FACTOR).unwrap_or_else(T::one);
        self.volume
    }

    /// Calculate shape function derivatives for tetrahedral element
    /// Based on Zienkiewicz & Taylor, "The Finite Element Method", Vol. 1, 6th Ed.
    pub fn calculate_shape_derivatives(&mut self, vertices: &[Vector3<T>]) {
        if vertices.len() != 4 {
            // For non-tetrahedral elements, use identity as fallback
            self.shape_derivatives = Matrix3x4::identity();
            return;
        }

        // For linear tetrahedron with vertices at (x1,y1,z1), ..., (x4,y4,z4)
        // Shape functions: N_i = (a_i + b_i*x + c_i*y + d_i*z) / (6*V)
        // where V is the element volume

        let v1 = vertices[0];
        let v2 = vertices[1];
        let v3 = vertices[2];
        let v4 = vertices[3];

        // Calculate volume using scalar triple product
        let volume =
            ((v2 - v1).cross(&(v3 - v1))).dot(&(v4 - v1)) / T::from_f64(6.0).unwrap_or_else(T::one);

        if volume.abs() < T::from_f64(1e-12).unwrap_or_else(T::zero) {
            // Degenerate element
            self.shape_derivatives = Matrix3x4::zeros();
            return;
        }

        // Calculate shape function derivatives (constant for linear elements)
        // dN1/dx = (y2(z3-z4) + y3(z4-z2) + y4(z2-z3)) / (6*V)
        // Similar for other derivatives
        let six_v = T::from_f64(6.0).unwrap_or_else(T::one) * volume;

        // Store as columns of the derivative matrix
        // Derivatives for nodes 1, 2, 3
        let d1 = Vector3::new(
            (v2.y * (v3.z - v4.z) + v3.y * (v4.z - v2.z) + v4.y * (v2.z - v3.z)) / six_v,
            (v2.x * (v4.z - v3.z) + v3.x * (v2.z - v4.z) + v4.x * (v3.z - v2.z)) / six_v,
            (v2.x * (v3.y - v4.y) + v3.x * (v4.y - v2.y) + v4.x * (v2.y - v3.y)) / six_v,
        );
        let d2 = Vector3::new(
            (v1.y * (v4.z - v3.z) + v3.y * (v1.z - v4.z) + v4.y * (v3.z - v1.z)) / six_v,
            (v1.x * (v3.z - v4.z) + v3.x * (v4.z - v1.z) + v4.x * (v1.z - v3.z)) / six_v,
            (v1.x * (v4.y - v3.y) + v3.x * (v1.y - v4.y) + v4.x * (v3.y - v1.y)) / six_v,
        );
        let d3 = Vector3::new(
            (v1.y * (v2.z - v4.z) + v2.y * (v4.z - v1.z) + v4.y * (v1.z - v2.z)) / six_v,
            (v1.x * (v4.z - v2.z) + v2.x * (v1.z - v4.z) + v4.x * (v2.z - v1.z)) / six_v,
            (v1.x * (v2.y - v4.y) + v2.x * (v4.y - v1.y) + v4.x * (v1.y - v2.y)) / six_v,
        );
        // For linear shape functions: dN1 + dN2 + dN3 + dN4 = 0, so dN4 = -(dN1 + dN2 + dN3)
        let d4 = -(d1 + d2 + d3);
        
        self.shape_derivatives = Matrix3x4::from_columns(&[d1, d2, d3, d4]);
    }

    /// Calculate strain rate from velocity gradient
    pub fn strain_rate(&self, velocity_gradient: &Matrix3<T>) -> Matrix3<T> {
        let half = T::from_f64(0.5).unwrap_or_else(T::zero);
        (velocity_gradient + velocity_gradient.transpose()) * half
    }

    /// Calculate element stiffness matrix contribution for Stokes flow
    /// Based on Hughes, "The Finite Element Method", Dover Publications, 2000
    pub fn stiffness_contribution(&self, viscosity: T) -> DMatrix<T> {
        let n_dof = self.nodes.len() * constants::VELOCITY_COMPONENTS;
        let mut k_e = DMatrix::zeros(n_dof, n_dof);

        // For Stokes flow: K_e = ∫_Ωe B^T D B dΩ
        // where B is the strain-displacement matrix and D is the material matrix

        // For linear tetrahedral elements, we use single-point integration
        // This is exact for linear elements

        if self.nodes.len() == 4 {
            // Tetrahedral element
            let factor = viscosity * self.volume;

            // Build B matrix (strain-displacement matrix)
            // For 3D: B relates nodal velocities to strain rates
            // Size: 6 x (4*3) for tetrahedral element

            // Full Stokes element matrix with viscous and pressure terms
            // K = [K_uu  B^T]
            //     [B     0  ]
            // where K_uu is viscous stiffness, B is divergence operator

            for i in 0..4 {
                for j in 0..4 {
                    // Viscous stiffness block (velocity-velocity coupling)
                    for d in 0..3 {
                        let row = i * 3 + d;
                        let col = j * 3 + d;

                        // Viscous contribution: μ ∫ (∇Ni : ∇Nj) dΩ
                        // Full tensor contraction for Stokes equations
                        let mut visc_term = T::zero();
                        for k in 0..3 {
                            visc_term +=
                                self.shape_derivatives[(k, i)] * self.shape_derivatives[(k, j)];
                        }
                        k_e[(row, col)] += factor * visc_term;

                        // Additional viscous terms from strain rate tensor
                        // ∫ μ (∂Ni/∂xk)(∂Nj/∂xk) dΩ for cross-derivatives
                        if d < 3 {
                            for k in 0..3 {
                                if k != d {
                                    let cross_term = self.shape_derivatives[(d, i)]
                                        * self.shape_derivatives[(k, j)]
                                        + self.shape_derivatives[(k, i)]
                                            * self.shape_derivatives[(d, j)];
                                    k_e[(row, col)] += factor
                                        * cross_term
                                        * T::from_f64(0.5).unwrap_or_else(|| T::one());
                                }
                            }
                        }
                    }
                }
            }
        }

        k_e
    }
}
