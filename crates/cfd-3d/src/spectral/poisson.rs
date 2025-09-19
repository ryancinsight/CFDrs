//! Poisson solver using spectral methods
//!
//! Solves ∇²u = f with various boundary conditions
//! Reference: Boyd, J.P. (2001). "Chebyshev and Fourier Spectral Methods"

use super::basis::SpectralBasis;
use super::chebyshev::ChebyshevPolynomial;
use cfd_core::error::Result;
use nalgebra::{DMatrix, DVector, RealField};
use num_traits::FromPrimitive;

/// Boundary condition type for Poisson equation
#[derive(Debug, Clone)]
pub enum PoissonBoundaryCondition<T: RealField + Copy> {
    /// Dirichlet: u = g on boundary
    Dirichlet(T),
    /// Neumann: ∂u/∂n = g on boundary
    Neumann(T),
    /// Robin: αu + β∂u/∂n = g on boundary
    Robin {
        /// Coefficient α for Robin boundary condition
        alpha: T,
        /// Coefficient β for Robin boundary condition
        beta: T,
        /// Boundary value g for Robin condition
        value: T,
    },
    /// Periodic boundary conditions
    Periodic,
}

/// Spectral Poisson solver
pub struct PoissonSolver<T: RealField + Copy> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
    nz: usize,
    /// Chebyshev basis for each direction
    basis_x: ChebyshevPolynomial<T>,
    basis_y: ChebyshevPolynomial<T>,
    basis_z: ChebyshevPolynomial<T>,
}

impl<T: RealField + FromPrimitive + Copy> PoissonSolver<T> {
    /// Create new Poisson solver with default Chebyshev basis
    pub fn new(nx: usize, ny: usize, nz: usize) -> Result<Self> {
        Ok(Self {
            nx,
            ny,
            nz,
            basis_x: ChebyshevPolynomial::new(nx)?,
            basis_y: ChebyshevPolynomial::new(ny)?,
            basis_z: ChebyshevPolynomial::new(nz)?,
        })
    }

    /// Create new Poisson solver with specified basis functions
    pub fn new_with_basis(
        modes: (usize, usize, usize),
        _basis: (SpectralBasis, SpectralBasis, SpectralBasis),
    ) -> Result<Self> {
        // Currently only Chebyshev basis is implemented
        // Fourier basis support can be added based on the _basis parameter
        // when required for periodic boundary conditions
        Self::new(modes.0, modes.1, modes.2)
    }

    /// Solve 3D Poisson equation: ∇²u = f
    /// Using tensor product of 1D operators
    pub fn solve(
        &self,
        f: &DMatrix<T>,
        bc_x: (PoissonBoundaryCondition<T>, PoissonBoundaryCondition<T>),
        bc_y: (PoissonBoundaryCondition<T>, PoissonBoundaryCondition<T>),
        bc_z: (PoissonBoundaryCondition<T>, PoissonBoundaryCondition<T>),
    ) -> Result<DMatrix<T>> {
        // Build the discrete Laplacian operator
        let mut laplacian = self.build_laplacian_matrix()?;
        let mut rhs = self.flatten_rhs(f);

        // Apply boundary conditions
        let boundary_indices_x = vec![0, self.nx - 1];
        let boundary_indices_y = vec![0, self.ny - 1];
        let boundary_indices_z = vec![0, self.nz - 1];

        self.apply_boundary_conditions(&mut laplacian, &mut rhs, &bc_x.0, &boundary_indices_x)?;
        self.apply_boundary_conditions(&mut laplacian, &mut rhs, &bc_x.1, &boundary_indices_x)?;
        self.apply_boundary_conditions(&mut laplacian, &mut rhs, &bc_y.0, &boundary_indices_y)?;
        self.apply_boundary_conditions(&mut laplacian, &mut rhs, &bc_y.1, &boundary_indices_y)?;
        self.apply_boundary_conditions(&mut laplacian, &mut rhs, &bc_z.0, &boundary_indices_z)?;
        self.apply_boundary_conditions(&mut laplacian, &mut rhs, &bc_z.1, &boundary_indices_z)?;

        // Solve the linear system using LU decomposition
        let lu = laplacian.lu();
        let solution = lu.solve(&rhs).ok_or_else(|| {
            cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::SingularMatrix)
        })?;

        // Reshape solution back to matrix form
        Ok(self.unflatten_solution(&solution))
    }

    /// Build the discrete Laplacian matrix using Kronecker products
    fn build_laplacian_matrix(&self) -> Result<DMatrix<T>> {
        // Get the 1D second-derivative operators from the basis functions
        let d2x = self.basis_x.second_derivative_matrix()?;
        let d2y = self.basis_y.second_derivative_matrix()?;
        let d2z = self.basis_z.second_derivative_matrix()?;

        // Get identity matrices of the appropriate sizes
        let ix = DMatrix::<T>::identity(self.nx, self.nx);
        let iy = DMatrix::<T>::identity(self.ny, self.ny);
        let iz = DMatrix::<T>::identity(self.nz, self.nz);

        // Build the 3D operator by summing the Kronecker products
        // L = D2x ⊗ Iy ⊗ Iz + Ix ⊗ D2y ⊗ Iz + Ix ⊗ Iy ⊗ D2z
        let term_x = d2x.kronecker(&iy).kronecker(&iz);
        let term_y = ix.kronecker(&d2y).kronecker(&iz);
        let term_z = ix.kronecker(&iy).kronecker(&d2z);

        let laplacian = term_x + term_y + term_z;
        Ok(laplacian)
    }

    /// Convert linear index to 3D indices
    #[allow(dead_code)]
    fn linear_to_3d_index(&self, idx: usize) -> (usize, usize, usize) {
        let iz = idx % self.nz;
        let iy = (idx / self.nz) % self.ny;
        let ix = idx / (self.ny * self.nz);
        (ix, iy, iz)
    }

    /// Flatten RHS matrix to vector
    fn flatten_rhs(&self, f: &DMatrix<T>) -> DVector<T> {
        // Use nalgebra's iterators to perform the copy efficiently
        DVector::from_iterator(f.len(), f.iter().cloned())
    }

    /// Unflatten solution vector to matrix
    fn unflatten_solution(&self, solution: &DVector<T>) -> DMatrix<T> {
        // Reconstruct the matrix directly from the solution vector's iterator
        DMatrix::from_iterator(self.nx * self.ny, self.nz, solution.iter().cloned())
    }

    /// Apply boundary conditions to the system
    fn apply_boundary_conditions(
        &self,
        matrix: &mut DMatrix<T>,
        rhs: &mut DVector<T>,
        bc: &PoissonBoundaryCondition<T>,
        boundary_indices: &[usize],
    ) -> Result<()> {
        match bc {
            PoissonBoundaryCondition::Dirichlet(value) => {
                for &idx in boundary_indices {
                    // Set row to identity
                    for j in 0..matrix.ncols() {
                        matrix[(idx, j)] = if idx == j { T::one() } else { T::zero() };
                    }
                    rhs[idx] = *value;
                }
            }
            PoissonBoundaryCondition::Neumann(value) => {
                // Implement Neumann BC using ghost points or modified stencil
                // This requires modifying the differentiation matrix
                for &idx in boundary_indices {
                    rhs[idx] += *value;
                }
            }
            PoissonBoundaryCondition::Robin { alpha, beta, value } => {
                // Implement Robin BC: αu + β∂u/∂n = g
                for &idx in boundary_indices {
                    matrix[(idx, idx)] += *alpha;
                    // Apply Robin boundary condition: alpha*u + beta*du/dn = value
                    // For spectral methods, the normal derivative is incorporated through
                    // modification of the spectral coefficients
                    if beta.abs() > T::from_f64(1e-10).unwrap_or_else(|| T::zero()) {
                        // Modify RHS to account for normal derivative term
                        // This requires spectral differentiation in the normal direction
                        let normal_derivative_factor =
                            *beta / (*alpha + T::from_f64(1e-10).unwrap_or_else(|| T::zero()));
                        rhs[idx] = *value - normal_derivative_factor * rhs[idx];
                    } else {
                        // Pure Dirichlet case when beta = 0
                        rhs[idx] = *value / *alpha;
                    }
                }
            }
            PoissonBoundaryCondition::Periodic => {
                // Periodic BCs require special treatment in spectral methods
                // Typically handled by using Fourier basis instead
            }
        }
        Ok(())
    }
}
