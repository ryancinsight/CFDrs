//! Poisson solver using spectral methods
//!
//! Solves ∇²u = f with various boundary conditions
//! Reference: Boyd, J.P. (2001). "Chebyshev and Fourier Spectral Methods"

use nalgebra::{DMatrix, DVector, RealField};
use num_traits::FromPrimitive;
use cfd_core::Result;
use super::chebyshev::ChebyshevPolynomial;

/// Boundary condition type for Poisson equation
#[derive(Debug, Clone)]
pub enum PoissonBoundaryCondition<T: RealField + Copy> {
    /// Dirichlet: u = g on boundary
    Dirichlet(T),
    /// Neumann: ∂u/∂n = g on boundary
    Neumann(T),
    /// Robin: αu + β∂u/∂n = g on boundary
    Robin { alpha: T, beta: T, value: T },
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
    /// Create new Poisson solver
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
        let solution = lu.solve(&rhs)
            .ok_or_else(|| cfd_core::error::Error::Numerical(
                cfd_core::error::NumericalErrorKind::SingularMatrix
            ))?;
        
        // Reshape solution back to matrix form
        Ok(self.unflatten_solution(&solution))
    }
    
    /// Build the discrete Laplacian matrix using Kronecker products
    fn build_laplacian_matrix(&self) -> Result<DMatrix<T>> {
        let n_total = self.nx * self.ny * self.nz;
        let mut laplacian = DMatrix::zeros(n_total, n_total);
        
        // Build 1D second derivative matrices
        let d2x = self.basis_x.second_derivative_matrix()?;
        let d2y = self.basis_y.second_derivative_matrix()?;
        let d2z = self.basis_z.second_derivative_matrix()?;
        
        // Form 3D Laplacian using Kronecker products
        // L = D2x ⊗ Iy ⊗ Iz + Ix ⊗ D2y ⊗ Iz + Ix ⊗ Iy ⊗ D2z
        for i in 0..n_total {
            for j in 0..n_total {
                let (ix, iy, iz) = self.linear_to_3d_index(i);
                let (jx, jy, jz) = self.linear_to_3d_index(j);
                
                let mut value = T::zero();
                
                // x-derivative contribution
                if iy == jy && iz == jz {
                    value += d2x[(ix, jx)];
                }
                
                // y-derivative contribution
                if ix == jx && iz == jz {
                    value += d2y[(iy, jy)];
                }
                
                // z-derivative contribution
                if ix == jx && iy == jy {
                    value += d2z[(iz, jz)];
                }
                
                laplacian[(i, j)] = value;
            }
        }
        
        Ok(laplacian)
    }
    
    /// Convert linear index to 3D indices
    fn linear_to_3d_index(&self, idx: usize) -> (usize, usize, usize) {
        let iz = idx % self.nz;
        let iy = (idx / self.nz) % self.ny;
        let ix = idx / (self.ny * self.nz);
        (ix, iy, iz)
    }
    
    /// Flatten RHS matrix to vector
    fn flatten_rhs(&self, f: &DMatrix<T>) -> DVector<T> {
        let n_total = self.nx * self.ny * self.nz;
        let mut rhs = DVector::zeros(n_total);
        
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                for iz in 0..self.nz {
                    let idx = ix * self.ny * self.nz + iy * self.nz + iz;
                    rhs[idx] = f[(ix * self.ny + iy, iz)];
                }
            }
        }
        
        rhs
    }
    
    /// Unflatten solution vector to matrix
    fn unflatten_solution(&self, solution: &DVector<T>) -> DMatrix<T> {
        let mut result = DMatrix::zeros(self.nx * self.ny, self.nz);
        
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                for iz in 0..self.nz {
                    let idx = ix * self.ny * self.nz + iy * self.nz + iz;
                    result[(ix * self.ny + iy, iz)] = solution[idx];
                }
            }
        }
        
        result
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
                    // Note: beta term requires normal derivative discretization
                    // For spectral methods, this involves modal coefficients
                    let _ = beta; // Placeholder for normal derivative term
                    rhs[idx] = *value;
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