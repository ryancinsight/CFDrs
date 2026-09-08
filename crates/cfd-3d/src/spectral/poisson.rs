//! Poisson solver using spectral methods
//!
//! Solves ∇²u = f with various boundary conditions
//!
//! # Theorem — Spectral Poisson Convergence (Boyd 2001)
//!
//! For an analytic right-hand side $f$, the Chebyshev–tau solution $u_N$ of
//! $\nabla^2 u = f$ converges spectrally:
//!
//! ```text
//! ‖u − u_N‖_∞ = O(e^{−cN})
//! ```
//!
//! For $f \in H^s$, the convergence is algebraic: $O(N^{2-s})$ in the $H^1$ norm.
//!
//! Reference: Boyd, J.P. (2001). "Chebyshev and Fourier Spectral Methods"

use super::basis::SpectralBasis;
use super::chebyshev::ChebyshevPolynomial;
use crate::scalar;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use eunomia::RealField;
use leto::{Array1, Array2};
use leto_ops::{MatrixProduct, MatrixSolve, RealScalar};

/// Boundary condition type for Poisson equation
#[derive(Debug, Clone)]
pub enum PoissonBoundaryCondition<T: cfd_mesh::domain::core::Scalar + RealField + RealScalar + Copy>
{
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
pub struct PoissonSolver<T: cfd_mesh::domain::core::Scalar + RealField + RealScalar + Copy> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
    nz: usize,
    /// Chebyshev basis for each direction
    basis_x: ChebyshevPolynomial<T>,
    basis_y: ChebyshevPolynomial<T>,
    basis_z: ChebyshevPolynomial<T>,
}

impl<T> PoissonSolver<T>
where
    T: cfd_mesh::domain::core::Scalar + RealField + RealScalar + Copy,
{
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

    /// Create new Poisson solver with specified basis functions.
    ///
    /// All three directions use the Chebyshev–Gauss–Lobatto basis
    /// regardless of the `_basis` parameter.  The parameter is
    /// retained for API forward-compatibility; Fourier and Legendre
    /// bases require different quadrature and boundary treatment.
    pub fn new_with_basis(
        modes: (usize, usize, usize),
        _basis: (SpectralBasis, SpectralBasis, SpectralBasis),
    ) -> Result<Self> {
        Self::new(modes.0, modes.1, modes.2)
    }

    /// Solve 3D Poisson equation: ∇²u = f
    /// Using tensor product of 1D operators
    pub fn solve(
        &self,
        f: &Array1<T>,
        bc_x: &(PoissonBoundaryCondition<T>, PoissonBoundaryCondition<T>),
        bc_y: &(PoissonBoundaryCondition<T>, PoissonBoundaryCondition<T>),
        bc_z: &(PoissonBoundaryCondition<T>, PoissonBoundaryCondition<T>),
    ) -> Result<Array1<T>> {
        let expected_len = self.nx * self.ny * self.nz;
        if f.size() != expected_len {
            return Err(Error::DimensionMismatch {
                expected: expected_len,
                actual: f.size(),
            });
        }

        // Build the discrete Laplacian operator.
        let mut laplacian = self.build_laplacian_matrix()?;
        let mut rhs = Array1::from_vec([expected_len], f.iter().copied().collect())
            .map_err(|error| Error::InvalidConfiguration(format!("Poisson RHS shape: {error}")))?;

        // Apply boundary conditions to all nodes on the boundary surfaces
        // X-boundaries
        self.apply_boundary_conditions(
            &mut laplacian,
            &mut rhs,
            &bc_x.0,
            &self.get_indices_x(0),
            &self.basis_x,
            0,
            self.ny * self.nz,
        )?;
        self.apply_boundary_conditions(
            &mut laplacian,
            &mut rhs,
            &bc_x.1,
            &self.get_indices_x(self.nx - 1),
            &self.basis_x,
            self.nx - 1,
            self.ny * self.nz,
        )?;

        // Y-boundaries
        self.apply_boundary_conditions(
            &mut laplacian,
            &mut rhs,
            &bc_y.0,
            &self.get_indices_y(0),
            &self.basis_y,
            0,
            self.nz,
        )?;
        self.apply_boundary_conditions(
            &mut laplacian,
            &mut rhs,
            &bc_y.1,
            &self.get_indices_y(self.ny - 1),
            &self.basis_y,
            self.ny - 1,
            self.nz,
        )?;

        // Z-boundaries
        self.apply_boundary_conditions(
            &mut laplacian,
            &mut rhs,
            &bc_z.0,
            &self.get_indices_z(0),
            &self.basis_z,
            0,
            1,
        )?;
        self.apply_boundary_conditions(
            &mut laplacian,
            &mut rhs,
            &bc_z.1,
            &self.get_indices_z(self.nz - 1),
            &self.basis_z,
            self.nz - 1,
            1,
        )?;

        // Solve the linear system using Leto-ops LU decomposition.
        laplacian
            .solve(&rhs.view())
            .map_err(|error| Self::map_leto_error("PoissonSolver::solve", error))
    }

    fn get_indices_x(&self, i: usize) -> Vec<usize> {
        let mut indices = Vec::with_capacity(self.ny * self.nz);
        for j in 0..self.ny {
            for k in 0..self.nz {
                indices.push(i * self.ny * self.nz + j * self.nz + k);
            }
        }
        indices
    }

    fn get_indices_y(&self, j: usize) -> Vec<usize> {
        let mut indices = Vec::with_capacity(self.nx * self.nz);
        for i in 0..self.nx {
            for k in 0..self.nz {
                indices.push(i * self.ny * self.nz + j * self.nz + k);
            }
        }
        indices
    }

    fn get_indices_z(&self, k: usize) -> Vec<usize> {
        let mut indices = Vec::with_capacity(self.nx * self.ny);
        for i in 0..self.nx {
            for j in 0..self.ny {
                indices.push(i * self.ny * self.nz + j * self.nz + k);
            }
        }
        indices
    }

    /// Build the discrete Laplacian matrix using Kronecker products
    fn build_laplacian_matrix(&self) -> Result<Array2<T>> {
        // Get the 1D second-derivative operators from the basis functions
        let d2x = self.basis_x.second_derivative_matrix()?;
        let d2y = self.basis_y.second_derivative_matrix()?;
        let d2z = self.basis_z.second_derivative_matrix()?;

        // Get identity matrices of the appropriate sizes
        let ix = Self::identity_matrix(self.nx);
        let iy = Self::identity_matrix(self.ny);
        let iz = Self::identity_matrix(self.nz);

        // Build the 3D operator by summing the Kronecker products
        // L = D2x ⊗ Iy ⊗ Iz + Ix ⊗ D2y ⊗ Iz + Ix ⊗ Iy ⊗ D2z
        let term_x = d2x
            .kron(&iy)
            .and_then(|term| term.kron(&iz))
            .map_err(|error| Self::map_leto_error("Poisson x-Laplacian kron", error))?;
        let term_y = ix
            .kron(&d2y)
            .and_then(|term| term.kron(&iz))
            .map_err(|error| Self::map_leto_error("Poisson y-Laplacian kron", error))?;
        let term_z = ix
            .kron(&iy)
            .and_then(|term| term.kron(&d2z))
            .map_err(|error| Self::map_leto_error("Poisson z-Laplacian kron", error))?;

        Self::add_three_matrices(&term_x, &term_y, &term_z)
    }

    fn identity_matrix(n: usize) -> Array2<T> {
        Array2::from_shape_fn([n, n], |[row, col]| {
            if row == col {
                scalar::one::<T>()
            } else {
                scalar::zero::<T>()
            }
        })
    }

    fn add_three_matrices(
        left: &Array2<T>,
        middle: &Array2<T>,
        right: &Array2<T>,
    ) -> Result<Array2<T>> {
        let shape = left.shape();
        if middle.shape() != shape {
            return Err(Error::DimensionMismatch {
                expected: left.size(),
                actual: middle.size(),
            });
        }
        if right.shape() != shape {
            return Err(Error::DimensionMismatch {
                expected: left.size(),
                actual: right.size(),
            });
        }

        Array2::from_vec(
            shape,
            left.iter()
                .zip(middle.iter())
                .zip(right.iter())
                .map(|((&x, &y), &z)| x + y + z)
                .collect(),
        )
        .map_err(|error| Error::InvalidConfiguration(format!("Poisson Laplacian sum: {error}")))
    }

    fn map_leto_error(context: &str, error: leto::LetoError) -> Error {
        if matches!(
            &error,
            leto::LetoError::StorageError { reason } if reason.contains("singular")
        ) {
            Error::Numerical(NumericalErrorKind::SingularMatrix)
        } else {
            Error::InvalidConfiguration(format!("{context}: {error}"))
        }
    }

    /// Apply boundary conditions to the system
    fn apply_boundary_conditions(
        &self,
        matrix: &mut Array2<T>,
        rhs: &mut Array1<T>,
        bc: &PoissonBoundaryCondition<T>,
        boundary_indices: &[usize],
        basis: &ChebyshevPolynomial<T>,
        local_idx: usize,
        stride: usize,
    ) -> Result<()> {
        let n = basis.num_points();
        let diff_matrix = basis.diff_matrix();

        // Normal direction factor:
        // Index 0 is x=1, normal is +direction
        // Index n-1 is x=-1, normal is -direction
        let normal_factor = if local_idx == 0 {
            scalar::one::<T>()
        } else {
            -scalar::one::<T>()
        };

        match bc {
            PoissonBoundaryCondition::Dirichlet(value) => {
                for &idx in boundary_indices {
                    // Set row to identity
                    for j in 0..matrix.shape()[1] {
                        matrix[[idx, j]] = if idx == j {
                            scalar::one::<T>()
                        } else {
                            scalar::zero::<T>()
                        };
                    }
                    rhs[idx] = *value;
                }
            }
            PoissonBoundaryCondition::Neumann(value) => {
                for &idx in boundary_indices {
                    // Clear row
                    for j in 0..matrix.shape()[1] {
                        matrix[[idx, j]] = scalar::zero::<T>();
                    }

                    // ∂u/∂n = normal_factor * ∂u/∂x = value
                    // Row idx should be: normal_factor * sum_m D[local_idx, m] * u(m, pencil)
                    for m in 0..n {
                        let pencil_idx = if m >= local_idx {
                            idx + (m - local_idx) * stride
                        } else {
                            idx - (local_idx - m) * stride
                        };
                        matrix[[idx, pencil_idx]] = normal_factor * diff_matrix[[local_idx, m]];
                    }
                    rhs[idx] = *value;
                }
            }
            PoissonBoundaryCondition::Robin { alpha, beta, value } => {
                for &idx in boundary_indices {
                    // Clear row
                    for j in 0..matrix.shape()[1] {
                        matrix[[idx, j]] = scalar::zero::<T>();
                    }

                    // αu + β∂u/∂n = value
                    // Row idx: α * u(local_idx) + β * normal_factor * sum_m D[local_idx, m] * u(m, pencil)
                    for m in 0..n {
                        let pencil_idx = if m >= local_idx {
                            idx + (m - local_idx) * stride
                        } else {
                            idx - (local_idx - m) * stride
                        };

                        let mut coeff = *beta * normal_factor * diff_matrix[[local_idx, m]];
                        if m == local_idx {
                            coeff += *alpha;
                        }
                        matrix[[idx, pencil_idx]] = coeff;
                    }
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
