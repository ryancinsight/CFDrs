//! # Spectral Method Solver for 3D PDEs
//!
//! This module implements high-accuracy spectral methods using Fourier and
//! Chebyshev expansions for solving partial differential equations.
//!
//! ## Mathematical Foundation
//!
//! ### Spectral Expansion
//! Solutions are represented as truncated series:
//!
//! ```math
//! u(x,y,z) ≈ ∑_{k=0}^N ∑_{m=0}^M ∑_{n=0}^P û_{kmn} φ_k(x) ψ_m(y) χ_n(z)
//! ```
//!
//! ### Fourier Spectral Method
//! - **Basis Functions**: φ_k(x) = e^{i k x} (periodic domains)
//! - **Transform**: FFT provides O(N log N) evaluation
//! - **Accuracy**: Exponential convergence for smooth, periodic functions
//!
//! ### Chebyshev Spectral Method
//! - **Basis Functions**: φ_k(x) = T_k(x) Chebyshev polynomials
//! - **Collocation**: Solutions evaluated at Gauss-Lobatto points
//! - **Accuracy**: Exponential convergence for smooth functions with boundaries
//!
//! ## Convergence Theory
//!
//! **Theorem (Spectral Convergence)**: For analytic solutions,
//! spectral methods achieve exponential convergence:
//!
//! ```math
//! ||u - u_N|| ≤ C e^{-c N}
//! ```
//!
//! where N is the number of modes and c depends on solution regularity.
//!
//! **Theorem (Aliasing Error)**: For nonlinear terms, dealiasing is required
//! to maintain spectral accuracy. The 3/2-rule ensures aliasing-free computation.
//!
//! ## Algorithm Implementation
//!
//! ### Fourier Transform Acceleration
//! - **Cooley-Tukey FFT**: O(N log N) complexity
//! - **Bit-reversal Permutation**: Optimal memory access pattern
//! - **Butterfly Operations**: Danielson-Lanczos lemma
//!
//! ### Boundary Conditions
//! - **Fourier**: Periodic boundaries only
//! - **Chebyshev**: Non-periodic boundaries via basis recombination
//! - **Mixed**: Fourier in periodic directions, Chebyshev in bounded directions
//!
//! ### Nonlinear Terms
//! - **Convolution Sum**: Direct evaluation (expensive)
//! - **Transform Method**: FFT-based convolution (efficient)
//! - **Dealiasing**: 3/2-rule to prevent aliasing errors
//!
//! ## Applications
//!
//! - **Direct Numerical Simulation (DNS)**: Turbulent flow simulation
//! - **Boundary Layer Flows**: High accuracy near walls
//! - **Wave Propagation**: Exact dispersion relation preservation
//! - **Stability Analysis**: Linear eigenvalue problems
//!
//! ## Performance Characteristics
//!
//! - **Memory**: O(N³) for 3D transforms (N modes per direction)
//! - **Time**: O(N³ log N) for nonlinear terms with FFT
//! - **Accuracy**: Machine precision for smooth solutions
//! - **Scalability**: Excellent parallel efficiency
//!
//! ## References
//!
//! - Boyd, J.P. (2001). *Chebyshev and Fourier Spectral Methods*
//! - Canuto, C. et al. (2006). *Spectral Methods: Fundamentals in Single Domains*
//! - Trefethen, L.N. (2000). *Spectral Methods in MATLAB*
//! - Cooley, J.W. & Tukey, J.W. (1965). "An algorithm for the machine calculation of complex Fourier series"

use super::basis::SpectralBasis;
use super::poisson::{PoissonBoundaryCondition, PoissonSolver};
use cfd_core::error::Result;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Problem definition for Poisson equation
#[derive(Debug, Clone)]
pub struct PoissonProblem<T: RealField + Copy> {
    /// Source term field (right-hand side of Poisson equation)
    /// Flattened 3D field of size nx * ny * nz
    pub source_term: DVector<T>,
    /// Boundary conditions in x-direction (min, max)
    pub bc_x: (PoissonBoundaryCondition<T>, PoissonBoundaryCondition<T>),
    /// Boundary conditions in y-direction (min, max)
    pub bc_y: (PoissonBoundaryCondition<T>, PoissonBoundaryCondition<T>),
    /// Boundary conditions in z-direction (min, max)
    pub bc_z: (PoissonBoundaryCondition<T>, PoissonBoundaryCondition<T>),
}

/// Spectral method configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectralConfig<T: RealField + Copy> {
    /// Base solver configuration
    pub base: cfd_core::solver::SolverConfig<T>,
    /// Number of modes in x direction
    pub nx_modes: usize,
    /// Number of modes in y direction
    pub ny_modes: usize,
    /// Number of modes in z direction
    pub nz_modes: usize,
    /// Basis type for each direction
    pub basis_x: SpectralBasis,
    /// Basis function type for y-direction (Fourier, Chebyshev, etc.)
    pub basis_y: SpectralBasis,
    /// Basis function type for z-direction (Fourier, Chebyshev, etc.)
    pub basis_z: SpectralBasis,
    /// Time step for time-dependent problems
    pub dt: Option<T>,
}

impl<T: RealField + FromPrimitive + Copy> SpectralConfig<T> {
    /// Create new configuration with validation
    pub fn new(nx: usize, ny: usize, nz: usize) -> Result<Self> {
        Ok(Self {
            base: cfd_core::solver::SolverConfig::builder()
                .tolerance(T::from_f64(1e-10).ok_or_else(|| {
                    cfd_core::error::Error::InvalidConfiguration("Cannot convert tolerance".into())
                })?)
                .max_iterations(100)
                .build(),
            nx_modes: nx,
            ny_modes: ny,
            nz_modes: nz,
            basis_x: SpectralBasis::Chebyshev,
            basis_y: SpectralBasis::Chebyshev,
            basis_z: SpectralBasis::Chebyshev,
            dt: None,
        })
    }
}

/// Spectral solver for 3D problems
pub struct SpectralSolver<T: RealField + Copy> {
    config: SpectralConfig<T>,
    poisson_solver: PoissonSolver<T>,
}

impl<T: RealField + FromPrimitive + Copy> SpectralSolver<T> {
    /// Create new spectral solver
    pub fn new(config: SpectralConfig<T>) -> Result<Self> {
        // Pass the full configuration details to the sub-solver
        let poisson_solver = PoissonSolver::new_with_basis(
            (config.nx_modes, config.ny_modes, config.nz_modes),
            (config.basis_x, config.basis_y, config.basis_z),
        )?;

        Ok(Self {
            config,
            poisson_solver,
        })
    }

    /// Solve a problem using spectral methods
    pub fn solve(&mut self, problem: &PoissonProblem<T>) -> Result<SpectralSolution<T>> {
        // Initialize solution with proper dimensions
        let mut solution = SpectralSolution::new(
            self.config.nx_modes,
            self.config.ny_modes,
            self.config.nz_modes,
        );

        // Use the source and BCs from the user-provided problem definition (zero-copy: pass by reference)
        let potential = self.poisson_solver.solve(
            &problem.source_term,
            &problem.bc_x,
            &problem.bc_y,
            &problem.bc_z,
        )?;

        // Store the result in the solution
        solution.u = potential;

        Ok(solution)
    }
}

/// Solution from spectral solver
pub struct SpectralSolution<T: RealField + Copy> {
    /// Solution field (flattened 3D)
    pub u: DVector<T>,
    /// Grid dimensions
    pub nx: usize,
    /// Number of grid points in y-direction
    pub ny: usize,
    /// Number of grid points in z-direction  
    pub nz: usize,
}

impl<T: RealField + Copy> SpectralSolution<T> {
    /// Create new solution
    #[must_use]
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            u: DVector::zeros(nx * ny * nz),
            nx,
            ny,
            nz,
        }
    }

    /// Get solution at a point
    #[must_use]
    pub fn at(&self, i: usize, j: usize, k: usize) -> T
    where
        T: Clone,
    {
        let idx = i * self.ny * self.nz + j * self.nz + k;
        self.u[idx].clone()
    }
}
