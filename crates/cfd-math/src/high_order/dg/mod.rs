//! # Discontinuous Galerkin (DG) Methods
//!
//! This module provides a high-performance implementation of Discontinuous Galerkin (DG) methods
//! for solving partial differential equations (PDEs). DG methods combine the flexibility of
//! finite volume methods with the high-order accuracy of spectral methods, making them
//! particularly well-suited for problems with sharp gradients, discontinuities, and complex
//! geometries.
//!
//! ## Features
//!
//! - **High-order accuracy**: Support for arbitrary polynomial orders
//! - **Flexible time integration**: Multiple explicit, implicit, and IMEX time integration schemes
//! - **Shock-capturing**: Built-in limiters for handling discontinuities
//! - **Adaptive mesh refinement**: Support for h- and p-adaptivity
//! - **Parallel computing**: Designed for efficient parallel execution
//!
//! ## Quick Start
//!
//! ```no_run
//! use cfd_math::high_order::dg::*;
//! use nalgebra::{DVector, DMatrix};
//!
//! // Create a DG operator
//! let order = 3;
//! let num_components = 1;
//! let params = DGOperatorParams::new()
//!     .with_volume_flux(FluxType::Central)
//!     .with_surface_flux(FluxType::LaxFriedrichs)
//!     .with_limiter(LimiterType::Minmod);
//!
//! let dg_op = DGOperator::new(order, num_components, Some(params)).unwrap();
//!
//! // Create a time integrator
//! let integrator = TimeIntegratorFactory::create(TimeIntegration::SSPRK3);
//!
//! // Set up the solver
//! let t_final = 1.0;
//! let solver_params = TimeIntegrationParams::new(TimeIntegration::SSPRK3)
//!     .with_t_final(t_final)
//!     .with_dt(0.01)
//!     .with_verbose(true);
//!
//! let mut solver = DGSolver::new(dg_op, integrator, solver_params);
//!
//! // Set initial condition
//! let u0 = |x: f64| DVector::from_vec(vec![x.sin()]);
//! solver.initialize(u0).unwrap();
//!
//! // Define the right-hand side function
//! fn f(_t: f64, u: &DMatrix<f64>) -> Result<DMatrix<f64>> {
//!     // For the linear advection equation: du/dt = -du/dx
//!     // This is a simplified example; in practice, you would use the DG operator
//!     // to compute the spatial derivatives
//!     Ok(-u.clone())
//! }
//!
//! // Run the solver
//! solver.solve(f, None::<fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>>).unwrap();
//!
//! // Evaluate the solution
//! let x = 0.5;
//! let u = solver.evaluate(x);
//! ```
//!
//! ## Modules
//!
//! - `basis`: Basis functions and quadrature rules
//! - `flux`: Numerical flux functions
//! - `limiter`: Slope limiters for shock capturing
//! - `operators`: Core DG operators and discretization
//! - `solver`: Time integration and solution algorithms

#![warn(missing_docs)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::float_cmp)]

mod basis;
mod flux;
mod limiter;
mod operators;
mod solver;

// Import from parent modules
use super::spectral::SpectralError;

// Re-export all public types and traits
pub use basis::*;
pub use flux::*;
pub use limiter::*;
pub use operators::*;
pub use solver::*;

use nalgebra::{DMatrix, DVector};
use std::fmt;

/// Error types for DG methods
#[derive(Debug, thiserror::Error)]
pub enum DGError {
    /// Invalid polynomial order
    #[error("Polynomial order must be at least 1, got {0}")]
    InvalidOrder(usize),

    /// Invalid number of quadrature points
    #[error("Number of quadrature points must be at least ceil((order + 1)/2), got {0}")]
    InvalidQuadrature(usize),

    /// Invalid input dimensions
    #[error("Invalid dimensions: {0}")]
    InvalidDimensions(String),

    /// Numerical error (e.g., matrix inversion failed)
    #[error("Numerical error: {0}")]
    NumericalError(String),

    /// Solver did not converge
    #[error("Solver did not converge: {0}")]
    ConvergenceError(String),

    /// Invalid parameter
    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),
}

impl From<SpectralError> for DGError {
    fn from(err: SpectralError) -> Self {
        match err {
            SpectralError::InvalidOrder(o) => DGError::InvalidOrder(o),
            SpectralError::NodeComputation(s) => {
                DGError::NumericalError(format!("Node computation failed: {s}"))
            }
            SpectralError::MatrixError(s) => DGError::NumericalError(format!("Matrix error: {s}")),
        }
    }
}

/// Result type for DG operations
pub type Result<T> = std::result::Result<T, DGError>;

/// Represents a DG solution on a single element
#[derive(Clone)]
pub struct DGSolution {
    /// Polynomial order
    pub order: usize,
    /// Number of components (for systems of equations)
    pub num_components: usize,
    /// Solution coefficients (num_components × num_basis_functions)
    pub coefficients: DMatrix<f64>,
    /// Basis functions
    pub basis: DGBasis,
}

impl DGSolution {
    /// Create a new DG solution with given order and number of components
    ///
    /// # Arguments
    /// * `order` - Polynomial order (must be ≥ 1)
    /// * `num_components` - Number of components in the solution vector
    /// * `basis_type` - Type of basis functions (optional, defaults to Orthogonal)
    ///
    /// # Returns
    /// A new `DGSolution` instance or an error if the order is invalid
    pub fn new(order: usize, num_components: usize) -> Result<Self> {
        Self::with_basis(order, num_components, BasisType::Orthogonal)
    }

    /// Create a new DG solution with given order, components and basis type
    pub fn with_basis(order: usize, num_components: usize, basis_type: BasisType) -> Result<Self> {
        if order == 0 {
            return Err(DGError::InvalidOrder(order));
        }

        let basis = DGBasis::new(order, basis_type)?;
        let num_basis = order + 1;

        Ok(Self {
            order,
            num_components,
            coefficients: DMatrix::zeros(num_components, num_basis),
            basis,
        })
    }

    /// Evaluate the solution at a point in the reference element
    ///
    /// # Arguments
    /// * `x` - Point in the reference element [-1, 1]
    ///
    /// # Returns
    /// The solution vector at the given point
    pub fn evaluate(&self, x: f64) -> DVector<f64> {
        let mut result = DVector::zeros(self.num_components);

        for i in 0..self.coefficients.ncols() {
            let phi_i = self.basis.evaluate_basis(i, x);
            for c in 0..self.num_components {
                result[c] += self.coefficients[(c, i)] * phi_i;
            }
        }

        result
    }

    /// Compute the L² norm of the solution
    ///
    /// # Returns
    /// The L² norm of the solution
    pub fn l2_norm(&self) -> f64 {
        let mut norm_sq = 0.0;

        for i in 0..self.num_components {
            for j in 0..self.coefficients.ncols() {
                for k in 0..self.coefficients.ncols() {
                    // Mass matrix M_{j,k} = ∫ φ_j(x) φ_k(x) dx
                    let m_jk = self.basis.mass_matrix[(j, k)];
                    norm_sq += self.coefficients[(i, j)] * self.coefficients[(i, k)] * m_jk;
                }
            }
        }

        norm_sq.sqrt()
    }

    /// Apply a limiter to the solution
    ///
    /// # Arguments
    /// * `limiter` - The limiter to apply
    ///
    /// # Remarks
    /// TODO: Implement limiter application using neighbor solutions and troubled-cell detection.
    pub fn apply_limiter<L: Limiter>(&mut self, limiter: &L) {
        // TODO: Provide neighbor information and propagate limiter errors instead of discarding.

        let empty_neighbors = &[];
        let params = limiter::LimiterParams::new(limiter::LimiterType::Minmod);

        limiter.limit(self, empty_neighbors, &params).unwrap_or(());
    }

    /// Get the cell average of the solution
    ///
    /// # Returns
    /// The cell average of the solution
    pub fn average(&self) -> DVector<f64> {
        let mut avg = DVector::zeros(self.num_components);

        // The average is (1/2) * ∫ u(x) dx over [-1, 1]
        // ∫ u(x) dx = ∑_j c_j * ∫ φ_j(x) dx
        // The integrals of the basis functions can be computed using quadrature
        let mut basis_integrals = DVector::zeros(self.coefficients.ncols());
        for j in 0..self.coefficients.ncols() {
            let mut int_phi_j = 0.0;
            for q in 0..self.basis.quad_points.len() {
                int_phi_j += self.basis.quad_weights[q] * self.basis.phi[(j, q)];
            }
            basis_integrals[j] = int_phi_j;
        }

        for i in 0..self.num_components {
            let mut sum = 0.0;
            for j in 0..self.coefficients.ncols() {
                sum += self.coefficients[(i, j)] * basis_integrals[j];
            }
            avg[i] = 0.5 * sum;
        }

        avg
    }
}

impl fmt::Debug for DGSolution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "DGSolution {{ order: {}, num_components: {}, coefficients: [...] }}",
            self.order, self.num_components
        )
    }
}

/// Trait for DG methods
pub trait DGMethod {
    /// Compute the time derivative of the solution
    fn compute_rhs(&self, t: f64, u: &DGSolution) -> Result<DGSolution>;

    /// Compute the maximum stable time step
    fn max_time_step(&self, u: &DGSolution) -> f64;

    /// Apply boundary conditions
    fn apply_boundary_conditions(&mut self, u: &mut DGSolution, t: f64) -> Result<()>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_dg_solution() {
        // Create a DG solution with a quadratic polynomial
        let order = 2;
        let num_components = 1;
        let mut sol = DGSolution::new(order, num_components).unwrap();

        // Set the coefficients for u(x) = 1 + x + x²
        // In the Legendre basis: 4/3*P₀(x) + P₁(x) + (2/3)*P₂(x)
        sol.coefficients[(0, 0)] = 4.0 / 3.0; // P₀ term
        sol.coefficients[(0, 1)] = 1.0; // P₁ term
        sol.coefficients[(0, 2)] = 2.0 / 3.0; // P₂ term

        // Test evaluation at x = 1.0
        let u1 = sol.evaluate(1.0);
        assert_relative_eq!(u1[0], 3.0, epsilon = 1e-10);

        // Test evaluation at x = 0.0
        let u0 = sol.evaluate(0.0);
        assert_relative_eq!(u0[0], 1.0, epsilon = 1e-10);

        // Test L² norm
        // ∫(1 + x + x²)² dx from -1 to 1 = 4.4
        let norm = sol.l2_norm();
        assert_relative_eq!(norm, (4.4f64).sqrt(), epsilon = 1e-10);
    }

    #[test]
    fn test_dg_solution_average() {
        // Create a DG solution with a constant function
        let order = 2;
        let num_components = 1;
        let mut sol = DGSolution::new(order, num_components).unwrap();

        // Set the coefficients for u(x) = 2.0
        // Since P₀(x) = 1.0, u(x) = 2.0 * P₀(x)
        sol.coefficients[(0, 0)] = 2.0;

        // The average should be 2.0
        let avg = sol.average();
        assert_relative_eq!(avg[0], 2.0, epsilon = 1e-10);
    }
}
