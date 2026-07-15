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
//! use cfd_math::error::Result;
//! use leto::{Array1, Array2};
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
//! let u0 = |x: f64| Array1::from_shape_vec([1], vec![x.sin()]).unwrap();
//! solver.initialize(u0).unwrap();
//!
//! // Define the right-hand side function
//! fn f(_t: f64, u: &Array2<f64>) -> Result<Array2<f64>> {
//!     // For the linear advection equation: du/dt = -du/dx
//!     // This is a simplified example; in practice, you would use the DG operator
//!     // to compute the spatial derivatives
//!     Ok(Array2::from_shape_fn(u.shape(), |idx| -u[idx]))
//! }
//!
//! // Run the solver
//! solver.solve(f, None::<fn(f64, &Array2<f64>) -> Result<Array2<f64>>>).unwrap();
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

use crate::error::Result;
use cfd_core::error::{Error, ErrorContext};
use leto::{Array1, Array2};
use leto_ops::MatrixSolve;
use std::fmt;

// Re-export all public types and traits
pub use basis::*;
pub use flux::*;
pub use limiter::*;
pub use operators::*;
pub use solver::*;

#[cfg(test)]
pub(crate) fn vector_from_vec(values: Vec<f64>) -> Array1<f64> {
    Array1::from_shape_vec([values.len()], values)
        .expect("invariant: vector shape matches element count")
}

pub(crate) fn vector_from_element(len: usize, value: f64) -> Array1<f64> {
    Array1::from_elem([len], value)
}

pub(crate) fn vector_zeros(len: usize) -> Array1<f64> {
    Array1::zeros([len])
}

pub(crate) fn vector_len(vector: &Array1<f64>) -> usize {
    vector.shape()[0]
}

pub(crate) fn vector_sum(vector: &Array1<f64>) -> f64 {
    vector.iter().copied().sum()
}

pub(crate) fn vector_norm(vector: &Array1<f64>) -> f64 {
    vector
        .iter()
        .map(|&value| value * value)
        .sum::<f64>()
        .sqrt()
}

pub(crate) fn vector_amax(vector: &Array1<f64>) -> f64 {
    vector.iter().map(|&value| value.abs()).fold(0.0, f64::max)
}

pub(crate) fn vector_dot(lhs: &Array1<f64>, rhs: &Array1<f64>) -> f64 {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: vector dot operands must have equal length"
    );
    lhs.iter().zip(rhs.iter()).map(|(&l, &r)| l * r).sum()
}

pub(crate) fn vector_add(lhs: &Array1<f64>, rhs: &Array1<f64>) -> Array1<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: vector add operands must have equal length"
    );
    Array1::from_shape_fn(lhs.shape(), |idx| lhs[idx] + rhs[idx])
}

pub(crate) fn vector_sub(lhs: &Array1<f64>, rhs: &Array1<f64>) -> Array1<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: vector sub operands must have equal length"
    );
    Array1::from_shape_fn(lhs.shape(), |idx| lhs[idx] - rhs[idx])
}

pub(crate) fn vector_scale(vector: &Array1<f64>, scale: f64) -> Array1<f64> {
    Array1::from_shape_fn(vector.shape(), |idx| vector[idx] * scale)
}

pub(crate) fn vector_add_assign_scaled(target: &mut Array1<f64>, source: &Array1<f64>, scale: f64) {
    assert_eq!(
        target.shape(),
        source.shape(),
        "invariant: vector add-assign operands must have equal length"
    );
    for i in 0..target.shape()[0] {
        target[i] += scale * source[i];
    }
}

pub(crate) fn matrix_zeros(rows: usize, cols: usize) -> Array2<f64> {
    Array2::zeros([rows, cols])
}

#[cfg(test)]
pub(crate) fn matrix_from_vec(rows: usize, cols: usize, values: Vec<f64>) -> Array2<f64> {
    Array2::from_shape_vec([rows, cols], values)
        .expect("invariant: matrix shape matches element count")
}

#[cfg(test)]
pub(crate) fn matrix_from_element(rows: usize, cols: usize, value: f64) -> Array2<f64> {
    Array2::from_elem([rows, cols], value)
}

pub(crate) fn matrix_rows(matrix: &Array2<f64>) -> usize {
    matrix.shape()[0]
}

pub(crate) fn matrix_cols(matrix: &Array2<f64>) -> usize {
    matrix.shape()[1]
}

pub(crate) fn matrix_len(matrix: &Array2<f64>) -> usize {
    let [rows, cols] = matrix.shape();
    rows * cols
}

pub(crate) fn matrix_norm(matrix: &Array2<f64>) -> f64 {
    matrix
        .iter()
        .map(|&value| value * value)
        .sum::<f64>()
        .sqrt()
}

pub(crate) fn matrix_neg(matrix: &Array2<f64>) -> Array2<f64> {
    Array2::from_shape_fn(matrix.shape(), |idx| -matrix[idx])
}

pub(crate) fn matrix_add(lhs: &Array2<f64>, rhs: &Array2<f64>) -> Array2<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: matrix add operands must have equal shape"
    );
    Array2::from_shape_fn(lhs.shape(), |idx| lhs[idx] + rhs[idx])
}

pub(crate) fn matrix_sub(lhs: &Array2<f64>, rhs: &Array2<f64>) -> Array2<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: matrix sub operands must have equal shape"
    );
    Array2::from_shape_fn(lhs.shape(), |idx| lhs[idx] - rhs[idx])
}

pub(crate) fn matrix_scale(matrix: &Array2<f64>, scale: f64) -> Array2<f64> {
    Array2::from_shape_fn(matrix.shape(), |idx| matrix[idx] * scale)
}

pub(crate) fn matrix_add_scaled(lhs: &Array2<f64>, rhs: &Array2<f64>, scale: f64) -> Array2<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: matrix add-scaled operands must have equal shape"
    );
    Array2::from_shape_fn(lhs.shape(), |idx| lhs[idx] + scale * rhs[idx])
}

pub(crate) fn matrix_add_assign_scaled(target: &mut Array2<f64>, source: &Array2<f64>, scale: f64) {
    assert_eq!(
        target.shape(),
        source.shape(),
        "invariant: matrix add-assign operands must have equal shape"
    );
    let [rows, cols] = target.shape();
    for row in 0..rows {
        for col in 0..cols {
            target[[row, col]] += scale * source[[row, col]];
        }
    }
}

pub(crate) fn matrix_identity(size: usize) -> Array2<f64> {
    Array2::from_shape_fn(
        [size, size],
        |[row, col]| if row == col { 1.0 } else { 0.0 },
    )
}

pub(crate) fn matrix_flatten(matrix: &Array2<f64>) -> Array1<f64> {
    let [rows, cols] = matrix.shape();
    Array1::from_shape_fn([rows * cols], |[idx]| matrix[[idx / cols, idx % cols]])
}

pub(crate) fn matrix_from_flat(rows: usize, cols: usize, values: &Array1<f64>) -> Array2<f64> {
    assert_eq!(
        values.shape(),
        [rows * cols],
        "invariant: flattened vector length must match matrix shape"
    );
    Array2::from_shape_fn([rows, cols], |[row, col]| values[row * cols + col])
}

pub(crate) fn matrix_solve(matrix: &Array2<f64>, rhs: &Array1<f64>) -> Result<Array1<f64>> {
    matrix.solve(&rhs.view()).map_err(|err| {
        Error::Solver(format!(
            "Leto dense solve failed for {}x{} matrix and RHS length {}: {err}",
            matrix_rows(matrix),
            matrix_cols(matrix),
            vector_len(rhs)
        ))
    })
}

pub(crate) fn matrix_transpose_vector_mul(
    matrix: &Array2<f64>,
    vector: &Array1<f64>,
) -> Array1<f64> {
    let [rows, cols] = matrix.shape();
    assert_eq!(
        vector.shape(),
        [rows],
        "invariant: transposed matrix-vector dimensions must match"
    );
    Array1::from_shape_fn([cols], |[col]| {
        (0..rows).map(|row| matrix[[row, col]] * vector[row]).sum()
    })
}

pub(crate) fn row_vector(matrix: &Array2<f64>, row: usize) -> Array1<f64> {
    let [rows, cols] = matrix.shape();
    assert!(row < rows, "invariant: requested row is in bounds");
    Array1::from_shape_fn([cols], |[col]| matrix[[row, col]])
}

pub(crate) fn set_row(matrix: &mut Array2<f64>, row: usize, values: &Array1<f64>) {
    let [rows, cols] = matrix.shape();
    assert!(row < rows, "invariant: target row is in bounds");
    assert_eq!(
        values.shape(),
        [cols],
        "invariant: row value length must match matrix columns"
    );
    for col in 0..cols {
        matrix[[row, col]] = values[col];
    }
}

pub(crate) fn set_column(matrix: &mut Array2<f64>, col: usize, values: &Array1<f64>) {
    let [rows, cols] = matrix.shape();
    assert!(col < cols, "invariant: target column is in bounds");
    assert_eq!(
        values.shape(),
        [rows],
        "invariant: column value length must match matrix rows"
    );
    for row in 0..rows {
        matrix[[row, col]] = values[row];
    }
}

pub(crate) fn column_vector(matrix: &Array2<f64>, col: usize) -> Array1<f64> {
    let [rows, cols] = matrix.shape();
    assert!(col < cols, "invariant: requested column is in bounds");
    Array1::from_shape_fn([rows], |[row]| matrix[[row, col]])
}

/// Represents a DG solution on a single element
#[derive(Clone)]
pub struct DGSolution {
    /// Polynomial order
    pub order: usize,
    /// Number of components (for systems of equations)
    pub num_components: usize,
    /// Solution coefficients (num_components × num_basis_functions)
    pub coefficients: Array2<f64>,
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
            return Err(Error::InvalidInput(format!(
                "Polynomial order must be at least 1, got {order}"
            )));
        }

        let basis = DGBasis::new(order, basis_type).context("constructing DG basis functions")?;
        let num_basis = order + 1;

        Ok(Self {
            order,
            num_components,
            coefficients: matrix_zeros(num_components, num_basis),
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
    pub fn evaluate(&self, x: f64) -> Array1<f64> {
        let mut result = vector_zeros(self.num_components);

        for i in 0..matrix_cols(&self.coefficients) {
            let phi_i = self.basis.evaluate_basis(i, x);
            for c in 0..self.num_components {
                result[c] += self.coefficients[[c, i]] * phi_i;
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
            for j in 0..matrix_cols(&self.coefficients) {
                for k in 0..matrix_cols(&self.coefficients) {
                    // Mass matrix M_{j,k} = ∫ φ_j(x) φ_k(x) dx
                    let m_jk = self.basis.mass_matrix[[j, k]];
                    norm_sq += self.coefficients[[i, j]] * self.coefficients[[i, k]] * m_jk;
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
    pub fn apply_limiter<L: Limiter>(
        &mut self,
        limiter: &L,
        neighbors: &[DGSolution],
        params: &LimiterParams,
    ) -> Result<()> {
        if !params.adaptive || limiter.is_troubled_cell(self, neighbors, params) {
            limiter
                .limit(self, neighbors, params)
                .context("applying slope limiter to DG solution")?;
        }

        Ok(())
    }

    /// Get the cell average of the solution
    ///
    /// # Returns
    /// The cell average of the solution
    pub fn average(&self) -> Array1<f64> {
        let mut avg = vector_zeros(self.num_components);

        // The average is (1/2) * ∫ u(x) dx over [-1, 1]
        // ∫ u(x) dx = ∑_j c_j * ∫ φ_j(x) dx
        // The integrals of the basis functions can be computed using quadrature
        let mut basis_integrals = vector_zeros(matrix_cols(&self.coefficients));
        for j in 0..matrix_cols(&self.coefficients) {
            let mut int_phi_j = 0.0;
            for q in 0..vector_len(&self.basis.quad_points) {
                int_phi_j += self.basis.quad_weights[q] * self.basis.phi[[j, q]];
            }
            basis_integrals[j] = int_phi_j;
        }

        for i in 0..self.num_components {
            let mut sum = 0.0;
            for j in 0..matrix_cols(&self.coefficients) {
                sum += self.coefficients[[i, j]] * basis_integrals[j];
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
        sol.coefficients[[0, 0]] = 4.0 / 3.0; // P₀ term
        sol.coefficients[[0, 1]] = 1.0; // P₁ term
        sol.coefficients[[0, 2]] = 2.0 / 3.0; // P₂ term

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
        sol.coefficients[[0, 0]] = 2.0;

        // The average should be 2.0
        let avg = sol.average();
        assert_relative_eq!(avg[0], 2.0, epsilon = 1e-10);
    }
}
