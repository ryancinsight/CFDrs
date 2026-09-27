use super::basis::{BasisType, DGBasis};
use super::{Limiter, LimiterParams, matrix_cols, matrix_zeros, vector_len, vector_zeros};
use crate::error::Result;
use cfd_core::error::{Error, ErrorContext};
use leto::{Array1, Array2};
use std::fmt;

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
