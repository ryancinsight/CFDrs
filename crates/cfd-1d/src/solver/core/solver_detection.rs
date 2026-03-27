//! Solver detection and validation utilities.
//!
//! Provides SPD detection for assembled conductance matrices and
//! validation of linear systems before solve.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use super::{LinearSolverMethod, NetworkSolver};
use crate::domain::network::Network;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use cfd_core::physics::fluid::FluidTrait;

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + Float, F: FluidTrait<T> + Clone>
    NetworkSolver<T, F>
{
    /// Detect whether the network has only linear (flow-independent) resistances.
    pub(super) fn is_linear_static_network(network: &Network<T, F>) -> bool {
        let eps = T::default_epsilon();
        let all_edges_linear = network
            .graph
            .edge_weights()
            .all(|edge| <T as Float>::abs(edge.quad_coeff) <= eps);
        let no_geometry_updates = network
            .properties
            .values()
            .all(|props| props.geometry.is_none());
        all_edges_linear && no_geometry_updates
    }

    /// Detect whether the assembled matrix is SPD via diagonal dominance check.
    ///
    /// The Laplacian sign structure (positive diagonal, non-positive off-diagonal)
    /// is topologically invariant, so this classification is stable across
    /// Picard iterations.
    pub(super) fn detect_solver_method(
        matrix: &nalgebra_sparse::CsrMatrix<T>,
    ) -> LinearSolverMethod {
        let mut is_spd = true;
        for i in 0..matrix.nrows() {
            let row = matrix.row(i);
            let mut diag = T::zero();
            let mut sum_off = T::zero();
            for (j, val) in row.col_indices().iter().zip(row.values()) {
                if *j == i {
                    diag = *val;
                } else if *val > T::zero() {
                    is_spd = false;
                    break;
                } else {
                    sum_off += <T as Float>::abs(*val);
                }
            }
            if !is_spd {
                break;
            }
            let is_identity_dirichlet = diag == T::one() && sum_off == T::zero();
            if (diag < sum_off || diag <= T::zero()) && !is_identity_dirichlet {
                is_spd = false;
                break;
            }
        }
        if is_spd {
            LinearSolverMethod::ConjugateGradient
        } else {
            LinearSolverMethod::BiCGSTAB
        }
    }

    /// Validate assembled linear system for well-formedness.
    pub(super) fn validate_linear_system(
        matrix: &nalgebra_sparse::CsrMatrix<T>,
        rhs: &nalgebra::DVector<T>,
    ) -> Result<()> {
        if matrix.nrows() == 0 || matrix.ncols() == 0 {
            return Err(Error::InvalidConfiguration(
                "Assembled network matrix is empty".to_string(),
            ));
        }
        for row_idx in 0..matrix.nrows() {
            let row = matrix.row(row_idx);
            for value in row.values() {
                if !value.is_finite() {
                    return Err(Error::Numerical(NumericalErrorKind::InvalidValue {
                        value: format!("matrix[{row_idx}] is non-finite"),
                    }));
                }
            }
        }
        if rhs.iter().any(|value| !value.is_finite()) {
            return Err(Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "RHS contains non-finite entries".to_string(),
            }));
        }
        Ok(())
    }

    /// Compute the L2 norm of the linear-system residual ||Ax - b||₂.
    pub(super) fn compute_residual_norm(
        matrix: &nalgebra_sparse::CsrMatrix<T>,
        solution: &nalgebra::DVector<T>,
        rhs: &nalgebra::DVector<T>,
        n: usize,
    ) -> T {
        let mut norm = T::zero();
        for i in 0..n {
            let row = matrix.row(i);
            let mut ax_i = T::zero();
            for (j, val) in row.col_indices().iter().zip(row.values()) {
                ax_i += *val * solution[*j];
            }
            let r_i = ax_i - rhs[i];
            norm += r_i * r_i;
        }
        <T as Float>::sqrt(norm)
    }

    /// Check all entries in a vector are finite.
    pub(super) fn vector_is_finite(values: &nalgebra::DVector<T>) -> bool {
        values.iter().all(|value| value.is_finite())
    }

    /// Convert a scalar T to f64 for diagnostics.
    pub(super) fn scalar_to_f64(value: T) -> Option<f64> {
        value.to_f64()
    }
}
