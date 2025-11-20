//! Algebraic multigrid preconditioner

use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;

// Multigrid constants
const COARSENING_THRESHOLD: f64 = 0.25;
const MIN_COARSE_SIZE: usize = 10;
const MAX_LEVELS: usize = 10;

/// Algebraic Multigrid (AMG) preconditioner
///
/// Implements a V-cycle multigrid method for sparse linear systems.
pub struct AlgebraicMultigrid<T: RealField + Copy> {
    /// Hierarchy of matrices at different levels
    levels: Vec<MultigridLevel<T>>,
    /// Number of pre-smoothing iterations
    pre_smooth: usize,
    /// Number of post-smoothing iterations
    post_smooth: usize,
}

/// Single level in the multigrid hierarchy
struct MultigridLevel<T: RealField + Copy> {
    /// System matrix at this level
    matrix: CsrMatrix<T>,
    /// Restriction operator (fine to coarse)
    restriction: Option<CsrMatrix<T>>,
    /// Prolongation operator (coarse to fine)
    prolongation: Option<CsrMatrix<T>>,
}

impl<T: RealField + Copy + FromPrimitive> AlgebraicMultigrid<T> {
    /// Create AMG preconditioner
    pub fn new(matrix: CsrMatrix<T>) -> Result<Self> {
        let levels = Self::build_hierarchy(matrix)?;

        Ok(Self {
            levels,
            pre_smooth: 2,
            post_smooth: 2,
        })
    }

    /// Build the multigrid hierarchy
    fn build_hierarchy(matrix: CsrMatrix<T>) -> Result<Vec<MultigridLevel<T>>> {
        let mut levels = Vec::new();
        let mut current_matrix = matrix;

        for _level in 0..MAX_LEVELS {
            let n = current_matrix.nrows();

            // Stop if matrix is small enough
            if n <= MIN_COARSE_SIZE {
                levels.push(MultigridLevel {
                    matrix: current_matrix,
                    restriction: None,
                    prolongation: None,
                });
                break;
            }

            // Build coarsening operators
            let (restriction, prolongation) = Self::build_transfer_operators(&current_matrix)?;

            // Compute coarse matrix: A_c = R * A * P
            let coarse_matrix =
                Self::galerkin_product(&restriction, &current_matrix, &prolongation)?;

            levels.push(MultigridLevel {
                matrix: current_matrix,
                restriction: Some(restriction.clone()),
                prolongation: Some(prolongation.clone()),
            });

            current_matrix = coarse_matrix;
        }

        if levels.is_empty() {
            return Err(Error::InvalidInput(
                "Failed to build multigrid hierarchy".to_string(),
            ));
        }

        Ok(levels)
    }

    /// Build restriction and prolongation operators
    fn build_transfer_operators(matrix: &CsrMatrix<T>) -> Result<(CsrMatrix<T>, CsrMatrix<T>)> {
        let n = matrix.nrows();
        let threshold = T::from_f64(COARSENING_THRESHOLD).unwrap_or_else(T::zero);

        // Simple coarsening: aggregate strongly connected nodes
        let mut coarse_map = vec![None; n];
        let mut coarse_count = 0;

        for i in 0..n {
            if coarse_map[i].is_none() {
                // Create new coarse node
                coarse_map[i] = Some(coarse_count);

                // Find strongly connected neighbors
                let row_start = matrix.row_offsets()[i];
                let row_end = matrix.row_offsets()[i + 1];
                let diag = Self::get_diagonal_value(matrix, i);

                for idx in row_start..row_end {
                    let j = matrix.col_indices()[idx];
                    if j != i && coarse_map[j].is_none() {
                        let val = matrix.values()[idx].abs();
                        if val > threshold * diag.abs() {
                            coarse_map[j] = Some(coarse_count);
                        }
                    }
                }

                coarse_count += 1;
            }
        }

        // Build prolongation matrix
        let mut p_vals = Vec::new();
        let mut p_indices = Vec::new();
        let mut p_offsets = vec![0];

        for i in 0..n {
            if let Some(c) = coarse_map[i] {
                p_indices.push(c);
                p_vals.push(T::one());
            }
            p_offsets.push(p_indices.len());
        }

        let prolongation =
            CsrMatrix::try_from_csr_data(n, coarse_count, p_offsets, p_indices, p_vals).map_err(
                |_| Error::InvalidInput("Failed to create prolongation matrix".to_string()),
            )?;

        // Restriction is transpose of prolongation (for simplicity)
        let restriction = prolongation.transpose();

        Ok((restriction, prolongation))
    }

    /// Compute Galerkin product: R * A * P
    fn galerkin_product(
        r: &CsrMatrix<T>,
        a: &CsrMatrix<T>,
        p: &CsrMatrix<T>,
    ) -> Result<CsrMatrix<T>> {
        // First compute A * P
        let ap = a * p;
        // Then compute R * (A * P)
        let rap = r * ap;
        Ok(rap)
    }

    /// Get diagonal value of matrix
    fn get_diagonal_value(matrix: &CsrMatrix<T>, row: usize) -> T {
        let row_start = matrix.row_offsets()[row];
        let row_end = matrix.row_offsets()[row + 1];

        for idx in row_start..row_end {
            if matrix.col_indices()[idx] == row {
                return matrix.values()[idx];
            }
        }

        T::one()
    }

    /// Smooth using weighted Jacobi
    fn smooth(&self, level: usize, b: &DVector<T>, x: &mut DVector<T>, iterations: usize) {
        let matrix = &self.levels[level].matrix;
        let n = matrix.nrows();
        let omega = T::from_f64(0.67).unwrap_or_else(T::one); // Standard Jacobi weight

        for _ in 0..iterations {
            let x_old = x.clone();

            for i in 0..n {
                let mut sum = b[i];
                let mut diag = T::one();

                let row_start = matrix.row_offsets()[i];
                let row_end = matrix.row_offsets()[i + 1];

                for idx in row_start..row_end {
                    let j = matrix.col_indices()[idx];
                    let val = matrix.values()[idx];

                    if j == i {
                        diag = val;
                    } else {
                        sum -= val * x_old[j];
                    }
                }

                x[i] = omega * sum / diag + (T::one() - omega) * x_old[i];
            }
        }
    }

    /// V-cycle recursive implementation
    fn v_cycle(&self, level: usize, b: &DVector<T>, x: &mut DVector<T>) {
        if level == self.levels.len() - 1 {
            // Coarsest level - solve directly
            self.smooth(level, b, x, 10);
        } else {
            // Pre-smooth
            self.smooth(level, b, x, self.pre_smooth);

            // Compute residual
            let matrix = &self.levels[level].matrix;
            let residual = b - matrix * &*x;

            // Restrict residual to coarse grid
            let restriction = self.levels[level].restriction.as_ref().unwrap();
            let coarse_b = restriction * &residual;

            // Solve on coarse grid
            let mut coarse_x = DVector::zeros(coarse_b.len());
            self.v_cycle(level + 1, &coarse_b, &mut coarse_x);

            // Prolongate correction to fine grid
            let prolongation = self.levels[level].prolongation.as_ref().unwrap();
            let correction = prolongation * coarse_x;

            // Apply correction
            *x += correction;

            // Post-smooth
            self.smooth(level, b, x, self.post_smooth);
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> Preconditioner<T> for AlgebraicMultigrid<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        // Initialize z to zero
        z.fill(T::zero());

        // Apply one V-cycle
        self.v_cycle(0, r, z);

        Ok(())
    }
}
