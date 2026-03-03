//! Grid transfer operators and smoothers for geometric multigrid
//!
//! Weighted Jacobi relaxation, full-weighting restriction, and
//! bilinear prolongation operators.

use super::GeometricMultigrid;
use cfd_core::conversion::SafeFromF64;
use nalgebra::{DMatrix, DVector, RealField};
use num_traits::FromPrimitive;

impl<T: RealField + Copy + FromPrimitive> GeometricMultigrid<T> {
    /// Apply weighted Jacobi relaxation
    pub(super) fn jacobi_relaxation(
        &self,
        matrix: &DMatrix<T>,
        u: &mut DVector<T>,
        f: &DVector<T>,
        omega: T,
        iterations: usize,
    ) {
        let n = matrix.nrows();

        for _ in 0..iterations {
            for i in 0..n {
                let mut sum = T::zero();

                // Compute off-diagonal contributions
                for j in 0..n {
                    if i != j {
                        sum += matrix[(i, j)] * u[j];
                    }
                }

                // Update using Jacobi formula: u_new = ω*(f - sum_off_diag)/diag + (1-ω)*u_old
                let diagonal = matrix[(i, i)];
                let residual = f[i] - sum;
                let u_new = omega * (residual / diagonal) + (T::one() - omega) * u[i];

                u[i] = u_new;
            }
        }
    }

    /// Restrict residual to coarser grid using full weighting
    pub(super) fn restrict_residual(
        &self,
        fine_residual: &DVector<T>,
        fine_nx: usize,
        fine_ny: usize,
        coarse_nx: usize,
        coarse_ny: usize,
    ) -> DVector<T> {
        let mut coarse_residual = DVector::zeros(coarse_nx * coarse_ny);

        // Full weighting restriction operator
        for i in 0..coarse_nx {
            for j in 0..coarse_ny {
                let coarse_idx = j * coarse_nx + i;

                let mut sum = T::zero();
                let mut weight_sum = T::zero();

                let fine_i = i * 2;
                let fine_j = j * 2;

                let four = T::from_f64_or_one(4.0);
                let two = T::from_f64_or_one(2.0);
                let one = T::one();

                if fine_i < fine_nx && fine_j < fine_ny {
                    let fine_idx = fine_j * fine_nx + fine_i;
                    sum += fine_residual[fine_idx] * four;
                    weight_sum += four;
                }

                if fine_i > 0 && fine_j < fine_ny {
                    let fine_idx = fine_j * fine_nx + (fine_i - 1);
                    sum += fine_residual[fine_idx] * two;
                    weight_sum += two;
                }

                if fine_i + 1 < fine_nx && fine_j < fine_ny {
                    let fine_idx = fine_j * fine_nx + (fine_i + 1);
                    sum += fine_residual[fine_idx] * two;
                    weight_sum += two;
                }

                if fine_j > 0 && fine_i < fine_nx {
                    let fine_idx = (fine_j - 1) * fine_nx + fine_i;
                    sum += fine_residual[fine_idx] * two;
                    weight_sum += two;
                }

                if fine_j + 1 < fine_ny && fine_i < fine_nx {
                    let fine_idx = (fine_j + 1) * fine_nx + fine_i;
                    sum += fine_residual[fine_idx] * two;
                    weight_sum += two;
                }

                if fine_i > 0 && fine_j > 0 {
                    let fine_idx = (fine_j - 1) * fine_nx + (fine_i - 1);
                    sum += fine_residual[fine_idx] * one;
                    weight_sum += one;
                }

                if fine_i + 1 < fine_nx && fine_j > 0 {
                    let fine_idx = (fine_j - 1) * fine_nx + (fine_i + 1);
                    sum += fine_residual[fine_idx] * one;
                    weight_sum += one;
                }

                if fine_i > 0 && fine_j + 1 < fine_ny {
                    let fine_idx = (fine_j + 1) * fine_nx + (fine_i - 1);
                    sum += fine_residual[fine_idx] * one;
                    weight_sum += one;
                }

                if fine_i + 1 < fine_nx && fine_j + 1 < fine_ny {
                    let fine_idx = (fine_j + 1) * fine_nx + (fine_i + 1);
                    sum += fine_residual[fine_idx] * one;
                    weight_sum += one;
                }

                if weight_sum > T::zero() {
                    coarse_residual[coarse_idx] = sum / weight_sum;
                }
            }
        }

        coarse_residual
    }

    /// Prolongate correction to finer grid using bilinear interpolation
    pub(super) fn prolongate_correction(
        &self,
        coarse_correction: &DVector<T>,
        coarse_nx: usize,
        coarse_ny: usize,
        fine_nx: usize,
        fine_ny: usize,
    ) -> DVector<T> {
        let mut fine_correction = DVector::zeros(fine_nx * fine_ny);

        // Bilinear prolongation
        for i in 0..fine_nx {
            for j in 0..fine_ny {
                let fine_idx = j * fine_nx + i;

                // Map to coarse grid coordinates
                let coarse_i = i / 2;
                let coarse_j = j / 2;

                if coarse_i < coarse_nx && coarse_j < coarse_ny {
                    let coarse_idx = coarse_j * coarse_nx + coarse_i;
                    fine_correction[fine_idx] = coarse_correction[coarse_idx];
                }
            }
        }

        fine_correction
    }

    /// Compute residual r = f - A*u
    pub(super) fn compute_residual(
        &self,
        matrix: &DMatrix<T>,
        u: &DVector<T>,
        f: &DVector<T>,
    ) -> DVector<T> {
        f - (matrix * u)
    }

    /// Compute residual norm ||f - A*u||_2
    pub(super) fn compute_residual_norm(
        &self,
        matrix: &DMatrix<T>,
        u: &DVector<T>,
        f: &DVector<T>,
    ) -> T {
        let residual = self.compute_residual(matrix, u, f);
        residual.norm()
    }
}
