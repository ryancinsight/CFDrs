//! Diagnostic and residual computation for SIMPLEC/PIMPLE solvers
//!
//! Provides velocity residual, continuity residual, and convergence monitoring
//! utilities used by both the SIMPLEC and PIMPLE algorithm implementations.
//!
//! # Theorem (Discrete Continuity Residual)
//!
//! For incompressible flows on a collocated grid the discrete continuity constraint
//! is `∇·u = 0`. The continuity residual `‖∇·u‖_∞` measures the maximum cell-wise
//! divergence and must converge to zero (within solver tolerance) for the velocity
//! field to be physically admissible.
//!
//! **Proof sketch**: By the discrete divergence theorem, the net flux through all faces
//! of a control volume must vanish. The infinity-norm of the cell divergence therefore
//! equals the worst-case mass imbalance. When the pressure-correction equation is solved
//! to residual `ε`, the continuity residual is bounded by `O(ε / Δx²)` via the
//! Laplacian scaling of the Poisson operator.

use super::solver::SimplecPimpleSolver;
use crate::fields::SimulationFields;
use nalgebra::{RealField, Vector2};
use num_traits::{FromPrimitive, ToPrimitive};

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::LowerExp>
    SimplecPimpleSolver<T>
{
    /// Calculate velocity residual between two velocity fields (L∞ norm)
    pub(super) fn calculate_velocity_residual_from_vectors(
        &self,
        old_velocity: &[Vec<Vector2<T>>],
        new_velocity: &[Vec<Vector2<T>>],
    ) -> T {
        let mut max_residual = T::zero();
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                let residual = (new_velocity[i][j] - old_velocity[i][j]).norm();
                if residual > max_residual {
                    max_residual = residual;
                }
            }
        }
        max_residual
    }

    /// Calculate continuity residual `‖∇·u‖_∞` from cell-centred fields
    ///
    /// Uses central-difference approximation of the divergence operator.
    pub(super) fn calculate_continuity_residual(&self, fields: &SimulationFields<T>) -> T {
        let mut max_divergence = T::zero();
        let two = T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);

        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                let du_dx = (fields.u.at(i + 1, j) - fields.u.at(i - 1, j)) / (two * self.grid.dx);
                let dv_dy = (fields.v.at(i, j + 1) - fields.v.at(i, j - 1)) / (two * self.grid.dy);
                let abs_div = (du_dx + dv_dy).abs();
                if abs_div > max_divergence {
                    max_divergence = abs_div;
                }
            }
        }
        max_divergence
    }

    /// Calculate continuity residual from Rhie-Chow consistent face velocities
    ///
    /// More accurate than cell-centred divergence when Rhie-Chow interpolation
    /// is active, because it uses the same flux balance that the pressure
    /// correction equation solves.
    pub(super) fn calculate_continuity_residual_from_faces(
        &self,
        face_velocity: &[Vec<Vector2<T>>],
    ) -> T {
        let mut max_divergence = T::zero();

        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                let du_dx = (face_velocity[i][j].x - face_velocity[i - 1][j].x) / self.grid.dx;
                let dv_dy = (face_velocity[i][j].y - face_velocity[i][j - 1].y) / self.grid.dy;
                let abs_div = (du_dx + dv_dy).abs();
                if abs_div > max_divergence {
                    max_divergence = abs_div;
                }
            }
        }
        max_divergence
    }

    /// Compute continuity residual with Rhie-Chow awareness
    ///
    /// Selects the appropriate residual computation method based on whether
    /// Rhie-Chow interpolation is active.
    pub(super) fn compute_continuity_residual(
        &self,
        fields: &SimulationFields<T>,
        dt: Option<T>,
    ) -> T {
        if let Some(ref rhie_chow) = self.rhie_chow {
            let u_consistent = self.interpolate_consistent_velocity(rhie_chow, fields, dt);
            self.calculate_continuity_residual_from_faces(&u_consistent)
        } else {
            self.calculate_continuity_residual(fields)
        }
    }

    /// Extract velocity field from simulation fields as 2D vector array
    pub(super) fn extract_velocity_field(
        &self,
        fields: &SimulationFields<T>,
    ) -> Vec<Vec<Vector2<T>>> {
        let mut velocity = vec![vec![Vector2::zeros(); self.grid.ny]; self.grid.nx];
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                velocity[i][j] = Vector2::new(fields.u.at(i, j), fields.v.at(i, j));
            }
        }
        velocity
    }

    /// Convert `Field2D<T>` to `Vec<Vec<T>>`
    pub(super) fn field2d_to_vec2d(&self, field: &crate::fields::Field2D<T>) -> Vec<Vec<T>> {
        let mut result = vec![vec![T::zero(); self.grid.ny]; self.grid.nx];
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                result[i][j] = field[(i, j)];
            }
        }
        result
    }

    /// Convert `Vec<Vec<T>>` to `Field2D<T>`
    pub(super) fn vec2d_to_field2d(&self, field: &mut crate::fields::Field2D<T>, vec: &[Vec<T>]) {
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                field[(i, j)] = vec[i][j];
            }
        }
    }
}
