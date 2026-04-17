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
use crate::grid::array2d::Array2D;
use crate::fields::SimulationFields;
use crate::solvers::continuity::{max_central_continuity_residual, max_face_continuity_residual};
use nalgebra::{RealField, Vector2};
use num_traits::{FromPrimitive, ToPrimitive};

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::LowerExp>
    SimplecPimpleSolver<T>
{
    /// Calculate velocity residual between two velocity fields (L∞ norm)
    pub(super) fn calculate_velocity_residual_from_vectors(
        &self,
        old_velocity: &Array2D<Vector2<T>>,
        new_velocity: &Array2D<Vector2<T>>,
    ) -> T {
        let mut max_residual = T::zero();
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                let residual = (new_velocity[(i, j)] - old_velocity[(i, j)]).norm();
                if residual > max_residual {
                    max_residual = residual;
                }
            }
        }
        max_residual
    }

    /// Calculate continuity residual `‖∇·u‖_∞` from cell-centred fields
    ///
    /// Uses the standard central-difference discretization of the divergence operator.
    pub(super) fn calculate_continuity_residual(&self, fields: &SimulationFields<T>) -> T {
        max_central_continuity_residual(
            self.grid.nx,
            self.grid.ny,
            self.grid.dx,
            self.grid.dy,
            |i, j| fields.u.at(i, j),
            |i, j| fields.v.at(i, j),
        )
    }

    /// Calculate continuity residual from Rhie-Chow consistent face velocities
    ///
    /// More accurate than cell-centred divergence when Rhie-Chow interpolation
    /// is active, because it uses the same flux balance that the pressure
    /// correction equation solves.
    pub(super) fn calculate_continuity_residual_from_faces(
        &self,
        face_velocity: &Array2D<Vector2<T>>,
    ) -> T {
        max_face_continuity_residual(
            self.grid.nx,
            self.grid.ny,
            self.grid.dx,
            self.grid.dy,
            |i, j| face_velocity[(i, j)].x,
            |i, j| face_velocity[(i, j)].y,
        )
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
            let mut vfc = self._vel_field_cache.borrow_mut();
            if vfc.as_ref().is_none_or(|v| {
                let (nx, ny) = v.dimensions();
                nx != self.grid.nx || ny != self.grid.ny
            }) {
                *vfc = Some(crate::fields::Field2D::new(self.grid.nx, self.grid.ny, Vector2::zeros()));
            }
            let velocity_field = vfc.as_mut().unwrap();

            let mut cvc = self._cons_vel_cache.borrow_mut();
            if cvc
                .as_ref()
                .is_none_or(|v| v.rows() != self.grid.nx || v.cols() != self.grid.ny)
            {
                *cvc = Some(Array2D::new(self.grid.nx, self.grid.ny, Vector2::zeros()));
            }
            let consistent_velocity = cvc.as_mut().unwrap();

            self.interpolate_consistent_velocity(rhie_chow, fields, dt, velocity_field, consistent_velocity);
            self.calculate_continuity_residual_from_faces(consistent_velocity)
        } else {
            self.calculate_continuity_residual(fields)
        }
    }

    /// Extract velocity field from simulation fields as 2D vector array
    pub(super) fn extract_velocity_field(
        &self,
        fields: &SimulationFields<T>,
    ) -> Array2D<Vector2<T>> {
        let mut velocity = Array2D::new(self.grid.nx, self.grid.ny, Vector2::zeros());
        velocity.as_mut_slice()
            .iter_mut()
            .zip(fields.u.as_slice().iter())
            .zip(fields.v.as_slice().iter())
            .for_each(|((vel, &u), &v)| {
                *vel = Vector2::new(u, v);
            });
        velocity
    }
}
