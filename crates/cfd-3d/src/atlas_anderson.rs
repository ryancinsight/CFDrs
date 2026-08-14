//! Leto boundary helpers for Anderson acceleration over FEM solution storage.

use crate::fem::FemDofVector;
use cfd_math::nonlinear_solver::AndersonAccelerator;

/// Run Anderson acceleration through the Leto-backed cfd-math boundary.
///
/// The FEM solution type stores velocity in Leto-backed DOF vectors; this
/// function keeps Anderson acceleration on the canonical Leto array boundary.
pub(crate) fn accelerate_velocity(
    accelerator: &mut AndersonAccelerator<f64>,
    previous_velocity: &FemDofVector<f64>,
    current_velocity: &FemDofVector<f64>,
) -> FemDofVector<f64> {
    let accelerated =
        accelerator.compute_next(previous_velocity.as_array(), current_velocity.as_array());
    FemDofVector::from_vec(accelerated.into_vec())
}
