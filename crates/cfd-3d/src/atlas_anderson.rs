//! Leto boundary helpers for Anderson acceleration over FEM solution storage.

use crate::fem::FemDofVector;
use cfd_math::nonlinear_solver::AndersonAccelerator;
use leto::Array1;

/// Run Anderson acceleration through the Leto-backed cfd-math boundary.
///
/// The FEM solution type stores velocity in Leto-backed DOF vectors; this
/// function keeps Anderson acceleration on the canonical Leto array boundary.
pub(crate) fn accelerate_velocity(
    accelerator: &mut AndersonAccelerator<f64>,
    previous_velocity: &FemDofVector<f64>,
    current_velocity: &FemDofVector<f64>,
) -> FemDofVector<f64> {
    let previous = Array1::from_shape_vec(
        [previous_velocity.len()],
        previous_velocity.as_slice().to_vec(),
    )
    .expect("invariant: previous velocity vector length matches Leto shape");
    let current = Array1::from_shape_vec(
        [current_velocity.len()],
        current_velocity.as_slice().to_vec(),
    )
    .expect("invariant: current velocity vector length matches Leto shape");
    let accelerated = accelerator.compute_next(&previous, &current);
    FemDofVector::from_vec(accelerated.into_vec())
}
