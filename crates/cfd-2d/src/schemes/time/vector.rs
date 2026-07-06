//! Leto-backed state vector helpers for time integration.

use crate::scalar::Cfd2dScalar;
use leto::Array1;
use leto_ops::norm_l2;

/// Owned one-dimensional time-integration state vector.
pub type StateVector<T> = Array1<T>;

#[inline]
pub(super) fn l2_norm<T: Cfd2dScalar>(state: &StateVector<T>) -> T {
    norm_l2(&state.view()).expect("invariant: Leto state vector has a valid layout")
}
