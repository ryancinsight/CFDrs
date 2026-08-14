//! Leto-backed state vector helpers for time integration.

use cfd_core::CfdScalar;
use leto::Array1;
use leto_ops::norm_l2;

/// Owned one-dimensional time-integration state vector.
pub type StateVector<T> = Array1<T>;

#[inline]
pub(super) fn l2_norm<T: CfdScalar>(state: &StateVector<T>) -> T {
    norm_l2(&state.view()).expect("invariant: Leto state vector has a valid layout")
}
