//! Internal primitive array ops — retained for potential future use.
#![allow(dead_code)]

use cfd_core::error::{Error, Result};
use leto::Array1;
use leto_ops::Scalar;

#[inline]
pub(super) fn vector_len<T>(vector: &Array1<T>) -> usize {
    vector.shape()[0]
}

pub(super) fn validate_vector_len<T>(
    name: &str,
    vector: &Array1<T>,
    expected: usize,
) -> Result<()> {
    let actual = vector_len(vector);
    if actual != expected {
        return Err(Error::InvalidConfiguration(format!(
            "{name} length mismatch: expected {expected}, got {actual}"
        )));
    }
    Ok(())
}

pub(super) fn dot<T: Scalar>(lhs: &Array1<T>, rhs: &Array1<T>) -> T {
    T::dot_slice(lhs.as_slice().unwrap(), rhs.as_slice().unwrap())
}

pub(super) fn norm<T: Scalar>(vector: &Array1<T>) -> T {
    eunomia::NumericElement::sqrt(dot(vector, vector))
}

pub(super) fn copy_array<T: Scalar>(src: &Array1<T>, dst: &mut Array1<T>) {
    dst.as_slice_mut()
        .unwrap()
        .copy_from_slice(src.as_slice().unwrap());
}

pub(super) fn assign_residual<T: Scalar>(residual: &mut Array1<T>, rhs: &Array1<T>, ax: &Array1<T>) {
    T::sub_slice(
        rhs.as_slice().unwrap(),
        ax.as_slice().unwrap(),
        residual.as_slice_mut().unwrap(),
    );
}

pub(super) fn axpy<T: Scalar>(x: &mut Array1<T>, alpha: T, y: &Array1<T>) {
    T::axpy_slice(
        alpha,
        y.as_slice().unwrap(),
        x.as_slice_mut().unwrap(),
    );
}

pub(super) fn scale_add<T: Scalar>(x: &mut Array1<T>, scale: T, y: &Array1<T>) {
    let x_slice = x.as_slice_mut().unwrap();
    let y_slice = y.as_slice().unwrap();
    for i in 0..x_slice.len() {
        x_slice[i] = x_slice[i] * scale + y_slice[i];
    }
}
