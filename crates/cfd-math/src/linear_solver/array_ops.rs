use cfd_core::error::{Error, Result};
use eunomia::{NumericElement, RealField};
use leto::Array1;

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

pub(super) fn dot<T: RealField + Copy + NumericElement>(lhs: &Array1<T>, rhs: &Array1<T>) -> T {
    let mut sum = <T as NumericElement>::ZERO;
    for idx in 0..vector_len(lhs) {
        sum += lhs[idx] * rhs[idx];
    }
    sum
}

pub(super) fn norm<T: RealField + Copy + NumericElement>(vector: &Array1<T>) -> T {
    NumericElement::sqrt(dot(vector, vector))
}

pub(super) fn copy_array<T: Copy>(src: &Array1<T>, dst: &mut Array1<T>) {
    for idx in 0..vector_len(src) {
        dst[idx] = src[idx];
    }
}

pub(super) fn assign_residual<T: RealField + Copy + NumericElement>(
    residual: &mut Array1<T>,
    rhs: &Array1<T>,
    ax: &Array1<T>,
) {
    for idx in 0..vector_len(rhs) {
        residual[idx] = rhs[idx] - ax[idx];
    }
}

pub(super) fn axpy<T: RealField + Copy>(x: &mut Array1<T>, alpha: T, y: &Array1<T>) {
    for idx in 0..vector_len(x) {
        x[idx] += alpha * y[idx];
    }
}

pub(super) fn scale_add<T: RealField + Copy>(x: &mut Array1<T>, scale: T, y: &Array1<T>) {
    for idx in 0..vector_len(x) {
        x[idx] = x[idx] * scale + y[idx];
    }
}
