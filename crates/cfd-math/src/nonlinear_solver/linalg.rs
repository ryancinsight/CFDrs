use eunomia::NumericElement;
use leto::{Array1, Array2};

#[inline]
pub(super) fn vector_from_vec<T>(values: Vec<T>) -> Array1<T> {
    Array1::from_shape_vec([values.len()], values)
        .expect("invariant: vector shape matches produced element count")
}

#[inline]
pub(super) fn vector_zeros<T: NumericElement>(len: usize) -> Array1<T> {
    Array1::from_elem([len], T::ZERO)
}

#[inline]
pub(super) fn matrix_zeros<T: NumericElement>(rows: usize, cols: usize) -> Array2<T> {
    Array2::from_elem([rows, cols], T::ZERO)
}

#[inline]
pub(super) fn vector_len<T>(vector: &Array1<T>) -> usize {
    vector.shape()[0]
}

#[inline]
pub(super) fn dot<T: NumericElement>(lhs: &Array1<T>, rhs: &Array1<T>) -> T {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: vector dot requires equal lengths"
    );
    (0..vector_len(lhs)).fold(T::ZERO, |acc, idx| acc + lhs[idx] * rhs[idx])
}

#[inline]
pub(super) fn norm<T: NumericElement>(vector: &Array1<T>) -> T {
    NumericElement::sqrt(dot(vector, vector))
}

#[inline]
pub(super) fn add<T: NumericElement>(lhs: &Array1<T>, rhs: &Array1<T>) -> Array1<T> {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: vector addition requires equal lengths"
    );
    vector_from_vec(
        (0..vector_len(lhs))
            .map(|idx| lhs[idx] + rhs[idx])
            .collect(),
    )
}

#[inline]
pub(super) fn sub<T: NumericElement>(lhs: &Array1<T>, rhs: &Array1<T>) -> Array1<T> {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: vector subtraction requires equal lengths"
    );
    vector_from_vec(
        (0..vector_len(lhs))
            .map(|idx| lhs[idx] - rhs[idx])
            .collect(),
    )
}

#[inline]
pub(super) fn add_scaled<T: NumericElement>(
    lhs: &Array1<T>,
    rhs: &Array1<T>,
    scale: T,
) -> Array1<T> {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: scaled vector addition requires equal lengths"
    );
    vector_from_vec(
        (0..vector_len(lhs))
            .map(|idx| lhs[idx] + rhs[idx] * scale)
            .collect(),
    )
}

#[inline]
pub(super) fn add_scaled_in_place<T: NumericElement>(
    lhs: &mut Array1<T>,
    rhs: &Array1<T>,
    scale: T,
) {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: scaled vector update requires equal lengths"
    );
    for idx in 0..vector_len(lhs) {
        lhs[idx] += rhs[idx] * scale;
    }
}

#[inline]
pub(super) fn sub_in_place_scaled<T: NumericElement>(
    lhs: &mut Array1<T>,
    rhs: &Array1<T>,
    scale: T,
) {
    add_scaled_in_place(lhs, rhs, T::ZERO - scale);
}

#[inline]
pub(super) fn scale<T: NumericElement>(vector: &Array1<T>, factor: T) -> Array1<T> {
    vector_from_vec(
        (0..vector_len(vector))
            .map(|idx| vector[idx] * factor)
            .collect(),
    )
}

#[inline]
pub(super) fn neg<T: NumericElement>(vector: &Array1<T>) -> Array1<T> {
    scale(vector, T::ZERO - T::ONE)
}
