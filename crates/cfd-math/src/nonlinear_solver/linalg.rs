use eunomia::NumericElement;
use leto::Array1;
use leto_ops::Scalar;

#[inline]
pub(super) fn vector_zeros<T: Scalar>(len: usize) -> Array1<T> {
    Array1::from_elem([len], T::ZERO)
}

#[inline]
pub(super) fn vector_len<T>(vector: &Array1<T>) -> usize {
    vector.shape()[0]
}

#[inline]
pub(super) fn dot<T: Scalar>(lhs: &Array1<T>, rhs: &Array1<T>) -> T {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: vector dot requires equal lengths"
    );
    T::dot_slice(lhs.as_slice().unwrap(), rhs.as_slice().unwrap())
}

#[inline]
pub(super) fn norm<T: Scalar>(vector: &Array1<T>) -> T {
    NumericElement::sqrt(dot(vector, vector))
}

#[inline]
pub(super) fn add<T: Scalar>(lhs: &Array1<T>, rhs: &Array1<T>) -> Array1<T> {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: vector addition requires equal lengths"
    );
    let n = vector_len(lhs);
    let mut out = Array1::from_elem([n], T::ZERO);
    T::add_slice(
        lhs.as_slice().unwrap(),
        rhs.as_slice().unwrap(),
        out.as_slice_mut().unwrap(),
    );
    out
}

#[inline]
pub(super) fn sub<T: Scalar>(lhs: &Array1<T>, rhs: &Array1<T>) -> Array1<T> {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: vector subtraction requires equal lengths"
    );
    let n = vector_len(lhs);
    let mut out = Array1::from_elem([n], T::ZERO);
    T::sub_slice(
        lhs.as_slice().unwrap(),
        rhs.as_slice().unwrap(),
        out.as_slice_mut().unwrap(),
    );
    out
}

#[inline]
pub(super) fn add_scaled<T: Scalar>(lhs: &Array1<T>, rhs: &Array1<T>, scale: T) -> Array1<T> {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: scaled vector addition requires equal lengths"
    );
    let mut out = lhs.clone();
    T::axpy_slice(scale, rhs.as_slice().unwrap(), out.as_slice_mut().unwrap());
    out
}

#[inline]
pub(super) fn add_scaled_in_place<T: Scalar>(lhs: &mut Array1<T>, rhs: &Array1<T>, scale: T) {
    assert_eq!(
        vector_len(lhs),
        vector_len(rhs),
        "invariant: scaled vector update requires equal lengths"
    );
    T::axpy_slice(
        scale,
        rhs.as_slice().unwrap(),
        lhs.as_slice_mut().unwrap(),
    );
}

#[inline]
pub(super) fn sub_in_place_scaled<T: Scalar>(lhs: &mut Array1<T>, rhs: &Array1<T>, scale: T) {
    add_scaled_in_place(lhs, rhs, T::ZERO - scale);
}

#[inline]
pub(super) fn scale<T: Scalar>(vector: &Array1<T>, factor: T) -> Array1<T> {
    let n = vector_len(vector);
    let mut out = Array1::from_elem([n], T::ZERO);
    T::axpy_slice(
        factor,
        vector.as_slice().unwrap(),
        out.as_slice_mut().unwrap(),
    );
    out
}

#[inline]
pub(super) fn neg<T: Scalar>(vector: &Array1<T>) -> Array1<T> {
    scale(vector, T::ZERO - T::ONE)
}
