use super::{GmgMatrix, GmgVector};
use eunomia::{FloatElement, NumericElement, RealField};

#[inline]
pub(super) fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
pub(super) fn from_usize<T: FloatElement>(value: usize) -> T {
    let value_u64 = u64::try_from(value).expect("invariant: grid dimension fits into u64");
    from_f64(<u64 as NumericElement>::to_f64(value_u64))
}

pub(super) fn l2_norm<T: RealField>(vector: &GmgVector<T>) -> T {
    let mut sum = <T as NumericElement>::ZERO;
    for i in 0..vector.shape()[0] {
        sum += vector[i] * vector[i];
    }
    <T as NumericElement>::sqrt(sum)
}

pub(super) fn matrix_vector_product<T: RealField>(
    matrix: &GmgMatrix<T>,
    vector: &GmgVector<T>,
) -> GmgVector<T> {
    let [rows, cols] = matrix.shape();
    debug_assert_eq!(cols, vector.shape()[0]);
    let mut output = GmgVector::zeros([rows]);
    for row in 0..rows {
        let mut sum = <T as NumericElement>::ZERO;
        for col in 0..cols {
            sum += matrix[[row, col]] * vector[col];
        }
        output[row] = sum;
    }
    output
}

pub(super) fn vector_add<T: RealField>(lhs: &GmgVector<T>, rhs: &GmgVector<T>) -> GmgVector<T> {
    debug_assert_eq!(lhs.shape(), rhs.shape());
    let mut output = GmgVector::zeros([lhs.shape()[0]]);
    for i in 0..lhs.shape()[0] {
        output[i] = lhs[i] + rhs[i];
    }
    output
}

pub(super) fn vector_sub<T: RealField>(lhs: &GmgVector<T>, rhs: &GmgVector<T>) -> GmgVector<T> {
    debug_assert_eq!(lhs.shape(), rhs.shape());
    let mut output = GmgVector::zeros([lhs.shape()[0]]);
    for i in 0..lhs.shape()[0] {
        output[i] = lhs[i] - rhs[i];
    }
    output
}

pub(super) fn vector_add_assign<T: RealField>(lhs: &mut GmgVector<T>, rhs: &GmgVector<T>) {
    debug_assert_eq!(lhs.shape(), rhs.shape());
    for i in 0..lhs.shape()[0] {
        lhs[i] += rhs[i];
    }
}
