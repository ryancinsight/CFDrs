use crate::{scalar, scalar::Cfd3dScalar};
use leto::{Array1, Array2, FixedMatrix};
use leto::geometry::Vector3;

pub(crate) type Matrix3<T> = FixedMatrix<T, 3, 3>;
pub(crate) type Matrix3x4<T> = FixedMatrix<T, 3, 4>;

#[inline]
pub(crate) fn array1_len<T>(array: &Array1<T>) -> usize {
    array.shape()[0]
}

#[inline]
pub(crate) fn array1_subarray<T: Copy>(array: &Array1<T>, start: usize, len: usize) -> Array1<T> {
    Array1::from_shape_fn([len], |[i]| array[start + i])
}

#[inline]
pub(crate) fn array1_copy<T: Copy>(source: &Array1<T>, target: &mut Array1<T>) {
    let source = source
        .as_slice()
        .expect("invariant: solver arrays are dense one-dimensional Leto arrays");
    let target = target
        .as_slice_mut()
        .expect("invariant: solver arrays are dense one-dimensional Leto arrays");
    target.copy_from_slice(source);
}

#[inline]
pub(crate) fn array1_l2_norm<T: Cfd3dScalar>(array: &Array1<T>) -> T {
    eunomia::NumericElement::sqrt(
        array
            .iter()
            .fold(scalar::zero::<T>(), |acc, &value| acc + value * value),
    )
}

#[inline]
pub(crate) fn vector3_from_indexed<T: Copy, V>(value: &V) -> Vector3<T>
where
    V: core::ops::Index<usize, Output = T>,
{
    Vector3::new(value[0], value[1], value[2])
}

#[inline]
pub(crate) fn matrix3_from_columns<T: Cfd3dScalar>(
    c0: Vector3<T>,
    c1: Vector3<T>,
    c2: Vector3<T>,
) -> Matrix3<T> {
    FixedMatrix::from_rows([
        [c0[0], c1[0], c2[0]],
        [c0[1], c1[1], c2[1]],
        [c0[2], c1[2], c2[2]],
    ])
}

#[inline]
pub(crate) fn matrix3x4_from_columns<T: Cfd3dScalar>(columns: [Vector3<T>; 4]) -> Matrix3x4<T> {
    FixedMatrix::from_rows([
        [columns[0][0], columns[1][0], columns[2][0], columns[3][0]],
        [columns[0][1], columns[1][1], columns[2][1], columns[3][1]],
        [columns[0][2], columns[1][2], columns[2][2], columns[3][2]],
    ])
}

#[inline]
pub(crate) fn reference_tet_gradients<T: Cfd3dScalar>() -> Matrix3x4<T> {
    FixedMatrix::from_rows([
        [
            -scalar::one::<T>(),
            scalar::one::<T>(),
            scalar::zero::<T>(),
            scalar::zero::<T>(),
        ],
        [
            -scalar::one::<T>(),
            scalar::zero::<T>(),
            scalar::one::<T>(),
            scalar::zero::<T>(),
        ],
        [
            -scalar::one::<T>(),
            scalar::zero::<T>(),
            scalar::zero::<T>(),
            scalar::one::<T>(),
        ],
    ])
}

#[inline]
pub(crate) fn matrix3x4_column<T: Copy>(matrix: &Matrix3x4<T>, column: usize) -> Vector3<T> {
    Vector3::new(matrix[(0, column)], matrix[(1, column)], matrix[(2, column)])
}

#[inline]
pub(crate) fn array2_column3<T: Copy>(matrix: &Array2<T>, column: usize) -> Vector3<T> {
    Vector3::new(
        matrix[[0, column]],
        matrix[[1, column]],
        matrix[[2, column]],
    )
}

#[inline]
pub(crate) fn array2_set_column3<T: Copy>(
    matrix: &mut Array2<T>,
    column: usize,
    vector: Vector3<T>,
) {
    matrix[[0, column]] = vector[0];
    matrix[[1, column]] = vector[1];
    matrix[[2, column]] = vector[2];
}

#[inline]
pub(crate) fn matrix3_determinant<T: Cfd3dScalar>(matrix: &Matrix3<T>) -> T {
    let a = matrix[(0, 0)];
    let b = matrix[(0, 1)];
    let c = matrix[(0, 2)];
    let d = matrix[(1, 0)];
    let e = matrix[(1, 1)];
    let f = matrix[(1, 2)];
    let g = matrix[(2, 0)];
    let h = matrix[(2, 1)];
    let i = matrix[(2, 2)];

    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

#[inline]
pub(crate) fn matrix3_try_inverse<T: Cfd3dScalar>(matrix: &Matrix3<T>) -> Option<Matrix3<T>> {
    let a = matrix[(0, 0)];
    let b = matrix[(0, 1)];
    let c = matrix[(0, 2)];
    let d = matrix[(1, 0)];
    let e = matrix[(1, 1)];
    let f = matrix[(1, 2)];
    let g = matrix[(2, 0)];
    let h = matrix[(2, 1)];
    let i = matrix[(2, 2)];

    let det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if det == scalar::zero::<T>() {
        return None;
    }

    let inv_det = scalar::one::<T>() / det;
    Some(FixedMatrix::from_rows([
        [
            (e * i - f * h) * inv_det,
            (c * h - b * i) * inv_det,
            (b * f - c * e) * inv_det,
        ],
        [
            (f * g - d * i) * inv_det,
            (a * i - c * g) * inv_det,
            (c * d - a * f) * inv_det,
        ],
        [
            (d * h - e * g) * inv_det,
            (b * g - a * h) * inv_det,
            (a * e - b * d) * inv_det,
        ],
    ]))
}

#[inline]
pub(crate) fn matrix3_scale<T>(matrix: &Matrix3<T>, scale: T) -> Matrix3<T>
where
    T: Copy + core::ops::Mul<Output = T>,
{
    FixedMatrix::from_rows(std::array::from_fn(|row| {
        std::array::from_fn(|col| matrix[(row, col)] * scale)
    }))
}

#[inline]
pub(crate) fn symmetric_part<T>(matrix: &Matrix3<T>) -> Matrix3<T>
where
    T: eunomia::FloatElement + Copy + Default + core::ops::Add<Output = T> + core::ops::Mul<Output = T>,
{
    let half = <T as eunomia::FloatElement>::from_f64(0.5);
    let transpose = matrix.transpose();
    FixedMatrix::from_rows(std::array::from_fn(|row| {
        std::array::from_fn(|col| (matrix[(row, col)] + transpose[(row, col)]) * half)
    }))
}
