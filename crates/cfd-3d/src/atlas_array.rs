//! Private Leto-owned array helpers for Apollo FFT and NUFFT call sites.
//!
//! CFDrs keeps dense spectral storage in Leto arrays and calls Apollo's typed
//! Leto transform surfaces directly.

use apollo_fft::{
    fft_1d_array_typed as apollo_fft_1d_array, fft_3d_array_typed as apollo_fft_3d_array,
    ifft_1d_array_typed as apollo_ifft_1d_array, ifft_3d_array_typed as apollo_ifft_3d_array,
    Complex64,
};
use apollo_nufft::nufft_type2_3d as apollo_nufft_type2_3d;
use apollo_nufft::{nufft_type1_3d as apollo_nufft_type1_3d, UniformGrid3D};
use cfd_core::error::{Error, Result};
use leto::{Array, Array1, Array3, VecStorage};

pub(crate) fn value<T: Copy, const N: usize>(
    array: &Array<T, VecStorage<T>, N>,
    index: [usize; N],
) -> T {
    *array
        .get(index)
        .expect("invariant: array index is bounded by its validated shape")
}

pub(crate) fn set<T, const N: usize>(
    array: &mut Array<T, VecStorage<T>, N>,
    index: [usize; N],
    value: T,
) {
    *array
        .get_mut(index)
        .expect("invariant: array index is bounded by its validated shape") = value;
}

pub(crate) fn values<T: Clone, const N: usize>(array: &Array<T, VecStorage<T>, N>) -> Vec<T> {
    array.clone().into_vec()
}

pub(crate) fn map<T: Copy, U, const N: usize>(
    array: &Array<T, VecStorage<T>, N>,
    mut f: impl FnMut(T) -> U,
) -> Array<U, VecStorage<U>, N> {
    Array::from_shape_fn(array.shape(), |index| f(value(array, index)))
}

pub(crate) fn fft_1d_array(field: &Array1<f64>) -> Result<Array1<Complex64>> {
    Ok(apollo_fft_1d_array(field))
}

pub(crate) fn ifft_1d_array(field_hat: &Array1<Complex64>) -> Result<Array1<f64>> {
    Ok(apollo_ifft_1d_array(field_hat))
}

pub(crate) fn fft_3d_array(field: &Array3<f64>) -> Result<Array3<Complex64>> {
    Ok(apollo_fft_3d_array(field))
}

pub(crate) fn ifft_3d_array(field_hat: &Array3<Complex64>) -> Result<Array3<f64>> {
    Ok(apollo_ifft_3d_array(field_hat))
}

pub(crate) fn nufft_type1_3d(
    positions: &[(f64, f64, f64)],
    values: &[Complex64],
    grid: UniformGrid3D,
) -> Result<Array3<Complex64>> {
    if positions.len() != values.len() {
        return Err(Error::DimensionMismatch {
            expected: positions.len(),
            actual: values.len(),
        });
    }
    Ok(apollo_nufft_type1_3d(positions, values, grid))
}

pub(crate) fn nufft_type2_3d(
    positions: &[(f64, f64, f64)],
    modes: &Array3<Complex64>,
    grid: UniformGrid3D,
) -> Result<Vec<Complex64>> {
    if modes.shape() != [grid.nx, grid.ny, grid.nz] {
        return Err(Error::DimensionMismatch {
            expected: grid.nx * grid.ny * grid.nz,
            actual: modes.size(),
        });
    }
    Ok(apollo_nufft_type2_3d(positions, modes, grid))
}
