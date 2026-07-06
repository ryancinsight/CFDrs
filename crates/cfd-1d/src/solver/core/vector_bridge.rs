use eunomia::NumericElement;
use leto::{Array1, Storage};
use nalgebra::DVector;

pub(super) fn array_from_dvector<T: nalgebra::Scalar + Copy>(vector: &DVector<T>) -> Array1<T> {
    Array1::from_shape_vec([vector.len()], vector.as_slice().to_vec())
        .expect("invariant: DVector length matches produced Leto vector shape")
}

pub(super) fn dvector_from_array<T: nalgebra::Scalar + Clone>(array: &Array1<T>) -> DVector<T> {
    DVector::from_vec(array.storage().as_slice().to_vec())
}

pub(super) fn dvector_from_owned_array<T: nalgebra::Scalar>(array: Array1<T>) -> DVector<T> {
    DVector::from_vec(array.into_vec())
}

pub(super) fn array_l2_norm<T: NumericElement + Copy>(values: &Array1<T>) -> T {
    let norm_sq = values
        .storage()
        .as_slice()
        .iter()
        .copied()
        .fold(<T as NumericElement>::ZERO, |acc, value| {
            acc + value * value
        });
    <T as NumericElement>::sqrt(norm_sq)
}

pub(super) fn copy_array<T: Copy>(source: &Array1<T>, target: &mut Array1<T>) {
    let source = source.storage().as_slice();
    let target = target
        .as_slice_mut()
        .expect("invariant: solver workspace arrays are dense C-order");
    assert_eq!(
        source.len(),
        target.len(),
        "invariant: solver workspace copy preserves vector length"
    );
    target.copy_from_slice(source);
}
