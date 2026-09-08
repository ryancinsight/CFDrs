use eunomia::NumericElement;
use leto::{Array1, Storage};

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
