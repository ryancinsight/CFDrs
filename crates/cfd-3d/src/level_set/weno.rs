//! WENO-Z reconstruction helpers for the 3D level-set solver.

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Compute the upwind WENO5-Z spatial derivative `dφ/dx` at point `i`
/// via flux differencing.
///
/// Takes a 7-point stencil `v = [φ_{i-3}, …, φ_{i+3}]`, cell spacing `h`,
/// and the local velocity component `u`.
///
/// # Algorithm
///
/// The derivative is assembled from left-biased (`D⁻`) and right-biased (`D⁺`)
/// flux-difference reconstructions:
///
/// ```text
/// D⁻φ_i = (φ̂⁻_{i+½} − φ̂⁻_{i−½}) / h   (u > 0: upwind from left)
/// D⁺φ_i = (φ̂⁺_{i+½} − φ̂⁺_{i−½}) / h   (u < 0: upwind from right)
/// ```
///
/// where each `φ̂` is a WENO5 reconstruction at a cell face from 5 point values.
pub(super) fn weno5_derivative<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    v: [T; 7],
    h: T,
    u: T,
) -> T {
    let half =
        <T as FromPrimitive>::from_f64(0.5).expect("0.5 is exactly representable in IEEE 754");

    let fl_right = weno5_reconstruct_left([v[1], v[2], v[3], v[4], v[5]]);
    let fl_left = weno5_reconstruct_left([v[0], v[1], v[2], v[3], v[4]]);
    let dm = (fl_right - fl_left) / h;

    let fr_right = weno5_reconstruct_right([v[2], v[3], v[4], v[5], v[6]]);
    let fr_left = weno5_reconstruct_right([v[1], v[2], v[3], v[4], v[5]]);
    let dp = (fr_right - fr_left) / h;

    if u > T::zero() {
        dm
    } else if u < T::zero() {
        dp
    } else {
        (dm + dp) * half
    }
}

/// Left-biased WENO5 reconstruction of `φ` at the right face of the middle cell.
///
/// Given `[φ_{i-2}, φ_{i-1}, φ_i, φ_{i+1}, φ_{i+2}]`, reconstructs
/// `φ̂⁻_{i+½}` using three overlapping 3rd-order sub-stencils.
pub(super) fn weno5_reconstruct_left<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    v: [T; 5],
) -> T {
    // Keep the regularization representable for both f32 and f64 so the
    // fallback path remains well-defined for constant fields.
    let eps = <T as FromPrimitive>::from_f64(1e-12)
        .expect("1e-12 is an IEEE 754 representable f64 constant");
    let one = T::one();
    let two = one + one;
    let three = two + one;
    let five = three + two;
    let six = three + three;
    let seven = three + three + one;
    let eleven = five + six;

    let q0 = (two * v[0] - seven * v[1] + eleven * v[2]) / six;
    let q1 = (-v[1] + five * v[2] + two * v[3]) / six;
    let q2 = (two * v[2] + five * v[3] - v[4]) / six;

    let b0 = smoothness_indicator(v[0], v[1], v[2]);
    let b1 = smoothness_indicator(v[1], v[2], v[3]);
    let b2 = smoothness_indicator(v[2], v[3], v[4]);

    let (w0, w1, w2) = nonlinear_weights(
        b0,
        b1,
        b2,
        <T as FromPrimitive>::from_f64(0.1).expect("0.1 is an IEEE 754 representable f64 constant"),
        <T as FromPrimitive>::from_f64(0.6).expect("0.6 is an IEEE 754 representable f64 constant"),
        <T as FromPrimitive>::from_f64(0.3).expect("0.3 is an IEEE 754 representable f64 constant"),
        eps,
    );

    w0 * q0 + w1 * q1 + w2 * q2
}

/// Right-biased WENO5 reconstruction of `φ` at the left face of the middle cell,
/// obtained by mirroring the left-biased stencil.
pub(super) fn weno5_reconstruct_right<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    v: [T; 5],
) -> T {
    weno5_reconstruct_left([v[4], v[3], v[2], v[1], v[0]])
}

/// WENO5 smoothness indicator for a 3-point sub-stencil (Jiang & Shu 1996).
#[inline]
pub(super) fn smoothness_indicator<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    v0: T,
    v1: T,
    v2: T,
) -> T {
    let thirteen_over_twelve = <T as FromPrimitive>::from_f64(13.0 / 12.0)
        .expect("13/12 is an IEEE 754 representable f64 constant");
    let quarter =
        <T as FromPrimitive>::from_f64(0.25).expect("0.25 is exactly representable in IEEE 754");
    let two = T::one() + T::one();

    let diff1 = v0 - two * v1 + v2;
    let diff2 = v0 - v2;
    thirteen_over_twelve * diff1 * diff1 + quarter * diff2 * diff2
}

/// Compute normalized WENO5-Z nonlinear weights from smoothness indicators.
#[inline]
pub(super) fn nonlinear_weights<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    b0: T,
    b1: T,
    b2: T,
    d0: T,
    d1: T,
    d2: T,
    eps: T,
) -> (T, T, T) {
    let tau5 = global_smoothness_indicator(b0, b1, b2);
    let one = T::one();
    let a0 = d0 * (one + num_traits::Float::powi(tau5 / (eps + b0), 2));
    let a1 = d1 * (one + num_traits::Float::powi(tau5 / (eps + b1), 2));
    let a2 = d2 * (one + num_traits::Float::powi(tau5 / (eps + b2), 2));
    let sum = a0 + a1 + a2;
    if sum < eps {
        return (
            d0 / (d0 + d1 + d2),
            d1 / (d0 + d1 + d2),
            d2 / (d0 + d1 + d2),
        );
    }
    (a0 / sum, a1 / sum, a2 / sum)
}

#[inline]
fn global_smoothness_indicator<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    b0: T,
    _b1: T,
    b2: T,
) -> T {
    num_traits::Float::abs(b0 - b2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_smoothness_indicator_non_negative(
            v0 in -10.0..10.0f64,
            v1 in -10.0..10.0f64,
            v2 in -10.0..10.0f64,
        ) {
            let beta = smoothness_indicator(v0, v1, v2);
            assert!(beta >= 0.0);
        }

        #[test]
        fn test_nonlinear_weights_sum_to_one(
            b0 in 0.0..10.0f64,
            b1 in 0.0..10.0f64,
            b2 in 0.0..10.0f64,
        ) {
            let d0 = 0.1;
            let d1 = 0.6;
            let d2 = 0.3;
            let eps = 1e-12;

            let (w0, w1, w2) = nonlinear_weights(b0, b1, b2, d0, d1, d2, eps);
            let sum = w0 + w1 + w2;

            assert!((sum - 1.0).abs() < 1e-14);
            assert!(w0 >= 0.0);
            assert!(w1 >= 0.0);
            assert!(w2 >= 0.0);
        }
    }

    #[test]
    fn weno5_constant_derivative_is_zero() {
        let c = 3.7_f64;
        let h = 0.1_f64;
        let v = [c; 7];

        let d_pos = weno5_derivative(v, h, 1.0);
        assert!((d_pos).abs() < 1e-14, "d_pos = {d_pos}");

        let d_neg = weno5_derivative(v, h, -1.0);
        assert!((d_neg).abs() < 1e-14, "d_neg = {d_neg}");

        let d_zero = weno5_derivative(v, h, 0.0);
        assert!((d_zero).abs() < 1e-14, "d_zero = {d_zero}");
    }

    #[test]
    fn weno5_linear_exactness() {
        let h = 0.1_f64;
        let slope = 2.5_f64;
        let phi_center = 1.0_f64;
        let v: [f64; 7] = std::array::from_fn(|k| phi_center + slope * ((k as f64 - 3.0) * h));

        let d_pos = weno5_derivative(v, h, 1.0);
        assert!(
            (d_pos - slope).abs() < 1e-12,
            "linear d_pos = {d_pos}, expected {slope}"
        );

        let d_neg = weno5_derivative(v, h, -1.0);
        assert!(
            (d_neg - slope).abs() < 1e-12,
            "linear d_neg = {d_neg}, expected {slope}"
        );
    }

    #[test]
    fn weno5_constant_derivative_is_zero_f32() {
        let c = 1.25_f32;
        let h = 0.125_f32;
        let v = [c; 7];

        let d_pos = weno5_derivative(v, h, 1.0);
        let d_neg = weno5_derivative(v, h, -1.0);
        let d_zero = weno5_derivative(v, h, 0.0);

        assert!(d_pos.abs() < 1e-6, "d_pos = {d_pos}");
        assert!(d_neg.abs() < 1e-6, "d_neg = {d_neg}");
        assert!(d_zero.abs() < 1e-6, "d_zero = {d_zero}");
    }
}
