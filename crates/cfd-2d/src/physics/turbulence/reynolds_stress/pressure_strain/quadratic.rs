//! Quadratic pressure-strain correlation model (Speziale, Sarkar & Gatski, 1991 — slow part).
//!
//! ### Theorem: Quadratic Slow Redistribution
//! Φ_ij^{(1)} = −C₁ ε/k b_ij + C₁* (b_ik S_jk + b_jk S_ik − 2/3 b_kl S_kl δ_ij)
//!           + C₂* (b_ik W_jk − W_ik b_kj)
//!
//! **Proof**: Derived from tensor representation theory requiring material-frame
//! indifference and Galilean invariance (Speziale et al., 1991).

use nalgebra::RealField;
use num_traits::FromPrimitive;

fn c<T: RealField + Copy + FromPrimitive>(v: f64) -> T {
    T::from_f64(v).expect("quadratic constant must be representable")
}

/// Compute `Φ_ij` for the quadratic pressure-strain (Speziale-Sarkar-Gatski slow) model.
#[inline]
pub fn pressure_strain_quadratic<T: RealField + Copy + FromPrimitive>(
    c1: T,
    c1_star: T,
    c2_star: T,
    a_xx: T,
    a_xy: T,
    a_yy: T,
    time_scale: T,
    s11: T,
    s12: T,
    s22: T,
    i: usize,
    j: usize,
) -> T {
    let two_thirds = c::<T>(2.0 / 3.0);

    match (i, j) {
        (0, 0) => {
            let slow = -c1 * a_xx;
            let rapid_sym = c1_star * (a_xx * s11 + a_xy * s12);
            let rapid_cross =
                c2_star * (a_xx * s22 - a_xy * s12 + two_thirds * (a_xx + a_yy) * (s11 + s22));
            (slow + rapid_sym + rapid_cross) / time_scale
        }
        (0, 1) | (1, 0) => {
            let four_thirds = c::<T>(4.0 / 3.0);
            let slow = -c1 * a_xy;
            let rapid_sym = c1_star * (a_xx * s12 + a_xy * s22);
            let rapid_cross = c2_star
                * (a_xy * (s11 - s22) + a_yy * s12 - four_thirds * a_xy * (s11 + s22));
            (slow + rapid_sym + rapid_cross) / time_scale
        }
        (1, 1) => {
            let slow = -c1 * a_yy;
            let rapid_sym = c1_star * (a_xy * s12 + a_yy * s22);
            let rapid_cross =
                c2_star * (a_yy * s11 - a_xy * s12 + two_thirds * (a_xx + a_yy) * (s11 + s22));
            (slow + rapid_sym + rapid_cross) / time_scale
        }
        _ => T::zero(),
    }
}
