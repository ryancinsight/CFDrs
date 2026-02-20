//! Streamline curvature correction for the pressure-strain correlation.
//!
//! ### Theorem: Suga-Craft Curvature Correction (2003)
//!
//! Convex curvature stabilises turbulence (reduces Reynolds stress anisotropy),
//! while concave curvature destabilises it. The correction reads:
//!
//! ```text
//! Φ_ij^curv = C_curv (k/ε) K A_ij(S)
//! ```
//!
//! where K is the curvature parameter and C_curv > 0 for convex, < 0 for concave
//! (Suga & Craft, 2003).

use nalgebra::RealField;
use num_traits::FromPrimitive;

fn c<T: RealField + Copy + FromPrimitive>(v: f64) -> T {
    T::from_f64(v).expect("curvature constant must be representable")
}

/// Compute the Suga-Craft (2003) curvature correction term.
///
/// # Arguments
/// * `a_xx, a_xy, a_yy` — anisotropy tensor components
/// * `time_scale` — turbulence time scale k/ε
/// * `s11, s12, s22` — mean strain rate tensor components
/// * `w12, w21` — rotation rate tensor components
/// * `i, j` — tensor indices
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn curvature_correction<T: RealField + Copy + FromPrimitive>(
    a_xx: T,
    a_xy: T,
    a_yy: T,
    time_scale: T,
    s11: T,
    s12: T,
    s22: T,
    w12: T,
    i: usize,
    j: usize,
) -> T {
    let curvature_param = calculate_curvature_parameter(s11, s12, s22, w12);

    let curvature_strength = curvature_param.abs();
    if curvature_strength <= c::<T>(1e-6) {
        return T::zero();
    }

    let c_convex = c::<T>(0.15);
    let c_concave = c::<T>(-0.1);
    let curvature_factor = if curvature_param >= T::zero() { c_convex } else { c_concave };

    let correction = match (i, j) {
        (0, 0) => curvature_factor * curvature_strength
            * (a_xx * s11 + c::<T>(2.0) * a_xy * s12 - c::<T>(2.0 / 3.0) * (a_xx + a_yy) * s11),
        (0, 1) | (1, 0) => curvature_factor * curvature_strength
            * (a_xx * s12 + a_xy * s22 + a_yy * s12 - c::<T>(2.0 / 3.0) * (a_xx + a_yy) * s12),
        (1, 1) => curvature_factor * curvature_strength
            * (a_xy * s12 + a_yy * s22 - c::<T>(2.0 / 3.0) * (a_xx + a_yy) * s22),
        _ => T::zero(),
    };

    correction / time_scale
}

/// Compute dimensionless curvature parameter K.
///
/// K = (|S|² − |W|²) / (|S|² + |W|²)
/// - K > 0: convex (strain-dominated)
/// - K < 0: concave (rotation-dominated)
#[inline]
pub fn calculate_curvature_parameter<T: RealField + Copy + FromPrimitive>(
    s11: T,
    s12: T,
    s22: T,
    w12: T,
) -> T {
    let strain_sq = s11 * s11 + c::<T>(2.0) * s12 * s12 + s22 * s22;
    let rotation_sq = c::<T>(2.0) * w12 * w12;
    let denom = strain_sq + rotation_sq;
    if denom > c::<T>(1e-12) {
        (strain_sq - rotation_sq) / denom
    } else {
        T::zero()
    }
}
