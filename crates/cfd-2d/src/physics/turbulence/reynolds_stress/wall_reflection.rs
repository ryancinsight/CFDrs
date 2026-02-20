//! Wall-reflection correction for the pressure-strain model.
//!
//! ### Theorem: Gibson-Launder Wall-Reflection Hypothesis (1978)
//!
//! Near a solid wall, image sources in the pressure-Poisson equation create an
//! additional redistribution contribution:
//!
//! ```text
//! Φ_ij^wall = −C_w (k/ε) [a_in a_jn + a_ip a_jp − (2/3) δ_ij (a_nn + a_pp)] f(y+)
//! ```
//!
//! where `n` is the wall-normal direction and `f(y+)` is a damping function
//! that vanishes in the log-layer (Gibson & Launder, 1978).

use nalgebra::RealField;
use num_traits::FromPrimitive;

fn c<T: RealField + Copy + FromPrimitive>(v: f64) -> T {
    T::from_f64(v).expect("wall-reflection constant must be representable")
}

/// Compute the Gibson-Launder (1978) wall-reflection correction to `Φ_ij`.
///
/// # Arguments
/// * `a_nn`, `a_pp`, `a_ps` — anisotropy components in wall-normal / parallel directions
/// * `k`, `epsilon`, `y` — turbulence variables + wall distance `y`
/// * `nu` — kinematic viscosity (for y⁺ calculation)
/// * `wall_normal_index` — grid axis perpendicular to the wall (0 = x, 1 = y)
/// * `i`, `j` — tensor indices
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn wall_reflection_correction<T: RealField + Copy + FromPrimitive>(
    a_nn: T,
    a_pp: T,
    a_ps: T,
    k: T,
    epsilon: T,
    y: T,
    nu: T,
    wall_normal_index: usize,
    i: usize,
    j: usize,
) -> T {
    if y <= T::zero() || epsilon <= T::zero() || k <= T::zero() {
        return T::zero();
    }

    let time_scale = k / epsilon;

    // Wall-reflection constants (Gibson & Launder, 1978)
    let c_w1 = c::<T>(0.5);
    let c_w2 = c::<T>(0.3);

    // y⁺ damping function: (1 − exp(−y⁺/A))/ y⁺  (Van Driest form)
    let u_tau = (epsilon * c::<T>(0.09).sqrt() * k).sqrt(); // u_τ ≈ C_μ^{1/4} k^{1/2}
    let y_plus = if nu > T::zero() { y * u_tau / nu } else { T::zero() };

    let damping_factor = if y_plus > T::zero() {
        let a_damping = c::<T>(25.0);
        (T::one() - (-y_plus / a_damping).exp()) / y_plus
    } else {
        T::zero()
    };

    let reflection_term = match (i, j) {
        (0, 0) => {
            let delta = if wall_normal_index == 0 { T::one() } else { T::zero() };
            -c_w1 * damping_factor * (a_nn * a_nn + a_pp * a_pp - c::<T>(2.0 / 3.0) * (a_nn + a_pp) * delta)
        }
        (0, 1) | (1, 0) => -c_w2 * damping_factor * (a_nn * a_ps + a_pp * a_ps),
        (1, 1) => {
            let delta = if wall_normal_index == 1 { T::one() } else { T::zero() };
            -c_w1 * damping_factor * (a_nn * a_nn + a_pp * a_pp - c::<T>(2.0 / 3.0) * (a_nn + a_pp) * delta)
        }
        _ => T::zero(),
    };

    reflection_term / time_scale
}
