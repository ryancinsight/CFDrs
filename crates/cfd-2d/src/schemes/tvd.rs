//! Total Variation Diminishing (TVD) schemes

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Flux limiter for TVD schemes
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum FluxLimiter {
    /// No limiter (unlimited)
    None,
    /// Van Leer limiter
    VanLeer,
    /// Van Albada limiter
    VanAlbada,
    /// Superbee limiter
    Superbee,
    /// MC (Monotonized Central) limiter
    MC,
    /// Minmod limiter
    Minmod,
}

impl FluxLimiter {
    /// Apply flux limiter function
    pub fn apply<T: RealField + Copy + FromPrimitive + Copy>(&self, r: T) -> T {
        match self {
            FluxLimiter::None => T::one(),
            FluxLimiter::VanLeer => {
                if r > T::zero() {
                    let two = T::from_f64(2.0).unwrap_or_else(T::zero);
                    (r + r.abs()) / (T::one() + r.abs())
                } else {
                    T::zero()
                }
            }
            FluxLimiter::VanAlbada => {
                if r > T::zero() {
                    (r * r + r) / (r * r + T::one())
                } else {
                    T::zero()
                }
            }
            FluxLimiter::Superbee => {
                let two = T::from_f64(2.0).unwrap_or_else(T::zero);
                // Superbee limiter: φ(r) = max(0, min(1, 2r), min(2, r))
                T::zero().max(T::one().min(two * r).max(two.min(r)))
            }
            FluxLimiter::MC => {
                let two = T::from_f64(2.0).unwrap_or_else(T::zero);
                T::zero().max(((T::one() + r) / two).min(two.min(two * r)))
            }
            FluxLimiter::Minmod => {
                if r > T::zero() {
                    T::one().min(r)
                } else {
                    T::zero()
                }
            }
        }
    }
}

/// MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws)
pub struct MUSCLScheme<T: RealField + Copy> {
    limiter: FluxLimiter,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> MUSCLScheme<T> {
    /// Create new MUSCL scheme with specified limiter
    #[must_use]
    pub fn new(limiter: FluxLimiter) -> Self {
        Self {
            limiter,
            _phantom: std::marker::PhantomData,
        }
    }
}

/// QUICK (Quadratic Upstream Interpolation for Convective Kinematics)
///
/// Third-order accurate scheme from Leonard (1979) "A stable and accurate convective modelling procedure"
/// based on quadratic upstream interpolation. Combines second-order central differencing with
/// third-order upwind-biased interpolation to minimize numerical diffusion.
pub struct QUICKScheme<T: RealField + Copy> {
    /// Courant number for stability analysis
    courant_max: T,
}

impl<T: RealField + Copy> Default for QUICKScheme<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> QUICKScheme<T> {
    /// Create new QUICK scheme with default Courant limit of 0.8
    #[must_use]
    pub fn new() -> Self {
        Self {
            courant_max: T::from_f64(0.8).unwrap_or_else(T::one),
        }
    }

    /// Compute face value using QUICK interpolation
    ///
    /// For uniform grid with flow from left to right:
    /// φ_face = 6/8 * φ_U + 3/8 * φ_C - 1/8 * φ_UU
    /// where U = upstream, C = central, UU = far upstream
    pub fn interpolate_face(
        &self,
        phi_uu: T, // Far upstream value
        phi_u: T,  // Upstream value
        phi_c: T,  // Central value
        phi_d: T,  // Downstream value
        velocity: T,
    ) -> T {
        if velocity > T::zero() {
            // Flow from left to right
            let six_eighths = T::from_f64(0.75).unwrap_or_else(T::one);
            let three_eighths = T::from_f64(0.375).unwrap_or_else(T::one);
            let one_eighth = T::from_f64(0.125).unwrap_or_else(T::one);

            six_eighths * phi_u + three_eighths * phi_c - one_eighth * phi_uu
        } else {
            // Flow from right to left (use downstream values)
            let six_eighths = T::from_f64(0.75).unwrap_or_else(T::one);
            let three_eighths = T::from_f64(0.375).unwrap_or_else(T::one);
            let one_eighth = T::from_f64(0.125).unwrap_or_else(T::one);

            six_eighths * phi_c + three_eighths * phi_u - one_eighth * phi_d
        }
    }

    /// Check Courant number for stability
    pub fn is_stable(&self, courant: T) -> bool {
        courant <= self.courant_max
    }
}
