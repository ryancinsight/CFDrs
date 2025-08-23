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
            },
            FluxLimiter::VanAlbada => {
                if r > T::zero() {
                    (r * r + r) / (r * r + T::one())
                } else {
                    T::zero()
                }
            },
            FluxLimiter::Superbee => {
                let two = T::from_f64(2.0).unwrap_or_else(T::zero);
                T::zero().max(T::one().min(two * r)).max(two.min(r))
            },
            FluxLimiter::MC => {
                let two = T::from_f64(2.0).unwrap_or_else(T::zero);
                T::zero().max(((T::one() + r) / two).min(two.min(two * r)))
            },
            FluxLimiter::Minmod => {
                if r > T::zero() {
                    T::one().min(r)
                } else {
                    T::zero()
                }
            },
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
    #[must_use] pub fn new(limiter: FluxLimiter) -> Self {
        Self {
            limiter,
            _phantom: std::marker::PhantomData,
        }
    }
}

/// QUICK (Quadratic Upstream Interpolation for Convective Kinematics)
pub struct QUICKScheme<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for QUICKScheme<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> QUICKScheme<T> {
    /// Create new QUICK scheme
    #[must_use] pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}
