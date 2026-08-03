use aequitas::systems::si::quantities::{Pressure, Time};

use super::{cavitation_amplified_hi, giersiepen_hi};

/// Composite haemolysis exposure combining shear stress, duration, and
/// cavitation potential.
///
/// Use [`HemolysisExposure::compute_index`] to obtain the amplified HI for
/// a single exposure event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HemolysisExposure {
    /// Wall shear stress at the exposure site \[Pa].
    pub shear: Pressure,
    /// Duration of the exposure event \[s].
    pub duration: Time,
    /// Local cavitation potential ∈ [0, 1]; 0 for non-venturi regions.
    pub cavitation_potential: f64,
}

impl HemolysisExposure {
    /// Create a new exposure event.
    #[must_use]
    pub fn new(shear: Pressure, duration: Time, cavitation_potential: f64) -> Self {
        Self {
            shear,
            duration,
            cavitation_potential,
        }
    }

    /// Create an exposure event with no cavitation contribution.
    #[must_use]
    pub fn shear_only(shear: Pressure, duration: Time) -> Self {
        Self::new(shear, duration, 0.0)
    }

    /// Compute the amplified haemolysis index for this exposure event.
    ///
    /// Combines [`giersiepen_hi`] with [`cavitation_amplified_hi`].
    #[must_use]
    pub fn compute_index(&self) -> f64 {
        cavitation_amplified_hi(
            giersiepen_hi(self.shear, self.duration),
            self.cavitation_potential,
        )
    }
}

/// Empirical collapse amplification coefficient α ≈ 3.0.
pub const RAYLEIGH_ALPHA: f64 = 3.0;

/// Standard atmospheric pressure \[Pa\] used as reference for amplification.
pub const P_REF_ATMOSPHERIC: f64 = 101_325.0;

/// Rayleigh inertial collapse time for a cavitation bubble \[s\].
///
/// ## Theorem — Rayleigh Collapse (Rayleigh 1917)
///
/// For a spherical bubble collapsing from radius R_max under external
/// pressure p_inf (neglecting gas content and surface tension):
///
/// ```text
/// t_collapse = 0.915 · R_max · √(ρ / p_inf)
/// ```
///
/// **Reference**: Rayleigh, Lord (1917). "On the pressure developed in a
/// liquid during the collapse of a spherical cavity",
/// *Phil. Mag.* 34:94-98.
///
/// # Arguments
///
/// * `r_max` — maximum bubble radius \[m\]
/// * `rho` — liquid density \[kg/m³\]
/// * `p_inf` — external (driving) pressure \[Pa\]
///
/// # Example
///
/// ```
/// use cfd_1d::rayleigh_collapse_time;
/// let t = rayleigh_collapse_time(50e-6, 1060.0, 101325.0);
/// assert!(t > 0.0);
/// ```
#[inline]
#[must_use]
pub fn rayleigh_collapse_time(r_max: f64, rho: f64, p_inf: f64) -> f64 {
    if r_max <= 0.0 || rho <= 0.0 || p_inf <= 0.0 {
        return 0.0;
    }

    0.915 * r_max * (rho / p_inf).sqrt()
}

/// Micro-jet velocity from Rayleigh bubble collapse \[m/s\].
///
/// The collapse generates a micro-jet with velocity:
///
/// ```text
/// v_jet ≈ √(2 · p_inf / ρ)
/// ```
///
/// # Arguments
///
/// * `p_inf` — external pressure \[Pa\]
/// * `rho` — liquid density \[kg/m³\]
///
/// # Example
///
/// ```
/// use cfd_1d::collapse_jet_velocity;
/// let v = collapse_jet_velocity(101325.0, 1060.0);
/// assert!(v > 10.0, "Jet velocity at atmospheric pressure should exceed 10 m/s");
/// ```
#[inline]
#[must_use]
pub fn collapse_jet_velocity(p_inf: f64, rho: f64) -> f64 {
    if p_inf <= 0.0 || rho <= 0.0 {
        return 0.0;
    }

    (2.0 * p_inf / rho).sqrt()
}

/// Hemolysis amplification factor from cavitation bubble collapse.
///
/// The amplification scales with collapse intensity:
///
/// ```text
/// A_collapse = 1 + α · (R_max / R_0)² · (p_inf / p_ref)
/// ```
///
/// where α ≈ 3.0 (empirical), R_0 is equilibrium radius,
/// p_ref = 101 325 Pa (1 atm).
///
/// # Arguments
///
/// * `r_max` — maximum bubble radius \[m\]
/// * `r_0` — equilibrium bubble radius \[m\]
/// * `p_inf` — external pressure \[Pa\]
///
/// # Example
///
/// ```
/// use cfd_1d::cavitation_hemolysis_amplification;
/// let a = cavitation_hemolysis_amplification(100e-6, 50e-6, 101325.0);
/// assert!(a > 1.0, "Amplification should exceed unity");
/// ```
#[inline]
#[must_use]
pub fn cavitation_hemolysis_amplification(r_max: f64, r_0: f64, p_inf: f64) -> f64 {
    if r_max <= 0.0 || r_0 <= 0.0 || p_inf <= 0.0 {
        return 1.0;
    }

    1.0 + RAYLEIGH_ALPHA * (r_max / r_0).powi(2) * (p_inf / P_REF_ATMOSPHERIC)
}
