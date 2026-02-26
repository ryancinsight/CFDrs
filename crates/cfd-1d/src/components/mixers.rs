//! Mixer components for microfluidic networks
//!
//! # Micromixer Resistance Theorem
//!
//! A passive micromixer is modelled as a straight channel section
//! (Hagen-Poiseuille base resistance) augmented by minor head-loss due to
//! flow-turning geometry:
//!
//! ```text
//! R_mixer = R_straight + R_minor
//!
//! R_straight = 128 μ L / (π D⁴)           [Hagen-Poiseuille, Pa·s/m³]
//! R_minor    = K_loss · μ / (A · Dh)      [linearised Idelchik minor loss]
//! ```
//!
//! where:
//! - `μ`       = dynamic viscosity [Pa·s]
//! - `L`       = total channel length [m]
//! - `D` / `Dh`= (hydraulic) diameter [m]
//! - `A`       = cross-sectional area `= π/4 · D²` [m²]
//! - `K_loss`  = dimensionless minor-loss coefficient per bend/geometry
//!
//! **Minor-loss coefficients by mixer type** (laminar regime):
//!
//! | Mixer type  | `K_loss`                    | Source                          |
//! |-------------|-----------------------------|---------------------------------|
//! | T-junction  | 1.5                         | Idelchik (1986) §7.23           |
//! | Y-junction  | 0.9                         | Idelchik (1986) §7.27           |
//! | Serpentine  | 0.3 × `n_bends`             | White (2003) Ch. 8              |
//! | Herringbone | 0.5 × `n_bends`             | Stroock et al. (2002) *Science* |
//!
//! **Linearisation**: The Idelchik minor-loss formula is originally expressed in
//! kinetic pressure `ΔP = K ½ρv²`. In the low-Re laminar regime of microfluidic
//! devices (`Re < 100`) the dynamic pressure term is negligible compared with the
//! viscous term and the minor loss can be expressed in resistance form as above
//! without degrading accuracy.
//!
//! # Mixing Efficiency
//!
//! `efficiency ∈ [0, 1]` is a user-set calibration parameter representing the
//! degree of mixing achieved (`1.0` = fully mixed). It has no effect on the
//! hydraulic resistance calculation; it is exposed for downstream composition
//! solvers that need to know the homogenisation degree.
//!
//! # Invariants
//!
//! - `hydraulic_diameter > 0`
//! - `length > 0`
//! - `n_bends ≥ 1` for serpentine/herringbone
//! - `efficiency ∈ [0, 1]` (clamped on set)

use super::Component;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::Fluid;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Mixer type enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MixerType {
    /// T-junction mixer (K_loss = 1.5, Idelchik §7.23)
    TJunction,
    /// Y-junction mixer (K_loss = 0.9, Idelchik §7.27)
    YJunction,
    /// Serpentine mixer (K_loss = 0.3 × n_bends, White 2003)
    Serpentine,
    /// Herringbone (staggered-groove) mixer (K_loss = 0.5 × n_bends, Stroock 2002)
    Herringbone,
}

impl MixerType {
    /// Return the dimensionless Idelchik minor-loss coefficient `K_loss`.
    ///
    /// For multi-bend geometries (`Serpentine`, `Herringbone`) the coefficient
    /// is multiplied externally by `n_bends` before use.
    #[must_use]
    pub fn base_k_loss(self) -> f64 {
        match self {
            MixerType::TJunction => 1.5,
            MixerType::YJunction => 0.9,
            MixerType::Serpentine => 0.3,
            MixerType::Herringbone => 0.5,
        }
    }

    /// Whether this geometry involves repeated bends.
    #[must_use]
    pub fn is_multi_bend(self) -> bool {
        matches!(self, MixerType::Serpentine | MixerType::Herringbone)
    }
}

/// Passive micromixer component with first-principles hydraulic resistance.
///
/// The total resistance is computed from geometry and fluid viscosity
/// via the Hagen-Poiseuille base term plus Idelchik minor losses.
/// See module documentation for the governing equations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Micromixer<T: RealField + Copy> {
    /// Mixer geometry type
    pub mixer_type: MixerType,
    /// Hydraulic diameter of the mixer channel [m] (must be > 0)
    pub hydraulic_diameter: T,
    /// Total channel path length [m] (must be > 0)
    pub length: T,
    /// Number of bends (used for Serpentine and Herringbone mixers; ≥ 1)
    pub n_bends: usize,
    /// Mixing efficiency in [0, 1] (calibration parameter, no effect on resistance)
    pub efficiency: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + Copy + FromPrimitive> Micromixer<T> {
    /// Create a new micromixer with validated geometry.
    ///
    /// # Errors
    ///
    /// Returns `Error::InvalidConfiguration` if `hydraulic_diameter ≤ 0` or `length ≤ 0`.
    pub fn new(
        mixer_type: MixerType,
        hydraulic_diameter: T,
        length: T,
        n_bends: usize,
    ) -> Result<Self> {
        if hydraulic_diameter <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Micromixer hydraulic_diameter must be > 0".into(),
            ));
        }
        if length <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Micromixer length must be > 0".into(),
            ));
        }
        let efficiency = T::from_f64(0.95).unwrap_or_else(T::one);
        Ok(Self {
            mixer_type,
            hydraulic_diameter,
            length,
            n_bends: n_bends.max(1),
            efficiency,
            parameters: HashMap::new(),
        })
    }

    /// Effective dimensionless Idelchik minor-loss coefficient for this mixer.
    ///
    /// ```text
    /// K_eff = K_base × n_bends   (if multi-bend)
    /// K_eff = K_base              (otherwise)
    /// ```
    #[must_use]
    pub fn effective_k_loss(&self) -> T {
        let k_base = T::from_f64(self.mixer_type.base_k_loss()).unwrap_or_else(T::zero);
        if self.mixer_type.is_multi_bend() {
            let n = T::from_usize(self.n_bends).unwrap_or_else(T::one);
            k_base * n
        } else {
            k_base
        }
    }

    /// Cross-sectional area from hydraulic diameter (circular equivalent).
    #[must_use]
    fn channel_area(&self) -> T {
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::one);
        let r = self.hydraulic_diameter / T::from_f64(2.0).unwrap_or_else(T::one);
        pi * r * r
    }
}

impl<T: RealField + Copy + FromPrimitive> Component<T> for Micromixer<T> {
    /// Compute the hydraulic resistance of the micromixer.
    ///
    /// ## Theorem: Micromixer Resistance
    ///
    /// ```text
    /// R = R_straight + R_minor
    ///
    /// R_straight = 128 μ L / (π D⁴)
    /// R_minor    = K_eff · μ / (A · Dh)
    /// ```
    ///
    /// **References**:
    /// - Idelchik, I. E. (1986). *Handbook of Hydraulic Resistance*. §7.23, §7.27.
    /// - White, F. M. (2003). *Fluid Mechanics* (5th ed.). Ch. 8.
    /// - Stroock, A. D. et al. (2002). Chaotic mixer for microchannels.
    ///   *Science*, 295(5555), 647–651.
    fn resistance(&self, fluid: &Fluid<T>) -> T {
        let mu = fluid.viscosity;

        let d = self.hydraulic_diameter;
        let l = self.length;
        let a = self.channel_area();

        let c128 = T::from_f64(128.0).unwrap_or_else(T::zero);
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::zero);

        // Guard against degenerate geometry
        let d4 = d * d * d * d;
        if d4 <= T::zero() || a <= T::zero() || d <= T::zero() {
            return T::from_f64(1e15).unwrap_or_else(T::zero);
        }

        // R_straight = 128 μ L / (π D⁴)
        let r_straight = c128 * mu * l / (pi * d4);

        // R_minor = K_eff · μ / (A · Dh)
        let k_eff = self.effective_k_loss();
        let r_minor = k_eff * mu / (a * d);

        r_straight + r_minor
    }

    fn component_type(&self) -> &'static str {
        "Micromixer"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "efficiency" => {
                self.efficiency = if value < T::zero() {
                    T::zero()
                } else if value > T::one() {
                    T::one()
                } else {
                    value
                };
            }
            "hydraulic_diameter" => {
                if value <= T::zero() {
                    return Err(Error::InvalidConfiguration(
                        "hydraulic_diameter must be > 0".into(),
                    ));
                }
                self.hydraulic_diameter = value;
            }
            "length" => {
                if value <= T::zero() {
                    return Err(Error::InvalidConfiguration("length must be > 0".into()));
                }
                self.length = value;
            }
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }
}
