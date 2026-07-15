//! Optimal bifurcation geometry computation.
//!
//! Calculates optimal diameters, angles, and flow distribution for a
//! bifurcation following Murray's Law principles.

use super::law::MurraysLaw;
use crate::scalar::Cfd1dScalar;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Complete optimal bifurcation geometry
///
/// Calculates optimal diameters, angles, and flow distribution for a
/// bifurcation following Murray's Law principles.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct OptimalBifurcation<T: Cfd1dScalar + Copy> {
    /// Parent vessel diameter \[m]
    pub parent_diameter: T,
    /// Major daughter diameter \[m]
    pub daughter1_diameter: T,
    /// Minor daughter diameter \[m]
    pub daughter2_diameter: T,
    /// Half-angle of major daughter branch from parent axis \[rad]
    pub angle1: T,
    /// Half-angle of minor daughter branch from parent axis \[rad]
    pub angle2: T,
    /// Flow rate in parent [m³/s]
    pub parent_flow: T,
    /// Flow rate in daughter 1 [m³/s]
    pub daughter1_flow: T,
    /// Flow rate in daughter 2 [m³/s]
    pub daughter2_flow: T,
}

impl<T: Cfd1dScalar + FloatElement + Copy> OptimalBifurcation<T> {
    /// Create symmetric bifurcation following Murray's Law
    ///
    /// # Theorem: Optimal Symmetric Branching Angle
    ///
    /// For a symmetric bifurcation (r₁ = r₂ = r₀ / 2^(1/3)), the optimal
    /// half-angle θ that minimises total power dissipation satisfies:
    ///
    /// ```text
    /// cos(θ) = (r₀⁴ + r₁⁴ − r₂⁴) / (2·r₀²·r₁²)
    /// ```
    ///
    /// Since r₁ = r₂ this simplifies to cos(θ) = r₀² / (2·r₁²) = 2^(−1/3).
    ///
    /// **Proof sketch**: substitute r₁ = r₂ into the general angle formula,
    /// cancelling the r₂⁴ terms, yielding cos(θ) = 2^(2/3) / 2 = 2^(−1/3).
    pub fn symmetric(parent_diameter: T, parent_flow: T) -> Self {
        let murray = MurraysLaw::<T>::new();
        let d_daughter = murray.symmetric_daughter_diameter(parent_diameter);

        let two = T::one() + T::one();
        let daughter_flow = parent_flow / two;

        // cos(θ) = 2^(−1/3)
        let one_third = T::one() / (T::one() + T::one() + T::one());
        let cos_theta = <T as FloatElement>::powf(two, -one_third);
        let angle = <T as FloatElement>::acos(cos_theta);

        Self {
            parent_diameter,
            daughter1_diameter: d_daughter,
            daughter2_diameter: d_daughter,
            angle1: angle,
            angle2: angle,
            parent_flow,
            daughter1_flow: daughter_flow,
            daughter2_flow: daughter_flow,
        }
    }

    /// Create asymmetric bifurcation with specified flow split ratio
    ///
    /// # Theorem: Optimal Asymmetric Branching Angle (Zamir 1978)
    ///
    /// For an asymmetric bifurcation with daughter radii r₁, r₂, the optimal
    /// half-angles θ₁, θ₂ that minimise total pumping power satisfy:
    ///
    /// ```text
    /// cos(θ₁) = (r₀⁴ + r₁⁴ − r₂⁴) / (2·r₀²·r₁²)
    /// cos(θ₂) = (r₀⁴ + r₂⁴ − r₁⁴) / (2·r₀²·r₂²)
    /// ```
    pub fn asymmetric(parent_diameter: T, parent_flow: T, flow_ratio: T) -> Self {
        let one = T::one();

        let q1 = parent_flow * flow_ratio;
        let q2 = parent_flow * (one - flow_ratio);

        // D₁ = D₀ · (Q₁/Q₀)^(1/3), D₂ = D₀ · (Q₂/Q₀)^(1/3)
        let one_third = one / (T::one() + T::one() + T::one());
        let d1 = parent_diameter * <T as FloatElement>::powf(flow_ratio, one_third);
        let d2 = parent_diameter * <T as FloatElement>::powf(one - flow_ratio, one_third);

        // Zamir (1978) branching angle formula
        let d0_sq = parent_diameter * parent_diameter;
        let d1_sq = d1 * d1;
        let d2_sq = d2 * d2;
        let d0_4 = d0_sq * d0_sq;
        let d1_4 = d1_sq * d1_sq;
        let d2_4 = d2_sq * d2_sq;
        let two = T::one() + T::one();

        let cos_theta1 = <T as eunomia::RealField>::clamp(
            (d0_4 + d1_4 - d2_4) / (two * d0_sq * d1_sq),
            -one,
            one,
        );
        let cos_theta2 = <T as eunomia::RealField>::clamp(
            (d0_4 + d2_4 - d1_4) / (two * d0_sq * d2_sq),
            -one,
            one,
        );

        Self {
            parent_diameter,
            daughter1_diameter: d1,
            daughter2_diameter: d2,
            angle1: <T as FloatElement>::acos(cos_theta1),
            angle2: <T as FloatElement>::acos(cos_theta2),
            parent_flow,
            daughter1_flow: q1,
            daughter2_flow: q2,
        }
    }

    /// Validate that this bifurcation satisfies Murray's Law
    pub fn is_murray_compliant(&self, tolerance: T) -> bool {
        let murray = MurraysLaw::<T>::new();
        murray.is_valid(
            self.parent_diameter,
            self.daughter1_diameter,
            self.daughter2_diameter,
            tolerance,
        )
    }

    /// Calculate Murray's Law deviation
    pub fn murray_deviation(&self) -> T {
        let murray = MurraysLaw::<T>::new();
        murray.deviation(
            self.parent_diameter,
            self.daughter1_diameter,
            self.daughter2_diameter,
        )
    }

    /// Calculate total daughter area to parent area ratio
    pub fn area_ratio(&self) -> T {
        let a0 = self.parent_diameter * self.parent_diameter;
        let a1 = self.daughter1_diameter * self.daughter1_diameter;
        let a2 = self.daughter2_diameter * self.daughter2_diameter;
        (a1 + a2) / a0
    }

    /// Check mass conservation: Q₀ = Q₁ + Q₂
    pub fn mass_conservation_error(&self) -> T {
        let sum = self.daughter1_flow + self.daughter2_flow;
        <T as NumericElement>::abs(self.parent_flow - sum) / self.parent_flow
    }

    /// Calculate pressure drop through major daughter branch
    ///
    /// Uses Poiseuille resistance: ΔP = 8μLQ/(πR⁴)
    pub fn pressure_drop_daughter1(&self, viscosity: T, length: T) -> T {
        let pi = T::pi();
        let eight = <T as FloatElement>::from_f64(8.0);
        let r = self.daughter1_diameter / (T::one() + T::one());
        eight * viscosity * length * self.daughter1_flow / (pi * <T as FloatElement>::powi(r, 4))
    }
}
