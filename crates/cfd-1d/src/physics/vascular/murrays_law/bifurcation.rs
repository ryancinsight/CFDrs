//! Optimal bifurcation geometry computation.
//!
//! Calculates optimal diameters, angles, and flow distribution for a
//! bifurcation following Murray's Law principles.

use super::law::MurraysLaw;
use aequitas::systems::si::quantities::{
    Angle, Dimensionless, DynamicViscosity, Length, Pressure, VolumetricFlowRate,
};
use cfd_core::CfdScalar;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Complete optimal bifurcation geometry
///
/// Calculates optimal diameters, angles, and flow distribution for a
/// bifurcation following Murray's Law principles.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct OptimalBifurcation<T: CfdScalar + Copy> {
    /// Parent vessel diameter \[m]
    pub parent_diameter: Length<T>,
    /// Major daughter diameter \[m]
    pub daughter1_diameter: Length<T>,
    /// Minor daughter diameter \[m]
    pub daughter2_diameter: Length<T>,
    /// Half-angle of major daughter branch from parent axis \[rad]
    pub angle1: Angle<T>,
    /// Half-angle of minor daughter branch from parent axis \[rad]
    pub angle2: Angle<T>,
    /// Flow rate in parent [m³/s]
    pub parent_flow: VolumetricFlowRate<T>,
    /// Flow rate in daughter 1 [m³/s]
    pub daughter1_flow: VolumetricFlowRate<T>,
    /// Flow rate in daughter 2 [m³/s]
    pub daughter2_flow: VolumetricFlowRate<T>,
}

impl<T: CfdScalar + FloatElement + Copy> OptimalBifurcation<T> {
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
    pub fn symmetric(parent_diameter: Length<T>, parent_flow: VolumetricFlowRate<T>) -> Self {
        let murray = MurraysLaw::<T>::new();
        let d_daughter = murray.symmetric_daughter_diameter(parent_diameter.into_base());

        let two = T::ONE + T::ONE;
        let daughter_flow = VolumetricFlowRate::from_base(parent_flow.into_base() / two);

        // cos(θ) = 2^(−1/3)
        let one_third = T::ONE / (T::ONE + T::ONE + T::ONE);
        let cos_theta = <T as FloatElement>::powf(two, -one_third);
        let angle = Angle::from_base(<T as FloatElement>::acos(cos_theta));

        Self {
            parent_diameter,
            daughter1_diameter: Length::from_base(d_daughter),
            daughter2_diameter: Length::from_base(d_daughter),
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
    pub fn asymmetric(
        parent_diameter: Length<T>,
        parent_flow: VolumetricFlowRate<T>,
        flow_ratio: T,
    ) -> Self {
        let one = T::ONE;
        let parent_diameter_base = parent_diameter.into_base();
        let parent_flow_base = parent_flow.into_base();

        let q1 = VolumetricFlowRate::from_base(parent_flow_base * flow_ratio);
        let q2 = VolumetricFlowRate::from_base(parent_flow_base * (one - flow_ratio));

        // D₁ = D₀ · (Q₁/Q₀)^(1/3), D₂ = D₀ · (Q₂/Q₀)^(1/3)
        let one_third = one / (T::ONE + T::ONE + T::ONE);
        let d1 = parent_diameter_base * <T as FloatElement>::powf(flow_ratio, one_third);
        let d2 = parent_diameter_base * <T as FloatElement>::powf(one - flow_ratio, one_third);

        // Zamir (1978) branching angle formula
        let d0_sq = parent_diameter_base * parent_diameter_base;
        let d1_sq = d1 * d1;
        let d2_sq = d2 * d2;
        let d0_4 = d0_sq * d0_sq;
        let d1_4 = d1_sq * d1_sq;
        let d2_4 = d2_sq * d2_sq;
        let two = T::ONE + T::ONE;

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
            daughter1_diameter: Length::from_base(d1),
            daughter2_diameter: Length::from_base(d2),
            angle1: Angle::from_base(<T as FloatElement>::acos(cos_theta1)),
            angle2: Angle::from_base(<T as FloatElement>::acos(cos_theta2)),
            parent_flow,
            daughter1_flow: q1,
            daughter2_flow: q2,
        }
    }

    /// Validate that this bifurcation satisfies Murray's Law
    pub fn is_murray_compliant(&self, tolerance: T) -> bool {
        let murray = MurraysLaw::<T>::new();
        murray.is_valid(
            self.parent_diameter.into_base(),
            self.daughter1_diameter.into_base(),
            self.daughter2_diameter.into_base(),
            tolerance,
        )
    }

    /// Calculate Murray's Law deviation
    pub fn murray_deviation(&self) -> Dimensionless<T> {
        let murray = MurraysLaw::<T>::new();
        Dimensionless::from_base(murray.deviation(
            self.parent_diameter.into_base(),
            self.daughter1_diameter.into_base(),
            self.daughter2_diameter.into_base(),
        ))
    }

    /// Calculate total daughter area to parent area ratio
    pub fn area_ratio(&self) -> Dimensionless<T> {
        let a0 = self.parent_diameter.into_base() * self.parent_diameter.into_base();
        let a1 = self.daughter1_diameter.into_base() * self.daughter1_diameter.into_base();
        let a2 = self.daughter2_diameter.into_base() * self.daughter2_diameter.into_base();
        Dimensionless::from_base((a1 + a2) / a0)
    }

    /// Check mass conservation: Q₀ = Q₁ + Q₂
    pub fn mass_conservation_error(&self) -> Dimensionless<T> {
        let parent_flow = self.parent_flow.into_base();
        let sum = self.daughter1_flow.into_base() + self.daughter2_flow.into_base();
        Dimensionless::from_base(<T as NumericElement>::abs(parent_flow - sum) / parent_flow)
    }

    /// Calculate pressure drop through major daughter branch
    ///
    /// Uses Poiseuille resistance: ΔP = 8μLQ/(πR⁴)
    pub fn pressure_drop_daughter1(
        &self,
        viscosity: DynamicViscosity<T>,
        length: Length<T>,
    ) -> Pressure<T> {
        let pi = T::pi();
        let eight = <T as FloatElement>::from_f64(8.0);
        let r = self.daughter1_diameter.into_base() / (T::ONE + T::ONE);
        Pressure::from_base(
            eight * viscosity.into_base() * length.into_base() * self.daughter1_flow.into_base()
                / (pi * <T as FloatElement>::powi(r, 4)),
        )
    }
}
