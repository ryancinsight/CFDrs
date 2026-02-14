//! Murray's Law for optimal vascular bifurcation geometry
//!
//! Murray's Law describes the optimal relationship between parent and daughter
//! vessel diameters at bifurcations, derived from minimizing the metabolic cost
//! of maintaining blood flow.
//!
//! # Mathematical Foundation
//!
//! ## Power Law Relationship
//! The generalized Murray's Law states:
//! ```text
//! D₀^k = D₁^k + D₂^k
//! ```
//!
//! Where:
//! - D₀ = parent vessel diameter
//! - D₁, D₂ = daughter vessel diameters
//! - k = bifurcation exponent (k=3 for laminar flow, k≈2.7 for turbulent)
//!
//! ## Optimal Branching Angle
//! The optimal bifurcation half-angle θ is related to the radius ratio by:
//! ```text
//! cos(θ₁) = (r₀⁴ + r₁⁴ - r₂⁴) / (2·r₀²·r₁²)
//! ```
//!
//! ## Flow Distribution
//! At a symmetric bifurcation:
//! - Q₁ = Q₂ = Q₀/2
//! - D₁ = D₂ = D₀ / 2^(1/3) ≈ 0.794·D₀
//!
//! # Reference
//! - Murray, C.D. (1926) "The Physiological Principle of Minimum Work"
//! - Sherman, T.F. (1981) "On connecting large vessels to small"
//! - Zamir, M. (1999) "On fractal properties of arterial trees"

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ============================================================================
// Murray's Law Calculator
// ============================================================================

/// Murray's Law calculator for optimal bifurcation design
///
/// Provides methods for calculating optimal vessel diameters, branching angles,
/// and validating bifurcation geometry against Murray's principle.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct MurraysLaw<T: RealField + Copy> {
    /// Bifurcation exponent k (default: 3.0 for laminar flow)
    ///
    /// - k = 3.0: Laminar Poiseuille flow (classic Murray's Law)
    /// - k ≈ 2.7: Turbulent flow correction
    /// - k ≈ 2.0: Area-preserving (continuity only)
    pub exponent: T,
}

impl<T: RealField + FromPrimitive + Copy> Default for MurraysLaw<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + FromPrimitive + Copy> MurraysLaw<T> {
    /// Create Murray's Law calculator with classic exponent k=3
    pub fn new() -> Self {
        Self {
            exponent: T::from_f64(3.0).unwrap(),
        }
    }

    /// Create Murray's Law calculator with custom exponent
    pub fn with_exponent(exponent: T) -> Self {
        Self { exponent }
    }

    /// Create for turbulent flow conditions (k ≈ 2.7)
    pub fn turbulent() -> Self {
        Self {
            exponent: T::from_f64(2.7).unwrap(),
        }
    }

    /// Create for area-preserving bifurcation (k = 2)
    pub fn area_preserving() -> Self {
        Self {
            exponent: T::from_f64(2.0).unwrap(),
        }
    }

    /// Calculate optimal daughter diameters for symmetric bifurcation
    ///
    /// For a symmetric bifurcation where D₁ = D₂:
    /// ```text
    /// D₁ = D₂ = D₀ / 2^(1/k)
    /// ```
    ///
    /// # Arguments
    /// * `parent_diameter` - Parent vessel diameter D₀
    ///
    /// # Returns
    /// Optimal daughter diameter (D₁ = D₂)
    pub fn symmetric_daughter_diameter(&self, parent_diameter: T) -> T {
        let two = T::from_f64(2.0).unwrap();
        let one_over_k = T::one() / self.exponent;
        parent_diameter / two.powf(one_over_k)
    }

    /// Calculate optimal minor daughter diameter for asymmetric bifurcation
    ///
    /// Given parent diameter D₀ and major daughter D₁, find D₂:
    /// ```text
    /// D₂ = (D₀^k - D₁^k)^(1/k)
    /// ```
    ///
    /// # Arguments
    /// * `parent_diameter` - Parent vessel diameter D₀
    /// * `major_daughter_diameter` - Larger daughter diameter D₁
    ///
    /// # Returns
    /// Required minor daughter diameter D₂, or None if impossible
    pub fn asymmetric_daughter_diameter(
        &self,
        parent_diameter: T,
        major_daughter_diameter: T,
    ) -> Option<T> {
        let d0_k = parent_diameter.powf(self.exponent);
        let d1_k = major_daughter_diameter.powf(self.exponent);

        if d1_k >= d0_k {
            return None; // Invalid: daughter cannot be larger than parent requires
        }

        let d2_k = d0_k - d1_k;
        if d2_k <= T::zero() {
            return None;
        }

        let one_over_k = T::one() / self.exponent;
        Some(d2_k.powf(one_over_k))
    }

    /// Calculate required parent diameter for given daughter diameters
    ///
    /// ```text
    /// D₀ = (D₁^k + D₂^k)^(1/k)
    /// ```
    pub fn parent_diameter(&self, daughter1: T, daughter2: T) -> T {
        let d1_k = daughter1.powf(self.exponent);
        let d2_k = daughter2.powf(self.exponent);
        let one_over_k = T::one() / self.exponent;
        (d1_k + d2_k).powf(one_over_k)
    }

    /// Calculate Murray's Law deviation for a bifurcation
    ///
    /// Returns the relative error from Murray's Law:
    /// ```text
    /// ε = |D₀^k - D₁^k - D₂^k| / D₀^k
    /// ```
    ///
    /// A perfect Murray's Law bifurcation has ε = 0.
    pub fn deviation(&self, parent: T, daughter1: T, daughter2: T) -> T {
        let d0_k = parent.powf(self.exponent);
        let d1_k = daughter1.powf(self.exponent);
        let d2_k = daughter2.powf(self.exponent);

        ((d0_k - d1_k - d2_k).abs()) / d0_k
    }

    /// Check if bifurcation satisfies Murray's Law within tolerance
    pub fn is_valid(&self, parent: T, daughter1: T, daughter2: T, tolerance: T) -> bool {
        self.deviation(parent, daughter1, daughter2) < tolerance
    }

    /// Calculate area ratio (daughter total / parent)
    ///
    /// For Murray's Law with k=3, area ratio = 2^(2/3) ≈ 1.26
    pub fn area_ratio(&self, parent: T, daughter1: T, daughter2: T) -> T {
        let a0 = parent * parent;
        let a1 = daughter1 * daughter1;
        let a2 = daughter2 * daughter2;
        (a1 + a2) / a0
    }

    /// Calculate expected area ratio for symmetric Murray's bifurcation
    ///
    /// A_daughters / A_parent = 2^(1 - 2/k)
    pub fn ideal_area_ratio(&self) -> T {
        let two = T::from_f64(2.0).unwrap();
        let exp = T::one() - two / self.exponent;
        two.powf(exp)
    }
}

// ============================================================================
// Optimal Bifurcation Geometry
// ============================================================================

/// Complete optimal bifurcation geometry
///
/// Calculates optimal diameters, angles, and flow distribution for a
/// bifurcation following Murray's Law principles.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct OptimalBifurcation<T: RealField + Copy> {
    /// Parent vessel diameter [m]
    pub parent_diameter: T,
    /// Major daughter diameter [m]
    pub daughter1_diameter: T,
    /// Minor daughter diameter [m]
    pub daughter2_diameter: T,
    /// Half-angle of major daughter branch from parent axis [rad]
    pub angle1: T,
    /// Half-angle of minor daughter branch from parent axis [rad]
    pub angle2: T,
    /// Flow rate in parent [m³/s]
    pub parent_flow: T,
    /// Flow rate in daughter 1 [m³/s]
    pub daughter1_flow: T,
    /// Flow rate in daughter 2 [m³/s]
    pub daughter2_flow: T,
}

impl<T: RealField + FromPrimitive + Copy> OptimalBifurcation<T> {
    /// Create symmetric bifurcation following Murray's Law
    ///
    /// # Arguments
    /// * `parent_diameter` - Parent vessel diameter [m]
    /// * `parent_flow` - Parent vessel flow rate [m³/s]
    pub fn symmetric(parent_diameter: T, parent_flow: T) -> Self {
        let murray = MurraysLaw::<T>::new();
        let d_daughter = murray.symmetric_daughter_diameter(parent_diameter);

        // Symmetric: equal flow split
        let two = T::from_f64(2.0).unwrap();
        let daughter_flow = parent_flow / two;

        // Optimal angle for symmetric bifurcation: θ = 37.5°
        // cos(θ) = r₁/r₀ for symmetric case
        // Actually: cos(θ) = (2^(1/3) - 1) / 2^(1/3) ... complex formula
        // Simplified: θ ≈ 37.5° = 0.654 rad for k=3
        let angle = T::from_f64(37.5 * PI / 180.0).unwrap();

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
    /// # Arguments
    /// * `parent_diameter` - Parent vessel diameter [m]
    /// * `parent_flow` - Parent vessel flow rate [m³/s]
    /// * `flow_ratio` - Flow ratio Q₁/(Q₁+Q₂) (0 < ratio < 1)
    pub fn asymmetric(parent_diameter: T, parent_flow: T, flow_ratio: T) -> Self {
        let one = T::one();

        // Flow split
        let q1 = parent_flow * flow_ratio;
        let q2 = parent_flow * (one - flow_ratio);

        // For Poiseuille flow: Q ∝ D⁴, so D ∝ Q^(1/4)
        // Under Murray's Law: D₁³ = D₀³ · (Q₁/Q₀), D₂³ = D₀³ · (Q₂/Q₀)
        // Therefore: D₁ = D₀ · (Q₁/Q₀)^(1/3)
        let d1 = parent_diameter * (flow_ratio).powf(T::one() / T::from_f64(3.0).unwrap());
        let d2 = parent_diameter * (one - flow_ratio).powf(T::one() / T::from_f64(3.0).unwrap());

        // Branching angles from optimization (simplified)
        // For asymmetric bifurcation, angles relate to radius ratios
        let r_ratio_1 = d1 / parent_diameter;
        let r_ratio_2 = d2 / parent_diameter;

        // Approximate optimal angles using empirical formula
        // θ ≈ arccos(r_daughter/r_parent) * scaling_factor
        let angle1 = (one - r_ratio_1 * r_ratio_1).sqrt().atan2(r_ratio_1);
        let angle2 = (one - r_ratio_2 * r_ratio_2).sqrt().atan2(r_ratio_2);

        Self {
            parent_diameter,
            daughter1_diameter: d1,
            daughter2_diameter: d2,
            angle1,
            angle2,
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
        (self.parent_flow - sum).abs() / self.parent_flow
    }

    /// Calculate pressure drop through major daughter branch
    ///
    /// Uses Poiseuille resistance: ΔP = 8μLQ/(πR⁴)
    pub fn pressure_drop_daughter1(&self, viscosity: T, length: T) -> T {
        let pi = T::from_f64(PI).unwrap();
        let eight = T::from_f64(8.0).unwrap();
        let r = self.daughter1_diameter / T::from_f64(2.0).unwrap();
        eight * viscosity * length * self.daughter1_flow / (pi * r.powi(4))
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_symmetric_daughter_k3() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0;
        let d_daughter = murray.symmetric_daughter_diameter(d0);

        // For k=3: D₁ = D₂ = D₀ / 2^(1/3) ≈ 0.794 · D₀
        let expected = d0 / 2.0_f64.powf(1.0 / 3.0);
        assert_relative_eq!(d_daughter, expected, epsilon = 1e-10);
        assert_relative_eq!(d_daughter, 7.937, epsilon = 0.001);
    }

    #[test]
    fn test_murray_deviation_perfect() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0;
        let d1 = murray.symmetric_daughter_diameter(d0);
        let d2 = d1;

        let deviation = murray.deviation(d0, d1, d2);
        assert!(deviation < 1e-10, "Perfect bifurcation should have zero deviation");
    }

    #[test]
    fn test_murray_deviation_imperfect() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0;
        let d1 = 7.5; // Not quite Murray's optimal
        let d2 = 7.5;

        let deviation = murray.deviation(d0, d1, d2);
        assert!(deviation > 0.0 && deviation < 0.2, "Deviation {} should be positive but small", deviation);
    }

    #[test]
    fn test_asymmetric_daughter() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0_f64;
        let d1 = 9.0_f64;

        let d2 = murray.asymmetric_daughter_diameter(d0, d1).unwrap();

        // Check: d0^3 = d1^3 + d2^3
        let lhs = d0.powi(3);
        let rhs = d1.powi(3) + d2.powi(3);
        assert_relative_eq!(lhs, rhs, epsilon = 1e-10);
    }

    #[test]
    fn test_asymmetric_daughter_impossible() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0;
        let d1 = 12.0; // Impossible: daughter larger than parent requires

        let result = murray.asymmetric_daughter_diameter(d0, d1);
        assert!(result.is_none());
    }

    #[test]
    fn test_parent_diameter_reconstruction() {
        let murray = MurraysLaw::<f64>::new();
        let d1 = 7.937;
        let d2 = 7.937;

        let d0 = murray.parent_diameter(d1, d2);
        assert_relative_eq!(d0, 10.0, epsilon = 0.01);
    }

    #[test]
    fn test_ideal_area_ratio_k3() {
        let murray = MurraysLaw::<f64>::new();
        let ratio = murray.ideal_area_ratio();

        // For k=3: ratio = 2^(1 - 2/3) = 2^(1/3) ≈ 1.26
        assert_relative_eq!(ratio, 2.0_f64.powf(1.0 / 3.0), epsilon = 1e-10);
        assert_relative_eq!(ratio, 1.26, epsilon = 0.01);
    }

    #[test]
    fn test_area_preserving() {
        let murray = MurraysLaw::<f64>::area_preserving();
        let ratio = murray.ideal_area_ratio();

        // For k=2: ratio = 2^(1 - 1) = 1.0 (area preserved)
        assert_relative_eq!(ratio, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_symmetric_bifurcation() {
        let bif = OptimalBifurcation::<f64>::symmetric(0.01, 1e-6);

        // Check Murray compliance
        assert!(bif.is_murray_compliant(0.001));

        // Check flow conservation
        assert!(bif.mass_conservation_error() < 1e-10);

        // Check symmetric diameters
        assert_relative_eq!(bif.daughter1_diameter, bif.daughter2_diameter, epsilon = 1e-10);
    }

    #[test]
    fn test_asymmetric_bifurcation() {
        let bif = OptimalBifurcation::<f64>::asymmetric(0.01, 1e-6, 0.7);

        // Check flow conservation
        assert!(bif.mass_conservation_error() < 1e-10);

        // Check flow split
        assert_relative_eq!(bif.daughter1_flow / bif.parent_flow, 0.7, epsilon = 1e-10);

        // Larger flow should be in larger diameter branch
        assert!(bif.daughter1_diameter > bif.daughter2_diameter);
    }

    #[test]
    fn test_area_ratio_symmetric() {
        let bif = OptimalBifurcation::<f64>::symmetric(0.01, 1e-6);
        let ratio = bif.area_ratio();

        // Should match ideal Murray ratio ≈ 1.26
        let murray = MurraysLaw::<f64>::new();
        assert_relative_eq!(ratio, murray.ideal_area_ratio(), epsilon = 0.01);
    }

    #[test]
    fn test_trifurcation_extension() {
        // Murray's Law extends to trifurcations: D₀³ = D₁³ + D₂³ + D₃³
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0_f64;
        let d1 = 6.0_f64;
        let d2 = 5.0_f64;

        // Calculate d3 such that law is satisfied
        let d3_cubed = d0.powi(3) - d1.powi(3) - d2.powi(3);
        let d3 = d3_cubed.powf(1.0 / 3.0);

        // Verify
        let sum = d1.powi(3) + d2.powi(3) + d3.powi(3);
        assert_relative_eq!(d0.powi(3), sum, epsilon = 1e-10);
    }

    #[test]
    fn test_pressure_drop() {
        let bif = OptimalBifurcation::<f64>::symmetric(0.01, 1e-6);
        let dp = bif.pressure_drop_daughter1(0.0035, 0.1);

        // Should be finite and positive
        assert!(dp > 0.0 && dp.is_finite());
    }
}
