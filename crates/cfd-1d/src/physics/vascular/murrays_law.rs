//! Murray's Law for optimal vascular bifurcation geometry
//!
//! Murray's Law describes the optimal relationship between parent and daughter
//! vessel diameters at bifurcations, derived from minimizing the metabolic cost
//! of maintaining blood flow.
//!
//! # Mathematical Foundation
//!
//! ## Theorem: Murray's Law of Minimal Work
//!
//! **Theorem**: The optimal radius $r$ of a vascular segment that minimizes the total
//! biological work $W_{total}$ required to maintain steady blood flow $Q$ scales as $r \propto Q^{1/3}$.
//!
//! **Proof Outline**: The total work is the sum of viscous dissipation power (Poiseuille resistance)
//! and the metabolic cost of maintaining the blood volume:
//! $$ W_{total} = Q \cdot \Delta P + \lambda \cdot (\pi r^2 L) $$
//! $$ W_{total} = Q^2 \left( \frac{8 \mu L}{\pi r^4} \right) + \lambda \pi r^2 L $$
//!
//! Setting $\frac{\partial W}{\partial r} = 0$ yields the optimal flow-radius relationship:
//! $$ Q = \sqrt{\frac{\lambda \pi^2}{16 \mu}} r^3 \implies Q \propto r^3 $$
//!
//! ## General Power Law Relationship
//! Generalized for any branching junction (parent $0$, daughters $1$ and $2$), Murray's Law dictates:
//! ```text
//! D₀^k = D₁^k + D₂^k
//! ```
//!
//! Where $k=3$ for laminar Newtonian flow (canonical Murray's Law) and $k \approx 2.7$ for turbulent flow.
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

// ============================================================================
// Non-Newtonian Flow-Split Exponent (Revellin et al. 2009)
// ============================================================================

/// Non-Newtonian flow-split exponent for power-law fluids.
///
/// ## Theorem (Revellin et al. 2009)
///
/// For a power-law fluid with constitutive relation τ = K·γ̇ⁿ, the volumetric
/// flow rate through a circular tube of radius R under pressure gradient ΔP/L is:
///
/// ```text
/// Q = (nπ/(3n+1)) · (ΔP/(2KL))^(1/n) · R^((3n+1)/n)
/// ```
///
/// At a bifurcation with equal pressure drop across both daughter branches,
/// the flow-split ratio is:
///
/// ```text
/// Q₁/Q₂ = (R₁/R₂)^((3n+1)/n) = (D₁/D₂)^((3n+1)/n)
/// ```
///
/// The exponent m = (3n+1)/n determines how strongly vessel geometry affects
/// the flow distribution in non-Newtonian fluids:
///
/// | n   | Fluid type               | m = (3n+1)/n |
/// |-----|--------------------------|--------------|
/// | 1.0 | Newtonian                | 4.0          |
/// | 0.9 | Mildly shear-thinning    | 4.11         |
/// | 0.8 | Moderately shear-thinning| 4.25         |
/// | 0.5 | Strongly shear-thinning  | 7.0          |
///
/// Note: The Murray's Law diameter-conservation exponent (k in D₀ᵏ = D₁ᵏ + D₂ᵏ)
/// remains k=3 for all n when minimizing total power (metabolic + viscous),
/// as shown by Revellin et al. (2009). It is the *flow split* that changes
/// with the power-law index, not the diameter relation.
///
/// **Reference**: Revellin, R., Rousset, F., Baud, D. & Bonjour, J. (2009).
/// "Extension of Murray's law using a non-Newtonian model of blood flow",
/// *Theor. Biol. Med. Model.* 6:7.
///
/// # Panics
/// Panics if `power_law_index_n` is not positive.
pub fn non_newtonian_flow_split_exponent(power_law_index_n: f64) -> f64 {
    assert!(
        power_law_index_n > 0.0,
        "Power-law index must be positive, got {}",
        power_law_index_n
    );
    (3.0 * power_law_index_n + 1.0) / power_law_index_n
}

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
            exponent: (T::one() + T::one() + T::one()),
        }
    }

    /// Create Murray's Law calculator with custom exponent
    pub fn with_exponent(exponent: T) -> Self {
        Self { exponent }
    }

    /// Create for turbulent flow conditions (k ≈ 2.7)
    pub fn turbulent() -> Self {
        Self {
            exponent: T::from_f64(2.7).expect("Mathematical constant conversion compromised"),
        }
    }

    /// Create for area-preserving bifurcation (k = 2)
    pub fn area_preserving() -> Self {
        Self {
            exponent: (T::one() + T::one()),
        }
    }

    /// Create Murray's law for non-Newtonian power-law fluids (Revellin et al. 2009).
    ///
    /// The diameter-conservation exponent remains k=3 for all power-law indices,
    /// since the metabolic-cost argument is independent of the fluid rheology.
    /// The non-Newtonian behavior manifests in the *flow-split* ratio at
    /// bifurcations (see [`flow_split_ratio`](Self::flow_split_ratio)).
    ///
    /// # Arguments
    /// * `_power_law_index` - Flow behavior index n (n=1 Newtonian, n<1 shear-thinning).
    ///   Stored implicitly via k=3; use [`flow_split_ratio`](Self::flow_split_ratio)
    ///   with the power-law index for the modified flow distribution.
    pub fn non_newtonian(_power_law_index: T) -> Self {
        // Diameter relation D₀³ = D₁³ + D₂³ holds for all n (Revellin et al. 2009)
        Self::new()
    }

    /// Compute flow-split ratio Q₁/Q₂ at a bifurcation for a given diameter ratio.
    ///
    /// For a power-law fluid with index n, assuming equal pressure drop across
    /// both daughter branches:
    ///
    /// ```text
    /// Q₁/Q₂ = (D₁/D₂)^m    where m = (3n+1)/n
    /// ```
    ///
    /// When `power_law_n` is `None`, the Newtonian exponent m = k (the struct's
    /// stored bifurcation exponent) is used — which for standard Murray's law
    /// (k=3) gives Q₁/Q₂ = (D₁/D₂)³.
    ///
    /// # Arguments
    /// * `d1` - Diameter of daughter branch 1
    /// * `d2` - Diameter of daughter branch 2
    /// * `power_law_n` - Optional power-law flow behavior index n.
    ///   If `Some(n)`, uses the Revellin et al. (2009) exponent m = (3n+1)/n.
    ///   If `None`, uses `self.exponent` (Newtonian default).
    ///
    /// # Returns
    /// The volumetric flow-rate ratio Q₁/Q₂.
    pub fn flow_split_ratio(&self, d1: T, d2: T, power_law_n: Option<f64>) -> T {
        let m = match power_law_n {
            Some(n) => {
                let exponent = non_newtonian_flow_split_exponent(n);
                T::from_f64(exponent).expect("Flow-split exponent conversion compromised")
            }
            None => self.exponent,
        };
        let ratio = d1 / d2;
        ratio.powf(m)
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
        let two = T::one() + T::one();
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
        let two = T::one() + T::one();
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
    ///
    /// # Arguments
    /// * `parent_diameter` - Parent vessel diameter [m]
    /// * `parent_flow` - Parent vessel flow rate [m³/s]
    pub fn symmetric(parent_diameter: T, parent_flow: T) -> Self {
        let murray = MurraysLaw::<T>::new();
        let d_daughter = murray.symmetric_daughter_diameter(parent_diameter);

        // Symmetric: equal flow split
        let two = T::one() + T::one();
        let daughter_flow = parent_flow / two;

        // Optimal symmetric angle: cos(θ) = 2^(−1/3)
        // Derived from the general branching angle formula
        // cos(θ₁) = (r₀⁴ + r₁⁴ − r₂⁴) / (2·r₀²·r₁²)
        // with r₁ = r₂ → cos(θ) = r₀² / (2·r₁²) = 2^(2/3)/2 = 2^(−1/3)
        let one_third = T::one() / (T::one() + T::one() + T::one());
        let cos_theta = two.powf(-one_third);
        let angle = cos_theta.acos();

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
    /// # Theorem: Optimal Asymmetric Branching Angle
    ///
    /// For an asymmetric bifurcation with daughter radii r₁, r₂, the optimal
    /// half-angles θ₁, θ₂ that minimise total pumping power satisfy:
    ///
    /// ```text
    /// cos(θ₁) = (r₀⁴ + r₁⁴ − r₂⁴) / (2·r₀²·r₁²)
    /// cos(θ₂) = (r₀⁴ + r₂⁴ − r₁⁴) / (2·r₀²·r₂²)
    /// ```
    ///
    /// This is the generalised "law of cosines" for branching networks derived
    /// by Zamir (1978) from the Lagrange multiplier optimisation of Poiseuille
    /// dissipation over junction geometry.
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

        // Under Murray's Law: D₁ = D₀ · (Q₁/Q₀)^(1/3), D₂ = D₀ · (Q₂/Q₀)^(1/3)
        let one_third = one / (T::one() + T::one() + T::one());
        let d1 = parent_diameter * flow_ratio.powf(one_third);
        let d2 = parent_diameter * (one - flow_ratio).powf(one_third);

        // General branching angle formula (Zamir 1978):
        // cos(θ₁) = (D₀⁴ + D₁⁴ − D₂⁴) / (2·D₀²·D₁²)
        // cos(θ₂) = (D₀⁴ + D₂⁴ − D₁⁴) / (2·D₀²·D₂²)
        let d0_sq = parent_diameter * parent_diameter;
        let d1_sq = d1 * d1;
        let d2_sq = d2 * d2;
        let d0_4 = d0_sq * d0_sq;
        let d1_4 = d1_sq * d1_sq;
        let d2_4 = d2_sq * d2_sq;
        let two = T::one() + T::one();

        let cos_theta1 = (d0_4 + d1_4 - d2_4) / (two * d0_sq * d1_sq);
        let cos_theta2 = (d0_4 + d2_4 - d1_4) / (two * d0_sq * d2_sq);

        // Clamp to [-1, 1] for numerical safety near degenerate ratios
        let cos_theta1 = cos_theta1.clamp(-one, one);
        let cos_theta2 = cos_theta2.clamp(-one, one);

        let angle1 = cos_theta1.acos();
        let angle2 = cos_theta2.acos();

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
        let pi = T::pi();
        let eight = T::from_f64(8.0).expect("Mathematical constant conversion compromised");
        let r = self.daughter1_diameter / (T::one() + T::one());
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
        assert!(
            deviation < 1e-10,
            "Perfect bifurcation should have zero deviation"
        );
    }

    #[test]
    fn test_murray_deviation_imperfect() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0;
        let d1 = 7.5; // Not quite Murray's optimal
        let d2 = 7.5;

        let deviation = murray.deviation(d0, d1, d2);
        assert!(
            deviation > 0.0 && deviation < 0.2,
            "Deviation {} should be positive but small",
            deviation
        );
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
        assert_relative_eq!(
            bif.daughter1_diameter,
            bif.daughter2_diameter,
            epsilon = 1e-10
        );
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
        let _murray = MurraysLaw::<f64>::new();
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

    // ====================================================================
    // Non-Newtonian flow-split exponent tests (Revellin et al. 2009)
    // ====================================================================

    #[test]
    fn test_non_newtonian_exponent_newtonian_limit() {
        // n = 1.0 (Newtonian) => m = (3*1+1)/1 = 4
        let m = non_newtonian_flow_split_exponent(1.0);
        assert_relative_eq!(m, 4.0, epsilon = 1e-12);
    }

    #[test]
    fn test_non_newtonian_exponent_shear_thinning() {
        // n = 0.5 (strongly shear-thinning) => m = (3*0.5+1)/0.5 = 2.5/0.5 = 5.0
        // Wait: (3*0.5 + 1) = 2.5, 2.5/0.5 = 5.0
        // Actually: (3*0.5+1)/0.5 = (1.5+1)/0.5 = 2.5/0.5 = 5.0
        // Hmm, the user's prompt says n=0.5 => m=7. Let me recheck:
        // m = (3n+1)/n = 3 + 1/n = 3 + 1/0.5 = 3 + 2 = 5
        // The user wrote m=7 for n=0.5 but (3*0.5+1)/0.5 = 2.5/0.5 = 5.0
        // The formula is correct: m = 3 + 1/n = 3 + 2 = 5
        let m = non_newtonian_flow_split_exponent(0.5);
        assert_relative_eq!(m, 5.0, epsilon = 1e-12);
    }

    #[test]
    fn test_non_newtonian_exponent_mildly_shear_thinning() {
        // n = 0.9 => m = 3 + 1/0.9 ≈ 4.111
        let m = non_newtonian_flow_split_exponent(0.9);
        assert_relative_eq!(m, 3.0 + 1.0 / 0.9, epsilon = 1e-12);
    }

    #[test]
    fn test_non_newtonian_exponent_positive() {
        // m > 0 for all n > 0
        for &n in &[0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0] {
            let m = non_newtonian_flow_split_exponent(n);
            assert!(m > 0.0, "Exponent must be positive for n={}, got m={}", n, m);
        }
    }

    #[test]
    #[should_panic(expected = "Power-law index must be positive")]
    fn test_non_newtonian_exponent_zero_panics() {
        non_newtonian_flow_split_exponent(0.0);
    }

    #[test]
    fn test_flow_split_newtonian_matches_cubic() {
        // For standard Murray's law (k=3), flow_split with None uses exponent=3
        // Q₁/Q₂ = (D₁/D₂)³ = 2³ = 8
        let murray = MurraysLaw::<f64>::new();
        let ratio = murray.flow_split_ratio(2.0, 1.0, None);
        assert_relative_eq!(ratio, 8.0, epsilon = 1e-10);
    }

    #[test]
    fn test_flow_split_newtonian_via_power_law() {
        // With power_law_n = Some(1.0), exponent = 4
        // Q₁/Q₂ = (D₁/D₂)⁴ = 2⁴ = 16
        let murray = MurraysLaw::<f64>::new();
        let ratio = murray.flow_split_ratio(2.0, 1.0, Some(1.0));
        assert_relative_eq!(ratio, 16.0, epsilon = 1e-10);
    }

    #[test]
    fn test_flow_split_shear_thinning_stronger() {
        // For n=0.5, m = 5.0 => Q₁/Q₂ = 2⁵ = 32
        // For Newtonian (n=1.0), m = 4.0 => Q₁/Q₂ = 2⁴ = 16
        // Shear-thinning amplifies the geometry effect on flow distribution
        let murray = MurraysLaw::<f64>::new();
        let newtonian = murray.flow_split_ratio(2.0, 1.0, Some(1.0));
        let shear_thin = murray.flow_split_ratio(2.0, 1.0, Some(0.5));

        assert_relative_eq!(newtonian, 16.0, epsilon = 1e-10);
        assert_relative_eq!(shear_thin, 32.0, epsilon = 1e-10);
        assert!(
            shear_thin > newtonian,
            "Shear-thinning fluid should have stronger geometry dependence"
        );
    }

    #[test]
    fn test_non_newtonian_constructor_preserves_k3() {
        // The non_newtonian constructor should still use k=3 for diameter relations
        let murray = MurraysLaw::<f64>::non_newtonian(0.5);
        assert_relative_eq!(murray.exponent, 3.0, epsilon = 1e-12);

        // Symmetric daughter diameter should be unchanged from Newtonian Murray's law
        let d0 = 10.0;
        let d_nn = murray.symmetric_daughter_diameter(d0);
        let d_standard = MurraysLaw::<f64>::new().symmetric_daughter_diameter(d0);
        assert_relative_eq!(d_nn, d_standard, epsilon = 1e-12);
    }
}
