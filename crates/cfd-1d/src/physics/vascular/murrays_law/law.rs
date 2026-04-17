//! Murray's Law calculator for optimal bifurcation design.
//!
//! # Theorem — Murray's Cube Law (Murray 1926)
//!
//! **Theorem**: For a bifurcating vascular network transporting a Newtonian
//! fluid under Poiseuille flow, the parent and daughter diameters satisfy:
//!
//! ```text
//! D₀ᵏ = D₁ᵏ + D₂ᵏ,   k = 3  (Newtonian, laminar)
//! ```
//!
//! **Proof sketch** (metabolic cost minimization):
//!
//! The total cost of maintaining a vessel segment of diameter D and length L is
//! the sum of viscous power dissipation and metabolic maintenance cost:
//!
//! ```text
//! W(D) = (128 μ L Q²) / (π D⁴)  +  b · (π D² L / 4)
//! ```
//!
//! where b is the metabolic cost per unit vessel volume. Setting ∂W/∂D = 0
//! yields Q ∝ D³. At a bifurcation node, conservation of flow Q₀ = Q₁ + Q₂
//! combined with Q ∝ D³ gives D₀³ = D₁³ + D₂³.  ∎
//!
//! ## Corollaries
//!
//! - **Symmetric bifurcation**: D₁ = D₂ = D₀ / 2^(1/k)
//! - **Area ratio**: A_daughters / A_parent = 2^(1 − 2/k). For k=3: 2^(1/3) ≈ 1.2599.
//! - **Deviation**: ε = |D₀ᵏ − D₁ᵏ − D₂ᵏ| / D₀ᵏ
//!
//! ## Non-Newtonian Extension (Revellin et al. 2009)
//!
//! For power-law fluids with index n, the diameter relation D₀³ = D₁³ + D₂³
//! holds for all n (metabolic-cost argument is rheology-independent). The
//! non-Newtonian effect appears in the flow-split exponent:
//!
//! ```text
//! Q₁/Q₂ = (D₁/D₂)^m,   m = (3n + 1)/n
//! ```
//!
//! For Newtonian (n = 1): m = 4. For shear-thinning blood (n ≈ 0.7): m ≈ 4.86.
//!
//! # References
//!
//! - Murray, C. D. (1926). "The Physiological Principle of Minimum Work:
//!   I. The Vascular System and the Cost of Blood Volume."
//!   *Proc. Natl. Acad. Sci.* 12(3):207–214.
//! - Revellin, R. et al. (2009). "Extension of Murray's law using a
//!   non-Newtonian model of blood flow." *Theor. Biol. Med. Model.* 6:7.

use super::non_newtonian_flow_split_exponent;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

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
    /// When `power_law_n` is `None`, the Newtonian exponent m = k is used.
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Symmetric daughter diameter: D₁ = D₂ = D₀ / 2^(1/3).
    /// Analytical: D₀ = 1.0 mm → D_d = 1.0 / 2^(1/3) ≈ 0.7937 mm.
    #[test]
    fn symmetric_daughter_diameter_analytical() {
        let murray = MurraysLaw::<f64>::new();
        let d0: f64 = 1.0e-3;
        let dd = murray.symmetric_daughter_diameter(d0);
        let expected = d0 / 2.0_f64.powf(1.0 / 3.0);
        assert_relative_eq!(dd, expected, max_relative = 1e-12);
    }

    /// Asymmetric daughter: given D₀ and D₁, D₂ = (D₀³ − D₁³)^(1/3).
    /// D₀ = 1.0 mm, D₁ = 0.8 mm → D₂ = (1.0 − 0.512)^(1/3) = 0.488^(1/3) ≈ 0.7874 mm.
    #[test]
    fn asymmetric_daughter_diameter_analytical() {
        let murray = MurraysLaw::<f64>::new();
        let d0: f64 = 1.0e-3;
        let d1: f64 = 0.8e-3;
        let d2 = murray
            .asymmetric_daughter_diameter(d0, d1)
            .expect("valid bifurcation");
        let expected = (d0.powi(3) - d1.powi(3)).powf(1.0 / 3.0);
        assert_relative_eq!(d2, expected, max_relative = 1e-12);
    }

    /// For a valid Murray bifurcation, deviation must be zero.
    /// D₀ = 1.0, D₁ = D₂ = D₀ / 2^(1/3) → ε = 0.
    #[test]
    fn valid_bifurcation_zero_deviation() {
        let murray = MurraysLaw::<f64>::new();
        let d0: f64 = 1.0e-3;
        let dd = murray.symmetric_daughter_diameter(d0);
        let eps = murray.deviation(d0, dd, dd);
        assert!(
            eps < 1e-12,
            "deviation ε = {eps} should be ~0 for valid Murray bifurcation"
        );
    }

    /// Area ratio for symmetric Murray bifurcation: 2^(1 − 2/3) = 2^(1/3) ≈ 1.2599.
    #[test]
    fn ideal_area_ratio_k3() {
        let murray = MurraysLaw::<f64>::new();
        let ar = murray.ideal_area_ratio();
        let expected = 2.0_f64.powf(1.0 / 3.0);
        assert_relative_eq!(ar, expected, max_relative = 1e-12);
    }

    /// Flow split ratio for Newtonian: (D₁/D₂)^3, but with power_law_n=None
    /// uses exponent k=3. D₁ = 2D₂ → ratio = 8.
    #[test]
    fn flow_split_newtonian_d_ratio_two() {
        let murray = MurraysLaw::<f64>::new();
        let ratio = murray.flow_split_ratio(2.0e-3, 1.0e-3, None);
        assert_relative_eq!(ratio, 8.0, max_relative = 1e-12);
    }

    /// Non-Newtonian flow split: m = (3n+1)/n. For n=0.7: m = 4.857...
    /// D₁/D₂ = 2 → ratio = 2^4.857 ≈ 29.0.
    #[test]
    fn flow_split_non_newtonian_power_law() {
        let murray = MurraysLaw::<f64>::new();
        let n = 0.7;
        let m_expected = (3.0 * n + 1.0) / n;
        let ratio = murray.flow_split_ratio(2.0e-3, 1.0e-3, Some(n));
        let expected = 2.0_f64.powf(m_expected);
        assert_relative_eq!(ratio, expected, max_relative = 1e-10);
    }

    /// Parent diameter round-trip: parent(d1, d2) → reconstruct (d1, d2).
    /// D₁ = 0.8 mm, D₂ = 0.6 mm → D₀ = (0.512 + 0.216)^(1/3) = 0.728^(1/3).
    #[test]
    fn parent_diameter_round_trip() {
        let murray = MurraysLaw::<f64>::new();
        let d1: f64 = 0.8e-3;
        let d2: f64 = 0.6e-3;
        let d0 = murray.parent_diameter(d1, d2);
        let d2_recon = murray.asymmetric_daughter_diameter(d0, d1).expect("valid");
        assert_relative_eq!(d2_recon, d2, max_relative = 1e-12);
    }

    /// Daughter larger than parent must return None.
    #[test]
    fn daughter_exceeds_parent_returns_none() {
        let murray = MurraysLaw::<f64>::new();
        let result = murray.asymmetric_daughter_diameter(1.0e-3, 1.5e-3);
        assert!(result.is_none(), "D₁ > D₀ should yield None");
    }
}
