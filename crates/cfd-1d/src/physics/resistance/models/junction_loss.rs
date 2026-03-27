//! Junction loss resistance models for T-junctions and cross-junctions.
//!
//! When channels intersect in a planar millifluidic device, additional
//! minor losses arise from flow splitting, merging, and impinging at the
//! junction.  These K-factor losses supplement the straight-channel
//! (Hagen-Poiseuille) resistance.
//!
//! # Theorem — Minor Loss Coefficient
//!
//! The pressure drop across a junction is:
//!
//! $$\Delta P = K \cdot \frac{\rho V^2}{2}$$
//!
//! where $K$ is the loss coefficient and $V$ is the velocity in the
//! reference branch (typically the upstream/combined branch).
//!
//! **Proof sketch**: Energy dissipation at a junction arises from
//! recirculation zones, flow separation, and mixing.  The $K$-factor
//! captures these losses as a fraction of the dynamic head.  Values
//! are determined experimentally (Idelchik 2007, Miller 1990) or via
//! validated CFD studies.
//!
//! # Cross-Junction K-Factors
//!
//! For a planar cross-junction (two channels crossing at 90°):
//!
//! | Configuration | K (combining) | K (dividing) |
//! |---|---|---|
//! | Equal-flow 90° cross | 1.2 | 1.5 |
//! | Low-Re millifluidic  | 0.8–1.0 | 1.0–1.3 |
//!
//! # T-Junction K-Factors (Idelchik 2007)
//!
//! | Configuration | K (branch) | K (run) |
//! |---|---|---|
//! | Dividing, equal area | 1.0 | 0.05 |
//! | Combining, equal area | 1.5 | 0.04 |
//! | 90° sharp-edge | 1.2 | 0.1 |
//!
//! # References
//!
//! - Idelchik, I. E. (2007). *Handbook of Hydraulic Resistance* (4th ed.).
//!   Begell House. Ch. 7: Resistance in merging and dividing flow branches.
//! - Miller, D. S. (1990). *Internal Flow Systems* (2nd ed.). BHRA.
//!   Ch. 13: Junction losses.
//! - Khodabandeh, E. et al. (2020). "Pressure drop and heat transfer in
//!   T-junction microchannels." *Int. J. Heat Mass Transfer*, 154, 119689.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Junction type classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum JunctionType {
    /// T-junction: one branch meets a straight run at 90°.
    Tee,
    /// Cross-junction: two channels cross at 90°.
    Cross,
}

/// Flow direction through the junction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum JunctionFlowDirection {
    /// Flow divides at the junction (one inlet, multiple outlets).
    Dividing,
    /// Flow combines at the junction (multiple inlets, one outlet).
    Combining,
}

/// Junction loss model for T-junction and cross-junction channels.
///
/// This model computes minor-loss pressure drops at channel intersections
/// using experimentally validated K-factor correlations.
///
/// # Theorem
///
/// The junction resistance is expressed as an equivalent hydraulic
/// resistance: $R_J = K \rho V / (2 A^2)$, which adds to the
/// Hagen-Poiseuille resistance of the connecting channels.
///
/// **Proof sketch**: Starting from $\Delta P = K \rho V^2 / 2$ and
/// $Q = V A$, we get $\Delta P = K \rho Q^2 / (2 A^2)$.  This is a
/// quadratic loss term; in linearised form (about a reference flow rate
/// $Q_0$) the effective resistance is $R_J = K \rho Q_0 / A^2$.
///
/// # Diameter-ratio correction (Idelchik 2007, §7.12)
///
/// When the branch and run channels have different diameters, the K-factor
/// is modulated by the area ratio $A_b / A_r$:
///
/// $$K_{\text{eff}} = K_{\text{equal}} \cdot \left(\frac{A_b}{A_r}\right)^{0.5}$$
///
/// A narrow branch meeting a wide run (ratio < 1) reduces K (less
/// recirculation); a wide branch into a narrow run increases K.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JunctionLossModel {
    /// Type of junction (T or cross).
    pub junction_type: JunctionType,
    /// Flow direction (dividing or combining).
    pub flow_direction: JunctionFlowDirection,
    /// Hydraulic diameter of the branch channel [m].
    pub branch_diameter_m: f64,
    /// Cross-sectional area of the branch channel [m²].
    pub branch_area_m2: f64,
    /// Length attributed to the junction (typically ≈ D_h) [m].
    pub junction_length_m: f64,
    /// Fluid density [kg/m³].
    pub density_kg_m3: f64,
    /// Optional cross-sectional area of the run (main) channel [m²].
    ///
    /// When provided, the K-factor is corrected for diameter mismatch
    /// using the Idelchik area-ratio correlation.  When `None`, the
    /// equal-area K-factor is used (backward-compatible).
    #[serde(default)]
    pub run_area_m2: Option<f64>,
}

impl JunctionLossModel {
    /// Retrieve the K-factor for this junction configuration.
    ///
    /// When [`run_area_m2`] is set, the equal-area K is corrected by the
    /// area ratio `(A_branch / A_run)^0.5` (Idelchik 2007, §7.12).
    /// This captures the additional recirculation loss when a narrow
    /// branch intersects a wide run channel, or the reduced loss when
    /// a wide branch meets a narrow run.
    #[must_use]
    pub fn loss_coefficient(&self) -> f64 {
        let k_equal = match (self.junction_type, self.flow_direction) {
            (JunctionType::Tee, JunctionFlowDirection::Dividing) => 1.0,
            (JunctionType::Tee, JunctionFlowDirection::Combining) => 1.5,
            (JunctionType::Cross, JunctionFlowDirection::Dividing) => 1.3,
            (JunctionType::Cross, JunctionFlowDirection::Combining) => 1.2,
        };

        if let Some(run_area) = self.run_area_m2 {
            if run_area > 1e-30 {
                let ratio = (self.branch_area_m2 / run_area).sqrt();
                return k_equal * ratio;
            }
        }
        k_equal
    }

    /// Compute the linearised junction resistance about a reference flow.
    ///
    /// Returns $R_J = K \rho Q_{\text{ref}} / A^2$ [Pa·s/m³].
    #[must_use]
    pub fn linearised_resistance(&self, reference_flow_m3_s: f64) -> f64 {
        let k = self.loss_coefficient();
        let a2 = self.branch_area_m2 * self.branch_area_m2;
        if a2 < 1e-30 {
            return 0.0;
        }
        k * self.density_kg_m3 * reference_flow_m3_s.abs() / a2
    }
}

impl<T: RealField + Copy + FromPrimitive> ResistanceModel<T> for JunctionLossModel {
    /// Calculate hydraulic resistance [Pa·s/m³].
    ///
    /// Uses the linear + quadratic coefficient model; the linear part
    /// accounts for viscous losses in the short junction length while the
    /// quadratic part captures the K-factor minor loss.
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let (r_lin, _k_quad) = self.calculate_coefficients(fluid, conditions)?;
        Ok(r_lin)
    }

    /// Return (R_linear, K_quadratic) for the junction.
    ///
    /// - R_linear: Hagen-Poiseuille resistance through the short junction
    ///   length (≈ one hydraulic diameter).
    /// - K_quadratic: Minor loss coefficient scaled to produce
    ///   $\Delta P = k Q |Q|$.
    fn calculate_coefficients<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        let state = fluid.properties_at(conditions.temperature, conditions.pressure)?;
        let density = state.density;

        // Estimate velocity from flow conditions for shear rate.
        let area = T::from_f64(self.branch_area_m2)
            .ok_or_else(|| Error::InvalidInput("branch_area_m2 conversion failed".into()))?;
        let v = if let Some(v) = conditions.velocity {
            v
        } else if let Some(q) = conditions.flow_rate {
            q / area
        } else {
            T::zero()
        };

        // Shear rate ≈ 8V/D for circular duct.
        let d = T::from_f64(self.branch_diameter_m)
            .ok_or_else(|| Error::InvalidInput("branch_diameter_m conversion failed".into()))?;
        let eight = T::from_f64(8.0).expect("Mathematical constant conversion compromised");
        let shear_rate = eight * v / d;

        let mu =
            fluid.viscosity_at_shear(shear_rate, conditions.temperature, conditions.pressure)?;

        let l = T::from_f64(self.junction_length_m)
            .ok_or_else(|| Error::InvalidInput("junction_length_m conversion failed".into()))?;

        // Linear component: Hagen-Poiseuille through the junction length.
        // R_lin = 128 μ L / (π D⁴)  [circular], or 12 μ L / (w h³) [rect]
        let pi = T::pi();
        let r_lin =
            T::from_f64(128.0).expect("Mathematical constant conversion compromised") * mu * l
                / (pi * d.powi(4));

        // Quadratic component: K ρ / (2 A²).
        let k_factor = T::from_f64(self.loss_coefficient()).unwrap_or_else(T::one);
        let two = T::one() + T::one();
        let k_quad = k_factor * density / (two * area * area);

        Ok((r_lin, k_quad))
    }

    fn model_name(&self) -> &str {
        match self.junction_type {
            JunctionType::Tee => "T-Junction Loss (Idelchik 2007)",
            JunctionType::Cross => "Cross-Junction Loss (Idelchik 2007)",
        }
    }

    fn reynolds_range(&self) -> (T, T) {
        // Applicable from creeping flow to low-turbulence millifluidic range.
        (
            T::zero(),
            T::from_f64(2300.0).expect("Mathematical constant conversion compromised"),
        )
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        self.validate_mach_number(fluid, conditions)?;
        if self.branch_diameter_m <= 0.0 {
            return Err(Error::PhysicsViolation("branch_diameter_m must be positive".into()));
        }
        if self.branch_area_m2 <= 0.0 {
            return Err(Error::PhysicsViolation("branch_area_m2 must be positive".into()));
        }
        if self.junction_length_m < 0.0 {
            return Err(Error::PhysicsViolation("junction_length_m must be non-negative".into()));
        }
        if let Some(run_area) = self.run_area_m2 {
            if run_area <= 0.0 {
                return Err(Error::PhysicsViolation("run_area_m2 must be positive if specified".into()));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::resistance::models::FlowConditions;

    fn make_model(jtype: JunctionType, direction: JunctionFlowDirection) -> JunctionLossModel {
        JunctionLossModel {
            junction_type: jtype,
            flow_direction: direction,
            branch_diameter_m: 500e-6,
            branch_area_m2: std::f64::consts::PI * (250e-6_f64).powi(2),
            junction_length_m: 500e-6,
            density_kg_m3: 1000.0,
            run_area_m2: None,
        }
    }

    #[test]
    fn tee_dividing_k_factor() {
        let m = make_model(JunctionType::Tee, JunctionFlowDirection::Dividing);
        assert!((m.loss_coefficient() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn cross_combining_k_factor() {
        let m = make_model(JunctionType::Cross, JunctionFlowDirection::Combining);
        assert!((m.loss_coefficient() - 1.2).abs() < 1e-10);
    }

    #[test]
    fn linearised_resistance_positive() {
        let m = make_model(JunctionType::Cross, JunctionFlowDirection::Dividing);
        let r = m.linearised_resistance(1e-6);
        assert!(r > 0.0, "linearised resistance must be positive");
    }

    #[test]
    fn resistance_model_trait_works() {
        let m = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let cond = FlowConditions::new(0.01);
        let r = m.calculate_resistance(&fluid, &cond).unwrap();
        assert!(r > 0.0, "resistance must be positive: got {r}");
    }

    #[test]
    fn coefficients_have_quadratic_term() {
        let m = make_model(JunctionType::Cross, JunctionFlowDirection::Dividing);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let cond = FlowConditions::new(0.01);
        let (r, k) = m.calculate_coefficients(&fluid, &cond).unwrap();
        assert!(r > 0.0, "linear resistance must be positive");
        assert!(k > 0.0, "quadratic loss coefficient must be positive");
    }

    #[test]
    fn combining_higher_loss_than_dividing_for_tee() {
        let div = make_model(JunctionType::Tee, JunctionFlowDirection::Dividing);
        let comb = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
        assert!(
            comb.loss_coefficient() > div.loss_coefficient(),
            "combining flow typically has higher losses than dividing"
        );
    }

    #[test]
    fn diameter_ratio_correction_narrows_branch_reduces_k() {
        let equal = make_model(JunctionType::Cross, JunctionFlowDirection::Dividing);
        let k_equal = equal.loss_coefficient();
        let narrow_branch = JunctionLossModel {
            run_area_m2: Some(equal.branch_area_m2 * 4.0), // run 4× larger than branch
            ..equal
        };
        let k_narrow = narrow_branch.loss_coefficient();
        assert!(
            k_narrow < k_equal,
            "narrow branch into wide run should reduce K: {k_narrow} vs {k_equal}"
        );
        // area ratio = 0.25, sqrt = 0.5 → K_narrow = K_equal × 0.5
        assert!((k_narrow - k_equal * 0.5).abs() < 1e-10);
    }

    #[test]
    fn diameter_ratio_none_gives_equal_area_k() {
        let m = make_model(JunctionType::Tee, JunctionFlowDirection::Dividing);
        assert!(m.run_area_m2.is_none());
        assert!((m.loss_coefficient() - 1.0).abs() < 1e-10);
    }
    
    #[test]
    fn invalid_geometry_rejected() {
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let cond = FlowConditions::new(0.01);
        
        let mut m = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
        m.branch_area_m2 = -1e-6;
        assert!(m.validate_invariants(&fluid, &cond).is_err());
        
        let mut m2 = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
        m2.run_area_m2 = Some(-1e-6);
        assert!(m2.validate_invariants(&fluid, &cond).is_err());
    }
}
