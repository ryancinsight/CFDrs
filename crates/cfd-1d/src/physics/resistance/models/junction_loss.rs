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
use aequitas::systems::si::quantities::{
    Area, Dimensionless, HydraulicResistance, Length, MassDensity, QuadraticHydraulicResistance,
    Velocity, VolumetricFlowRate,
};
use cfd_core::CfdScalar;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use eunomia::FloatElement;
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
/// resistance derived from $\Delta P = K \rho Q |Q| / (2 A^2)$, which adds to the
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
pub struct JunctionLossModel<T = f64> {
    /// Type of junction (T or cross).
    pub junction_type: JunctionType,
    /// Flow direction (dividing or combining).
    pub flow_direction: JunctionFlowDirection,
    /// Hydraulic diameter of the branch channel \[m].
    #[serde(rename = "branch_diameter_m")]
    pub branch_diameter: Length<T>,
    /// Cross-sectional area of the branch channel \[m²].
    #[serde(rename = "branch_area_m2")]
    pub branch_area: Area<T>,
    /// Length attributed to the junction (typically ≈ D_h) \[m].
    #[serde(rename = "junction_length_m")]
    pub junction_length: Length<T>,
    /// Fluid density [kg/m³].
    #[serde(rename = "density_kg_m3")]
    pub density: MassDensity<T>,
    /// Optional cross-sectional area of the run (main) channel \[m²].
    ///
    /// When provided, the K-factor is corrected for diameter mismatch
    /// using the Idelchik area-ratio correlation.  When `None`, the
    /// equal-area K-factor is used (backward-compatible).
    #[serde(default)]
    #[serde(rename = "run_area_m2")]
    pub run_area: Option<Area<T>>,
}

impl<T: CfdScalar> JunctionLossModel<T> {
    /// Retrieve the K-factor for this junction configuration.
    ///
    /// When `run_area` is set, the equal-area K is corrected by the
    /// area ratio `(A_branch / A_run)^0.5` (Idelchik 2007, §7.12).
    /// This captures the additional recirculation loss when a narrow
    /// branch intersects a wide run channel, or the reduced loss when
    /// a wide branch meets a narrow run.
    #[must_use]
    pub fn loss_coefficient(&self) -> Dimensionless<T> {
        let k_equal =
            <T as FloatElement>::from_f64(match (self.junction_type, self.flow_direction) {
                (JunctionType::Tee, JunctionFlowDirection::Dividing) => 1.0,
                (JunctionType::Tee, JunctionFlowDirection::Combining) => 1.5,
                (JunctionType::Cross, JunctionFlowDirection::Dividing) => 1.3,
                (JunctionType::Cross, JunctionFlowDirection::Combining) => 1.2,
            });

        let minimum_area = Area::from_base(<T as FloatElement>::from_f64(1e-30));
        let coefficient = self
            .run_area
            .filter(|area| *area > minimum_area)
            .map_or(k_equal, |area| {
                k_equal * (self.branch_area / area).into_base().sqrt()
            });
        Dimensionless::from_base(coefficient)
    }

    /// Compute the linearised resistance `K ρ |Q_ref| / A²` about a reference flow.
    #[must_use]
    pub fn linearised_resistance(
        &self,
        reference_flow: VolumetricFlowRate<T>,
    ) -> HydraulicResistance<T> {
        let area_squared = self.branch_area * self.branch_area;
        if area_squared.into_base() < <T as FloatElement>::from_f64(1e-30) {
            return HydraulicResistance::from_base(T::ZERO);
        }
        let gradient =
            minor_loss_coefficient(self.loss_coefficient(), self.density, self.branch_area);
        gradient * reference_flow.abs() * (T::ONE + T::ONE)
    }

    fn validate_geometry(&self) -> Result<()> {
        let required_positive = [
            ("branch diameter", self.branch_diameter.into_base()),
            ("branch area", self.branch_area.into_base()),
            ("reference density", self.density.into_base()),
        ];
        for (name, value) in required_positive {
            if !value.is_finite() || value <= T::ZERO {
                return Err(Error::PhysicsViolation(format!(
                    "Junction {name} must be finite and positive"
                )));
            }
        }
        let length = self.junction_length.into_base();
        if !length.is_finite() || length < T::ZERO {
            return Err(Error::PhysicsViolation(
                "Junction length must be finite and nonnegative".into(),
            ));
        }
        if let Some(area) = self.run_area
            && (!area.as_base().is_finite() || area.into_base() <= T::ZERO)
        {
            return Err(Error::PhysicsViolation(
                "Junction run area must be finite and positive".into(),
            ));
        }
        Ok(())
    }
}

/// The coefficient of `Q |Q|` in the signed pressure-drop law.
pub(crate) fn minor_loss_coefficient<T: CfdScalar>(
    loss: Dimensionless<T>,
    density: MassDensity<T>,
    area: Area<T>,
) -> QuadraticHydraulicResistance<T> {
    let two = T::ONE + T::ONE;
    loss * density / (area * two * area)
}

impl<T: CfdScalar> ResistanceModel<T> for JunctionLossModel<T> {
    /// Calculate hydraulic resistance [Pa·s/m³].
    ///
    /// Uses the linear + quadratic coefficient model; the linear part
    /// accounts for viscous losses in the short junction length while the
    /// quadratic part captures the K-factor minor loss.
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<HydraulicResistance<T>> {
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
    ) -> Result<(HydraulicResistance<T>, QuadraticHydraulicResistance<T>)> {
        self.validate_geometry()?;
        let state = fluid.properties_at(
            conditions.temperature.into_base(),
            conditions.pressure.into_base(),
        )?;
        let area = self.branch_area;
        let velocity = conditions
            .velocity
            .or_else(|| conditions.flow_rate.map(|flow| flow / area))
            .unwrap_or(Velocity::from_base(T::ZERO));

        let diameter = self.branch_diameter;
        let eight = <T as FloatElement>::from_f64(8.0);
        let shear_rate = if let Some(shear_rate) = conditions.shear_rate {
            if shear_rate.into_base() < T::ZERO {
                return Err(Error::PhysicsViolation(
                    "Junction wall shear rate must be nonnegative".to_string(),
                ));
            }
            shear_rate
        } else {
            velocity.abs() * eight / diameter
        };

        let viscosity = fluid.viscosity_at_shear(
            shear_rate.into_base(),
            conditions.temperature.into_base(),
            conditions.pressure.into_base(),
        )?;
        let diameter_squared = diameter * diameter;
        let resistance: HydraulicResistance<T> =
            viscosity * self.junction_length * <T as FloatElement>::from_f64(128.0)
                / (diameter_squared * diameter_squared * T::pi());
        let quadratic = minor_loss_coefficient(self.loss_coefficient(), state.density, area);
        Ok((resistance, quadratic))
    }

    fn model_name(&self) -> &str {
        match self.junction_type {
            JunctionType::Tee => "T-Junction Loss (Idelchik 2007)",
            JunctionType::Cross => "Cross-Junction Loss (Idelchik 2007)",
        }
    }

    fn reynolds_range(&self) -> (T, T) {
        // Applicable from creeping flow to low-turbulence millifluidic range.
        (T::ZERO, <T as FloatElement>::from_f64(2300.0))
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        self.validate_geometry()?;
        self.validate_mach_number(fluid, conditions)
    }
}

#[cfg(test)]
mod tests;
