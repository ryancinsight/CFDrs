//! Channel flow resistance computation.
//!
//! # Theorem: 1D Navier-Stokes Conservation Equations
//!
//! **Theorem**: For an incompressible Newtonian fluid in a rigid 1D channel,
//! the Navier-Stokes equations integrate over the cross-section `A` to yield the
//! 1D conservation laws:
//!
//! **Mass Conservation (Continuity)**:
//! ```text
//! ∂Q/∂x = 0
//! ```
//! This implies volumetric flow rate `Q` is constant along any unbranched channel segment.
//!
//! **Momentum Conservation**:
//! ```text
//! ρ ∂ū/∂t = −∂P/∂x − τ_w P_w / A
//! ```
//! where `τ_w` is the wall shear stress and `P_w` is the wetted perimeter.
//!
//! **Steady-State Pressure Drop (Darcy-Weisbach)**:
//! ```text
//! ΔP = f · (L / D_h) · (ρ ū² / 2),   R = ΔP / Q
//! ```
//!
//! # Theorem: Knudsen Slip Correction (Beskok-Karniadakis 1999)
//!
//! For `0.001 ≤ Kn ≤ 0.1`:
//! ```text
//! R_slip = R_laminar / (1 + 4 · Kn · σ_v)
//! ```
//! where `σ_v ≈ 1.0` is the tangential momentum accommodation coefficient (TMAC).
//!
//! # References
//! - White, F. M. (2011). *Fluid Mechanics*, 7th ed., §6.4.
//! - Shah, R. K. & London, A. L. (1978). *Laminar Flow Forced Convection in Ducts*.
//! - Haaland, S. E. (1983). J. Fluids Eng. 105(1), 89–90.
//! - Beskok, A. & Karniadakis, G. E. (1999). Microscale Thermophys. Eng. 3(1), 43–77.
//! - Schlichting, H. (1979). *Boundary-Layer Theory*, 7th ed.

use super::shape_factors::poiseuille_number;
use crate::domain::channel::flow::{Channel, FlowRegime, NumericalParameters, FlowState};
use crate::domain::channel::geometry::ChannelGeometry;
use crate::physics::resistance::models::{DarcyWeisbachModel, FlowConditions, ResistanceModel};
use cfd_core::error::Result;
use cfd_core::physics::constants::physics::dimensionless::reynolds::{
    PIPE_LAMINAR_MAX, PIPE_TURBULENT_MIN,
};
use cfd_core::physics::fluid::{ConstantFluid, ConstantPropertyFluid};
use nalgebra::RealField;
use num_traits::{cast::FromPrimitive, Float};

impl<T: RealField + Copy + FromPrimitive + Float> Channel<T> {
    /// Create a new channel with geometry.
    pub fn new(geometry: ChannelGeometry<T>) -> Self {
        Self {
            geometry,
            flow_state: FlowState {
                reynolds_number: None,
                knudsen_number: None,
                flow_regime: FlowRegime::Laminar,
                entrance_effects: false,
                secondary_flows: false,
            },
            numerical_params: NumericalParameters {
                discretization_points: 100,
                tolerance: T::from_f64(1e-6).expect("Mathematical constant conversion compromised"),
                entrance_effects: false,
                surface_tension_effects: false,
            },
        }
    }

    /// Calculate hydraulic resistance using physics models.
    ///
    /// Dispatches to the appropriate resistance model based on the current
    /// flow regime (Stokes, laminar, transitional, turbulent, or slip).
    ///
    /// # Errors
    /// Returns an error if flow state calculation or resistance computation fails.
    pub fn calculate_resistance(&mut self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        self.update_flow_state(fluid)?;

        match self.flow_state.flow_regime {
            FlowRegime::Stokes => self.stokes_resistance(fluid),
            FlowRegime::Laminar => self.laminar_resistance(fluid),
            FlowRegime::Transitional => self.transitional_resistance(fluid),
            FlowRegime::Turbulent => self.turbulent_resistance(fluid),
            FlowRegime::SlipFlow => self.slip_flow_resistance(fluid),
        }
    }

    /// Update flow state based on current conditions.
    ///
    /// Computes Knudsen number via Chapman-Enskog mean free path (White 2006, §1.7):
    /// ```text
    /// λ = μ · √(π/2) / (ρ · c)   [m]
    /// Kn = λ / D_h
    /// ```
    fn update_flow_state(&mut self, fluid: &ConstantPropertyFluid<T>) -> Result<()> {
        let dh = self.geometry.hydraulic_diameter();

        let kn_opt = if dh > T::zero() {
            let sqrt_half_pi =
                T::from_f64(std::f64::consts::FRAC_PI_2.sqrt()).unwrap_or_else(T::one);
            let lam =
                fluid.dynamic_viscosity() * sqrt_half_pi / (fluid.density * fluid.speed_of_sound);
            if lam > T::zero() {
                Some(lam / dh)
            } else {
                None
            }
        } else {
            None
        };
        self.flow_state.knudsen_number = kn_opt;

        if let Some(re) = self.flow_state.reynolds_number {
            self.flow_state.flow_regime = if let Some(kn) = kn_opt {
                FlowRegime::classify_with_knudsen(re, kn)
            } else {
                FlowRegime::from_reynolds_number(re)
            };

            let entrance_length = match self.flow_state.flow_regime {
                FlowRegime::Laminar | FlowRegime::Stokes => {
                    // Schlichting (1979): L_e / Dh = 0.06 Re
                    dh * T::from_f64(0.06).expect("Mathematical constant conversion compromised")
                        * re
                }
                FlowRegime::Transitional | FlowRegime::Turbulent => {
                    // L_e / Dh = 4.4 Re^(1/6)
                    let one_sixth = T::from_f64(1.0 / 6.0)
                        .expect("Mathematical constant conversion compromised");
                    dh * T::from_f64(4.4).expect("Mathematical constant conversion compromised")
                        * Float::powf(re, one_sixth)
                }
                FlowRegime::SlipFlow => {
                    dh * T::from_f64(0.06).expect("Mathematical constant conversion compromised")
                        * re
                }
            };
            self.flow_state.entrance_effects = self.geometry.length < entrance_length;
        }

        Ok(())
    }

    /// Stokes flow resistance (Re < 1).
    ///
    /// ## Theorem
    /// For Re < 1, `R = Po · μ · L / (2 · A · D_h²)`.
    fn stokes_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        self.laminar_resistance(fluid)
    }

    /// Laminar flow resistance (1 ≤ Re < 2300).
    ///
    /// ## Theorem: Hagen-Poiseuille Resistance
    /// `R = Po · μ · L / (2 · A · D_h²)`
    /// where `Po = f·Re` is the Poiseuille number (shape-dependent).
    pub fn laminar_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        let area = self.geometry.area();
        let dh = self.geometry.hydraulic_diameter();
        let po = poiseuille_number(&self.geometry.cross_section);
        let resistance =
            po * fluid.dynamic_viscosity() * self.geometry.length
                / ((T::one() + T::one()) * area * dh * dh);
        Ok(resistance)
    }

    /// Transitional flow resistance (2300 ≤ Re ≤ 4000).
    ///
    /// ## Theorem: Transitional Blending (White 2011, §6.6)
    /// `R_trans = R_L · (1 − w) + R_T · w`,  `w = (Re − 2300) / (4000 − 2300)`
    fn transitional_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        let re = self.flow_state.reynolds_number.ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration(
                "Reynolds number required for transitional resistance".to_string(),
            )
        })?;

        let re_l =
            T::from_f64(PIPE_LAMINAR_MAX).expect("Mathematical constant conversion compromised");
        let re_u =
            T::from_f64(PIPE_TURBULENT_MIN).expect("Mathematical constant conversion compromised");

        let r_l = self.laminar_resistance(fluid)?;
        let r_u = self.turbulent_resistance(fluid)?;

        let weight = (re - re_l) / (re_u - re_l);
        Ok(r_l * (T::one() - weight) + r_u * weight)
    }

    /// Turbulent flow resistance (Re > 4000).
    ///
    /// Delegates to the exact Darcy-Weisbach resistance model, which solves the
    /// Colebrook-White equation with Newton-Raphson instead of using an explicit
    /// friction-factor approximation.
    fn turbulent_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        let re = self.flow_state.reynolds_number.ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration(
                "Reynolds number required for turbulent resistance".to_string(),
            )
        })?;

        let area = self.geometry.area();
        let dh = self.geometry.hydraulic_diameter();
        let density = fluid.density;
        let viscosity = fluid.dynamic_viscosity();
        let velocity = (re * viscosity) / (density * dh);
        let model = DarcyWeisbachModel::new(
            dh,
            area,
            self.geometry.length,
            self.geometry.surface.roughness,
        );
        let mut conditions = FlowConditions::new(velocity);
        conditions.reynolds_number = Some(re);
        model.calculate_resistance(fluid, &conditions)
    }

    /// Slip flow resistance (Beskok-Karniadakis 1999).
    ///
    /// ## Theorem
    /// `R_slip = R_laminar / (1 + 4 · Kn · σ_v)`, `σ_v = 1`.
    fn slip_flow_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        let r_laminar = self.laminar_resistance(fluid)?;
        let kn = self.flow_state.knudsen_number.unwrap_or_else(T::zero);
        let four = T::one() + T::one() + T::one() + T::one();
        Ok(r_laminar / (T::one() + four * kn))
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::channel::geometry::ChannelGeometry;
    use crate::physics::resistance::models::ResistanceModel;
    use approx::assert_relative_eq;

    fn water() -> ConstantPropertyFluid<f64> {
        ConstantPropertyFluid::new(
            "water".to_string(),
            1000.0, // density [kg/m³]
            0.001,  // viscosity [Pa·s]
            4186.0, // specific heat
            0.598,  // thermal conductivity
            1480.0, // speed of sound
        )
    }

    fn circular_channel(diameter: f64, length: f64) -> Channel<f64> {
        let geom = ChannelGeometry::circular(length, diameter, 1e-6);
        Channel::new(geom)
    }

    /// Laminar resistance for circular channel: R = Po·μ·L / (2·A·Dh²).
    /// Po = 64 (circular), D = 1 mm, L = 10 mm, μ = 0.001.
    /// R = 64 × 0.001 × 0.01 / (2 × π/4 × 10⁻⁶ × 10⁻⁶)
    ///   = 6.4e-4 / (π/2 × 10⁻¹²) = 128/(π) × 10⁸ ≈ 4.074 × 10⁸ Pa·s/m³.
    /// This equals the Hagen-Poiseuille formula: R = 128μL / (πD⁴).
    #[test]
    fn laminar_resistance_matches_hagen_poiseuille() {
        let d: f64 = 1e-3;
        let l: f64 = 0.01;
        let mu: f64 = 0.001;
        let chan = circular_channel(d, l);
        let r = chan.laminar_resistance(&water()).unwrap();
        let r_hp = 128.0 * mu * l / (std::f64::consts::PI * d.powi(4));
        assert_relative_eq!(r, r_hp, max_relative = 1e-10);
    }

    /// Turbulent resistance should agree with the exact Darcy-Weisbach model.
    #[test]
    fn turbulent_resistance_matches_exact_darcy_weisbach_model() {
        use crate::physics::resistance::models::DarcyWeisbachModel;

        let mut chan = circular_channel(1e-3, 0.01);
        chan.flow_state.reynolds_number = Some(10_000.0);
        let fluid = water();

        let d = chan.geometry.hydraulic_diameter();
        let area = chan.geometry.area();
        let viscosity = fluid.dynamic_viscosity();
        let density = fluid.density;
        let velocity = (10_000.0 * viscosity) / (density * d);
        let mut conditions = FlowConditions::new(velocity);
        conditions.reynolds_number = Some(10_000.0);

        let expected = DarcyWeisbachModel::new(d, area, chan.geometry.length, 1e-6)
            .calculate_resistance(&fluid, &conditions)
            .unwrap();
        let actual = chan.turbulent_resistance(&fluid).unwrap();

        assert_relative_eq!(actual, expected, max_relative = 1e-12, epsilon = 1e-12);
    }

    /// Stokes resistance (Re < 1) routes to same formula as laminar.
    #[test]
    fn stokes_delegates_to_laminar() {
        let chan = circular_channel(1e-3, 0.01);
        let r_stokes = chan.stokes_resistance(&water()).unwrap();
        let r_laminar = chan.laminar_resistance(&water()).unwrap();
        assert_relative_eq!(r_stokes, r_laminar, max_relative = 1e-15);
    }

    /// Slip flow correction: R_slip = R_lam / (1 + 4·Kn).
    /// For Kn=0.01: R_slip = R_lam / 1.04.
    #[test]
    fn slip_flow_reduces_resistance() {
        let mut chan = circular_channel(1e-3, 0.01);
        chan.flow_state.knudsen_number = Some(0.01);
        let r_slip = chan.slip_flow_resistance(&water()).unwrap();
        let r_laminar = chan.laminar_resistance(&water()).unwrap();
        let expected = r_laminar / (1.0 + 4.0 * 0.01);
        assert_relative_eq!(r_slip, expected, max_relative = 1e-12);
    }

    /// Resistance scales linearly with viscosity.
    #[test]
    fn resistance_scales_with_viscosity() {
        let chan = circular_channel(1e-3, 0.01);
        let water1 = ConstantPropertyFluid::new(
            "w1".to_string(), 1000.0, 0.001, 4186.0, 0.598, 1480.0,
        );
        let water2 = ConstantPropertyFluid::new(
            "w2".to_string(), 1000.0, 0.003, 4186.0, 0.598, 1480.0,
        );
        let r1 = chan.laminar_resistance(&water1).unwrap();
        let r2 = chan.laminar_resistance(&water2).unwrap();
        assert_relative_eq!(r2 / r1, 3.0, max_relative = 1e-12);
    }
}
