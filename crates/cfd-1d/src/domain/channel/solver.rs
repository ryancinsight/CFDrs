//! Channel flow solver implementations
//!
//! ## Theorem: 1D Navier-Stokes Conservation Equations
//!
//! **Theorem**: For an incompressible Newtonian fluid in a rigid 1D channel,
//! the Navier-Stokes equations integrate over the cross-section $A$ to yield the
//! 1D conservation laws:
//!
//! **Mass Conservation (Continuity)**:
//! $$ \frac{\partial Q}{\partial x} = 0 $$
//! This implies volumetric flow rate $Q$ is constant along any unbranched channel segment.
//!
//! **Momentum Conservation**:
//! $$ \rho \frac{\partial \bar{u}}{\partial t} = -\frac{\partial P}{\partial x} - \frac{\tau_w \mathcal{P}}{A} $$
//! where $\tau_w$ is the wall shear stress and $\mathcal{P}$ is the wetted perimeter.
//!
//! **Steady-State Pressure Drop (Darcy-Weisbach)**:
//! $$ \Delta P = f \frac{L}{D_h} \frac{\rho \bar{u}^2}{2} \quad R = \frac{\Delta P}{Q} $$
//!
//! ## Theorem: Laminar Shape Factor (Shah-London 1978, Table 48)
//!
//! For a duct of arbitrary cross-section, the laminar friction factor satisfies:
//! $$ f \cdot Re_{D_h} = \text{Po} $$
//! The **Poiseuille number** Po depends only on cross-section geometry:
//! - Circular: `Po = 64`
//! - Square (AR = 1): `Po = 56.9`
//! - Infinite parallel plates (AR → ∞): `Po = 96`
//! - Rectangle AR = a/b ≥ 1: `Po = 96(1 − 1.3553/α + 1.9467/α² − 1.7012/α³ + 0.9564/α⁴ − 0.2537/α⁵)`
//!   where α = a/b ≥ 1 (this is the Shah-London 5-term polynomial).
//!
//! ## Theorem: Knudsen Slip Correction (Beskok-Karniadakis 1999)
//!
//! For `0.001 ≤ Kn ≤ 0.1`:
//! $$ R_{slip} = R_{laminar} / (1 + 4 \cdot Kn \cdot \sigma_v) $$
//! where $\sigma_v \approx 1.0$ is the tangential momentum accommodation coefficient (TMAC).
//! The factor `4·Kn` comes from the first-order slip boundary condition:
//! $u_{slip} = (2-\sigma_v)/\sigma_v \cdot Kn \cdot du/dn$   (Maxwell 1879).

use super::cross_section::CrossSection;
use super::flow::{Channel, FlowRegime, FlowState, NumericalParameters};
use super::geometry::ChannelGeometry;
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Result;
use cfd_core::physics::constants::mathematical::{numeric, PI};
use cfd_core::physics::constants::physics::dimensionless::reynolds::{
    PIPE_LAMINAR_MAX, PIPE_TURBULENT_MIN,
};
use cfd_core::physics::fluid::{ConstantFluid, ConstantPropertyFluid};
use nalgebra::RealField;
use num_traits::{cast::FromPrimitive, Float};

impl<T: RealField + Copy + FromPrimitive + Float> ChannelGeometry<T> {
    /// Create a rectangular channel geometry
    pub fn rectangular(length: T, width: T, height: T, roughness: T) -> Self {
        use super::surface::{SurfaceProperties, Wettability};
        Self {
            channel_type: super::geometry::ChannelType::Straight,
            length,
            cross_section: CrossSection::Rectangular { width, height },
            surface: SurfaceProperties {
                roughness,
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        }
    }

    /// Create a circular channel geometry
    pub fn circular(length: T, diameter: T, roughness: T) -> Self {
        use super::surface::{SurfaceProperties, Wettability};
        Self {
            channel_type: super::geometry::ChannelType::Straight,
            length,
            cross_section: CrossSection::Circular { diameter },
            surface: SurfaceProperties {
                roughness,
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        }
    }

    /// Get cross-sectional area
    pub fn area(&self) -> T {
        match &self.cross_section {
            CrossSection::Rectangular { width, height } => *width * *height,
            CrossSection::Circular { diameter } => {
                let pi = T::from_f64_or_zero(PI);
                let two = T::from_f64_or_zero(numeric::TWO);
                let radius = *diameter / two;
                pi * radius * radius
            }
            CrossSection::Elliptical {
                major_axis,
                minor_axis,
            } => {
                let pi = T::from_f64_or_zero(PI);
                let four = T::from_f64_or_zero(numeric::FOUR);
                pi * *major_axis * *minor_axis / four
            }
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => (*top_width + *bottom_width) * *height / (T::one() + T::one()),
            CrossSection::Custom { area, .. } => *area,
        }
    }

    /// Get hydraulic diameter
    pub fn hydraulic_diameter(&self) -> T {
        match &self.cross_section {
            CrossSection::Rectangular { width, height } => {
                let four = T::one() + T::one() + T::one() + T::one();
                four * self.area() / ((T::one() + T::one()) * (*width + *height))
            }
            CrossSection::Circular { diameter } => *diameter,
            CrossSection::Elliptical { .. } => {
                // Use definition Dh = 4 A / P with Ramanujan perimeter
                let four = T::one() + T::one() + T::one() + T::one();
                four * self.area() / self.wetted_perimeter()
            }
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => {
                let area = self.area();
                let hw = (*top_width - *bottom_width) / (T::one() + T::one());
                let side_length = Float::sqrt(Float::powi(*height, 2) + Float::powi(hw, 2));
                let perimeter = *top_width + *bottom_width + (T::one() + T::one()) * side_length;
                (T::one() + T::one() + T::one() + T::one()) * area / perimeter
            }
            CrossSection::Custom {
                hydraulic_diameter, ..
            } => *hydraulic_diameter,
        }
    }

    /// Get wetted perimeter
    pub fn wetted_perimeter(&self) -> T {
        match &self.cross_section {
            CrossSection::Rectangular { width, height } => {
                (T::one() + T::one()) * (*width + *height)
            }
            CrossSection::Circular { diameter } => {
                let pi = T::pi();
                pi * *diameter
            }
            CrossSection::Elliptical {
                major_axis,
                minor_axis,
            } => {
                // Exact formula using Arithmetic-Geometric Mean (AGM) for Complete Elliptic Integral of the Second Kind
                let pi = T::pi();
                let two = T::one() + T::one();
                let a_val = *major_axis / two;
                let b_val = *minor_axis / two;

                // Ensure a >= b for standard integral form
                let (a, b) = if a_val > b_val {
                    (a_val, b_val)
                } else {
                    (b_val, a_val)
                };

                if a == b || b == T::zero() {
                    return two * pi * a;
                }

                let m = T::one() - (b * b) / (a * a);

                let mut a_n = T::one();
                let mut b_n = Float::sqrt(T::one() - m);
                let mut c_n = Float::sqrt(m);

                let mut sum = c_n * c_n / two;
                let mut power = T::one();

                let tolerance =
                    T::from_f64(1e-14).expect("Mathematical constant conversion compromised");

                for _ in 0..20 {
                    let a_next = (a_n + b_n) / two;
                    let b_next = Float::sqrt(a_n * b_n);
                    let c_next = (a_n - b_n) / two;

                    a_n = a_next;
                    b_n = b_next;
                    c_n = c_next;

                    sum += power * c_n * c_n;
                    power *= two;

                    if c_n < tolerance || c_n == T::zero() {
                        break;
                    }
                }

                let e_m = (pi / (two * a_n)) * (T::one() - sum);
                let four = T::one() + T::one() + T::one() + T::one();
                four * a * e_m
            }
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => {
                let hw = (*top_width - *bottom_width) / (T::one() + T::one());
                let side_length = Float::sqrt(Float::powi(*height, 2) + Float::powi(hw, 2));
                *top_width + *bottom_width + (T::one() + T::one()) * side_length
            }
            CrossSection::Custom {
                area,
                hydraulic_diameter,
            } => (T::one() + T::one() + T::one() + T::one()) * *area / *hydraulic_diameter,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> Channel<T> {
    /// Create a new channel with geometry
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

    /// Calculate hydraulic resistance using physics models
    ///
    /// # Errors
    /// Returns an error if flow state calculation or resistance computation fails
    pub fn calculate_resistance(&mut self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        // Update flow state
        self.update_flow_state(fluid)?;

        // Calculate resistance based on flow regime
        match self.flow_state.flow_regime {
            FlowRegime::Stokes => self.calculate_stokes_resistance(fluid),
            FlowRegime::Laminar => self.calculate_laminar_resistance(fluid),
            FlowRegime::Transitional => self.calculate_transitional_resistance(fluid),
            FlowRegime::Turbulent => self.calculate_turbulent_resistance(fluid),
            FlowRegime::SlipFlow => self.calculate_slip_flow_resistance(fluid),
        }
    }

    /// Update flow state based on current conditions.
    ///
    /// Computes the Knudsen number from the fluid's speed of sound and thermodynamic
    /// state using the Chapman-Enskog mean free path equation (White 2006, §1.7):
    ///
    /// ```text
    /// λ = μ · √(π R_gas T / 2 M) / P   [m]
    /// Kn = λ / Dh
    /// ```
    ///
    /// For liquid-phase millifluidics, Kn << 0.001 (continuum limit).
    fn update_flow_state(&mut self, fluid: &ConstantPropertyFluid<T>) -> Result<()> {
        // Compute hydraulic diameter
        let dh = self.geometry.hydraulic_diameter();

        // Compute Knudsen number from fluid properties
        // For ideal gas: λ = μ / (ρ · √(2RT/πM)) using speed of sound c = √(γRT/M)
        // Simplified: λ ≈ μ / (ρ · c / √(π/2·γ))
        let kn_opt = if dh > T::zero() {
            // Chapman-Enskog: λ = μ √(π/2) / (ρ c)
            // c = speed of sound = √(γ·P/ρ) ≈ bulk speed_of_sound field in Fluid
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

        // Classify flow regime, giving Kn priority over Re for slip conditions
        if let Some(re) = self.flow_state.reynolds_number {
            self.flow_state.flow_regime = if let Some(kn) = kn_opt {
                FlowRegime::classify_with_knudsen(re, kn)
            } else {
                FlowRegime::from_reynolds_number(re)
            };

            // Check for entrance effects using Schlichting correlations
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
                    // Slip flow — entrance effects are weaker; use laminar correlation
                    dh * T::from_f64(0.06).expect("Mathematical constant conversion compromised")
                        * re
                }
            };
            self.flow_state.entrance_effects = self.geometry.length < entrance_length;
        }

        Ok(())
    }

    /// Calculate Stokes flow (creeping flow) resistance.
    ///
    /// ## Theorem: Stokes Flow Resistance (Re < 1)
    ///
    /// For Re < 1, inertial terms are negligible. The momentum equation reduces to:
    /// $$ \nabla P = \mu \nabla^2 \mathbf{u} $$
    /// For a duct of cross-section with Poiseuille number Po = f·Re:
    /// $$ R = \frac{\text{Po} \cdot \mu \cdot L}{2 \cdot A \cdot D_h^2} $$
    ///
    /// **Validity**: Re < 1 (Stokes regime), any cross-section shape.
    /// **Reference**: White, F. M. (2011). *Fluid Mechanics*, 7th ed., §6.4.
    fn calculate_stokes_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        let area = self.geometry.area();
        let dh = self.geometry.hydraulic_diameter();
        let length = self.geometry.length;
        let viscosity = fluid.dynamic_viscosity();
        let shape_factor = self.get_shape_factor();
        let resistance =
            shape_factor * viscosity * length / ((T::one() + T::one()) * area * dh * dh);
        Ok(resistance)
    }

    /// Calculate laminar flow resistance.
    ///
    /// ## Theorem: Hagen-Poiseuille Resistance (1 ≤ Re < 2300)
    ///
    /// For fully-developed laminar flow in a duct:
    /// $$ R = \frac{\text{Po} \cdot \mu \cdot L}{2 \cdot A \cdot D_h^2} $$
    /// where Po = f·Re is the Poiseuille number depending only on cross-section shape.
    ///
    /// For a circular pipe, Po = 64, giving the classical Hagen-Poiseuille formula:
    /// $$ R = \frac{128 \mu L}{\pi D^4} = \frac{64 \mu L}{2 A D_h^2} \quad (A = \pi D^2/4, D_h = D) $$
    ///
    /// **Validity**: 1 ≤ Re < 2300, fully-developed, incompressible, Newtonian.
    /// **Reference**: Shah, R. K. & London, A. L. (1978). *Laminar Flow Forced Convection*.
    fn calculate_laminar_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        let area = self.geometry.area();
        let dh = self.geometry.hydraulic_diameter();
        let length = self.geometry.length;
        let viscosity = fluid.dynamic_viscosity();
        let shape_factor = self.get_shape_factor();
        let resistance =
            shape_factor * viscosity * length / ((T::one() + T::one()) * area * dh * dh);
        Ok(resistance)
    }

    /// Calculate transitional flow resistance.
    ///
    /// ## Theorem: Transitional Blending (2300 ≤ Re ≤ 4000)
    ///
    /// No universal analytical closure exists in the transitional regime.
    /// We use linear blending between laminar and turbulent resistances:
    ///
    /// $$ R_{trans} = R_L \cdot (1 - w) + R_T \cdot w, \quad w = \frac{Re - Re_L}{Re_T - Re_L} $$
    ///
    /// where `Re_L = 2300` and `Re_T = 4000`.
    /// This is consistent with the Moody diagram interpolation used in
    /// engineering practice (White 2011, §6.6).
    ///
    /// **Validity**: 2300 ≤ Re ≤ 4000 (transitional region).
    fn calculate_transitional_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        let re = self.flow_state.reynolds_number.ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration(
                "Reynolds number required for transitional resistance".to_string(),
            )
        })?;

        // Linear interpolation between laminar and turbulent at transition bounds
        let re_l =
            T::from_f64(PIPE_LAMINAR_MAX).expect("Mathematical constant conversion compromised");
        let re_u =
            T::from_f64(PIPE_TURBULENT_MIN).expect("Mathematical constant conversion compromised");

        let r_l = self.calculate_laminar_resistance(fluid)?;
        let r_u = self.calculate_turbulent_resistance(fluid)?;

        let weight = (re - re_l) / (re_u - re_l);
        Ok(r_l * (T::one() - weight) + r_u * weight)
    }

    /// Calculate turbulent flow resistance.
    ///
    /// ## Theorem: Darcy-Weisbach Turbulent Resistance (Re > 4000)
    ///
    /// The Haaland (1983) approximation for the Darcy friction factor:
    ///
    /// $$ \frac{1}{\sqrt{f}} = -1.8 \log_{10}\!\left[\left(\frac{\varepsilon/D_h}{3.7}\right)^{1.11} + \frac{6.9}{Re}\right] $$
    ///
    /// Accuracy: within 2% of the Colebrook-White equation for `Re ≥ 4000`.
    ///
    /// The effective linearized resistance for the nonlinear ΔP = k·Q²:
    /// $$ R_{eff} = k \cdot |Q| = \frac{f \rho L |V|}{2 D_h A} $$
    ///
    /// **Reference**: Haaland, S. E. (1983). J. Fluids Eng. 105(1), 89–90.
    fn calculate_turbulent_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        let re = self.flow_state.reynolds_number.ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration(
                "Reynolds number required for turbulent resistance".to_string(),
            )
        })?;

        let area = self.geometry.area();
        let dh = self.geometry.hydraulic_diameter();
        let length = self.geometry.length;
        let density = fluid.density;
        let f = self.calculate_friction_factor(re);

        // R_eff = f·ρ·L·|V| / (2·Dh·A),  V = Re·μ / (ρ·Dh)
        let viscosity = fluid.dynamic_viscosity();
        let velocity = (re * viscosity) / (density * dh);
        let resistance = (f * density * length * velocity) / ((T::one() + T::one()) * dh * area);
        Ok(resistance)
    }

    /// Calculate slip flow resistance using the Beskok-Karniadakis correction.
    ///
    /// ## Theorem: Beskok-Karniadakis Slip Correction (1999)
    ///
    /// For `0.001 ≤ Kn < 0.1`, the velocity slip at the wall reduces resistance:
    ///
    /// $$ R_{slip} = \frac{R_{laminar}}{1 + 4 \cdot Kn \cdot \sigma_v} $$
    ///
    /// where $\sigma_v \approx 1.0$ is the tangential momentum accommodation coefficient (TMAC).
    ///
    /// The Knudsen number `Kn = λ / Dh` comes from `FlowState::knudsen_number`,
    /// which is computed in `update_flow_state` via Chapman-Enskog kinetic theory:
    ///
    /// $$ \lambda = \frac{\mu\,\sqrt{\pi/2}}{\rho \, c} $$
    ///
    /// where `c` is the fluid speed of sound. For liquid-phase millifluidic channels
    /// (Dh ≈ 0.1–1 mm), Kn << 0.001, so the slip correction is negligible (< 0.4%).
    ///
    /// **References**:
    /// - Beskok, A. & Karniadakis, G. E. (1999). Microscale Thermophys. Eng. 3(1), 43–77.
    /// - Maxwell, J. C. (1879). Phil. Trans. R. Soc. Lond. 170, 231–256.
    fn calculate_slip_flow_resistance(&self, fluid: &ConstantPropertyFluid<T>) -> Result<T> {
        let r_laminar = self.calculate_laminar_resistance(fluid)?;

        // Use Kn from flow state (computed by update_flow_state from mean free path).
        // If Kn is unavailable (pure liquid, no speed_of_sound data), fall back to
        // laminar resistance (Kn → 0 → slip correction = 1).
        let kn = self.flow_state.knudsen_number.unwrap_or_else(T::zero);

        // Beskok-Karniadakis first-order slip: σ_v = 1 (air-like TMAC)
        // R_slip = R_lam / (1 + 4·Kn·σ_v)
        let four = T::one() + T::one() + T::one() + T::one();
        Ok(r_laminar / (T::one() + four * kn))
    }

    fn calculate_friction_factor(&self, reynolds: T) -> T {
        let re_val = reynolds.to_f64().unwrap_or(0.0);
        let roughness = self.geometry.surface.roughness;
        let dh = self.geometry.hydraulic_diameter();
        let roughness_ratio = (roughness / dh).to_f64().unwrap_or(0.0);

        if re_val < 0.1 {
            return T::from_f64(64.0 / 0.1).expect("Mathematical constant conversion compromised");
        }

        // Haaland approximation for Darcy friction factor
        let factor = -1.8 * f64::log10(f64::powf(roughness_ratio / 3.7, 1.11) + 6.9 / re_val);
        let f = 1.0 / (factor * factor);
        T::from_f64(f).expect("Mathematical constant conversion compromised")
    }

    fn get_shape_factor(&self) -> T {
        // Poiseuille number Po = f·Re for laminar flow in various cross-sections.
        // Source: Shah & London (1978), Table 48.
        match &self.geometry.cross_section {
            CrossSection::Circular { .. } => {
                // Po = 64 (exact, Hagen-Poiseuille)
                T::from_f64(64.0).expect("Mathematical constant conversion compromised")
            }
            CrossSection::Rectangular { width, height } => {
                // Shah-London 5-term polynomial (Shah & London 1978, Table 48):
                // Po = 96(1 - 1.3553/α + 1.9467/α² - 1.7012/α³ + 0.9564/α⁴ - 0.2537/α⁵)
                // where α = a/b ≥ 1, a = longer side, b = shorter side
                let (a, b) = if *width >= *height {
                    (*width, *height)
                } else {
                    (*height, *width)
                };
                // Guard against degenerate aspect ratio
                if b <= T::zero() {
                    return T::from_f64(64.0)
                        .expect("Mathematical constant conversion compromised");
                }
                let alpha = a / b;
                let c1 = T::from_f64(1.3553).expect("Mathematical constant conversion compromised");
                let c2 = T::from_f64(1.9467).expect("Mathematical constant conversion compromised");
                let c3 = T::from_f64(1.7012).expect("Mathematical constant conversion compromised");
                let c4 = T::from_f64(0.9564).expect("Mathematical constant conversion compromised");
                let c5 = T::from_f64(0.2537).expect("Mathematical constant conversion compromised");
                let inv_a = T::one() / alpha;
                let inv_a2 = inv_a * inv_a;
                let inv_a3 = inv_a2 * inv_a;
                let inv_a4 = inv_a3 * inv_a;
                let inv_a5 = inv_a4 * inv_a;

                T::from_f64(96.0).expect("Mathematical constant conversion compromised")
                    * (T::one() - c1 * inv_a + c2 * inv_a2 - c3 * inv_a3 + c4 * inv_a4
                        - c5 * inv_a5)
            }
            CrossSection::Elliptical {
                major_axis,
                minor_axis,
            } => {
                // For an ellipse with semi-axes a, b (a ≥ b):
                // Exact Stokes solution: Po = 2π(a² + b²) / (a·b)
                // (Dryden, Murnaghan & Bateman, 1932)
                let two_pi = T::from_f64(2.0 * std::f64::consts::PI)
                    .expect("Mathematical constant conversion compromised");
                let a = *major_axis / (T::one() + T::one());
                let b = *minor_axis / (T::one() + T::one());
                if a <= T::zero() || b <= T::zero() {
                    return T::from_f64(64.0)
                        .expect("Mathematical constant conversion compromised");
                }
                two_pi * (a * a + b * b) / (a * b)
            }
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => {
                // For trapezoidal channels, use the equivalent rectangular approach:
                // The Shah-London polynomial applied to the effective aspect ratio
                // computed from the trapezoid's area and wetted perimeter.
                //
                // Effective half-width: w_eff = A / h = (top + bottom)/2
                // Effective aspect ratio: α_eff = w_eff / height (if w_eff >= height)
                //
                // Reference: Muzychka, Y. S. & Yovanovich, M. M. (2002).
                // ASME J. Heat Transfer, 124(2), 260–267 — equivalent rectangle method.
                let avg_width = (*top_width + *bottom_width) / (T::one() + T::one());
                if *height <= T::zero() || avg_width <= T::zero() {
                    return T::from_f64(64.0)
                        .expect("Mathematical constant conversion compromised");
                }
                let (a, b) = if avg_width >= *height {
                    (avg_width, *height)
                } else {
                    (*height, avg_width)
                };
                let alpha = a / b;
                let c1 = T::from_f64(1.3553).expect("Mathematical constant conversion compromised");
                let c2 = T::from_f64(1.9467).expect("Mathematical constant conversion compromised");
                let c3 = T::from_f64(1.7012).expect("Mathematical constant conversion compromised");
                let c4 = T::from_f64(0.9564).expect("Mathematical constant conversion compromised");
                let c5 = T::from_f64(0.2537).expect("Mathematical constant conversion compromised");
                let inv_a = T::one() / alpha;
                let inv_a2 = inv_a * inv_a;
                let inv_a3 = inv_a2 * inv_a;
                let inv_a4 = inv_a3 * inv_a;
                let inv_a5 = inv_a4 * inv_a;
                T::from_f64(96.0).expect("Mathematical constant conversion compromised")
                    * (T::one() - c1 * inv_a + c2 * inv_a2 - c3 * inv_a3 + c4 * inv_a4
                        - c5 * inv_a5)
            }
            CrossSection::Custom { .. } => {
                // For custom cross-sections, assume circular (Po = 64).
                // Users should supply a more accurate shape factor via a custom
                // ResistanceModel if the geometry differs significantly.
                T::from_f64(64.0).expect("Mathematical constant conversion compromised")
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_circular_shape_factor_is_64() {
        let geom = ChannelGeometry::circular(0.01_f64, 0.001, 0.0);
        let channel = Channel::new(geom);
        // get_shape_factor is private, so test via Stokes/laminar resistance ratio.
        // For circular: Po = 64. We verify indirectly by checking resistance formula:
        // R = Po * mu * L / (2 * A * Dh^2)
        // For a circular pipe: A = pi*D^2/4, Dh = D
        // R = 64 * mu * L / (2 * pi*D^2/4 * D^2) = 128 * mu * L / (pi * D^4)
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let mut ch = channel;
        ch.flow_state.reynolds_number = Some(100.0); // laminar regime
        ch.flow_state.flow_regime = FlowRegime::Laminar;
        let r = ch.calculate_laminar_resistance(&fluid).unwrap();

        let d = 0.001_f64;
        let l = 0.01_f64;
        let mu = fluid.dynamic_viscosity();
        let expected = 128.0 * mu * l / (std::f64::consts::PI * d.powi(4));
        assert_relative_eq!(r, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_flow_regime_stokes_below_1() {
        let regime = FlowRegime::from_reynolds_number(0.5_f64);
        assert_eq!(regime, FlowRegime::Stokes);
    }

    #[test]
    fn test_flow_regime_laminar_below_2300() {
        let regime = FlowRegime::from_reynolds_number(500.0_f64);
        assert_eq!(regime, FlowRegime::Laminar);
    }

    #[test]
    fn test_flow_regime_transitional() {
        let regime = FlowRegime::from_reynolds_number(3000.0_f64);
        assert_eq!(regime, FlowRegime::Transitional);
    }

    #[test]
    fn test_flow_regime_turbulent_above_4000() {
        let regime = FlowRegime::from_reynolds_number(5000.0_f64);
        assert_eq!(regime, FlowRegime::Turbulent);
    }

    #[test]
    fn test_square_channel_shape_factor_approx_56_9() {
        // Square channel (AR=1): Po should be ~56.9 per Shah-London
        let geom = ChannelGeometry::rectangular(0.01_f64, 0.001, 0.001, 0.0);
        let channel = Channel::new(geom);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let mut ch = channel;
        ch.flow_state.reynolds_number = Some(100.0);
        ch.flow_state.flow_regime = FlowRegime::Laminar;
        let r = ch.calculate_laminar_resistance(&fluid).unwrap();

        // R = Po * mu * L / (2 * A * Dh^2)
        let w = 0.001_f64;
        let h = 0.001_f64;
        let a = w * h;
        let dh = 4.0 * a / (2.0 * (w + h));
        let mu = fluid.dynamic_viscosity();
        let l = 0.01_f64;
        // Extract Po from the measured resistance
        let po = r * 2.0 * a * dh * dh / (mu * l);
        assert_relative_eq!(po, 56.908, epsilon = 0.1);
    }
}
