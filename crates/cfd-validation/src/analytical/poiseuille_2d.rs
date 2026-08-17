//! Non-Newtonian Poiseuille flow in 2D channels
//!
//! Analytical and semi-analytical solutions for non-Newtonian fluid flow
//! between parallel plates, including power-law and Casson models.
//!
//! # Mathematical Foundation
//!
//! ## Power-Law Fluid
//! For a power-law fluid: τ = K·γ̇ⁿ in a 2D channel (width 2H):
//! ```text
//! u(y) = u_c · [1 - |y/H|^((n+1)/n)]^(n/(n+1))
//! ```
//! where:
//! - u_c = centerline velocity
//! - n = power-law index (n<1: shear-thinning, n>1: shear-thickening, n=1: Newtonian)
//! - K = consistency index [Pa·sⁿ]
//!
//! ## Casson Blood Model
//! For Casson fluid: √τ = √τ_y + √(μ_∞·γ̇)
//! The velocity profile must be solved numerically, but satisfies:
//! ```text
//! γ̇(y) = -du/dy > 0 for y > 0
//! τ(y) = (dp/dx)·y (linear stress distribution)
//! ```
//!
//! # References
//! - Bird, R.B., Stewart, W.E., Lightfoot, E.N. (2002) "Transport Phenomena"
//! - Chhabra, R.P., Richardson, J.F. (2008) "Non-Newtonian Flow and Applied Rheology"
//! - Merrill, E.W. et al. (1969) "Rheology of blood"

use super::AnalyticalSolution;
use crate::scalar;
use aequitas::systems::si::quantities::{
    AreaPerTime, Dimensionless, DynamicViscosity, Length, MassDensity, Pressure, PressureGradient,
    ReciprocalTime, Velocity,
};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::CfdScalar;
use eunomia::FloatElement;
use eunomia::RealField;
use leto::geometry::Vector3;

// ============================================================================
// Rheological Model Trait
// ============================================================================

/// Trait for rheological models in Poiseuille flow
pub trait RheologicalModel<T: RealField + Copy> {
    /// Compute shear stress from shear rate
    fn shear_stress(&self, shear_rate: ReciprocalTime<T>) -> Pressure<T>;

    /// Compute viscosity from shear rate
    fn viscosity(&self, shear_rate: ReciprocalTime<T>) -> DynamicViscosity<T>;

    /// Name of the rheological model
    fn model_name(&self) -> &str;
}

// ============================================================================
// Power-Law Model
// ============================================================================

/// Power-law rheological model: τ = K·γ̇ⁿ
///
/// The consistency coefficient has exponent-dependent SI units `Pa·sⁿ`.
/// Since the exponent is runtime data, Aequitas cannot represent that unit as
/// one fixed dimension; the coefficient is therefore confined to the
/// formula-bound rheology parameter newtype below.
#[derive(Debug, Clone, Copy)]
pub struct PowerLawModel<T: RealField + Copy> {
    /// Consistency index K [Pa·sⁿ]
    pub consistency: PowerLawConsistency<T>,
    /// Power-law index n (dimensionless)
    pub index: Dimensionless<T>,
}

/// Exponent-dependent power-law consistency coefficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PowerLawConsistency<T: RealField + Copy>(T);

impl<T: RealField + Copy> PowerLawConsistency<T> {
    /// Construct a consistency coefficient from its formula-bound value.
    #[must_use]
    pub fn new(value: T) -> Self {
        Self(value)
    }

    fn into_base(self) -> T {
        self.0
    }
}

impl<T: RealField + FloatElement + Copy> PowerLawModel<T> {
    /// Create power-law model
    ///
    /// # Arguments
    /// * `consistency` - K [Pa·sⁿ]
    /// * `index` - n (n<1: shear-thinning, n=1: Newtonian, n>1: shear-thickening)
    pub fn new(consistency: PowerLawConsistency<T>, index: Dimensionless<T>) -> Self {
        Self { consistency, index }
    }

    /// Create Newtonian model (n=1)
    pub fn newtonian(viscosity: DynamicViscosity<T>) -> Self {
        Self {
            consistency: PowerLawConsistency::new(viscosity.into_base()),
            index: Dimensionless::from_base(scalar::one::<T>()),
        }
    }

    /// Create typical blood-like shear-thinning model (n=0.6)
    pub fn blood_like(consistency: PowerLawConsistency<T>) -> Self {
        Self {
            consistency,
            index: Dimensionless::from_base(<T as FloatElement>::from_f64(0.6)),
        }
    }
}

impl<T: RealField + FloatElement + Copy> RheologicalModel<T> for PowerLawModel<T> {
    fn shear_stress(&self, shear_rate: ReciprocalTime<T>) -> Pressure<T> {
        let shear_rate = shear_rate.into_base();
        let index = self.index.into_base();
        let magnitude = scalar::powf(scalar::abs(shear_rate), index);
        let sign = if shear_rate >= scalar::zero::<T>() {
            scalar::one::<T>()
        } else {
            -scalar::one::<T>()
        };
        Pressure::from_base(self.consistency.into_base() * magnitude * sign)
    }

    fn viscosity(&self, shear_rate: ReciprocalTime<T>) -> DynamicViscosity<T> {
        let shear_rate = shear_rate.into_base();
        if scalar::abs(shear_rate) < <T as FloatElement>::from_f64(1e-12) {
            // Avoid division by zero at centerline
            return DynamicViscosity::from_base(<T as FloatElement>::from_f64(1e12));
        }
        DynamicViscosity::from_base(
            self.consistency.into_base()
                * scalar::powf(
                    scalar::abs(shear_rate),
                    self.index.into_base() - scalar::one::<T>(),
                ),
        )
    }

    fn model_name(&self) -> &'static str {
        "Power-Law"
    }
}

// ============================================================================
// Non-Newtonian Poiseuille Flow - Power-Law
// ============================================================================

/// Non-Newtonian Poiseuille flow with power-law rheology
///
/// Analytical solution for flow between parallel plates with half-width H.
#[derive(Debug, Clone)]
pub struct PowerLawPoiseuille<T: RealField + Copy> {
    /// Power-law model
    pub model: PowerLawModel<T>,
    /// Channel half-width H \[m]
    pub half_width: Length<T>,
    /// Pressure gradient magnitude |dp/dx| [Pa/m] (positive value)
    pub pressure_gradient: PressureGradient<T>,
    /// Channel length for domain \[m]
    pub length: Length<T>,
}

impl<T: RealField + FloatElement + Copy> PowerLawPoiseuille<T> {
    /// Create power-law Poiseuille flow
    ///
    /// # Arguments
    /// * `consistency` - K [Pa·sⁿ]
    /// * `index` - n (power-law index)
    /// * `half_width` - H \[m]
    /// * `pressure_gradient` - |dp/dx| [Pa/m]
    /// * `length` - L \[m] (domain length)
    pub fn new(
        consistency: PowerLawConsistency<T>,
        index: Dimensionless<T>,
        half_width: Length<T>,
        pressure_gradient: PressureGradient<T>,
        length: Length<T>,
    ) -> Self {
        Self {
            model: PowerLawModel::new(consistency, index),
            half_width,
            pressure_gradient,
            length,
        }
    }

    /// Centerline velocity u_c
    ///
    /// For power-law fluid:
    /// ```text
    /// u_c = [(n/(n+1)) · (H^(n+1)/K) · (dp/dx)^(1/n)]
    /// ```
    pub fn centerline_velocity(&self) -> Velocity<T> {
        let n = self.model.index.into_base();
        let one = scalar::one::<T>();
        let h = self.half_width.into_base();
        let dp_dx = self.pressure_gradient.into_base();
        let k = self.model.consistency.into_base();

        // u_c = [n/(n+1)] · (H/K)^(1/n) · (dp/dx·H)^(1/n)
        let factor = n / (n + one);
        let term = (dp_dx * h) / k;

        Velocity::from_base(factor * h * scalar::powf(term, one / n))
    }

    /// Velocity profile u(y) for power-law fluid
    ///
    /// ```text
    /// u(y) = u_c · [1 - |y/H|^((n+1)/n)]
    /// ```
    pub fn velocity_at(&self, y: Length<T>) -> Velocity<T> {
        let n = self.model.index.into_base();
        let h = self.half_width.into_base();
        let u_c = self.centerline_velocity().into_base();

        let y_normalized = scalar::abs(y.into_base() / h);

        if y_normalized >= scalar::one::<T>() {
            // Outside channel
            return Velocity::from_base(scalar::zero::<T>());
        }

        // Exponent: (n+1)/n
        let exponent = (n + scalar::one::<T>()) / n;
        Velocity::from_base(u_c * (scalar::one::<T>() - scalar::powf(y_normalized, exponent)))
    }

    /// Flow rate per unit depth Q' [m²/s]
    ///
    /// ```text
    /// Q' = 2·∫[0 to H] u(y) dy
    /// ```
    pub fn flow_rate_per_depth(&self) -> AreaPerTime<T> {
        let n = self.model.index.into_base();
        let h = self.half_width.into_base();
        let u_c = self.centerline_velocity().into_base();

        // Q' = 2·u_c·H·[n/(2n+1)]
        let two = <T as FloatElement>::from_f64(2.0);
        let factor = n / (two * n + scalar::one::<T>());

        AreaPerTime::from_base(two * u_c * h * factor)
    }

    /// Wall shear stress τ_w \[Pa]
    ///
    /// ```text
    /// τ_w = H · dp/dx
    /// ```
    pub fn wall_shear_stress(&self) -> Pressure<T> {
        Pressure::from_base(self.half_width.into_base() * self.pressure_gradient.into_base())
    }

    /// Wall shear rate γ̇_w [1/s]
    ///
    /// From τ_w = K·γ̇_w^n:
    /// ```text
    /// γ̇_w = (τ_w / K)^(1/n)
    /// ```
    pub fn wall_shear_rate(&self) -> ReciprocalTime<T> {
        let tau_w = self.wall_shear_stress().into_base();
        let k = self.model.consistency.into_base();
        let n = self.model.index.into_base();

        ReciprocalTime::from_base(scalar::powf(tau_w / k, scalar::one::<T>() / n))
    }

    /// Reynolds number for power-law fluid
    ///
    /// Defined using generalized Reynolds number:
    /// ```text
    /// Re_PL = ρ·u_c^(2-n)·H^n / K
    /// ```
    pub fn reynolds_number(&self, density: MassDensity<T>) -> Dimensionless<T> {
        let u_c = self.centerline_velocity().into_base();
        let n = self.model.index.into_base();
        let h = self.half_width.into_base();
        let k = self.model.consistency.into_base();

        Dimensionless::from_base(
            density.into_base()
                * scalar::powf(u_c, <T as FloatElement>::from_f64(2.0) - n)
                * scalar::powf(h, n)
                / k,
        )
    }
}

impl<T: RealField + FloatElement + Copy> AnalyticalSolution<T> for PowerLawPoiseuille<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        let u = self.velocity_at(Length::from_base(y)).into_base();
        Vector3::new(u, scalar::zero::<T>(), scalar::zero::<T>())
    }

    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        // Linear pressure drop: p(x) = p0 - (dp/dx)·x
        -self.pressure_gradient.into_base() * x
    }

    fn name(&self) -> &'static str {
        "Power-Law Poiseuille (2D Channel)"
    }

    fn domain_bounds(&self) -> [T; 6] {
        [
            scalar::zero::<T>(),
            self.length.into_base(), // x: [0, L]
            -self.half_width.into_base(),
            self.half_width.into_base(), // y: [-H, H]
            scalar::zero::<T>(),
            scalar::zero::<T>(), // z: 0 (2D)
        ]
    }

    fn length_scale(&self) -> T {
        self.half_width.into_base()
    }

    fn velocity_scale(&self) -> T {
        self.centerline_velocity().into_base()
    }
}

// ============================================================================
// Casson Blood Model Wrapper
// ============================================================================

/// Wrapper for Casson blood rheology in Poiseuille flow
impl<T: CfdScalar> RheologicalModel<T> for CassonBlood<T> {
    fn shear_stress(&self, shear_rate: ReciprocalTime<T>) -> Pressure<T> {
        let shear_rate = shear_rate.into_base();
        let gamma_dot = scalar::abs(shear_rate);
        let sqrt_tau_y = scalar::sqrt(self.yield_stress.into_base());
        let sqrt_mu_inf = scalar::sqrt(self.infinite_shear_viscosity.into_base());
        let sqrt_gamma = scalar::sqrt(gamma_dot);

        let sqrt_tau = sqrt_tau_y + sqrt_mu_inf * sqrt_gamma;
        let sign = if shear_rate >= scalar::zero::<T>() {
            scalar::one::<T>()
        } else {
            -scalar::one::<T>()
        };
        Pressure::from_base(sqrt_tau * sqrt_tau * sign)
    }

    fn viscosity(&self, shear_rate: ReciprocalTime<T>) -> DynamicViscosity<T> {
        DynamicViscosity::from_base(self.apparent_viscosity(scalar::abs(shear_rate.into_base())))
    }

    fn model_name(&self) -> &'static str {
        "Casson"
    }
}

// ============================================================================
// Non-Newtonian Poiseuille - Casson Blood (Numerical Solution)
// ============================================================================

/// Casson blood Poiseuille flow (numerical solution required)
///
/// The velocity profile for Casson fluid must be solved numerically.
/// We provide numerical integration and validation methods.
#[derive(Debug, Clone)]
pub struct CassonPoiseuille<T: CfdScalar> {
    /// Casson blood model
    pub model: CassonBlood<T>,
    /// Channel half-width H \[m]
    pub half_width: Length<T>,
    /// Pressure gradient magnitude |dp/dx| [Pa/m]
    pub pressure_gradient: PressureGradient<T>,
    /// Channel length \[m]
    pub length: Length<T>,
    /// Plug flow radius (where τ < τ_y) \[m]
    pub plug_radius: Length<T>,
}

impl<T: CfdScalar> CassonPoiseuille<T> {
    /// Create Casson Poiseuille flow
    ///
    /// # Arguments
    /// * `model` - Casson blood rheology
    /// * `half_width` - H \[m]
    /// * `pressure_gradient` - |dp/dx| [Pa/m]
    /// * `length` - L \[m]
    pub fn new(
        model: CassonBlood<T>,
        half_width: Length<T>,
        pressure_gradient: PressureGradient<T>,
        length: Length<T>,
    ) -> Self {
        // Calculate plug radius: y_p where τ(y_p) = τ_y
        // τ(y) = (dp/dx)·y, so y_p = τ_y / (dp/dx)
        let tau_y = model.yield_stress.into_base();
        let y_p = tau_y / pressure_gradient.into_base();

        // Clamp to channel width
        let plug_radius = if y_p > half_width.into_base() {
            half_width
        } else {
            Length::from_base(y_p)
        };

        Self {
            model,
            half_width,
            pressure_gradient,
            length,
            plug_radius,
        }
    }

    /// Velocity in plug region (|y| ≤ y_p)
    ///
    /// In the plug region where τ < τ_y, the fluid moves as a rigid body
    /// with velocity equal to the boundary velocity at y = y_p.
    pub fn plug_velocity(&self) -> Velocity<T> {
        // Integrate from plug boundary to wall
        // This requires numerical integration of du/dy = γ̇(y)
        self.velocity_at_numerical(self.plug_radius)
    }

    /// Velocity profile u(y) - numerical integration
    ///
    /// Integrate du/dy from y to H:
    /// ```text
    /// u(y) = ∫[y to H] γ̇(η) dη
    /// ```
    /// where γ̇(η) is obtained from Casson equation inverted.
    pub fn velocity_at_numerical(&self, y: Length<T>) -> Velocity<T> {
        let y_abs = scalar::abs(y.into_base());
        let half_width = self.half_width.into_base();

        if y_abs >= half_width {
            return Velocity::from_base(scalar::zero::<T>());
        }

        if y_abs < self.plug_radius.into_base() {
            // Plug region: constant velocity
            return self.plug_velocity();
        }

        // Numerical integration using Simpson's rule
        let n_points = 100;
        let dy = (half_width - y_abs) / scalar::from_usize::<T>(n_points);
        let mut integral = scalar::zero::<T>();

        let two = <T as FloatElement>::from_f64(2.0);
        let four = <T as FloatElement>::from_f64(4.0);
        let three = <T as FloatElement>::from_f64(3.0);
        let sqrt_tau_y = scalar::sqrt(self.model.yield_stress.into_base());
        let sqrt_mu_inf = scalar::sqrt(self.model.infinite_shear_viscosity.into_base());

        for i in 0..=n_points {
            let eta = y_abs + scalar::from_usize::<T>(i) * dy;
            let tau = self.pressure_gradient.into_base() * eta;

            // Casson equation: √τ = √τ_y + √(μ_∞·γ̇)
            // Solve for γ̇: γ̇ = [(√τ - √τ_y)² / μ_∞]
            let sqrt_tau = scalar::sqrt(tau);

            let gamma_dot = if tau > self.model.yield_stress.into_base() {
                let diff = sqrt_tau - sqrt_tau_y;
                (diff / sqrt_mu_inf) * (diff / sqrt_mu_inf)
            } else {
                scalar::zero::<T>()
            };

            // Simpson's rule weights
            let weight = if i == 0 || i == n_points {
                scalar::one::<T>()
            } else if i % 2 == 0 {
                two
            } else {
                four
            };

            integral += weight * gamma_dot;
        }

        Velocity::from_base(integral * dy / three)
    }

    /// Wall shear stress \[Pa]
    pub fn wall_shear_stress(&self) -> Pressure<T> {
        Pressure::from_base(self.half_width.into_base() * self.pressure_gradient.into_base())
    }

    /// Flow rate per unit depth (numerical integration)
    pub fn flow_rate_per_depth(&self) -> AreaPerTime<T> {
        // Q' = 2·∫[0 to H] u(y) dy
        let n_points = 100;
        let dy = self.half_width.into_base() / scalar::from_usize::<T>(n_points);
        let mut integral = scalar::zero::<T>();

        let two = <T as FloatElement>::from_f64(2.0);
        let four = <T as FloatElement>::from_f64(4.0);
        let three = <T as FloatElement>::from_f64(3.0);

        for i in 0..=n_points {
            let y = scalar::from_usize::<T>(i) * dy;
            let u = self.velocity_at_numerical(Length::from_base(y)).into_base();

            let weight = if i == 0 || i == n_points {
                scalar::one::<T>()
            } else if i % 2 == 0 {
                two
            } else {
                four
            };

            integral += weight * u;
        }

        AreaPerTime::from_base(two * integral * dy / three)
    }
}

impl<T: CfdScalar> AnalyticalSolution<T> for CassonPoiseuille<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        let u = self.velocity_at_numerical(Length::from_base(y)).into_base();
        Vector3::new(u, scalar::zero::<T>(), scalar::zero::<T>())
    }

    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        -self.pressure_gradient.into_base() * x
    }

    fn name(&self) -> &'static str {
        "Casson Blood Poiseuille (2D Channel)"
    }

    fn domain_bounds(&self) -> [T; 6] {
        [
            scalar::zero::<T>(),
            self.length.into_base(),
            -self.half_width.into_base(),
            self.half_width.into_base(),
            scalar::zero::<T>(),
            scalar::zero::<T>(),
        ]
    }

    fn length_scale(&self) -> T {
        self.half_width.into_base()
    }

    fn velocity_scale(&self) -> T {
        self.plug_velocity().into_base()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn test_power_law_newtonian_limit() {
        // For n=1, power-law should match Newtonian Poiseuille
        let mu = 0.001_f64; // Pa·s
        let h = 0.001_f64; // 1 mm
        let dp_dx = 100.0_f64; // 100 Pa/m
        let l = 0.1_f64; // 10 cm

        let flow = PowerLawPoiseuille::new(
            PowerLawConsistency::new(mu),
            Dimensionless::from_base(1.0),
            Length::from_base(h),
            PressureGradient::from_base(dp_dx),
            Length::from_base(l),
        );

        // Newtonian: u_c = H²·dp/dx / (2μ)
        let u_c_expected = h * h * dp_dx / (2.0 * mu);
        let u_c_actual = flow.centerline_velocity().into_base();

        assert_relative_eq!(u_c_actual, u_c_expected, epsilon = 1e-10);

        // Check parabolic profile: u(y) = u_c·(1 - y²/H²)
        let y = 0.5 * h;
        let u_expected = u_c_expected * (1.0 - (y / h).powi(2));
        let u_actual = flow.velocity_at(Length::from_base(y)).into_base();

        assert_relative_eq!(u_actual, u_expected, epsilon = 1e-10);
    }

    #[test]
    fn test_power_law_shear_thinning() {
        // Shear-thinning (n=0.6) should have higher centerline velocity than Newtonian
        let k = 0.01_f64; // Pa·s^n
        let n = 0.6_f64;
        let h = 0.001_f64;
        let dp_dx = 100.0_f64;
        let l = 0.1_f64;

        let flow = PowerLawPoiseuille::<f64>::new(
            PowerLawConsistency::new(k),
            Dimensionless::from_base(n),
            Length::from_base(h),
            PressureGradient::from_base(dp_dx),
            Length::from_base(l),
        );
        let u_c = flow.centerline_velocity().into_base();

        // Should be finite and positive
        assert!(u_c > 0.0 && u_c.is_finite());

        // Velocity should decrease from center to wall
        let u_center = flow.velocity_at(Length::from_base(0.0)).into_base();
        let u_mid = flow.velocity_at(Length::from_base(0.5 * h)).into_base();
        let u_wall = flow.velocity_at(Length::from_base(h)).into_base();

        assert!(u_center > u_mid);
        assert!(u_mid > u_wall);
        assert_relative_eq!(u_wall, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_power_law_wall_shear() {
        let k = 0.01_f64;
        let n = 0.6_f64;
        let h = 0.001_f64;
        let dp_dx = 100.0_f64;
        let l = 0.1_f64;

        let flow = PowerLawPoiseuille::new(
            PowerLawConsistency::new(k),
            Dimensionless::from_base(n),
            Length::from_base(h),
            PressureGradient::from_base(dp_dx),
            Length::from_base(l),
        );
        let tau_w = flow.wall_shear_stress().into_base();

        // τ_w = H·dp/dx
        assert_relative_eq!(tau_w, h * dp_dx, epsilon = 1e-10);
    }

    #[test]
    fn test_power_law_typed_derived_metrics() {
        let viscosity = 0.001_f64;
        let h = 0.001_f64;
        let dp_dx = 100.0_f64;
        let flow = PowerLawPoiseuille::new(
            PowerLawConsistency::new(viscosity),
            Dimensionless::from_base(1.0),
            Length::from_base(h),
            PressureGradient::from_base(dp_dx),
            Length::from_base(0.1),
        );

        let centerline = flow.centerline_velocity().into_base();
        let expected_flow_rate = 2.0 * h * centerline / 3.0;
        assert_relative_eq!(
            flow.flow_rate_per_depth().into_base(),
            expected_flow_rate,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            flow.wall_shear_rate().into_base(),
            h * dp_dx / viscosity,
            epsilon = 1e-10
        );

        let reynolds = flow.reynolds_number(MassDensity::from_base(1000.0));
        assert!(reynolds.into_base().is_finite() && reynolds.into_base() > 0.0);

        let model = PowerLawModel::newtonian(DynamicViscosity::from_base(viscosity));
        assert_relative_eq!(
            model
                .shear_stress(ReciprocalTime::from_base(10.0))
                .into_base(),
            10.0 * viscosity,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            model.viscosity(ReciprocalTime::from_base(10.0)).into_base(),
            viscosity,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_casson_plug_flow() {
        // High yield stress should create large plug region
        let tau_y = 0.01_f64; // 10 mPa (yield stress)
        let mu_inf = 0.0035_f64; // 3.5 cP

        let model = CassonBlood::new(
            MassDensity::from_base(1060.0),
            Pressure::from_base(tau_y),
            DynamicViscosity::from_base(mu_inf),
            Dimensionless::from_base(0.45),
        );

        let h = 0.001_f64; // 1 mm
        let dp_dx = 100.0_f64; // 100 Pa/m
        let l = 0.1_f64;

        let flow = CassonPoiseuille::new(
            model,
            Length::from_base(h),
            PressureGradient::from_base(dp_dx),
            Length::from_base(l),
        );

        // Plug radius: y_p = τ_y / (dp/dx)
        let y_p_expected = tau_y / dp_dx;
        assert_relative_eq!(flow.plug_radius.into_base(), y_p_expected, epsilon = 1e-10);

        // In plug region, velocity should be constant
        let u_plug = flow.plug_velocity().into_base();
        let u_at_half_plug = flow
            .velocity_at_numerical(Length::from_base(0.5 * flow.plug_radius.into_base()))
            .into_base();

        // Should be approximately equal (within numerical error)
        assert_relative_eq!(u_at_half_plug, u_plug, epsilon = 0.01);
    }

    #[test]
    fn test_casson_wall_velocity_zero() {
        let model = CassonBlood::new(
            MassDensity::from_base(1060.0),
            Pressure::from_base(0.0056),
            DynamicViscosity::from_base(0.00345),
            Dimensionless::from_base(0.45),
        );
        let h = 0.001_f64;
        let dp_dx = 100.0_f64;
        let l = 0.1_f64;

        let flow = CassonPoiseuille::<f64>::new(
            model,
            Length::from_base(h),
            PressureGradient::from_base(dp_dx),
            Length::from_base(l),
        );

        // Wall velocity should be zero
        let u_wall = flow.velocity_at_numerical(Length::from_base(h)).into_base();
        assert!(u_wall.abs() < 1e-6, "Wall velocity {u_wall} should be ~0");
    }

    #[test]
    fn test_casson_flow_rate() {
        let model = CassonBlood::new(
            MassDensity::from_base(1060.0),
            Pressure::from_base(0.0056),
            DynamicViscosity::from_base(0.00345),
            Dimensionless::from_base(0.45),
        );
        let h = 0.001_f64;
        let dp_dx = 100.0_f64;
        let l = 0.1_f64;

        let flow = CassonPoiseuille::<f64>::new(
            model,
            Length::from_base(h),
            PressureGradient::from_base(dp_dx),
            Length::from_base(l),
        );
        let q = flow.flow_rate_per_depth().into_base();

        // Flow rate should be positive and finite
        assert!(q > 0.0 && q.is_finite(), "Flow rate {q} invalid");
    }
}
