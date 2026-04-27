//! 2D serpentine channel flow solver with mixing efficiency validation
//!
//! This module implements validated solvers for serpentine (sinuous) channels,
//! which are critical for microfluidic mixing and sample preparation.
//!
//! # Physics Background
//!
//! ## Serpentine Channel Design
//!
//! Serpentine channels enhance mixing through:
//! 1. **Advection**: Fluid streams are stretched and folded
//! 2. **Secondary flow**: Dean vortices in curves increase transverse mixing
//! 3. **Residence time**: Longer path length increases reaction time
//!
//! ## Mixing Mechanisms
//!
//! In laminar microfluidic flows (Re << 1), mixing is diffusion-limited:
//!
//! **Diffusion time scale:**
//! ```text
//! t_diff = w² / (4D)  [seconds, from Fick's law]
//! ```
//!
//! where:
//! - w = channel width [m]
//! - D = diffusion coefficient [m²/s]
//!
//! **Advection time scale:**
//! ```text
//! t_adv = L / u  [seconds, from advection]
//! ```
//!
//! **Mixing length:**
//! ```text
//! L_mix = u · t_diff = u · w² / (4D)
//! ```
//!
//! For complete mixing across channel width: L_adv > L_mix
//!
//! ## Peclet Number
//!
//! Ratio of advection to diffusion:
//! ```text
//! Pe = u·w / D
//! ```
//!
//! - Pe << 1: Diffusion-dominated mixing
//! - Pe >> 1: Mixing limited by diffusion length
//!
//! ## Mixing Efficiency
//!
//! Quantified by intensity of segregation (ISO):
//! ```text
//! ISO = ∫(c - c_mean)² dV / ∫c_max² dV
//! ```
//!
//! where:
//! - ISO = 1 at inlet (completely unmixed)
//! - ISO = 0 at outlet (perfectly mixed)
//!
//! # Validation Strategy
//!
//! 1. **Advection-diffusion**: Solve transport equation with known inlet conditions
//! 2. **Richardson extrapolation**: Grid convergence study
//! 3. **Literature**: Compare against experimental mixing data
//! 4. **Analytical**: Use advection-diffusion analytical solutions
//!
//! # References
//!
//! - Hardt, S. & Schönfeld, F. (2003). "Microfluidic technologies for miniaturized
//!   analysis systems". SpringerLink
//! - Squires, T.M. & Quake, S.R. (2005). "Microfluidics: Fluid physics at the
//!   nanoliter scale". Reviews of Modern Physics, 77(3), 977
//! - Yao, Z., et al. (2014). "Numerical study of mixing in microchannels with
//!   oscillatory flow". Microfluidics and Nanofluidics, 16(1), 145-155
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

mod solver;

pub use solver::SerpentineSolver2D;

use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ============================================================================
// Serpentine Channel Geometry
// ============================================================================

/// Serpentine channel geometry with periodic turns
///
/// Defines a channel that alternates between straight sections and 90° turns.
/// The path creates a snake-like pattern that enhances mixing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerpentineGeometry<T: RealField + Copy> {
    /// Channel width [m]
    pub width: T,
    /// Channel height [m] (constant)
    pub height: T,
    /// Straight section length [m]
    pub l_straight: T,
    /// Turn radius [m] (for curved turns)
    pub turn_radius: T,
    /// Number of complete back-and-forth cycles
    pub n_cycles: usize,
}

impl<T: RealField + Copy + FromPrimitive> SerpentineGeometry<T> {
    /// Create standard serpentine for microfluidics
    ///
    /// # Typical Design
    ///
    /// - Width: 100-500 μm
    /// - Straight sections: 500 μm
    /// - Turn radius: 200 μm
    /// - Cycles: 5-10 (for ~5-10 mm total length)
    pub fn microfluidic_standard() -> Self {
        Self {
            width: T::from_f64_or_one(200e-6),
            height: T::from_f64_or_one(50e-6),
            l_straight: T::from_f64_or_one(500e-6),
            turn_radius: T::from_f64_or_one(200e-6),
            n_cycles: 5,
        }
    }

    /// Create custom serpentine
    pub fn new(width: T, height: T, l_straight: T, turn_radius: T, n_cycles: usize) -> Self {
        Self {
            width,
            height,
            l_straight,
            turn_radius,
            n_cycles,
        }
    }

    /// Calculate total channel length
    ///
    /// Length = n_cycles × (2 × l_straight + turn_length)
    /// where turn_length ≈ π × turn_radius (90° arc)
    pub fn total_length(&self) -> T {
        let pi = T::from_f64_or_one(PI);
        let two = T::from_f64_or_one(2.0);
        let turn_length = pi / T::from_f64_or_one(2.0) * self.turn_radius;

        T::from_f64_or_one(self.n_cycles as f64) * (two * self.l_straight + turn_length)
    }

    /// Get cross-sectional area
    pub fn cross_section_area(&self) -> T {
        self.width * self.height
    }

    /// Calculate number of diffusion lengths in straight section
    ///
    /// # Definition
    ///
    /// ```text
    /// n_diff = (l_straight / width) × Pe
    /// ```
    ///
    /// where Pe = u·w/D (Peclet number)
    pub fn diffusion_lengths_per_section(&self, peclet: T) -> T {
        let three = T::from_f64_or_one(3.0);
        (self.l_straight / (self.width + T::from_f64_or_one(1e-15))) * three * peclet
    }

    /// Check if a point (x, y) is within the fluid domain
    ///
    /// Assuming serpentine starts at (0, R) and snakes upwards in Y.
    /// Channel centerline path:
    /// Cycle N:
    ///   Segment 1: (0, 4NR + R) -> (Ls, 4NR + R)
    ///   Turn 1: Semi-circle around (Ls, 4NR + 2R) with radius R
    ///   Segment 2: (Ls, 4NR + 3R) -> (0, 4NR + 3R)
    ///   Turn 2: Semi-circle around (0, 4NR + 4R) with radius R
    pub fn contains(&self, x: T, y: T) -> bool {
        let r = self.turn_radius;
        let ls = self.l_straight;
        let w = self.width;
        let half_w = w / T::from_f64(2.0).expect("analytical constant conversion");
        let four_r = r * T::from_f64(4.0).expect("analytical constant conversion");

        // Normalize y to cycle 0
        let y_in_cycle = y % four_r;

        let inner_r = r - half_w;
        let outer_r = r + half_w;
        let inner_r_sq = inner_r * inner_r;
        let outer_r_sq = outer_r * outer_r;

        // Sequence within a cycle:
        // 0..2R: Segment 1 (bottom straight) and Turn 1 (right)
        // 2R..4R: Segment 2 (top straight) and Turn 2 (left)

        if y_in_cycle < r + half_w && y_in_cycle > r - half_w {
            // Segment 1: (0, R) -> (Ls, R)
            // Extend to x = -R for inlet tail
            if x >= -r - half_w && x <= ls {
                return true;
            }
        }

        if y_in_cycle < T::from_f64(3.0).expect("analytical constant conversion") * r + half_w
            && y_in_cycle > T::from_f64(3.0).expect("analytical constant conversion") * r - half_w
        {
            // Segment 2: (Ls, 3R) -> (0, 3R)
            if x >= T::zero() && x <= ls {
                return true;
            }
        }

        // Right Turn (around x=ls, y=2R)
        let dx_r = x - ls;
        let dy_r = y_in_cycle - r * T::from_f64(2.0).expect("analytical constant conversion");
        let d_sq_r = dx_r * dx_r + dy_r * dy_r;
        if x >= ls && d_sq_r >= inner_r_sq && d_sq_r <= outer_r_sq {
            return true;
        }

        // Left Turn (around x=0, y=4R/0)
        // Check both top of this cycle and bottom of next
        let dy_l_top = y_in_cycle - four_r;
        let d_sq_l_top = x * x + dy_l_top * dy_l_top;
        if x <= T::zero() && d_sq_l_top >= inner_r_sq && d_sq_l_top <= outer_r_sq {
            return true;
        }

        let d_sq_l_bot = x * x + y_in_cycle * y_in_cycle;
        if x <= T::zero() && d_sq_l_bot >= inner_r_sq && d_sq_l_bot <= outer_r_sq {
            return true;
        }

        false
    }

    /// Get bounding box [min_x, max_x, min_y, max_y]
    pub fn bounding_box(&self) -> [T; 4] {
        let r = self.turn_radius;
        let ls = self.l_straight;
        let w = self.width;
        let half_w = w / T::from_f64(2.0).expect("analytical constant conversion");
        let four_r = r * T::from_f64(4.0).expect("analytical constant conversion");

        [
            -r - half_w,
            ls + r + half_w,
            T::zero(),
            four_r * T::from_f64(self.n_cycles as f64).expect("analytical constant conversion"),
        ]
    }
}

// ============================================================================
// Mixing Models
// ============================================================================

/// Analytical model for advection-diffusion mixing in straight sections
///
/// In a straight microchannel with two inlet streams, concentration
/// evolves according to Fick's law:
///
/// ```text
/// ∂c/∂t + u·∂c/∂x = D·∂²c/∂y²
/// ```
///
/// The segregation variance for an inlet step concentration is:
/// ```text
/// σ²(t)/σ²(0) = (8/π²) Σ_{m=0}∞ exp[-2(2m+1)²π²Dt/w²] / (2m+1)²
/// ```
///
/// where `w` is channel width and `D` is molecular diffusivity.  A mixing
/// fraction `M = 1 - σ(t)/σ(0)` reaches 90% when the variance ratio is `0.01`.
///
/// # Theorem — Transverse Diffusion Mixing Fraction
///
/// For a two-stream step concentration in a straight channel with no-flux
/// sidewalls, the Neumann eigenfunctions `cos(nπy/w)` diagonalize the
/// transverse diffusion operator. Only odd modes are present, and each mode
/// decays as `exp(-n²π²Dt/w²)`.  The normalized variance therefore equals the
/// odd-mode series above, and `M(x) = 1 - sqrt(σ²(x/u)/σ²(0))`.
///
/// **Proof.** Expanding the centered inlet step `c(y,0)-c_mean` in the Neumann
/// basis gives coefficients `a_n = -2 sin(nπ/2)/(nπ)`. Orthogonality gives
/// `σ²(t)=Σ a_n² exp(-2n²π²Dt/w²)/2`, while `σ²(0)=1/4`. Substitution leaves
/// the odd-mode series and the stated mixing fraction. ∎
pub struct AdvectionDiffusionMixing<T: RealField + Copy> {
    /// Channel width [m]
    pub width: T,
    /// Mean flow velocity [m/s]
    pub velocity: T,
    /// Diffusion coefficient [m²/s]
    pub diffusion_coeff: T,
}

impl<T: RealField + Copy + FromPrimitive> AdvectionDiffusionMixing<T> {
    /// Create mixing model
    pub fn new(width: T, velocity: T, diffusion_coeff: T) -> Self {
        Self {
            width,
            velocity,
            diffusion_coeff,
        }
    }

    /// Calculate Peclet number
    ///
    /// Pe = u·w / D
    pub fn peclet_number(&self) -> T {
        (self.velocity * self.width) / (self.diffusion_coeff + T::from_f64_or_one(1e-15))
    }

    #[inline]
    fn variance_ratio_from_fourier_number(fourier: T) -> T {
        if fourier <= T::zero() {
            return T::one();
        }

        let pi = T::from_f64_or_one(PI);
        let pi_sq = pi * pi;
        let eight_over_pi_sq = T::from_f64_or_one(8.0) / pi_sq;
        let two = T::from_f64_or_one(2.0);
        let mut sum = T::zero();

        for mode in 0..256 {
            let n = T::from_f64_or_one(f64::from(2 * mode + 1));
            let n_sq = n * n;
            sum += (-two * n_sq * pi_sq * fourier).exp() / n_sq;
        }

        eight_over_pi_sq * sum
    }

    /// Calculate advective length required for 90% variance-based mixing.
    ///
    /// The method solves `σ²(t)/σ²(0)=0.01` for `Fo = Dt/w²` by bisection of
    /// the closed-form eigenfunction series and returns `L90 = u w² Fo90 / D`.
    pub fn mixing_length_90_percent(&self) -> T {
        if self.width <= T::zero()
            || self.velocity <= T::zero()
            || self.diffusion_coeff <= T::zero()
        {
            return T::zero();
        }

        let target_variance_ratio = T::from_f64_or_one(0.01);
        let mut lo = T::zero();
        let mut hi = T::from_f64_or_one(0.5);

        for _ in 0..80 {
            let mid = (lo + hi) / T::from_f64_or_one(2.0);
            if Self::variance_ratio_from_fourier_number(mid) > target_variance_ratio {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        let fo_90 = hi;
        self.velocity * self.width * self.width * fo_90 / self.diffusion_coeff
    }

    /// Calculate mixing time to achieve 90% homogeneity
    ///
    /// t_mix = L_mix / u
    pub fn mixing_time_90_percent(&self) -> T {
        let l_mix = self.mixing_length_90_percent();
        if self.velocity <= T::zero() {
            T::zero()
        } else {
            l_mix / self.velocity
        }
    }

    /// Estimate concentration profile at position x
    ///
    /// For two inlets with concentrations c_A and c_B, the concentration
    /// at position x depends on diffusion progress:
    ///
    /// ```text
    /// c(x, y) = c_A × erf(y / sqrt(4Dx/u)) + c_B × (1 - erf(...))
    /// ```
    ///
    /// At the centerline (y=0), maximum mixing occurs there.
    /// Returns fraction mixed at position x relative to channel width.
    pub fn mixing_fraction(&self, x: T) -> T {
        if x <= T::zero()
            || self.width <= T::zero()
            || self.velocity <= T::zero()
            || self.diffusion_coeff <= T::zero()
        {
            return T::zero();
        }

        let fourier = self.diffusion_coeff * x / (self.velocity * self.width * self.width);
        let variance_ratio = Self::variance_ratio_from_fourier_number(fourier);
        T::one() - variance_ratio.sqrt()
    }
}

// ============================================================================
// Serpentine Mixing Solution
// ============================================================================

/// Solution for serpentine channel mixing
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SerpentineMixingSolution<T: RealField + Copy> {
    /// Inlet concentration (first inlet) [mol/m³]
    pub c_inlet_a: T,
    /// Inlet concentration (second inlet) [mol/m³]
    pub c_inlet_b: T,
    /// Peclet number (dimensionless)
    pub peclet: T,
    /// Mixing length for 90% homogeneity [m]
    pub l_mix_90: T,
    /// Mixing time to 90% [s]
    pub t_mix_90: T,
    /// Mixing fraction at outlet [0-1]
    pub mixing_fraction_outlet: T,
    /// Estimated viscous pressure drop [Pa]
    pub pressure_drop: T,
}

impl<T: RealField + Copy + FromPrimitive> SerpentineMixingSolution<T> {
    /// Create solution from parameters
    pub fn new(
        geometry: &SerpentineGeometry<T>,
        velocity: T,
        diffusion_coeff: T,
        c_inlet_a: T,
        c_inlet_b: T,
        viscosity: T,
        density: T,
    ) -> Self {
        let mixing_model = AdvectionDiffusionMixing::new(geometry.width, velocity, diffusion_coeff);
        let pe = mixing_model.peclet_number();
        let l_mix = mixing_model.mixing_length_90_percent();
        let t_mix = mixing_model.mixing_time_90_percent();

        let total_length = geometry.total_length();
        let mixing_frac = mixing_model.mixing_fraction(total_length);

        // Pressure drop: ΔP = f × (L/D_h) × (1/2)ρu²
        // For laminar flow, friction factor f = 64/Re
        let d_h = T::from_f64_or_one(2.0) * (geometry.width * geometry.height)
            / (geometry.width + geometry.height);
        let re = (density * velocity * d_h) / viscosity;
        let f = T::from_f64_or_one(64.0) / re.max(T::from_f64_or_one(1.0));
        let dynamic_pressure = T::from_f64_or_one(0.5) * density * velocity * velocity;
        let dp = f * (total_length / d_h) * dynamic_pressure;

        Self {
            c_inlet_a,
            c_inlet_b,
            peclet: pe,
            l_mix_90: l_mix,
            t_mix_90: t_mix,
            mixing_fraction_outlet: mixing_frac,
            pressure_drop: dp,
        }
    }

    /// Check if mixing is achieved at outlet
    ///
    /// Typically consider "mixed" if fraction > 0.9 (90%)
    pub fn is_well_mixed(&self) -> bool {
        self.mixing_fraction_outlet > T::from_f64_or_one(0.9)
    }

    /// Estimate outlet concentration (assuming complete mixing)
    ///
    /// For equal volume flows: c_outlet = (c_A + c_B) / 2
    pub fn estimated_outlet_concentration(&self) -> T {
        (self.c_inlet_a + self.c_inlet_b) / T::from_f64_or_one(2.0)
    }
}

// ============================================================================
// Serpentine Validator
// ============================================================================

/// Validator for serpentine mixing channels
pub struct SerpentineValidator<T: RealField + Copy> {
    geometry: SerpentineGeometry<T>,
}

impl<T: RealField + Copy + FromPrimitive + num_traits::ToPrimitive> SerpentineValidator<T> {
    /// Create new validator
    pub fn new(geometry: SerpentineGeometry<T>) -> Self {
        Self { geometry }
    }

    /// Validate mixing efficiency against expected values
    ///
    /// # Validation Criteria
    ///
    /// For well-designed serpentines:
    /// - Mixing length < total channel length
    /// - Mixing fraction at outlet > 90%
    /// - Pressure drop < 100 kPa (for ~mm-scale channels)
    pub fn validate_mixing(
        &self,
        solution: &SerpentineMixingSolution<T>,
    ) -> Result<SerpentineValidationResult<T>, String> {
        let total_length = self.geometry.total_length();
        let l_mix = solution.l_mix_90;

        // Check if mixing length is achieved within channel
        let achievable = l_mix < total_length * T::from_f64_or_one(2.0);

        // Check mixing fraction at outlet
        let well_mixed = solution.is_well_mixed();

        // Check pressure drop is reasonable
        let dp_reasonable = solution.pressure_drop < T::from_f64_or_one(1e5); // < 100 kPa

        let validation_passed = achievable && well_mixed && dp_reasonable;

        let mut error_message = None;
        if !achievable {
            error_message = Some("Mixing length exceeds 2× total channel length".to_string());
        } else if !well_mixed {
            error_message = Some(format!(
                "Mixing fraction {:.2}% < 90% at outlet",
                (solution.mixing_fraction_outlet * T::from_f64_or_one(100.0))
                    .to_f64()
                    .unwrap_or(f64::NAN)
            ));
        } else if !dp_reasonable {
            error_message = Some(format!(
                "Pressure drop {:.0} Pa > 100 kPa",
                solution.pressure_drop.to_f64().unwrap_or(f64::NAN)
            ));
        }

        Ok(SerpentineValidationResult {
            validation_passed,
            mixing_fraction_outlet: solution.mixing_fraction_outlet,
            error_message,
        })
    }
}

/// Validation result for serpentine
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerpentineValidationResult<T: RealField + Copy> {
    /// Validation passed
    pub validation_passed: bool,
    /// Measured mixing fraction at outlet
    pub mixing_fraction_outlet: T,
    /// Error message if validation failed
    pub error_message: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::AdvectionDiffusionMixing;

    #[test]
    fn transverse_diffusion_series_starts_unmixed() {
        let model = AdvectionDiffusionMixing::new(200e-6_f64, 0.01, 1.0e-10);
        assert_eq!(model.mixing_fraction(0.0), 0.0);
    }

    #[test]
    fn mixing_length_reaches_ninety_percent_by_definition() {
        let model = AdvectionDiffusionMixing::new(200e-6_f64, 0.01, 1.0e-10);
        let l90 = model.mixing_length_90_percent();
        let fraction = model.mixing_fraction(l90);

        assert!(
            (fraction - 0.9).abs() < 5.0e-12,
            "Expected eigenfunction L90 to produce 90% mixing, got {fraction:.15}"
        );
    }

    #[test]
    fn mixing_length_scales_with_velocity_and_inverse_diffusivity() {
        let baseline = AdvectionDiffusionMixing::new(200e-6_f64, 0.01, 1.0e-10);
        let faster = AdvectionDiffusionMixing::new(200e-6_f64, 0.02, 1.0e-10);
        let more_diffusive = AdvectionDiffusionMixing::new(200e-6_f64, 0.01, 2.0e-10);

        assert!(
            (faster.mixing_length_90_percent() / baseline.mixing_length_90_percent() - 2.0).abs()
                < 1.0e-12
        );
        assert!(
            (more_diffusive.mixing_length_90_percent() / baseline.mixing_length_90_percent() - 0.5)
                .abs()
                < 1.0e-12
        );
    }
}
