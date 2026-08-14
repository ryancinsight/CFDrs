use super::SerpentineGeometry;
use crate::scalar::{self};
use cfd_core::CfdScalar;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

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
/// where `w` is channel width and `D` is molecular diffusivity. A mixing
/// fraction `M = 1 - σ(t)/σ(0)` reaches 90% when the variance ratio is `0.01`.
///
/// # Theorem — Transverse Diffusion Mixing Fraction
///
/// For a two-stream step concentration in a straight channel with no-flux
/// sidewalls, the Neumann eigenfunctions `cos(nπy/w)` diagonalize the
/// transverse diffusion operator. Only odd modes are present, and each mode
/// decays as `exp(-n²π²Dt/w²)`. The normalized variance therefore equals the
/// odd-mode series above, and `M(x) = 1 - sqrt(σ²(x/u)/σ²(0))`.
///
/// **Proof.** Expanding the centered inlet step `c(y,0)-c_mean` in the Neumann
/// basis gives coefficients `a_n = -2 sin(nπ/2)/(nπ)`. Orthogonality gives
/// `σ²(t)=Σ a_n² exp(-2n²π²Dt/w²)/2`, while `σ²(0)=1/4`. Substitution leaves
/// the odd-mode series and the stated mixing fraction. ∎
pub struct AdvectionDiffusionMixing<T: CfdScalar + Copy> {
    /// Channel width \[m]
    pub width: T,
    /// Mean flow velocity \[m/s]
    pub velocity: T,
    /// Diffusion coefficient [m²/s]
    pub diffusion_coeff: T,
}

impl<T: CfdScalar + Copy + FloatElement> AdvectionDiffusionMixing<T> {
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
        (self.velocity * self.width) / (self.diffusion_coeff + <T as FloatElement>::from_f64(1e-15))
    }

    #[inline]
    fn variance_ratio_from_fourier_number(fourier: T) -> T {
        if fourier <= scalar::zero() {
            return scalar::one();
        }

        let pi = <T as FloatElement>::from_f64(PI);
        let pi_sq = pi * pi;
        let eight_over_pi_sq = <T as FloatElement>::from_f64(8.0) / pi_sq;
        let two = <T as FloatElement>::from_f64(2.0);
        let mut sum: T = scalar::zero();

        for mode in 0..256 {
            let n = scalar::from_usize::<T>(2 * mode + 1);
            let n_sq = n * n;
            sum += FloatElement::exp(-two * n_sq * pi_sq * fourier) / n_sq;
        }

        eight_over_pi_sq * sum
    }

    /// Calculate advective length required for 90% variance-based mixing.
    ///
    /// The method solves `σ²(t)/σ²(0)=0.01` for `Fo = Dt/w²` by bisection of
    /// the closed-form eigenfunction series and returns `L90 = u w² Fo90 / D`.
    pub fn mixing_length_90_percent(&self) -> T {
        if self.width <= scalar::zero()
            || self.velocity <= scalar::zero()
            || self.diffusion_coeff <= scalar::zero()
        {
            return scalar::zero();
        }

        let target_variance_ratio = <T as FloatElement>::from_f64(0.01);
        let mut lo: T = scalar::zero();
        let mut hi = <T as FloatElement>::from_f64(0.5);

        for _ in 0..80 {
            let mid = (lo + hi) / <T as FloatElement>::from_f64(2.0);
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
    /// `t_mix` = `L_mix` / u
    pub fn mixing_time_90_percent(&self) -> T {
        let l_mix = self.mixing_length_90_percent();
        if self.velocity <= scalar::zero() {
            scalar::zero()
        } else {
            l_mix / self.velocity
        }
    }

    /// Estimate concentration profile at position x
    ///
    /// For two inlets with concentrations `c_A` and `c_B`, the concentration
    /// at position x depends on diffusion progress:
    ///
    /// ```text
    /// c(x, y) = c_A × erf(y / sqrt(4Dx/u)) + c_B × (1 - erf(...))
    /// ```
    ///
    /// At the centerline (y=0), maximum mixing occurs there.
    /// Returns fraction mixed at position x relative to channel width.
    pub fn mixing_fraction(&self, x: T) -> T {
        if x <= scalar::zero()
            || self.width <= scalar::zero()
            || self.velocity <= scalar::zero()
            || self.diffusion_coeff <= scalar::zero()
        {
            return scalar::zero();
        }

        let fourier = self.diffusion_coeff * x / (self.velocity * self.width * self.width);
        let variance_ratio = Self::variance_ratio_from_fourier_number(fourier);
        scalar::one::<T>() - <T as NumericElement>::sqrt(variance_ratio)
    }
}

/// Solution for serpentine channel mixing
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SerpentineMixingSolution<T: CfdScalar + Copy> {
    /// Inlet concentration (first inlet) [mol/m³]
    pub c_inlet_a: T,
    /// Inlet concentration (second inlet) [mol/m³]
    pub c_inlet_b: T,
    /// Peclet number (dimensionless)
    pub peclet: T,
    /// Mixing length for 90% homogeneity \[m]
    pub l_mix_90: T,
    /// Mixing time to 90% \[s]
    pub t_mix_90: T,
    /// Mixing fraction at outlet [0-1]
    pub mixing_fraction_outlet: T,
    /// Estimated viscous pressure drop \[Pa]
    pub pressure_drop: T,
}

impl<T: CfdScalar + Copy + FloatElement> SerpentineMixingSolution<T> {
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

        let d_h = <T as FloatElement>::from_f64(2.0) * (geometry.width * geometry.height)
            / (geometry.width + geometry.height);
        let re = (density * velocity * d_h) / viscosity;
        let f = <T as FloatElement>::from_f64(64.0) / re.max(<T as FloatElement>::from_f64(1.0));
        let dynamic_pressure = <T as FloatElement>::from_f64(0.5) * density * velocity * velocity;
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
        self.mixing_fraction_outlet > <T as FloatElement>::from_f64(0.9)
    }

    /// Estimate outlet concentration (assuming complete mixing)
    ///
    /// For equal volume flows: `c_outlet` = (`c_A` + `c_B`) / 2
    pub fn estimated_outlet_concentration(&self) -> T {
        (self.c_inlet_a + self.c_inlet_b) / <T as FloatElement>::from_f64(2.0)
    }
}

#[cfg(test)]
mod tests {
    use super::AdvectionDiffusionMixing;
    use eunomia::NumericElement;

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
            <f64 as NumericElement>::abs(fraction - 0.9) < 5.0e-12,
            "Expected eigenfunction L90 to produce 90% mixing, got {fraction:.15}"
        );
    }

    #[test]
    fn mixing_length_scales_with_velocity_and_inverse_diffusivity() {
        let baseline = AdvectionDiffusionMixing::new(200e-6_f64, 0.01, 1.0e-10);
        let faster = AdvectionDiffusionMixing::new(200e-6_f64, 0.02, 1.0e-10);
        let more_diffusive = AdvectionDiffusionMixing::new(200e-6_f64, 0.01, 2.0e-10);

        assert!(
            <f64 as NumericElement>::abs(
                faster.mixing_length_90_percent() / baseline.mixing_length_90_percent() - 2.0,
            ) < 1.0e-12
        );
        assert!(
            <f64 as NumericElement>::abs(
                more_diffusive.mixing_length_90_percent() / baseline.mixing_length_90_percent()
                    - 0.5,
            ) < 1.0e-12
        );
    }
}
