use super::constants;
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};

/// Fåhræus-Lindqvist effect for microvascular blood flow
///
/// In microvessels (D < 300 μm), blood exhibits apparent viscosity reduction due to:
/// 1. Axial migration of RBCs (cell-free layer near wall)
/// 2. Fåhræus effect (reduction in tube hematocrit vs feed hematocrit)
///
/// # Physical Behavior
/// The apparent viscosity DECREASES as vessel diameter decreases from 300 μm down to
/// approximately 7 μm (minimum viscosity), then increases again for very small vessels
/// where the RBC size becomes comparable to vessel diameter.
///
/// # Empirical Correlations
/// - **Pries et al. (1992)**: Standard fit for in-vitro apparent relative viscosity
/// - **Secomb (2017)**: Improved fit with smoothed transitions for microvascular networks
///
/// # References
/// - Pries, A.R., Neuhaus, D., Gaehtgens, P. (1992) "Blood viscosity in tube flow:
///   dependence on diameter and hematocrit"
/// - Secomb, T.W. (2017) "Blood Flow in the Microcirculation"
#[derive(Debug, Clone, Copy)]
pub struct FahraeuasLindqvist<T: RealField + Copy> {
    /// Vessel diameter \[m]
    pub diameter: T,
    /// Hematocrit (volume fraction of RBCs) [-]
    pub hematocrit: T,
    /// Plasma viscosity [Pa·s]
    pub plasma_viscosity: T,
}

impl<T: RealField + FloatElement + Copy> FahraeuasLindqvist<T> {
    /// Create new Fåhræus-Lindqvist calculator
    pub fn new(diameter: T, hematocrit: T) -> Self {
        Self {
            diameter,
            hematocrit,
            plasma_viscosity: scalar::<T>(constants::PLASMA_VISCOSITY_37C),
        }
    }

    /// Check if Fåhræus-Lindqvist effect is significant
    ///
    /// Effect is significant for D < 300 μm
    pub fn is_significant(&self) -> bool {
        self.diameter < scalar::<T>(constants::FAHRAEUS_LINDQVIST_CRITICAL_DIAMETER)
    }

    /// Calculate relative apparent viscosity using exact Pries et al. (1992) formulation
    ///
    /// Returns μ_rel = μ_app / μ_plasma
    ///
    /// # Exact Physical Model (Pries et al. 1992)
    /// For a vessel of diameter D (in μm) and hematocrit H_t:
    /// - Relative viscosity of plasma: μ_p_rel ≈ 1.0
    /// - Exact formulation for apparent viscosity:
    ///   μ_rel = 1 + (μ_45 - 1) * ((1-H_t)^C - 1) / ((1-0.45)^C - 1)
    pub fn pries_relative_viscosity(&self) -> T {
        let (mu_45, d_um) = self.compute_mu_45();
        let one = <T as NumericElement>::ONE;
        let zero = <T as NumericElement>::ZERO;

        // Pries (1992) exponent C
        // C = 0.8 + exp(-0.075D) * (-1 + 1/(1 + 10^-11 * D^12))
        let exponent_c = scalar::<T>(0.8)
            + <T as FloatElement>::exp((zero - scalar::<T>(0.075)) * d_um)
                * ((zero - one)
                    + one
                        / (one
                            + scalar::<T>(1e-11)
                                * <T as FloatElement>::powf(d_um, scalar::<T>(12.0))));

        self.compute_final_relative_viscosity(mu_45, exponent_c)
    }

    /// Calculate relative apparent viscosity using Secomb (2017) corrected formulation
    ///
    /// Returns μ_rel = μ_app / μ_plasma
    ///
    /// # Exact Physical Model (Secomb 2017)
    /// Incorporates a smoothed transition for intermediate diameters:
    /// C = (0.8 + exp(-0.075D)) * (-1 + 1/(1 + 10^-11 * D^12)) + (1/(1 + 10^-11 * D^12))
    pub fn secomb_relative_viscosity(&self) -> T {
        let (mu_45, d_um) = self.compute_mu_45();
        let one = <T as NumericElement>::ONE;
        let zero = <T as NumericElement>::ZERO;

        let d12_term =
            one / (one + scalar::<T>(1e-11) * <T as FloatElement>::powf(d_um, scalar::<T>(12.0)));

        // Secomb (2017) exponent C
        let exponent_c = (scalar::<T>(0.8)
            + <T as FloatElement>::exp((zero - scalar::<T>(0.075)) * d_um))
            * ((zero - one) + d12_term)
            + d12_term;

        self.compute_final_relative_viscosity(mu_45, exponent_c)
    }

    /// Default relative-viscosity model using the Secomb 2017 formulation.
    pub fn relative_viscosity(&self) -> T {
        self.secomb_relative_viscosity()
    }

    /// Helper to compute relative viscosity at 45% hematocrit (μ_45), shared between
    /// Pries (1992) and Secomb (2017) parameterisations.
    /// Caps effective D at 300 μm.
    #[inline]
    fn compute_mu_45(&self) -> (T, T) {
        let zero = <T as NumericElement>::ZERO;
        let mut d_um = self.diameter * scalar::<T>(1e6); // Convert to μm
                                                         // Fahraeus-Lindqvist scales back to bulk viscosity above 300μm
        let d_max = scalar::<T>(300.0);
        if d_um > d_max {
            d_um = d_max;
        }

        // Shared μ_45 fit from Pries et al. (1992)
        // μ_45 = 6 * exp(-0.085D) + 3.2 - 2.44 * exp(-0.06D^0.645)
        let mu_45 = scalar::<T>(6.0) * <T as FloatElement>::exp((zero - scalar::<T>(0.085)) * d_um)
            + scalar::<T>(3.2)
            - scalar::<T>(2.44)
                * <T as FloatElement>::exp(
                    (zero - scalar::<T>(0.06))
                        * <T as FloatElement>::powf(d_um, scalar::<T>(0.645)),
                );

        (mu_45, d_um)
    }

    /// Calculate final μ_rel avoiding division by zero
    #[inline]
    fn compute_final_relative_viscosity(&self, mu_45: T, exponent_c: T) -> T {
        let one = <T as NumericElement>::ONE;
        let zero = <T as NumericElement>::ZERO;
        let ht_clamp = <T as NumericElement>::max_scalar(self.hematocrit, zero);

        let num_base = <T as NumericElement>::max_scalar(one - ht_clamp, scalar::<T>(0.001));
        let den_base =
            <T as NumericElement>::max_scalar(one - scalar::<T>(0.45), scalar::<T>(0.001));

        let numerator = <T as FloatElement>::powf(num_base, exponent_c) - one;
        let denominator = <T as FloatElement>::powf(den_base, exponent_c) - one;

        let mu_rel = if <T as NumericElement>::abs(denominator) < scalar::<T>(1e-10) {
            one
        } else {
            one + (mu_45 - one) * (numerator / denominator)
        };

        <T as NumericElement>::max_scalar(mu_rel, one)
    }

    /// Calculate apparent viscosity in microvessel [Pa·s]
    pub fn apparent_viscosity(&self) -> T {
        self.plasma_viscosity * self.relative_viscosity()
    }

    /// Calculate tube hematocrit from feed hematocrit (Fåhræus effect)
    ///
    /// H_tube / H_feed = empirical correlation
    pub fn tube_hematocrit(&self) -> T {
        let d_um = self.diameter * scalar::<T>(1e6);
        let one = <T as NumericElement>::ONE;
        let zero = <T as NumericElement>::ZERO;

        // Empirical correlation for tube hematocrit reduction
        // Valid for D > 10 μm
        if d_um < scalar::<T>(10.0) {
            return self.hematocrit * scalar::<T>(0.5);
            // Approximate minimum
        }

        // Pries et al. empirical fit
        let reduction_factor =
            one - scalar::<T>(1.7) * <T as FloatElement>::exp((zero - d_um) / scalar::<T>(40.0));
        self.hematocrit * reduction_factor
    }
}

#[inline]
fn scalar<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fahraeus_lindqvist_significance() {
        // Large vessel - effect not significant
        let large = FahraeuasLindqvist::<f64>::new(1e-3, 0.45); // 1 mm
        assert!(!large.is_significant());

        // Small vessel - effect significant
        let small = FahraeuasLindqvist::<f64>::new(50e-6, 0.45); // 50 μm
        assert!(small.is_significant());
    }

    #[test]
    fn test_fahraeus_lindqvist_viscosity_reduction() {
        // In smaller vessels, apparent viscosity should be lower
        let d_100 = FahraeuasLindqvist::<f64>::new(100e-6, 0.45);
        let d_50 = FahraeuasLindqvist::<f64>::new(50e-6, 0.45);

        let mu_100 = d_100.apparent_viscosity();
        let mu_50 = d_50.apparent_viscosity();

        // Smaller vessel should have lower apparent viscosity
        assert!(
            mu_50 < mu_100,
            "μ(50 μm) = {mu_50} should be < μ(100 μm) = {mu_100}"
        );
    }

    #[test]
    fn test_fahraeus_lindqvist_tube_hematocrit() {
        let fl = FahraeuasLindqvist::<f64>::new(50e-6, 0.45);
        let h_tube = fl.tube_hematocrit();

        // Tube hematocrit should be less than feed due to Fåhræus effect
        assert!(
            h_tube < 0.45,
            "Tube hematocrit {} should be less than feed {}",
            h_tube,
            0.45
        );
    }
}
