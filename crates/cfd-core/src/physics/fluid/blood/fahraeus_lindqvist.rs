use super::constants;
use nalgebra::RealField;
use num_traits::FromPrimitive;

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
    /// Vessel diameter [m]
    pub diameter: T,
    /// Hematocrit (volume fraction of RBCs) [-]
    pub hematocrit: T,
    /// Plasma viscosity [Pa·s]
    pub plasma_viscosity: T,
}

impl<T: RealField + FromPrimitive + Copy> FahraeuasLindqvist<T> {
    /// Create new Fåhræus-Lindqvist calculator
    pub fn new(diameter: T, hematocrit: T) -> Self {
        Self {
            diameter,
            hematocrit,
            plasma_viscosity: T::from_f64(constants::PLASMA_VISCOSITY_37C)
                .unwrap_or_else(num_traits::Zero::zero),
        }
    }

    /// Check if Fåhræus-Lindqvist effect is significant
    ///
    /// Effect is significant for D < 300 μm
    pub fn is_significant(&self) -> bool {
        self.diameter
            < T::from_f64(constants::FAHRAEUS_LINDQVIST_CRITICAL_DIAMETER)
                .unwrap_or_else(num_traits::Zero::zero)
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
        let one = T::one();

        // Pries (1992) exponent C
        // C = 0.8 + exp(-0.075D) * (-1 + 1/(1 + 10^-11 * D^12))
        let exponent_c = T::from_f64(0.8).unwrap_or_else(num_traits::Zero::zero)
            + (-T::from_f64(0.075).unwrap_or_else(num_traits::Zero::zero) * d_um).exp()
                * (-one
                    + one
                        / (one
                            + T::from_f64(1e-11).unwrap_or_else(num_traits::Zero::zero)
                                * d_um.powf(
                                    T::from_f64(12.0).unwrap_or_else(num_traits::Zero::zero),
                                )));

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
        let one = T::one();

        let d12_term = one
            / (one
                + T::from_f64(1e-11).unwrap_or_else(num_traits::Zero::zero)
                    * d_um.powf(T::from_f64(12.0).unwrap_or_else(num_traits::Zero::zero)));

        // Secomb (2017) exponent C
        let exponent_c = (T::from_f64(0.8).unwrap_or_else(num_traits::Zero::zero)
            + (-T::from_f64(0.075).unwrap_or_else(num_traits::Zero::zero) * d_um).exp())
            * (-one + d12_term)
            + d12_term;

        self.compute_final_relative_viscosity(mu_45, exponent_c)
    }

    /// Backwards compatibility wrapper default to Secomb 2017 (preferred over Pries)
    pub fn relative_viscosity(&self) -> T {
        self.secomb_relative_viscosity()
    }

    /// Helper to compute relative viscosity at 45% hematocrit (μ_45), shared between
    /// Pries (1992) and Secomb (2017) parameterisations.
    /// Caps effective D at 300 μm.
    #[inline]
    fn compute_mu_45(&self) -> (T, T) {
        let mut d_um = self.diameter * T::from_f64(1e6).unwrap_or_else(num_traits::Zero::zero); // Convert to μm
                                                                                                // Fahraeus-Lindqvist scales back to bulk viscosity above 300μm
        let d_max = T::from_f64(300.0).unwrap_or_else(num_traits::Zero::zero);
        if d_um > d_max {
            d_um = d_max;
        }

        // Shared μ_45 fit from Pries et al. (1992)
        // μ_45 = 6 * exp(-0.085D) + 3.2 - 2.44 * exp(-0.06D^0.645)
        let mu_45 = T::from_f64(6.0).unwrap_or_else(num_traits::Zero::zero)
            * (-T::from_f64(0.085).unwrap_or_else(num_traits::Zero::zero) * d_um).exp()
            + T::from_f64(3.2).unwrap_or_else(num_traits::Zero::zero)
            - T::from_f64(2.44).unwrap_or_else(num_traits::Zero::zero)
                * (-T::from_f64(0.06).unwrap_or_else(num_traits::Zero::zero)
                    * d_um.powf(T::from_f64(0.645).unwrap_or_else(num_traits::Zero::zero)))
                .exp();

        (mu_45, d_um)
    }

    /// Calculate final μ_rel avoiding division by zero
    #[inline]
    fn compute_final_relative_viscosity(&self, mu_45: T, exponent_c: T) -> T {
        let one = T::one();
        let ht_clamp = self.hematocrit.max(T::zero());

        let num_base =
            (one - ht_clamp).max(T::from_f64(0.001).unwrap_or_else(num_traits::Zero::zero));
        let den_base = (one - T::from_f64(0.45).unwrap_or_else(num_traits::Zero::zero))
            .max(T::from_f64(0.001).unwrap_or_else(num_traits::Zero::zero));

        let numerator = num_base.powf(exponent_c) - one;
        let denominator = den_base.powf(exponent_c) - one;

        let mu_rel =
            if denominator.abs() < T::from_f64(1e-10).unwrap_or_else(num_traits::Zero::zero) {
                one
            } else {
                one + (mu_45 - one) * (numerator / denominator)
            };

        if mu_rel < one {
            one
        } else {
            mu_rel
        }
    }

    /// Calculate apparent viscosity in microvessel [Pa·s]
    pub fn apparent_viscosity(&self) -> T {
        self.plasma_viscosity * self.relative_viscosity()
    }

    /// Calculate tube hematocrit from feed hematocrit (Fåhræus effect)
    ///
    /// H_tube / H_feed = empirical correlation
    pub fn tube_hematocrit(&self) -> T {
        let d_um = self.diameter * T::from_f64(1e6).unwrap_or_else(num_traits::Zero::zero);
        let one = T::one();

        // Empirical correlation for tube hematocrit reduction
        // Valid for D > 10 μm
        if d_um < T::from_f64(10.0).unwrap_or_else(num_traits::Zero::zero) {
            return self.hematocrit * T::from_f64(0.5).unwrap_or_else(num_traits::Zero::zero);
            // Approximate minimum
        }

        // Pries et al. empirical fit
        let reduction_factor = one
            - T::from_f64(1.7).unwrap_or_else(num_traits::Zero::zero)
                * (-d_um / T::from_f64(40.0).unwrap_or_else(num_traits::Zero::zero)).exp();
        self.hematocrit * reduction_factor
    }
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
