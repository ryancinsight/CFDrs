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
/// # Empirical Correlation
/// Uses simplified Pries et al. (1994) correlation:
/// ```text
/// μ_rel(D) = 1 + (μ_bulk - 1) · f(D)
///
/// where f(D) accounts for the diameter dependence of the Fåhræus-Lindqvist effect
/// ```
///
/// # Reference
/// Pries, A.R., Neuhaus, D., Gaehtgens, P. (1992) "Blood viscosity in tube flow:
/// dependence on diameter and hematocrit"
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
            plasma_viscosity: T::from_f64(constants::PLASMA_VISCOSITY_37C).unwrap_or_else(num_traits::Zero::zero),
        }
    }

    /// Check if Fåhræus-Lindqvist effect is significant
    ///
    /// Effect is significant for D < 300 μm
    pub fn is_significant(&self) -> bool {
        self.diameter < T::from_f64(constants::FAHRAEUS_LINDQVIST_CRITICAL_DIAMETER).unwrap_or_else(num_traits::Zero::zero)
    }

    /// Calculate relative apparent viscosity using exact Pries et al. (1992) formulation
    ///
    /// Returns μ_rel = μ_app / μ_plasma
    ///
    /// # Exact Physical Model (Pries et al. 1992)
    /// For a vessel of diameter D (in μm) and hematocrit H_t:
    /// - Relative viscosity of plasma: μ_p_rel ≈ 1.0
    /// - Exact formulation for apparent viscosity:
    ///   μ_rel = [1 + (μ_bulk - 1) * ((D/(D-1.1))^2) * (D/(D-1.1))]
    ///
    /// The full in vivo viscosity equation from Pries 1992 incorporates
    /// complex fitting parameters for the cell-free layer.
    pub fn relative_viscosity(&self) -> T {
        let one = T::one();
        let d_um = self.diameter * T::from_f64(1e6).unwrap_or_else(num_traits::Zero::zero); // Convert to μm

        // Exact Pries et al. (1992) in vitro formulation:
        // μ_45 = 6 * exp(-0.85D) + 3.2 - 2.44 * exp(-0.06D^0.645)
        // With hematocrit dependence:
        // μ_rel_vit = 1 + (μ_45 - 1) * (( (1-H_t)^C - 1 ) / ( (1-0.45)^C - 1 ))
        // where C = (0.8 + exp(-0.075D)) * (-1 + 1/(1 + 10^-11 * D^12)) + (1/(1 + 10^-11 * D^12))

        let mu_45 = T::from_f64(6.0).unwrap_or_else(num_traits::Zero::zero) * (-T::from_f64(0.85).unwrap_or_else(num_traits::Zero::zero) * d_um).exp()
                  + T::from_f64(3.2).unwrap_or_else(num_traits::Zero::zero)
                  - T::from_f64(2.44).unwrap_or_else(num_traits::Zero::zero) * (-T::from_f64(0.06).unwrap_or_else(num_traits::Zero::zero) * d_um.powf(T::from_f64(0.645).unwrap_or_else(num_traits::Zero::zero))).exp();

        let exponent_c = (T::from_f64(0.8).unwrap_or_else(num_traits::Zero::zero) + (-T::from_f64(0.075).unwrap_or_else(num_traits::Zero::zero) * d_um).exp())
                       * (-one + one / (one + T::from_f64(1e-11).unwrap_or_else(num_traits::Zero::zero) * d_um.powf(T::from_f64(12.0).unwrap_or_else(num_traits::Zero::zero))))
                       + (one / (one + T::from_f64(1e-11).unwrap_or_else(num_traits::Zero::zero) * d_um.powf(T::from_f64(12.0).unwrap_or_else(num_traits::Zero::zero))));

        let ht_factor = ((one - self.hematocrit).powf(exponent_c) - one)
                      / ((one - T::from_f64(0.45).unwrap_or_else(num_traits::Zero::zero)).powf(exponent_c) - one);

        let mu_rel = one + (mu_45 - one) * ht_factor;

        // Ensure we don't drop below plasma viscosity
        if mu_rel < one { one } else { mu_rel }
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
            return self.hematocrit * T::from_f64(0.5).unwrap_or_else(num_traits::Zero::zero); // Approximate minimum
        }

        // Pries et al. empirical fit
        let reduction_factor =
            one - T::from_f64(1.7).unwrap_or_else(num_traits::Zero::zero) * (-d_um / T::from_f64(40.0).unwrap_or_else(num_traits::Zero::zero)).exp();
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
            "μ(50 μm) = {} should be < μ(100 μm) = {}",
            mu_50,
            mu_100
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
