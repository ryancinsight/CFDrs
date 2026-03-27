//! Carreau-Yasuda non-Newtonian viscosity model.
//!
//! # Theorem — Monotonic Shear Thinning
//!
//! The Carreau-Yasuda model describes pseudoplastic (shear-thinning) fluids
//! without the unphysical singularity at zero shear rate present in the
//! simple power-law Ostwald-de Waele model.
//!
//! ```text
//! μ(γ̇) = μ_∞ + (μ_0 - μ_∞) / [1 + (λ γ̇)^a]^((1-n)/a)
//! ```
//!
//! where:
//! - $\mu(γ̇)$ is the apparent dynamic viscosity [Pa·s]
//! - $\mu_0$ is the zero-shear viscosity limit [Pa·s]
//! - $\mu_\infty$ is the infinite-shear viscosity limit [Pa·s]
//! - $\lambda$ is the characteristic relaxation time [s]
//! - $n$ is the flow behavior index (power law index for shear-thinning < 1)
//! - $a$ is the Yasuda index controlling the transition curve
//! - $γ̇ = \sqrt{2 \mathbf{S}:\mathbf{S}}$ is the shear rate magnitude (2nd invariant)
//!
//! **Proof of monotonicity**:
//! For shear-thinning fluids, $n < 1$. $\mu_0 > \mu_\infty$.
//! Let $x = (\lambda γ̇)^a \ge 0$. The denominator $D(x) = (1+x)^{(1-n)/a}$
//! is strictly increasing with $x$ because $(1-n)/a > 0$.
//! Thus, $1/D(x)$ is strictly decreasing. Consequently, the apparent
//! viscosity monotonically decreases from $\mu_0$ as $γ̇ \to 0$, towards
//! $\mu_\infty$ as $γ̇ \to \infty$, correctly resolving the bounded limits.
//!
//! # References
//! - Yasuda, K., Armstrong, R.C., & Cohen, R.E. (1981). Shear flow properties
//!   of concentrated solutions of linear and star branched polystyrenes.
//!   *Rheol. Acta* 20:163-178.
//! - Boyd, J. et al. (2007). A common-sense approach to blood rheology.
//!   *Biophys. J.* 92:1565.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Carreau-Yasuda model parameters.
#[derive(Debug, Clone, Copy)]
pub struct CarreauYasudaModel<T: RealField + Copy> {
    /// Zero-shear viscosity limit $\mu_0$ [Pa·s]
    pub mu_0: T,
    /// Infinite-shear viscosity limit $\mu_\infty$ [Pa·s]
    pub mu_inf: T,
    /// Relaxation time $\lambda$ [s]
    pub lambda: T,
    /// Yasuda transition parameter $a$ (dimensionless)
    pub a: T,
    /// Power-law flow behavior index $n$ (dimensionless). $n<1$ for shear-thinning.
    pub n: T,
}

impl<T: RealField + Copy + Float + FromPrimitive> CarreauYasudaModel<T> {
    /// Standard Carreau-Yasuda parameters for bulk human blood
    ///
    /// Derived from Leuprecht & Perktold (2001) for healthy blood at 37°C.
    #[must_use]
    pub fn typical_blood() -> Self {
        Self {
            mu_0: T::from_f64(0.022).unwrap_or_else(num_traits::Zero::zero),
            mu_inf: T::from_f64(0.0022).unwrap_or_else(num_traits::Zero::zero),
            lambda: T::from_f64(0.11).unwrap_or_else(num_traits::Zero::zero),
            a: T::from_f64(0.644).unwrap_or_else(num_traits::Zero::zero),
            n: T::from_f64(0.392).unwrap_or_else(num_traits::Zero::zero),
        }
    }

    /// Compute the apparent dynamic viscosity given a local shear rate magnitude.
    ///
    /// # Arguments
    /// * `shear_rate` - Shear rate magnitude $|\dot{\gamma}|$ [1/s].
    ///                  Must be $\ge 0$.
    ///
    /// # Returns
    /// Apparent dynamic viscosity $\mu_{app}$ [Pa·s]. bounded in $[\mu_\infty, \mu_0]$.
    #[inline]
    #[must_use]
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        let sr = Float::max(shear_rate, T::zero());
        let one = T::one();

        // [1 + (λ γ̇)^a]
        let base = one + Float::powf(self.lambda * sr, self.a);
        // Exponent: (1-n)/a
        let exp = (one - self.n) / self.a;

        self.mu_inf + (self.mu_0 - self.mu_inf) * Float::powf(base, -exp)
    }

    /// Compute the kinematic viscosity $\nu = \mu / \rho$.
    #[inline]
    #[must_use]
    pub fn apparent_kinematic_viscosity(&self, shear_rate: T, density: T) -> T {
        self.apparent_viscosity(shear_rate) / density
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Theorem: The Carreau-Yasuda model is monotonically shear-thinning
    /// for bounded inputs ($n < 1$).
    #[test]
    fn monotonicity_and_bounds() {
        let model = CarreauYasudaModel::typical_blood();
        let rates = [
            0.0_f64, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0,
        ];

        let mut prev_visc = model.apparent_viscosity(rates[0]);

        // At shear_rate = 0, viscosity is exactly mu_0
        assert!(
            (prev_visc - model.mu_0).abs() < 1e-10,
            "Zero shear viscosity failed bounds"
        );

        for &sr in &rates[1..] {
            let visc = model.apparent_viscosity(sr);

            // Monotonically decreasing
            assert!(
                visc <= prev_visc,
                "Non-monotonic at shear rate {sr}: prev={prev_visc}, curr={visc}"
            );

            // Bounded by limits
            assert!(
                visc >= model.mu_inf && visc <= model.mu_0,
                "Viscosity {visc} out of bounds [{inf}, {zero}]",
                inf = model.mu_inf,
                zero = model.mu_0
            );

            prev_visc = visc;
        }

        // At very high shear rates, viscosity approaches mu_inf
        let inf_visc = model.apparent_viscosity(1e6);
        assert!(
            (inf_visc - model.mu_inf).abs() < 1e-3,
            "Infinite shear viscosity failed bounds: got {inf_visc}, expected {}",
            model.mu_inf
        );
    }

    /// Newtonian limit: when n = 1, the exponent (1−n)/a = 0, so
    /// the denominator is (1+x)^0 = 1, and μ = μ_∞ + (μ_0 − μ_∞) = μ_0
    /// for all shear rates. The fluid behaves as Newtonian.
    #[test]
    fn newtonian_limit_n_equals_one() {
        use approx::assert_relative_eq;
        let model = CarreauYasudaModel::<f64> {
            mu_0: 0.022,
            mu_inf: 0.0022,
            lambda: 0.11,
            a: 0.644,
            n: 1.0, // Newtonian limit
        };
        for &sr in &[0.0, 0.1, 1.0, 10.0, 100.0, 1e6] {
            assert_relative_eq!(
                model.apparent_viscosity(sr),
                model.mu_0,
                epsilon = 1e-10,
            );
        }
    }

    /// Kinematic viscosity: ν = μ/ρ, verified for typical blood parameters.
    #[test]
    fn kinematic_viscosity_consistency() {
        use approx::assert_relative_eq;
        let model = CarreauYasudaModel::<f64>::typical_blood();
        let rho = 1060.0_f64; // blood density
        let sr = 100.0_f64;
        let mu = model.apparent_viscosity(sr);
        let nu = model.apparent_kinematic_viscosity(sr, rho);
        assert_relative_eq!(nu, mu / rho, epsilon = 1e-15);
    }
}
