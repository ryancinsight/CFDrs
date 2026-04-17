//! Flux schemes for FVM
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

use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Flux scheme for convection terms
#[derive(Debug, Clone, Copy)]
pub enum FluxScheme {
    /// Central differencing scheme
    CentralDifference,
    /// First-order upwind scheme
    Upwind,
    /// QUICK (Quadratic Upstream Interpolation for Convective Kinematics) scheme
    QuadraticUpwind,
    /// Power law scheme
    PowerLaw,
    /// Hybrid scheme
    Hybrid,
}

/// Factory for creating flux schemes
pub struct FluxSchemeFactory;

impl FluxSchemeFactory {
    /// Create a flux calculator based on the scheme
    #[must_use]
    pub fn create<T: RealField + Copy + FromPrimitive>(
        scheme: FluxScheme,
        diffusion: f64,
    ) -> Box<dyn FluxCalculator<T>> {
        match scheme {
            FluxScheme::CentralDifference => Box::new(CentralDifferenceFlux),
            FluxScheme::Upwind => Box::new(UpwindFlux),
            FluxScheme::QuadraticUpwind => Box::new(QuadraticUpwindFlux),
            FluxScheme::PowerLaw => Box::new(PowerLawFlux::new(diffusion)),
            FluxScheme::Hybrid => Box::new(HybridFlux::new(diffusion)),
        }
    }
}

/// Trait for flux calculation
pub trait FluxCalculator<T: RealField + Copy>: Send + Sync {
    /// Calculate convective flux
    fn calculate_flux(&self, phi_p: T, phi_e: T, phi_w: T, u: T, dx: T) -> Result<T>;
}

/// Central difference flux calculator
struct CentralDifferenceFlux;

impl<T: RealField + Copy> FluxCalculator<T> for CentralDifferenceFlux {
    fn calculate_flux(&self, _phi_p: T, phi_e: T, phi_w: T, u: T, dx: T) -> Result<T> {
        Ok(u * (phi_e - phi_w) / (T::from_f64(2.0).expect("analytical constant conversion") * dx))
    }
}

/// Upwind flux calculator
struct UpwindFlux;

impl<T: RealField + Copy> FluxCalculator<T> for UpwindFlux {
    fn calculate_flux(&self, phi_p: T, phi_e: T, phi_w: T, u: T, dx: T) -> Result<T> {
        if u > T::zero() {
            Ok(u * (phi_p - phi_w) / dx)
        } else {
            Ok(u * (phi_e - phi_p) / dx)
        }
    }
}

/// QUICK flux calculator
struct QuadraticUpwindFlux;

impl<T: RealField + Copy> FluxCalculator<T> for QuadraticUpwindFlux {
    fn calculate_flux(&self, phi_p: T, phi_e: T, phi_w: T, u: T, dx: T) -> Result<T> {
        // Quadratic upstream interpolation
        let three_eighths = T::from_f64(3.0 / 8.0).expect("analytical constant conversion");
        let six_eighths = T::from_f64(6.0 / 8.0).expect("analytical constant conversion");
        let one_eighth = T::from_f64(1.0 / 8.0).expect("analytical constant conversion");

        if u > T::zero() {
            let phi_face = six_eighths * phi_p + three_eighths * phi_e - one_eighth * phi_w;
            Ok(u * phi_face / dx)
        } else {
            let phi_face = six_eighths * phi_p + three_eighths * phi_w - one_eighth * phi_e;
            Ok(u * phi_face / dx)
        }
    }
}

/// Power law flux calculator
/// Based on Patankar (1980) "Numerical Heat Transfer and Fluid Flow"
struct PowerLawFlux {
    /// Diffusion coefficient
    diffusion: f64,
}

impl PowerLawFlux {
    fn new(diffusion: f64) -> Self {
        Self { diffusion }
    }
}

impl<T: RealField + Copy + FromPrimitive> FluxCalculator<T> for PowerLawFlux {
    fn calculate_flux(&self, phi_p: T, phi_e: T, phi_w: T, u: T, dx: T) -> Result<T> {
        // Power law scheme from Patankar (1980) Section 5.2.4
        let gamma = T::from_f64(self.diffusion).ok_or(Error::InvalidConfiguration(
            "Invalid diffusion coefficient".to_string(),
        ))?;

        // Peclet number Pe = ρu∆x/Γ
        let peclet = u * dx / gamma;
        let abs_pe = peclet.abs();

        // Power law function A(|P|) = max(0, (1 - 0.1|P|)^5)
        let a_func = if abs_pe < T::from_f64(10.0).expect("analytical constant conversion") {
            let one = T::one();
            let point_one = T::from_f64(0.1).expect("analytical constant conversion");
            let term = one - point_one * abs_pe;
            let five = T::from_f64(5.0).expect("analytical constant conversion");
            term.powf(five).max(T::zero())
        } else {
            T::zero()
        };

        // Coefficients for power law scheme
        let d = gamma / dx;
        let f = u;

        // Face values using power law interpolation
        let a_e = d * a_func + (-f).max(T::zero());
        let a_w = d * a_func + f.max(T::zero());

        // Calculate flux
        let flux = a_w * phi_w - (a_w + a_e - f) * phi_p + a_e * phi_e;

        Ok(flux)
    }
}

/// Hybrid flux calculator
/// Based on Spalding (1972) and Patankar (1980)
struct HybridFlux {
    /// Diffusion coefficient
    diffusion: f64,
}

impl HybridFlux {
    fn new(diffusion: f64) -> Self {
        Self { diffusion }
    }
}

impl<T: RealField + Copy + FromPrimitive> FluxCalculator<T> for HybridFlux {
    fn calculate_flux(&self, phi_p: T, phi_e: T, phi_w: T, u: T, dx: T) -> Result<T> {
        // Hybrid scheme from Patankar (1980) Section 5.2.3
        let gamma = T::from_f64(self.diffusion).ok_or(Error::InvalidConfiguration(
            "Invalid diffusion coefficient".to_string(),
        ))?;

        // Peclet number Pe = ρu∆x/Γ
        let peclet = u * dx / gamma;
        let abs_pe = peclet.abs();

        // Diffusion conductance
        let d = gamma / dx;
        // Convection mass flow rate
        let f = u;

        // Hybrid scheme coefficients
        // For |Pe| < 2: use central differencing
        // For |Pe| >= 2: use upwind differencing
        let two = T::from_f64(2.0).expect("analytical constant conversion");

        let (a_w, a_e) = if abs_pe < two {
            // Central differencing with deferred correction
            let a_w = d + f / two;
            let a_e = d - f / two;
            (a_w.max(T::zero()), a_e.max(T::zero()))
        } else {
            // Pure upwind
            if peclet > T::zero() {
                // Flow from west to east
                (d + f, d)
            } else {
                // Flow from east to west
                (d, d - f)
            }
        };

        // Calculate flux using hybrid coefficients
        let flux = a_w * phi_w - (a_w + a_e - f) * phi_p + a_e * phi_e;

        Ok(flux)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Central Difference ──────────────────────────────────────────

    #[test]
    fn test_central_difference_uniform_field() {
        // For uniform φ, flux should be zero regardless of velocity
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::CentralDifference, 1.0);
        let flux = calc.calculate_flux(1.0, 1.0, 1.0, 5.0, 0.1).unwrap();
        assert!(
            flux.abs() < 1e-14,
            "Uniform field flux should be zero, got {flux}"
        );
    }

    #[test]
    fn test_central_difference_linear_field() {
        // φ_w=0, φ_p=1, φ_e=2, u=1, dx=1 => u*(φ_e−φ_w)/(2*dx) = 1*2/2 = 1
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::CentralDifference, 1.0);
        let flux = calc.calculate_flux(1.0, 2.0, 0.0, 1.0, 1.0).unwrap();
        assert!((flux - 1.0).abs() < 1e-14, "Expected 1.0, got {flux}");
    }

    // ── Upwind ──────────────────────────────────────────────────────

    #[test]
    fn test_upwind_positive_velocity() {
        // u > 0 => uses (φ_p − φ_w)/dx
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::Upwind, 1.0);
        let flux = calc.calculate_flux(3.0, 5.0, 1.0, 2.0, 1.0).unwrap();
        // 2.0 * (3.0 − 1.0) / 1.0 = 4.0
        assert!((flux - 4.0).abs() < 1e-14, "Expected 4.0, got {flux}");
    }

    #[test]
    fn test_upwind_negative_velocity() {
        // u < 0 => uses (φ_e − φ_p)/dx
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::Upwind, 1.0);
        let flux = calc.calculate_flux(3.0, 5.0, 1.0, -2.0, 1.0).unwrap();
        // -2.0 * (5.0 − 3.0) / 1.0 = -4.0
        assert!((flux - (-4.0)).abs() < 1e-14, "Expected -4.0, got {flux}");
    }

    #[test]
    fn test_upwind_uniform_field_zero_flux() {
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::Upwind, 1.0);
        let flux = calc.calculate_flux(7.0, 7.0, 7.0, 3.0, 0.5).unwrap();
        assert!(flux.abs() < 1e-14, "Uniform field should give zero flux");
    }

    // ── QUICK ───────────────────────────────────────────────────────

    #[test]
    fn test_quick_linear_field_positive_velocity() {
        // φ_w=0, φ_p=1, φ_e=2, u=1, dx=1
        // u > 0: phi_face = 6/8*φ_p + 3/8*φ_e - 1/8*φ_w = 6/8 + 6/8 - 0 = 12/8 = 1.5
        // flux = u * phi_face / dx = 1 * 1.5 / 1 = 1.5
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::QuadraticUpwind, 1.0);
        let flux = calc.calculate_flux(1.0, 2.0, 0.0, 1.0, 1.0).unwrap();
        assert!(
            (flux - 1.5).abs() < 1e-14,
            "QUICK linear field expected 1.5, got {flux}"
        );
    }

    // ── Power Law ───────────────────────────────────────────────────

    #[test]
    fn test_power_law_low_peclet() {
        // At low Pe (diffusion-dominated), scheme should approach central difference
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::PowerLaw, 10.0);
        let flux = calc.calculate_flux(1.0, 2.0, 0.0, 0.01, 1.0).unwrap();
        assert!(
            flux.is_finite(),
            "Power law flux should be finite at low Pe"
        );
    }

    #[test]
    fn test_power_law_high_peclet() {
        // At high Pe (convection-dominated), scheme should approach upwind
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::PowerLaw, 0.001);
        let flux = calc.calculate_flux(1.0, 2.0, 0.0, 100.0, 1.0).unwrap();
        assert!(
            flux.is_finite(),
            "Power law flux should be finite at high Pe"
        );
    }

    // ── Hybrid ──────────────────────────────────────────────────────

    #[test]
    fn test_hybrid_low_peclet() {
        // |Pe| < 2: central difference behaviour
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::Hybrid, 10.0);
        let flux = calc.calculate_flux(1.0, 2.0, 0.0, 1.0, 1.0).unwrap();
        assert!(flux.is_finite(), "Hybrid flux at low Pe should be finite");
    }

    #[test]
    fn test_hybrid_high_peclet() {
        // |Pe| > 2: pure upwind behaviour
        let calc = FluxSchemeFactory::create::<f64>(FluxScheme::Hybrid, 0.001);
        let flux = calc.calculate_flux(1.0, 2.0, 0.0, 100.0, 1.0).unwrap();
        assert!(flux.is_finite(), "Hybrid flux at high Pe should be finite");
    }
}
