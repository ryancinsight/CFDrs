//! Flux schemes for FVM

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
    /// QUICK scheme
    Quick,
    /// Power law scheme
    PowerLaw,
    /// Hybrid scheme
    Hybrid,
}

/// Factory for creating flux schemes
pub struct FluxSchemeFactory;

impl FluxSchemeFactory {
    /// Create a flux calculator based on the scheme
    pub fn create<T: RealField + Copy>(scheme: FluxScheme) -> Box<dyn FluxCalculator<T>> {
        match scheme {
            FluxScheme::CentralDifference => Box::new(CentralDifferenceFlux),
            FluxScheme::Upwind => Box::new(UpwindFlux),
            FluxScheme::Quick => Box::new(QuickFlux),
            FluxScheme::PowerLaw => Box::new(PowerLawFlux),
            FluxScheme::Hybrid => Box::new(HybridFlux),
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
        Ok(u * (phi_e - phi_w) / (T::from_f64(2.0).unwrap_or_else(T::zero) * dx))
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
struct QuickFlux;

impl<T: RealField + Copy> FluxCalculator<T> for QuickFlux {
    fn calculate_flux(&self, phi_p: T, phi_e: T, phi_w: T, u: T, dx: T) -> Result<T> {
        // Quadratic upstream interpolation
        let three_eighths = T::from_f64(3.0 / 8.0).unwrap_or_else(T::zero);
        let six_eighths = T::from_f64(6.0 / 8.0).unwrap_or_else(T::zero);
        let one_eighth = T::from_f64(1.0 / 8.0).unwrap_or_else(T::zero);

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
        let a_func = if abs_pe < T::from_f64(10.0).unwrap_or_else(T::zero) {
            let one = T::one();
            let point_one = T::from_f64(0.1).unwrap_or_else(T::zero);
            let term = one - point_one * abs_pe;
            let five = T::from_f64(5.0).unwrap_or_else(T::zero);
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
        let two = T::from_f64(2.0).unwrap_or_else(T::zero);

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
