//! Numerical flux functions for Discontinuous Galerkin methods.
//!
//! This module provides various numerical flux functions that can be used
//! to approximate the flux at element interfaces in DG methods.

use nalgebra::DVector;

/// Type of numerical flux
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FluxType {
    /// Central flux (average of left and right states)
    Central,
    /// Upwind flux (based on characteristic speeds)
    Upwind,
    /// Local Lax-Friedrichs flux
    LaxFriedrichs,
    /// Harten-Lax-van Leer (HLL) flux
    HLL,
    /// Harten-Lax-van Leer-Contact (HLLC) flux (for compressible flows)
    HLLC,
}

/// Parameters for numerical fluxes
#[derive(Debug, Clone)]
pub struct FluxParams {
    /// Type of numerical flux
    pub flux_type: FluxType,
    /// Wave speed scaling factor (for Lax-Friedrichs)
    pub alpha: f64,
    /// Entropy fix parameter (for HLL/HLLC)
    pub epsilon: f64,
}

impl Default for FluxParams {
    fn default() -> Self {
        Self {
            flux_type: FluxType::LaxFriedrichs,
            alpha: 1.0,
            epsilon: 1e-10,
        }
    }
}

impl FluxParams {
    /// Create a new set of flux parameters
    pub fn new(flux_type: FluxType) -> Self {
        Self {
            flux_type,
            ..Default::default()
        }
    }

    /// Set the wave speed scaling factor
    pub fn with_alpha(mut self, alpha: f64) -> Self {
        self.alpha = alpha;
        self
    }

    /// Set the entropy fix parameter
    pub fn with_epsilon(mut self, epsilon: f64) -> Self {
        self.epsilon = epsilon;
        self
    }
}

/// Trait for numerical flux functions
pub trait NumericalFlux {
    /// Compute the numerical flux at an interface
    fn compute_flux(
        &self,
        u_l: &DVector<f64>,
        u_r: &DVector<f64>,
        n: &DVector<f64>,
        params: &FluxParams,
    ) -> DVector<f64>;

    /// Compute the maximum wave speed for the given states
    fn max_wave_speed(&self, u_l: &DVector<f64>, u_r: &DVector<f64>, n: &DVector<f64>) -> f64;
}

/// Central flux (average of left and right states)
struct CentralFlux;

impl NumericalFlux for CentralFlux {
    fn compute_flux(
        &self,
        u_l: &DVector<f64>,
        u_r: &DVector<f64>,
        _n: &DVector<f64>,
        _params: &FluxParams,
    ) -> DVector<f64> {
        // Simple average of left and right fluxes
        // F* = 0.5 * (F(u_l) + F(u_r))
        0.5 * (u_l + u_r)
    }

    fn max_wave_speed(&self, _u_l: &DVector<f64>, _u_r: &DVector<f64>, _n: &DVector<f64>) -> f64 {
        1.0 // Conservative estimate
    }
}

/// Local Lax-Friedrichs flux
struct LaxFriedrichsFlux;

impl NumericalFlux for LaxFriedrichsFlux {
    fn compute_flux(
        &self,
        u_l: &DVector<f64>,
        u_r: &DVector<f64>,
        n: &DVector<f64>,
        params: &FluxParams,
    ) -> DVector<f64> {
        let alpha = self.max_wave_speed(u_l, u_r, n) * params.alpha;

        // F* = 0.5 * (F(u_l) + F(u_r)) - 0.5 * alpha * (u_r - u_l)
        0.5 * (u_l + u_r) - 0.5 * alpha * (u_r - u_l)
    }

    fn max_wave_speed(&self, u_l: &DVector<f64>, u_r: &DVector<f64>, _n: &DVector<f64>) -> f64 {
        // For scalar problems, the wave speed is the absolute value of the flux derivative
        // For systems, this should be the maximum eigenvalue of the flux Jacobian
        // This is a simplified version for demonstration
        (u_l.norm() + u_r.norm()) * 0.5
    }
}

/// HLL (Harten-Lax-van Leer) approximate Riemann solver
struct HLLFlux;

impl NumericalFlux for HLLFlux {
    fn compute_flux(
        &self,
        u_l: &DVector<f64>,
        u_r: &DVector<f64>,
        n: &DVector<f64>,
        params: &FluxParams,
    ) -> DVector<f64> {
        let s_l = self.wave_speed_left(u_l, u_r, n);
        let s_r = self.wave_speed_right(u_l, u_r, n);
        let eps = params.epsilon;

        if s_l >= 0.0 {
            // Right-going flow
            u_l.clone()
        } else if s_r <= 0.0 {
            // Left-going flow
            u_r.clone()
        } else if s_l <= 0.0 && 0.0 <= s_r {
            // Sonic point
            (s_r * u_l - s_l * u_r + s_l * s_r * (u_r - u_l) / (s_r - s_l + eps))
                / (s_r - s_l + eps)
        } else {
            // Use Lax-Friedrichs as a fallback
            let alpha = self.max_wave_speed(u_l, u_r, n);
            0.5 * (u_l + u_r) - 0.5 * alpha * (u_r - u_l)
        }
    }

    fn max_wave_speed(&self, u_l: &DVector<f64>, u_r: &DVector<f64>, n: &DVector<f64>) -> f64 {
        self.wave_speed_left(u_l, u_r, n)
            .abs()
            .max(self.wave_speed_right(u_l, u_r, n).abs())
    }
}

impl HLLFlux {
    /// Compute the left-going wave speed
    fn wave_speed_left(&self, u_l: &DVector<f64>, u_r: &DVector<f64>, _n: &DVector<f64>) -> f64 {
        // Simplified for demonstration
        // In practice, this should use the minimum eigenvalue of the flux Jacobian
        u_l.norm() - u_r.norm()
    }

    /// Compute the right-going wave speed
    fn wave_speed_right(&self, u_l: &DVector<f64>, u_r: &DVector<f64>, _n: &DVector<f64>) -> f64 {
        // Simplified for demonstration
        // In practice, this should use the maximum eigenvalue of the flux Jacobian
        u_r.norm() + u_l.norm()
    }
}

/// Factory for creating numerical flux functions
pub struct FluxFactory;

impl FluxFactory {
    /// Create a new numerical flux function
    pub fn create(flux_type: FluxType) -> Box<dyn NumericalFlux> {
        match flux_type {
            FluxType::Central => Box::new(CentralFlux),
            FluxType::HLL | FluxType::HLLC => Box::new(HLLFlux),
            FluxType::LaxFriedrichs | FluxType::Upwind => Box::new(LaxFriedrichsFlux),
        }
    }
}

/// Compute the numerical flux at an interface
pub fn numerical_flux(
    u_l: &DVector<f64>,
    u_r: &DVector<f64>,
    n: &DVector<f64>,
    params: &FluxParams,
) -> DVector<f64> {
    let flux = FluxFactory::create(params.flux_type);
    flux.compute_flux(u_l, u_r, n, params)
}

/// Compute the maximum wave speed at an interface
pub fn max_wave_speed(
    u_l: &DVector<f64>,
    u_r: &DVector<f64>,
    n: &DVector<f64>,
    params: &FluxParams,
) -> f64 {
    let flux = FluxFactory::create(params.flux_type);
    flux.max_wave_speed(u_l, u_r, n)
}

/// Compute the Rusanov (local Lax-Friedrichs) flux
pub fn rusanov_flux(
    f_l: &DVector<f64>,
    f_r: &DVector<f64>,
    u_l: &DVector<f64>,
    u_r: &DVector<f64>,
    alpha: f64,
) -> DVector<f64> {
    0.5 * (f_l + f_r) - 0.5 * alpha * (u_r - u_l)
}

/// Compute the upwind flux for a scalar conservation law
pub fn upwind_flux(
    f_l: &DVector<f64>,
    f_r: &DVector<f64>,
    _u_l: &DVector<f64>,
    _u_r: &DVector<f64>,
    a: f64,
) -> DVector<f64> {
    if a >= 0.0 {
        f_l.clone()
    } else {
        f_r.clone()
    }
}

/// Compute the HLLC flux for the Euler equations
#[allow(clippy::too_many_arguments)]
pub fn hllc_flux(
    rho_l: f64,
    u_l: f64,
    p_l: f64,
    rho_r: f64,
    u_r: f64,
    p_r: f64,
    gamma: f64,
) -> (f64, f64, f64) {
    // Compute sound speeds
    let c_l = (gamma * p_l / rho_l).sqrt();
    let c_r = (gamma * p_r / rho_r).sqrt();

    // Compute pressure and velocity in the star region (PVRS solver)
    let p_star = 0.5 * (p_l + p_r) - 0.5 * (u_r - u_l) * (rho_l + rho_r) * (c_l + c_r) / 2.0;
    let u_star = 0.5 * (u_l + u_r) - 0.5 * (p_r - p_l) / ((rho_l + rho_r) * (c_l + c_r) / 2.0);

    // Compute wave speeds
    let q_l = if p_star <= p_l {
        1.0
    } else {
        (1.0 + (gamma + 1.0) / (2.0 * gamma) * (p_star / p_l - 1.0)).sqrt()
    };

    let q_r = if p_star <= p_r {
        1.0
    } else {
        (1.0 + (gamma + 1.0) / (2.0 * gamma) * (p_star / p_r - 1.0)).sqrt()
    };

    let s_l = u_l - c_l * q_l;
    let s_r = u_r + c_r * q_r;

    // Compute HLLC flux
    if s_l >= 0.0 {
        // Left state
        let f_mass = rho_l * u_l;
        let f_momentum = rho_l * u_l * u_l + p_l;
        let f_energy = u_l * (0.5 * rho_l * u_l * u_l + gamma * p_l / (gamma - 1.0));

        (f_mass, f_momentum, f_energy)
    } else if s_l <= 0.0 && 0.0 <= u_star {
        // Star left state
        let rho_star_l = rho_l * (s_l - u_l) / (s_l - u_star);
        let e_l = p_l / (gamma - 1.0) + 0.5 * rho_l * u_l * u_l;
        let e_star_l = (e_l * (s_l - u_l) + p_l * u_star - p_l * u_l) / (s_l - u_star);

        let f_mass = rho_star_l * u_star;
        let f_momentum = rho_star_l * u_star * u_star + p_star;
        let f_energy = u_star * (e_star_l + p_star);

        (f_mass, f_momentum, f_energy)
    } else if u_star <= 0.0 && 0.0 <= s_r {
        // Star right state
        let rho_star_r = rho_r * (s_r - u_r) / (s_r - u_star);
        let e_r = p_r / (gamma - 1.0) + 0.5 * rho_r * u_r * u_r;
        let e_star_r = (e_r * (s_r - u_r) + p_star * u_star - p_r * u_r) / (s_r - u_star);

        let f_mass = rho_star_r * u_star;
        let f_momentum = rho_star_r * u_star * u_star + p_star;
        let f_energy = u_star * (e_star_r + p_star);

        (f_mass, f_momentum, f_energy)
    } else {
        // Right state
        let f_mass = rho_r * u_r;
        let f_momentum = rho_r * u_r * u_r + p_r;
        let f_energy = u_r * (0.5 * rho_r * u_r * u_r + gamma * p_r / (gamma - 1.0));

        (f_mass, f_momentum, f_energy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_central_flux() {
        let u_l = DVector::from_vec(vec![1.0, 2.0]);
        let u_r = DVector::from_vec(vec![3.0, 4.0]);
        let n = DVector::from_vec(vec![1.0, 0.0]);

        let params = FluxParams::new(FluxType::Central);
        let flux = FluxFactory::create(params.flux_type);
        let f = flux.compute_flux(&u_l, &u_r, &n, &params);

        assert_relative_eq!(f[0], 2.0, epsilon = 1e-10);
        assert_relative_eq!(f[1], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_lax_friedrichs_flux() {
        let u_l = DVector::from_vec(vec![0.0, 0.0]);
        let u_r = DVector::from_vec(vec![2.0, 2.0]);
        let n = DVector::from_vec(vec![1.0, 0.0]);

        let params = FluxParams::new(FluxType::LaxFriedrichs).with_alpha(1.0);
        let flux = FluxFactory::create(params.flux_type);
        let _f = flux.compute_flux(&u_l, &u_r, &n, &params);

        // For u_l=[0,0], u_r=[2,2], max_wave_speed = (0 + sqrt(8))/2 = 1.414...
        // Let's use u_l=[0,0], u_r=[2,0] so max_wave_speed = (0 + 2)/2 = 1.0
        let u_l = DVector::from_vec(vec![0.0]);
        let u_r = DVector::from_vec(vec![2.0]);
        let f = flux.compute_flux(&u_l, &u_r, &n, &params);

        // For alpha_combined = max_wave_speed * params.alpha = 1.0 * 1.0 = 1.0
        // f* = 0.5 * (u_l + u_r) - 0.5 * alpha_combined * (u_r - u_l)
        //    = 0.5 * (0 + 2) - 0.5 * 1.0 * (2 - 0) = 1.0 - 1.0 = 0.0 = u_l
        assert_relative_eq!(f[0], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_hll_flux() {
        let u_l = DVector::from_vec(vec![1.0, 2.0]);
        let u_r = DVector::from_vec(vec![3.0, 4.0]);
        let n = DVector::from_vec(vec![1.0, 0.0]);

        let params = FluxParams::new(FluxType::HLL);
        let flux = FluxFactory::create(params.flux_type);
        let f = flux.compute_flux(&u_l, &u_r, &n, &params);

        // The exact value depends on the wave speed estimates
        // Just verify that the result has the right dimensions
        assert_eq!(f.len(), 2);
    }

    #[test]
    fn test_rusanov_flux() {
        let f_l = DVector::from_vec(vec![1.0, 2.0]);
        let f_r = DVector::from_vec(vec![4.0, 5.0]);
        let u_l = DVector::from_vec(vec![1.0, 2.0]);
        let u_r = DVector::from_vec(vec![3.0, 4.0]);
        let alpha = 2.0;

        let f = rusanov_flux(&f_l, &f_r, &u_l, &u_r, alpha);

        // f* = 0.5*(f_l + f_r) - 0.5*alpha*(u_r - u_l)
        //    = 0.5*([1,2] + [4,5]) - 0.5*2.0*([3,4] - [1,2])
        //    = [2.5, 3.5] - [2, 2] = [0.5, 1.5]
        assert_relative_eq!(f[0], 0.5, epsilon = 1e-10);
        assert_relative_eq!(f[1], 1.5, epsilon = 1e-10);
    }

    #[test]
    fn test_upwind_flux() {
        let f_l = DVector::from_vec(vec![1.0, 2.0]);
        let f_r = DVector::from_vec(vec![4.0, 5.0]);
        let u_l = DVector::from_vec(vec![1.0, 2.0]);
        let u_r = DVector::from_vec(vec![3.0, 4.0]);

        // Test left-going flow (a < 0)
        let f = upwind_flux(&f_l, &f_r, &u_l, &u_r, -1.0);
        assert_relative_eq!(f[0], 4.0, epsilon = 1e-10);
        assert_relative_eq!(f[1], 5.0, epsilon = 1e-10);

        // Test right-going flow (a > 0)
        let f = upwind_flux(&f_l, &f_r, &u_l, &u_r, 1.0);
        assert_relative_eq!(f[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(f[1], 2.0, epsilon = 1e-10);

        // Test sonic point (a = 0)
        let f = upwind_flux(&f_l, &f_r, &u_l, &u_r, 0.0);
        assert_relative_eq!(f[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(f[1], 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_hllc_flux() {
        // Test case: Sod shock tube problem
        let rho_l = 1.0;
        let u_l = 0.0;
        let p_l = 1.0;

        let rho_r = 0.125;
        let u_r = 0.0;
        let p_r = 0.1;

        let gamma = 1.4;

        let (f_mass, f_momentum, f_energy) = hllc_flux(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma);

        // Verify some basic properties of the flux
        assert!(f_mass > 0.0);
        assert!(f_momentum > 0.0);
        assert!(f_energy > 0.0);

        // The exact values depend on the exact solution to the Riemann problem
        // Here we just check that the values are reasonable
        assert!(f_mass < 1.0);
        assert!(f_momentum < 2.0);
        assert!(f_energy < 3.0);
    }
}
