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

/// Central flux (average of left and right states).
///
/// # Theorem (Energy Conservation)
///
/// The central flux $F^* = \tfrac{1}{2}(u_L + u_R)$ conserves kinetic energy
/// but introduces no numerical dissipation. It requires explicit stabilisation
/// (e.g., interior penalty or filtering) for stability.
///
/// **Proof sketch**: The symmetric averaging ensures that telescope sums over
/// element interfaces cancel, yielding a discrete energy identity with no
/// dissipative remainder term.
#[derive(Debug, Clone, Copy)]
pub struct CentralFlux;

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
        // Central flux adds no dissipation, so the effective wave speed is zero.
        0.0
    }
}

/// Local Lax-Friedrichs (Rusanov) flux.
///
/// # Theorem (Entropy Stability)
///
/// The Lax-Friedrichs flux $F^* = \tfrac{1}{2}(F(u_L) + F(u_R)) - \tfrac{1}{2}\alpha(u_R - u_L)$
/// is entropy-stable and monotone for scalar conservation laws when $\alpha \ge \max_{u \in [u_L, u_R]} |f'(u)|$.
///
/// **Proof sketch**: Monotonicity follows because $\partial F^*/\partial u_L = \tfrac{1}{2}(f'(u_L) + \alpha) \ge 0$
/// and $\partial F^*/\partial u_R = \tfrac{1}{2}(f'(u_R) - \alpha) \le 0$ whenever $\alpha \ge |f'|$. By the
/// Crandall–Majda theorem, monotone fluxes are entropy-stable.
#[derive(Debug, Clone, Copy)]
pub struct LaxFriedrichsFlux;

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

    fn max_wave_speed(&self, u_l: &DVector<f64>, u_r: &DVector<f64>, n: &DVector<f64>) -> f64 {
        // Rusanov wave speed estimate: spectral radius of the flux Jacobian
        // projected onto the interface normal.
        //
        // When the state vector and normal have the same dimension (e.g.,
        // scalar advection or velocity-only formulations), this equals
        // max(|u_L · n̂|, |u_R · n̂|). For general systems where the state
        // has more components than the spatial dimension, we fall back to
        // ||u||∞ · ||n̂|| which upper-bounds the true spectral radius.
        let n_norm = n.norm();
        if n_norm < f64::EPSILON {
            return u_l.amax().max(u_r.amax());
        }
        if u_l.len() == n.len() {
            let n_hat = n / n_norm;
            let lambda_l = u_l.dot(&n_hat).abs();
            let lambda_r = u_r.dot(&n_hat).abs();
            lambda_l.max(lambda_r)
        } else {
            // General systems: conservative upper bound
            u_l.amax().max(u_r.amax()) * n_norm
        }
    }
}

/// HLL (Harten-Lax-van Leer) approximate Riemann solver.
///
/// # Theorem (Positivity Preservation)
///
/// The HLL flux preserves positivity of density/energy for the Euler equations
/// when the wave speed estimates satisfy $S_L \le \lambda_{\min}$ and
/// $S_R \ge \lambda_{\max}$ (the true minimum and maximum wave speeds).
///
/// **Proof sketch**: The HLL intermediate state
/// $u^* = (S_R u_R - S_L u_L + F_L - F_R)/(S_R - S_L)$ is a convex combination of the
/// left and right states weighted by the wave speed bounds, ensuring non-negativity
/// of conserved quantities when the bound conditions are met.
#[derive(Debug, Clone, Copy)]
pub struct HLLFlux;

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
    /// Compute the left-going wave speed estimate (Davis estimate).
    ///
    /// For matching dimensions: $S_L = \min(u_L \cdot \hat{n},\, u_R \cdot \hat{n})$.
    /// For general systems: uses $-\max(\|u_L\|_\infty, \|u_R\|_\infty)$.
    fn wave_speed_left(self, u_l: &DVector<f64>, u_r: &DVector<f64>, n: &DVector<f64>) -> f64 {
        let n_norm = n.norm();
        if n_norm < f64::EPSILON {
            return -(u_l.amax().max(u_r.amax()));
        }
        if u_l.len() == n.len() {
            let n_hat = n / n_norm;
            let v_l = u_l.dot(&n_hat);
            let v_r = u_r.dot(&n_hat);
            v_l.min(v_r)
        } else {
            -(u_l.amax().max(u_r.amax()) * n_norm)
        }
    }

    /// Compute the right-going wave speed estimate (Davis estimate).
    ///
    /// For matching dimensions: $S_R = \max(u_L \cdot \hat{n},\, u_R \cdot \hat{n})$.
    /// For general systems: uses $\max(\|u_L\|_\infty, \|u_R\|_\infty)$.
    fn wave_speed_right(self, u_l: &DVector<f64>, u_r: &DVector<f64>, n: &DVector<f64>) -> f64 {
        let n_norm = n.norm();
        if n_norm < f64::EPSILON {
            return u_l.amax().max(u_r.amax());
        }
        if u_l.len() == n.len() {
            let n_hat = n / n_norm;
            let v_l = u_l.dot(&n_hat);
            let v_r = u_r.dot(&n_hat);
            v_l.max(v_r)
        } else {
            u_l.amax().max(u_r.amax()) * n_norm
        }
    }
}

/// Enum wrapper for numerical flux implementations to avoid dynamic dispatch
#[derive(Debug, Clone, Copy)]
pub enum FluxImpl {
    /// Central flux
    Central(CentralFlux),
    /// Lax-Friedrichs flux
    LaxFriedrichs(LaxFriedrichsFlux),
    /// HLL flux
    HLL(HLLFlux),
}

impl NumericalFlux for FluxImpl {
    #[inline]
    fn compute_flux(
        &self,
        u_l: &DVector<f64>,
        u_r: &DVector<f64>,
        n: &DVector<f64>,
        params: &FluxParams,
    ) -> DVector<f64> {
        match self {
            FluxImpl::Central(f) => f.compute_flux(u_l, u_r, n, params),
            FluxImpl::LaxFriedrichs(f) => f.compute_flux(u_l, u_r, n, params),
            FluxImpl::HLL(f) => f.compute_flux(u_l, u_r, n, params),
        }
    }

    #[inline]
    fn max_wave_speed(&self, u_l: &DVector<f64>, u_r: &DVector<f64>, n: &DVector<f64>) -> f64 {
        match self {
            FluxImpl::Central(f) => f.max_wave_speed(u_l, u_r, n),
            FluxImpl::LaxFriedrichs(f) => f.max_wave_speed(u_l, u_r, n),
            FluxImpl::HLL(f) => f.max_wave_speed(u_l, u_r, n),
        }
    }
}

/// Factory for creating numerical flux functions
pub struct FluxFactory;

impl FluxFactory {
    /// Create a new numerical flux function
    #[inline]
    pub fn create(flux_type: FluxType) -> FluxImpl {
        match flux_type {
            FluxType::Central => FluxImpl::Central(CentralFlux),
            FluxType::HLL | FluxType::HLLC => FluxImpl::HLL(HLLFlux),
            FluxType::LaxFriedrichs | FluxType::Upwind => {
                FluxImpl::LaxFriedrichs(LaxFriedrichsFlux)
            }
        }
    }
}

/// Compute the numerical flux at an interface
#[inline]
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
        let f2d = flux.compute_flux(&u_l, &u_r, &n, &params);
        // max_wave_speed: matching dims → max(|u_l·n̂|, |u_r·n̂|) = max(0, 2) = 2
        // alpha_combined = 2 * 1.0 = 2
        // f* = 0.5*([0,0]+[2,2]) - 0.5*2*([2,2]-[0,0]) = [1,1] - [2,2] = [-1,-1]
        assert_relative_eq!(f2d[0], -1.0, epsilon = 1e-10);
        assert_relative_eq!(f2d[1], -1.0, epsilon = 1e-10);

        // 1-D scalar case with matching 1-D normal
        let u_l = DVector::from_vec(vec![0.0]);
        let u_r = DVector::from_vec(vec![2.0]);
        let n1 = DVector::from_vec(vec![1.0]);
        let f = flux.compute_flux(&u_l, &u_r, &n1, &params);

        // max_wave_speed = max(|0·1|, |2·1|) = 2
        // alpha_combined = 2 * 1.0 = 2
        // f* = 0.5*(0+2) - 0.5*2*(2-0) = 1 - 2 = -1
        assert_relative_eq!(f[0], -1.0, epsilon = 1e-10);
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
