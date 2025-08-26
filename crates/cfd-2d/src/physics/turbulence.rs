//! Turbulence modeling for 2D CFD simulations
//!
//! Implements various turbulence models including:
//! - k-ε model
//! - k-ω SST model
//! - Wall functions for near-wall treatment

use cfd_core::error::Result;
use cfd_core::numeric;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
/// Turbulence model constants
pub mod constants {
    /// von Kármán constant
    pub const KAPPA: f64 = cfd_core::constants::fluid::VON_KARMAN;
    /// Roughness parameter for smooth walls
    pub const E_WALL_FUNCTION: f64 = cfd_core::constants::fluid::WALL_FUNCTION_E;
    /// k-ε model constant Cμ
    pub const C_MU: f64 = 0.09;
    /// k-ε model constant C1ε
    pub const C1_EPSILON: f64 = 1.44;
    /// k-ε model constant C2ε
    pub const C2_EPSILON: f64 = 1.92;
    /// Turbulent Prandtl number for k
    pub const SIGMA_K: f64 = 1.0;
    /// Turbulent Prandtl number for ε
    pub const SIGMA_EPSILON: f64 = 1.3;
    /// Small value for numerical stability
    pub const EPSILON_MIN: f64 = 1e-10;
    /// Y+ threshold for viscous sublayer
    pub const Y_PLUS_VISCOUS_SUBLAYER: f64 = 5.0;
    /// Y+ threshold for log-law region  
    pub const Y_PLUS_LOG_LAW: f64 = cfd_core::constants::fluid::Y_PLUS_LAMINAR;
    /// K-epsilon coefficient for viscous sublayer
    pub const K_VISC_COEFFICIENT: f64 = 11.0;
    /// SST model constant `beta_1`
    pub const SST_BETA_1: f64 = 0.075;
    /// Omega wall coefficient for viscous sublayer
    pub const OMEGA_WALL_COEFFICIENT: f64 = 60.0;
}
/// Wall function types
#[derive(Debug, Clone, Copy)]
pub enum WallFunction {
    /// Standard wall function (log-law)
    Standard,
    /// Blended wall treatment (all y+)
    Blended,
    /// Low-Reynolds number (resolve to wall)
    LowReynolds,
/// k-ε turbulence model
pub struct KEpsilonModel<T: RealField + Copy> {
    /// Turbulent kinetic energy
    pub k: Vec<Vec<T>>,
    /// Turbulent dissipation rate
    pub epsilon: Vec<Vec<T>>,
    /// Turbulent viscosity
    pub nu_t: Vec<Vec<T>>,
    /// Wall function type
    pub wall_function: WallFunction,
    /// Grid dimensions
    nx: usize,
    ny: usize,
impl<T: RealField + Copy + FromPrimitive + Copy> KEpsilonModel<T> {
    /// Create new k-ε model
    #[must_use]
    pub fn new(nx: usize, ny: usize, wall_function: WallFunction) -> Self {
        let k_init = cfd_core::numeric::from_f64(1e-4)?;
        let epsilon_init = cfd_core::numeric::from_f64(1e-6)?;
        Self {
            k: vec![vec![k_init; ny]; nx],
            epsilon: vec![vec![epsilon_init; ny]; nx],
            nu_t: vec![vec![T::zero(); ny]; nx],
            wall_function,
            nx,
            ny,
        }
    }
    /// Update turbulent viscosity
    pub fn update_turbulent_viscosity(&mut self) {
        let c_mu = cfd_core::numeric::from_f64(constants::C_MU)?;
        for i in 0..self.nx {
            for j in 0..self.ny {
                let k_val = self.k[i][j];
                let eps_val = self.epsilon[i][j]
                    .max(cfd_core::numeric::from_f64(constants::EPSILON_MIN)?);
                self.nu_t[i][j] = c_mu * k_val * k_val / eps_val;
            }
    /// Apply wall functions at boundary
    pub fn apply_wall_functions(
        &mut self,
        u_velocity: &[Vec<T>],
        wall_distance: &[Vec<T>],
        nu: T,
    ) -> Result<()> {
        match self.wall_function {
            WallFunction::Standard => {
                self.apply_standard_wall_function(u_velocity, wall_distance, nu)
            WallFunction::Blended => {
                // Collect wall boundary points
                let mut wall_boundaries = Vec::new();
                for i in 0..self.nx {
                    wall_boundaries.push((i, 0)); // Bottom wall
                }
                self.apply_menter_sst_wall_treatment(
                    u_velocity,
                    wall_distance,
                    nu,
                    &wall_boundaries,
                )
            WallFunction::LowReynolds => Ok(()), // No special treatment, resolve to wall
    /// Standard wall function implementation
    /// WARNING: This implementation is hardcoded for walls at j=0 boundary
    /// Note: Current implementation assumes structured rectangular mesh with wall boundaries at j=0 and j=ny-1
    fn apply_standard_wall_function(
        let kappa = cfd_core::numeric::from_f64(constants::KAPPA)?;
        let e_wall_function = cfd_core::numeric::from_f64(constants::E_WALL_FUNCTION)?;
        // Apply at first cell from wall
            // Bottom wall
            let y = wall_distance[i][1];
            let u_p = u_velocity[i][1];
            // Calculate friction velocity using log-law
            let u_tau = self.calculate_friction_velocity(u_p, y, nu)?;
            // Calculate y+
            let y_plus = y * u_tau / nu;
            if y_plus
                > T::from_f64(cfd_core::constants::fluid::Y_PLUS_LAMINAR)
                    .unwrap_or_else(|| T::zero())
            {
                // Log-law region
                // Set k and ε based on equilibrium assumptions
                self.k[i][0] = u_tau * u_tau / c_mu.sqrt();
                self.epsilon[i][0] = u_tau.powi(3) / (kappa * y);
                // Wall shear stress
                let tau_wall = u_tau * u_tau;
                // Set turbulent viscosity at wall
                self.nu_t[i][0] = kappa * u_tau * y - nu;
            } else {
                // Viscous sublayer
                self.k[i][0] = T::zero();
                let two = cfd_core::numeric::from_f64(2.0)?;
                self.epsilon[i][0] = two * nu * self.k[i][1] / (y * y);
                self.nu_t[i][0] = T::zero();
        Ok(())
    /// Apply Menter's k-omega SST wall treatment
    ///
    /// Based on: Menter, F.R. (1994). "Two-equation eddy-viscosity turbulence models
    /// for engineering applications." AIAA Journal, 32(8), 1598-1605.
    /// This implementation uses automatic wall treatment that works for all y+ values:
    /// - For y+ < 5: Viscous sublayer (linear profile)
    /// - For y+ > 30: Log-law region
    /// - For 5 < y+ < 30: Blending between viscous and log-law
    fn apply_menter_sst_wall_treatment(
        wall_boundaries: &[(usize, usize)],
        let beta_star = cfd_core::numeric::from_f64(constants::C_MU)?; // SST model constant
        // Process each wall boundary point
        for &(i_wall, j_wall) in wall_boundaries {
            // Find the nearest interior point
            let (i_near, j_near) = if j_wall == 0 {
                (i_wall, 1) // Wall at bottom
            } else if j_wall == self.ny - 1 {
                (i_wall, self.ny - 2) // Wall at top
            } else if i_wall == 0 {
                (1, j_wall) // Wall at left
            } else if i_wall == self.nx - 1 {
                (self.nx - 2, j_wall) // Wall at right
                continue; // Not a boundary wall
            };
            let y = wall_distance[i_near][j_near];
            let u_p = u_velocity[i_near][j_near];
            // Calculate friction velocity iteratively
            // Menter SST blending function for near-wall treatment
            let arg1 = (y_plus / cfd_core::numeric::from_f64(2.5)?).min(T::one());
            let f1 = arg1.powi(3);
            // k boundary condition (Menter 1994)
            if y_plus < cfd_core::numeric::from_f64(5.0)? {
                // Viscous sublayer: k = 0 at wall
                self.k[i_wall][j_wall] = T::zero();
                self.k[i_near][j_near] =
                    u_tau * u_tau * y_plus / cfd_core::numeric::from_f64(11.0)?;
                // Log-law region: k from equilibrium assumption
                let k_log = u_tau * u_tau / beta_star.sqrt();
                let k_visc = T::zero();
                self.k[i_wall][j_wall] = (T::one() - f1) * k_visc + f1 * k_log;
                self.k[i_near][j_near] = k_log;
            // omega boundary condition (specific dissipation rate)
            // omega_wall = 60 * nu / (beta_1 * y^2) for viscous sublayer
            // omega_wall = u_tau / (sqrt(beta_star) * kappa * y) for log layer
            let beta_1 = cfd_core::numeric::from_f64(0.075)?; // SST model constant
                let omega_visc =
                    cfd_core::numeric::from_f64(60.0)? * nu / (beta_1 * y * y);
                self.epsilon[i_wall][j_wall] = omega_visc * self.k[i_wall][j_wall];
                self.epsilon[i_near][j_near] = omega_visc * self.k[i_near][j_near];
                let omega_log = u_tau / (beta_star.sqrt() * kappa * y);
                self.epsilon[i_wall][j_wall] = omega_log * self.k[i_wall][j_wall];
                self.epsilon[i_near][j_near] = omega_log * self.k[i_near][j_near];
            // Turbulent viscosity with damping
            let rev = u_tau * y / nu; // Reynolds number based on v_tau and y
            let f_mu = T::one() - (-rev / cfd_core::numeric::from_f64(70.0)?).exp();
            self.nu_t[i_wall][j_wall] = T::zero(); // Zero at wall
            self.nu_t[i_near][j_near] = c_mu * self.k[i_near][j_near] * self.k[i_near][j_near]
                / self.epsilon[i_near][j_near]
                * f_mu;
    /// Calculate friction velocity using Currentton-Raphson iteration
    fn calculate_friction_velocity(&self, u_p: T, y: T, nu: T) -> Result<T> {
        let tolerance = cfd_core::numeric::from_f64(1e-6)?;
        let max_iter = 20;
        // Initial guess
        let mut u_tau = cfd_core::numeric::from_f64(0.1)? * u_p;
        for _ in 0..max_iter {
                // Log-law
                let u_plus = u_p / u_tau;
                let f = u_plus - (y_plus.ln() / kappa + e_wall_function.ln() * kappa);
                let df = -u_p / (u_tau * u_tau) - T::one() / (kappa * u_tau);
                let delta = f / df;
                u_tau -= delta;
                if delta.abs() < tolerance {
                    break;
                // Linear law
                u_tau = (u_p * nu / y).sqrt();
                break;
        Ok(u_tau)
    /// Solve k-ε transport equations
    pub fn solve_transport_equations(
        velocity: &[Vec<Vector2<T>>],
        dt: T,
        dx: T,
        dy: T,
        let c1_eps = cfd_core::numeric::from_f64(constants::C1_EPSILON)?;
        let c2_eps = cfd_core::numeric::from_f64(constants::C2_EPSILON)?;
        let sigma_k = cfd_core::numeric::from_f64(constants::SIGMA_K)?;
        let sigma_eps = cfd_core::numeric::from_f64(constants::SIGMA_EPSILON)?;
        let mut current_k = self.k.clone();
        let mut current_epsilon = self.epsilon.clone();
        // Interior points only
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                // Calculate production term
                let du_dx = (velocity[i + 1][j].x - velocity[i - 1][j].x)
                    / (cfd_core::numeric::from_f64(2.0)? * dx);
                let du_dy = (velocity[i][j + 1].x - velocity[i][j - 1].x)
                    / (cfd_core::numeric::from_f64(2.0)? * dy);
                let dv_dx = (velocity[i + 1][j].y - velocity[i - 1][j].y)
                let dv_dy = (velocity[i][j + 1].y - velocity[i][j - 1].y)
                let s11 = du_dx;
                let s12 = cfd_core::numeric::from_f64(0.5)? * (du_dy + dv_dx);
                let s22 = dv_dy;
                let production = cfd_core::numeric::from_f64(2.0)?
                    * self.nu_t[i][j]
                    * (s11 * s11
                        + cfd_core::numeric::from_f64(2.0)? * s12 * s12
                        + s22 * s22);
                // k equation
                let k_diffusion =
                    self.calculate_diffusion(&self.k, i, j, nu + self.nu_t[i][j] / sigma_k, dx, dy);
                current_k[i][j] =
                    self.k[i][j] + dt * (production - self.epsilon[i][j] + k_diffusion);
                // ε equation with semi-implicit treatment for stability
                // Treat destruction term implicitly to avoid singularity when k is small
                let eps_diffusion = self.calculate_diffusion(
                    &self.epsilon,
                    i,
                    j,
                    nu + self.nu_t[i][j] / sigma_eps,
                    dx,
                    dy,
                );
                // Semi-implicit formulation: ε_new = (ε_old + dt * source) / (1 + dt * destruction_coeff)
                let source_term =
                    c1_eps * production * self.epsilon[i][j] / self.k[i][j] + eps_diffusion;
                let destruction_coeff = c2_eps * self.epsilon[i][j]
                    / self.k[i][j]
                        .max(cfd_core::numeric::from_f64(constants::EPSILON_MIN)?);
                current_epsilon[i][j] =
                    (self.epsilon[i][j] + dt * source_term) / (T::one() + dt * destruction_coeff);
                // Ensure positive values
                current_k[i][j] = current_k[i][j]
                current_epsilon[i][j] = current_epsilon[i][j]
        self.k = current_k;
        self.epsilon = current_epsilon;
        self.update_turbulent_viscosity();
    /// Calculate diffusion term using central differences
    fn calculate_diffusion(
        &self,
        field: &[Vec<T>],
        i: usize,
        j: usize,
        diffusivity: T,
    ) -> T {
        let d2f_dx2 = (field[i + 1][j]
            - cfd_core::numeric::from_f64(2.0)? * field[i][j]
            + field[i - 1][j])
            / (dx * dx);
        let d2f_dy2 = (field[i][j + 1]
            + field[i][j - 1])
            / (dy * dy);
        diffusivity * (d2f_dx2 + d2f_dy2)
