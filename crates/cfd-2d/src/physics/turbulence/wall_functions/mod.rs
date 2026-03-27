//! Wall function implementations for turbulence models.
//!
//! ## Mathematical Foundation
//!
//! Wall functions provide boundary conditions for turbulence models in high-Reynolds
//! number flows where resolving the viscous sublayer is computationally expensive.
//!
//! ### Log-Law Theory (von Kármán 1930, Prandtl 1925)
//!
//! In equilibrium turbulent boundary layers:
//! ```text
//! u⁺ = (1/κ) ln(y⁺) + B   for y⁺ > 30
//! ```
//!
//! where:
//! - u⁺ = u / u_τ (velocity in wall units)
//! - y⁺ = y u_τ / ν (wall distance in wall units)
//! - κ ≈ 0.41 (von Kármán constant)
//! - B ≈ 5.5 (log-law intercept)
//!
//! ## Literature References
//!
//! - Schlichting (1979): *Boundary Layer Theory*
//! - Spalding, D.B. (1961). *J. Appl. Mech.* 28(3):455–458.
//! - Cebeci & Bradshaw (1977): *Momentum Transfer in Boundary Layers*

/// Wall roughness classification and roughness function computation.
pub mod roughness;

/// Spalding (1961) implicit universal law of the wall.
pub mod spalding;

use super::constants::{
    BLENDING_FACTOR, C_MU, EPSILON_MIN, E_WALL_FUNCTION, KAPPA, OMEGA_WALL_COEFFICIENT,
    Y_PLUS_LOG_LAW, Y_PLUS_VISCOUS_SUBLAYER,
};
use nalgebra::RealField;
use num_traits::FromPrimitive;

pub use roughness::{RoughnessType, WallRoughness};
pub use spalding::spalding_u_plus;

/// Wall function types with enhanced physics.
#[derive(Debug, Clone, Copy)]
pub enum WallFunction {
    /// Standard wall function (log-law) — Launder & Spalding (1974).
    Standard,
    /// Blended wall treatment (all y⁺) — Menter (1994).
    Blended,
    /// Low-Reynolds number (resolve to wall).
    LowReynolds,
    /// Rough wall function with equivalent sand-grain roughness.
    RoughWall,
    /// Werner-Wengle wall function for complex surfaces.
    WernerWengle,
    /// Automatic wall treatment (blends all approaches).
    Automatic,
}

/// Enhanced wall treatment for turbulence models with roughness support.
#[derive(Debug, Clone)]
pub struct WallTreatment<T: RealField + Copy> {
    wall_function: WallFunction,
    kappa: T,
    e_wall: T,
    roughness: WallRoughness<T>,
    /// Validation mode for log-law compliance.
    validation_mode: bool,
}

impl<T: RealField + FromPrimitive + Copy + num_traits::ToPrimitive> WallTreatment<T> {
    /// Create a new wall treatment with smooth wall.
    pub fn new(wall_function: WallFunction) -> Self {
        Self {
            wall_function,
            kappa: T::from_f64(KAPPA).unwrap_or_else(T::one),
            e_wall: T::from_f64(E_WALL_FUNCTION).unwrap_or_else(T::one),
            roughness: WallRoughness::smooth(),
            validation_mode: false,
        }
    }

    /// Create wall treatment with specified roughness.
    pub fn with_roughness(wall_function: WallFunction, roughness: WallRoughness<T>) -> Self {
        Self {
            wall_function,
            kappa: T::from_f64(KAPPA).unwrap_or_else(T::one),
            e_wall: T::from_f64(E_WALL_FUNCTION).unwrap_or_else(T::one),
            roughness,
            validation_mode: false,
        }
    }

    /// Enable log-law validation mode.
    pub fn with_validation(mut self) -> Self {
        self.validation_mode = true;
        self
    }

    /// Get current roughness parameters.
    pub fn roughness(&self) -> &WallRoughness<T> {
        &self.roughness
    }

    /// Update roughness parameters.
    pub fn set_roughness(&mut self, roughness: WallRoughness<T>) {
        self.roughness = roughness;
    }

    /// Calculate u⁺ from y⁺ using wall function with roughness effects.
    pub fn u_plus(&self, y_plus: T) -> T {
        let result = match self.wall_function {
            WallFunction::Standard => self.standard_wall_function(y_plus),
            WallFunction::Blended => self.blended_wall_function(y_plus),
            WallFunction::LowReynolds => y_plus,
            WallFunction::RoughWall => self.rough_wall_function(y_plus),
            WallFunction::WernerWengle => self.werner_wengle_wall_function(y_plus),
            WallFunction::Automatic => self.automatic_wall_function(y_plus),
        };

        if self.validation_mode {
            self.validate_log_law(y_plus, result);
        }

        result
    }

    /// Standard log-law wall function.
    fn standard_wall_function(&self, y_plus: T) -> T {
        let y_visc = T::from_f64(Y_PLUS_VISCOUS_SUBLAYER).unwrap_or_else(T::one);
        let y_log = T::from_f64(Y_PLUS_LOG_LAW).unwrap_or_else(T::one);

        if y_plus <= y_visc {
            y_plus
        } else if y_plus >= y_log {
            (y_plus.ln() / self.kappa) + T::from_f64(5.5).unwrap_or_else(T::one)
        } else {
            let u_visc = y_visc;
            let u_log = (y_log.ln() / self.kappa) + T::from_f64(5.5).unwrap_or_else(T::one);
            let blend = (y_plus - y_visc) / (y_log - y_visc);
            u_visc * (T::one() - blend) + u_log * blend
        }
    }

    /// Blended wall function for all y⁺ values (Reichardt blending).
    fn blended_wall_function(&self, y_plus: T) -> T {
        let blending = T::from_f64(BLENDING_FACTOR).unwrap_or_else(T::one);
        let exp_factor = (-y_plus * blending).exp();

        let u_visc = y_plus;
        let u_log = (T::one() / self.kappa) * ((self.e_wall * y_plus).ln());

        u_visc * exp_factor + u_log * (T::one() - exp_factor)
    }

    /// Rough wall function with equivalent sand-grain roughness.
    ///
    /// Based on Cebeci & Bradshaw (1977) formulation for rough surfaces.
    fn rough_wall_function(&self, y_plus: T) -> T {
        let k_s = self.roughness.equivalent_sand_grain;
        let delta_b = self.roughness.roughness_function;

        let k_s_plus = k_s;

        if k_s_plus <= T::from_f64(0.1).unwrap_or_else(T::one) {
            self.standard_wall_function(y_plus)
        } else {
            let y_visc = T::from_f64(Y_PLUS_VISCOUS_SUBLAYER).unwrap_or_else(T::one);

            if y_plus <= y_visc {
                y_plus
            } else {
                let b_smooth = T::from_f64(5.5).unwrap_or_else(T::one);
                let b_rough = b_smooth - delta_b;
                (T::one() / self.kappa) * ((y_plus + k_s_plus) / k_s_plus).ln() + b_rough
            }
        }
    }

    /// Werner-Wengle wall function for complex surfaces.
    fn werner_wengle_wall_function(&self, y_plus: T) -> T {
        let k_s_plus = self.roughness.equivalent_sand_grain;

        if k_s_plus <= T::from_f64(0.1).unwrap_or_else(T::one) {
            self.blended_wall_function(y_plus)
        } else {
            let a1 = T::from_f64(8.3).unwrap_or_else(T::one);
            let a2 = T::from_f64(1.0 / 7.0).unwrap_or_else(T::one);
            let b = T::from_f64(1.0 / 7.0).unwrap_or_else(T::one);

            let y_eff = T::one() / self.kappa * (y_plus + k_s_plus).ln() - T::one();
            let exp_term = (-a1 * y_eff * y_eff).exp();

            a1 * (T::one() - exp_term)
                + a2 * (T::one() - (-b * y_eff).exp()) * (T::one() - exp_term)
        }
    }

    /// Automatic wall treatment that blends all approaches based on y⁺.
    fn automatic_wall_function(&self, y_plus: T) -> T {
        let y_visc = T::from_f64(Y_PLUS_VISCOUS_SUBLAYER).unwrap_or_else(T::one);
        let y_log = T::from_f64(Y_PLUS_LOG_LAW).unwrap_or_else(T::one);

        if y_plus <= y_visc {
            y_plus
        } else if y_plus >= y_log {
            if self.roughness.equivalent_sand_grain > T::from_f64(0.1).unwrap_or_else(T::one) {
                self.rough_wall_function(y_plus)
            } else {
                self.standard_wall_function(y_plus)
            }
        } else {
            let viscous_part = self.blended_wall_function(y_plus.min(y_visc));
            let log_part =
                if self.roughness.equivalent_sand_grain > T::from_f64(0.1).unwrap_or_else(T::one) {
                    self.rough_wall_function(y_plus.max(y_log))
                } else {
                    self.standard_wall_function(y_plus.max(y_log))
                };

            let blend_factor = (y_plus - y_visc) / (y_log - y_visc);
            viscous_part * (T::one() - blend_factor) + log_part * blend_factor
        }
    }

    /// Validate log-law compliance and issue warnings if needed.
    fn validate_log_law(&self, y_plus: T, u_plus: T) {
        let y_log = T::from_f64(Y_PLUS_LOG_LAW).unwrap_or_else(T::one);

        if y_plus >= y_log {
            let log_law_value =
                (y_plus.ln() / self.kappa) + T::from_f64(5.5).unwrap_or_else(T::one);
            let relative_error = (u_plus - log_law_value).abs() / log_law_value;

            if relative_error > T::from_f64(0.1).unwrap_or_else(T::one) {
                tracing::debug!(
                    "Warning: Wall function deviates from log-law at y+ = {:.1}, error = {:.1}%",
                    y_plus.to_f64().unwrap_or(0.0),
                    relative_error.to_f64().unwrap_or(0.0) * 100.0
                );
            }
        }
    }

    /// Calculate wall shear stress.
    pub fn wall_shear_stress(
        &self,
        velocity_parallel: T,
        wall_distance: T,
        density: T,
        viscosity: T,
    ) -> T {
        let y_plus = self.calculate_y_plus(wall_distance, velocity_parallel, density, viscosity);
        let u_plus = self.u_plus(y_plus);

        if u_plus > T::zero() {
            let u_tau = velocity_parallel / u_plus;
            density * u_tau * u_tau
        } else {
            viscosity * velocity_parallel / wall_distance
        }
    }

    /// Calculate y⁺ dimensionless wall distance.
    pub fn calculate_y_plus(
        &self,
        wall_distance: T,
        velocity_parallel: T,
        density: T,
        viscosity: T,
    ) -> T {
        let u_tau_init = (viscosity * velocity_parallel.abs() / (density * wall_distance)).sqrt();
        let y_plus = density * u_tau_init * wall_distance / viscosity;
        y_plus.max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero))
    }

    /// Calculate turbulent kinetic energy at wall-adjacent cell.
    pub fn wall_k(&self, u_tau: T, _c_mu: T) -> T {
        let c_mu_val = T::from_f64(C_MU).unwrap_or_else(T::one);
        u_tau * u_tau / c_mu_val.sqrt()
    }

    /// Calculate epsilon at wall-adjacent cell.
    pub fn wall_epsilon(&self, u_tau: T, wall_distance: T) -> T {
        let kappa = self.kappa;
        u_tau.powi(3) / (kappa * wall_distance)
    }

    /// Calculate omega at wall.
    pub fn wall_omega(&self, wall_distance: T, viscosity: T, density: T) -> T {
        if let WallFunction::LowReynolds = self.wall_function {
            let coeff = T::from_f64(OMEGA_WALL_COEFFICIENT).unwrap_or_else(T::one);
            let nu = viscosity / density;
            coeff * nu / (wall_distance * wall_distance)
        } else {
            let u_tau = (viscosity / (density * wall_distance)).sqrt();
            let y_plus = density * u_tau * wall_distance / viscosity;

            if y_plus < T::from_f64(Y_PLUS_VISCOUS_SUBLAYER).unwrap_or_else(T::one) {
                let coeff = T::from_f64(OMEGA_WALL_COEFFICIENT).unwrap_or_else(T::one);
                let nu = viscosity / density;
                coeff * nu / (wall_distance * wall_distance)
            } else {
                u_tau / (self.kappa * wall_distance)
            }
        }
    }

    /// Get wall function type.
    pub fn wall_function_type(&self) -> WallFunction {
        self.wall_function
    }
}

#[cfg(test)]
mod tests;
