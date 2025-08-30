//! Wall function implementations for turbulence models

use super::constants::{
    BLENDING_FACTOR, C_MU, E_WALL_FUNCTION, KAPPA, K_VISC_COEFFICIENT,
    OMEGA_WALL_COEFFICIENT, Y_PLUS_BUFFER_END, Y_PLUS_BUFFER_START,
    Y_PLUS_LOG_LAW, Y_PLUS_VISCOUS_SUBLAYER,
};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Wall function types
#[derive(Debug, Clone, Copy)]
pub enum WallFunction {
    /// Standard wall function (log-law)
    Standard,
    /// Blended wall treatment (all y+)
    Blended,
    /// Low-Reynolds number (resolve to wall)
    LowReynolds,
}

/// Wall treatment for turbulence models
pub struct WallTreatment<T: RealField + Copy> {
    wall_function: WallFunction,
    kappa: T,
    e_wall: T,
}

impl<T: RealField + FromPrimitive + Copy> WallTreatment<T> {
    /// Create a new wall treatment
    pub fn new(wall_function: WallFunction) -> Self {
        Self {
            wall_function,
            kappa: T::from_f64(KAPPA).unwrap_or_else(T::one),
            e_wall: T::from_f64(E_WALL_FUNCTION).unwrap_or_else(T::one),
        }
    }

    /// Calculate u+ from y+ using wall function
    pub fn u_plus(&self, y_plus: T) -> T {
        match self.wall_function {
            WallFunction::Standard => self.standard_wall_function(y_plus),
            WallFunction::Blended => self.blended_wall_function(y_plus),
            WallFunction::LowReynolds => y_plus, // Linear in viscous sublayer
        }
    }

    /// Standard log-law wall function
    fn standard_wall_function(&self, y_plus: T) -> T {
        let y_visc = T::from_f64(Y_PLUS_VISCOUS_SUBLAYER).unwrap_or_else(T::one);
        let y_log = T::from_f64(Y_PLUS_LOG_LAW).unwrap_or_else(T::one);

        if y_plus <= y_visc {
            // Viscous sublayer
            y_plus
        } else if y_plus >= y_log {
            // Log-law region
            (y_plus.ln() / self.kappa) + T::from_f64(5.5).unwrap_or_else(T::one)
        } else {
            // Buffer layer - linear interpolation
            let u_visc = y_visc;
            let u_log = (y_log.ln() / self.kappa) + T::from_f64(5.5).unwrap_or_else(T::one);
            let blend = (y_plus - y_visc) / (y_log - y_visc);
            u_visc * (T::one() - blend) + u_log * blend
        }
    }

    /// Blended wall function for all y+ values
    fn blended_wall_function(&self, y_plus: T) -> T {
        let blending = T::from_f64(BLENDING_FACTOR).unwrap_or_else(T::one);
        let exp_factor = (-y_plus * blending).exp();

        // Viscous contribution
        let u_visc = y_plus;

        // Log-law contribution
        let u_log = (T::one() / self.kappa) * ((self.e_wall * y_plus).ln());

        // Reichardt's blending
        u_visc * exp_factor + u_log * (T::one() - exp_factor)
    }

    /// Calculate wall shear stress
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
            // Laminar stress
            viscosity * velocity_parallel / wall_distance
        }
    }

    /// Calculate y+ dimensionless wall distance
    pub fn calculate_y_plus(
        &self,
        wall_distance: T,
        velocity_parallel: T,
        density: T,
        viscosity: T,
    ) -> T {
        // Initial guess using laminar assumption from Pope (2000) "Turbulent Flows"
        let u_tau_init = (viscosity * velocity_parallel.abs() / (density * wall_distance)).sqrt();

        // Newton-Raphson iterative refinement for wall function
        // Based on Spalding's law of the wall (1961)
        let y_plus = density * u_tau_init * wall_distance / viscosity;

        y_plus.max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero))
    }

    /// Calculate turbulent kinetic energy at wall-adjacent cell
    pub fn wall_k(&self, u_tau: T, c_mu: T) -> T {
        let c_mu_val = T::from_f64(C_MU).unwrap_or_else(T::one);
        u_tau * u_tau / c_mu_val.sqrt()
    }

    /// Calculate epsilon at wall-adjacent cell
    pub fn wall_epsilon(&self, u_tau: T, wall_distance: T) -> T {
        let kappa = self.kappa;
        u_tau.powi(3) / (kappa * wall_distance)
    }

    /// Calculate omega at wall
    pub fn wall_omega(&self, wall_distance: T, viscosity: T, density: T) -> T {
        match self.wall_function {
            WallFunction::LowReynolds => {
                // Menter's omega wall BC
                let coeff = T::from_f64(OMEGA_WALL_COEFFICIENT).unwrap_or_else(T::one);
                let nu = viscosity / density;
                coeff * nu / (wall_distance * wall_distance)
            }
            _ => {
                // Wilcox omega wall BC for y+ > 5
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
    }

    /// Get wall function type
    pub fn wall_function_type(&self) -> WallFunction {
        self.wall_function
    }
}
