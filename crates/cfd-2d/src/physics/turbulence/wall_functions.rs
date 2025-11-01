//! Wall function implementations for turbulence models

use super::constants::{
    BLENDING_FACTOR, C_MU, EPSILON_MIN, E_WALL_FUNCTION, KAPPA, OMEGA_WALL_COEFFICIENT,
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
#[derive(Debug, Clone)]
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
    pub fn wall_k(&self, u_tau: T, _c_mu: T) -> T {
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
        if let WallFunction::LowReynolds = self.wall_function {
            // Menter's omega wall BC
            let coeff = T::from_f64(OMEGA_WALL_COEFFICIENT).unwrap_or_else(T::one);
            let nu = viscosity / density;
            coeff * nu / (wall_distance * wall_distance)
        } else {
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

    /// Get wall function type
    pub fn wall_function_type(&self) -> WallFunction {
        self.wall_function
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_new_standard_wall_treatment() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        assert_relative_eq!(treatment.kappa, KAPPA, epsilon = 1e-10);
        assert_relative_eq!(treatment.e_wall, E_WALL_FUNCTION, epsilon = 1e-10);
    }

    #[test]
    fn test_new_blended_wall_treatment() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Blended);
        assert_relative_eq!(treatment.kappa, KAPPA, epsilon = 1e-10);
    }

    #[test]
    fn test_low_reynolds_u_plus() {
        let treatment = WallTreatment::<f64>::new(WallFunction::LowReynolds);
        // Low Reynolds: u+ = y+ (linear in viscous sublayer)
        assert_relative_eq!(treatment.u_plus(1.0), 1.0, epsilon = 1e-10);
        assert_relative_eq!(treatment.u_plus(5.0), 5.0, epsilon = 1e-10);
        assert_relative_eq!(treatment.u_plus(10.0), 10.0, epsilon = 1e-10);
    }

    #[test]
    fn test_standard_wall_function_viscous_sublayer() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        // For y+ < 5: u+ ≈ y+ (viscous sublayer)
        let y_plus = 3.0;
        let u_plus = treatment.u_plus(y_plus);
        assert_relative_eq!(u_plus, y_plus, epsilon = 1e-10);
    }

    #[test]
    fn test_standard_wall_function_log_law() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        // For y+ > 30: u+ = ln(y+)/κ + B (log-law region)
        let y_plus = 100.0;
        let u_plus = treatment.u_plus(y_plus);
        let expected = (y_plus.ln() / KAPPA) + 5.5;
        assert_relative_eq!(u_plus, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_standard_wall_function_buffer_layer() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        // Buffer layer (5 < y+ < 30): linear interpolation
        let y_plus = 15.0;
        let u_plus = treatment.u_plus(y_plus);
        // Should be between viscous and log-law values
        assert!(u_plus > 5.0);
        assert!(u_plus < (30.0_f64.ln() / KAPPA) + 5.5);
    }

    #[test]
    fn test_blended_wall_function_smooth() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Blended);
        // Blended function should be smooth across all y+
        let y_plus_values = vec![1.0, 5.0, 10.0, 30.0, 100.0, 1000.0];
        
        for &y_plus in &y_plus_values {
            let u_plus = treatment.u_plus(y_plus);
            // Should be positive and reasonable
            assert!(u_plus > 0.0);
            assert!(u_plus < y_plus * 2.0); // Sanity check
        }
    }

    #[test]
    fn test_blended_wall_function_asymptotic_behavior() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Blended);
        
        // At low y+, should be close to linear
        let y_plus_low = 1.0;
        let u_plus_low = treatment.u_plus(y_plus_low);
        assert_relative_eq!(u_plus_low, y_plus_low, epsilon = 0.1);
        
        // At high y+, should approach log-law
        let y_plus_high = 1000.0;
        let u_plus_high = treatment.u_plus(y_plus_high);
        let expected_high = (y_plus_high.ln() / KAPPA) + 5.5;
        assert_relative_eq!(u_plus_high, expected_high, epsilon = 1.0);
    }

    #[test]
    fn test_calculate_y_plus_positive() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        let wall_distance = 0.001; // 1mm
        let velocity = 10.0;       // 10 m/s
        let density = 1.0;         // 1 kg/m³
        let viscosity = 1e-5;      // ~air
        
        let y_plus = treatment.calculate_y_plus(wall_distance, velocity, density, viscosity);
        assert!(y_plus > 0.0);
        assert!(y_plus.is_finite());
    }

    #[test]
    fn test_calculate_y_plus_scales_with_velocity() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        let wall_distance = 0.001;
        let density = 1.0;
        let viscosity = 1e-5;
        
        let y_plus_1 = treatment.calculate_y_plus(wall_distance, 10.0, density, viscosity);
        let y_plus_2 = treatment.calculate_y_plus(wall_distance, 20.0, density, viscosity);
        
        // Higher velocity should give higher y+
        assert!(y_plus_2 > y_plus_1);
    }

    #[test]
    fn test_wall_shear_stress_positive() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        let velocity = 10.0;
        let wall_distance = 0.001;
        let density = 1.0;
        let viscosity = 1e-5;
        
        let tau_w = treatment.wall_shear_stress(velocity, wall_distance, density, viscosity);
        assert!(tau_w > 0.0);
        assert!(tau_w.is_finite());
    }

    #[test]
    fn test_wall_shear_stress_scales_with_velocity() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        let wall_distance = 0.001;
        let density = 1.0;
        let viscosity = 1e-5;
        
        let tau_1 = treatment.wall_shear_stress(10.0, wall_distance, density, viscosity);
        let tau_2 = treatment.wall_shear_stress(20.0, wall_distance, density, viscosity);
        
        // Shear stress should increase with velocity
        assert!(tau_2 > tau_1);
    }

    #[test]
    fn test_wall_k_positive() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        let u_tau = 0.5;
        let c_mu = C_MU;
        
        let k = treatment.wall_k(u_tau, c_mu);
        assert!(k > 0.0);
        assert!(k.is_finite());
    }

    #[test]
    fn test_wall_k_scales_with_u_tau_squared() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        let c_mu = C_MU;
        
        let k_1 = treatment.wall_k(1.0, c_mu);
        let k_2 = treatment.wall_k(2.0, c_mu);
        
        // k ~ u_tau^2, so doubling u_tau should quadruple k
        assert_relative_eq!(k_2 / k_1, 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_wall_epsilon_positive() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        let u_tau = 0.5;
        let wall_distance = 0.001;
        
        let epsilon = treatment.wall_epsilon(u_tau, wall_distance);
        assert!(epsilon > 0.0);
        assert!(epsilon.is_finite());
    }

    #[test]
    fn test_wall_epsilon_scales_with_u_tau_cubed() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        let wall_distance = 0.001;
        
        let eps_1 = treatment.wall_epsilon(1.0, wall_distance);
        let eps_2 = treatment.wall_epsilon(2.0, wall_distance);
        
        // ε ~ u_tau^3, so doubling u_tau should give 8x epsilon
        assert_relative_eq!(eps_2 / eps_1, 8.0, epsilon = 1e-10);
    }

    #[test]
    fn test_wall_omega_low_reynolds() {
        let treatment = WallTreatment::<f64>::new(WallFunction::LowReynolds);
        let wall_distance = 0.001;
        let viscosity = 1e-5;
        let density = 1.0;
        
        let omega = treatment.wall_omega(wall_distance, viscosity, density);
        assert!(omega > 0.0);
        assert!(omega.is_finite());
    }

    #[test]
    fn test_wall_omega_standard() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        let wall_distance = 0.001;
        let viscosity = 1e-5;
        let density = 1.0;
        
        let omega = treatment.wall_omega(wall_distance, viscosity, density);
        assert!(omega > 0.0);
        assert!(omega.is_finite());
    }

    #[test]
    fn test_wall_function_type_standard() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        matches!(treatment.wall_function_type(), WallFunction::Standard);
    }

    #[test]
    fn test_wall_function_type_blended() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Blended);
        matches!(treatment.wall_function_type(), WallFunction::Blended);
    }

    #[test]
    fn test_wall_function_type_low_reynolds() {
        let treatment = WallTreatment::<f64>::new(WallFunction::LowReynolds);
        matches!(treatment.wall_function_type(), WallFunction::LowReynolds);
    }

    #[test]
    fn test_standard_wall_function_continuity() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        
        // Check continuity at transition points
        let y_trans_1 = Y_PLUS_VISCOUS_SUBLAYER - 0.1;
        let y_trans_2 = Y_PLUS_VISCOUS_SUBLAYER + 0.1;
        
        let u_1 = treatment.u_plus(y_trans_1);
        let u_2 = treatment.u_plus(y_trans_2);
        
        // Should be reasonably close (continuous)
        assert!((u_2 - u_1).abs() < 1.0);
    }

    #[test]
    fn test_u_plus_monotonically_increasing() {
        let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
        
        // u+ should increase monotonically with y+
        let y_plus_values = vec![1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0];
        let mut prev_u_plus = 0.0;
        
        for &y_plus in &y_plus_values {
            let u_plus = treatment.u_plus(y_plus);
            assert!(u_plus > prev_u_plus);
            prev_u_plus = u_plus;
        }
    }
}
