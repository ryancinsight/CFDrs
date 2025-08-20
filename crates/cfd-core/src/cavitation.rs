//! Cavitation physics models for CFD simulations
//!
//! This module implements hydrodynamic cavitation models based on established
//! literature including Brennen (1995), Franc & Michel (2004), and Rayleigh-Plesset
//! dynamics for bubble growth and collapse.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

/// Physical constants for cavitation
pub mod constants {
    /// Surface tension of water at 20°C (N/m)
    pub const SURFACE_TENSION_WATER: f64 = 0.0728;
    
    /// Vapor pressure of water at 20°C (Pa)
    pub const VAPOR_PRESSURE_WATER_20C: f64 = 2339.0;
    
    /// Saturation pressure ratio threshold for inception
    pub const CAVITATION_INCEPTION_THRESHOLD: f64 = 0.9;
    
    /// Blake critical radius coefficient
    pub const BLAKE_CRITICAL_COEFFICIENT: f64 = 0.85;
    
    /// Minimum bubble radius for numerical stability (m)
    pub const MIN_BUBBLE_RADIUS: f64 = 1e-9;
    
    /// Maximum void fraction before numerical issues
    pub const MAX_VOID_FRACTION: f64 = 0.999;
    
    /// Nucleation site density in clean water (#/m³)
    pub const NUCLEATION_DENSITY_CLEAN: f64 = 1e6;
    
    /// Nucleation site density in technical water (#/m³)
    pub const NUCLEATION_DENSITY_TECHNICAL: f64 = 1e8;
}

/// Cavitation number (dimensionless)
/// σ = (p - p_v) / (0.5 * ρ * v²)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CavitationNumber<T: RealField> {
    /// Reference pressure (Pa)
    pub reference_pressure: T,
    /// Vapor pressure (Pa)
    pub vapor_pressure: T,
    /// Reference density (kg/m³)
    pub density: T,
    /// Reference velocity (m/s)
    pub velocity: T,
}

impl<T: RealField + FromPrimitive + Copy> CavitationNumber<T> {
    /// Calculate cavitation number
    pub fn calculate(&self) -> T {
        let dynamic_pressure = T::from_f64(0.5).unwrap_or_else(|| T::one()) * 
            self.density * self.velocity * self.velocity;
        
        if dynamic_pressure > T::from_f64(1e-10).unwrap_or_else(|| T::one()) {
            (self.reference_pressure - self.vapor_pressure) / dynamic_pressure
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one()) // Large number for zero velocity
        }
    }
    
    /// Check if cavitation is likely (σ < threshold)
    pub fn is_cavitating(&self, threshold: T) -> bool {
        self.calculate() < threshold
    }
}

/// Rayleigh-Plesset equation for bubble dynamics
/// R*R̈ + (3/2)*Ṙ² = (p_B - p_∞)/ρ
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RayleighPlesset<T: RealField> {
    /// Bubble radius (m)
    pub radius: T,
    /// Bubble wall velocity (m/s)
    pub velocity: T,
    /// Liquid density (kg/m³)
    pub liquid_density: T,
    /// Surface tension (N/m)
    pub surface_tension: T,
    /// Liquid viscosity (Pa·s)
    pub viscosity: T,
    /// Vapor pressure (Pa)
    pub vapor_pressure: T,
    /// Initial radius (m)
    pub initial_radius: T,
}

impl<T: RealField + FromPrimitive + Copy> RayleighPlesset<T> {
    /// Create new Rayleigh-Plesset model
    pub fn new(
        initial_radius: T,
        liquid_density: T,
        surface_tension: T,
        viscosity: T,
        vapor_pressure: T,
    ) -> Self {
        Self {
            radius: initial_radius,
            velocity: T::zero(),
            liquid_density,
            surface_tension,
            viscosity,
            vapor_pressure,
            initial_radius,
        }
    }
    
    /// Calculate bubble pressure including surface tension
    pub fn bubble_pressure(&self, _ambient_pressure: T) -> T {
        // p_B = p_v - 2σ/R - 4μṘ/R
        let surface_term = T::from_f64(2.0).unwrap_or_else(|| T::one()) * self.surface_tension / self.radius;
        let viscous_term = T::from_f64(4.0).unwrap_or_else(|| T::one()) * self.viscosity * 
            self.velocity / self.radius;
        
        self.vapor_pressure - surface_term - viscous_term
    }
    
    /// Calculate bubble acceleration (R̈)
    pub fn acceleration(&self, ambient_pressure: T) -> T {
        let p_bubble = self.bubble_pressure(ambient_pressure);
        let pressure_diff = p_bubble - ambient_pressure;
        
        // R̈ = (1/R) * [(p_B - p_∞)/ρ - (3/2)Ṙ²]
        let pressure_term = pressure_diff / self.liquid_density;
        let kinetic_term = T::from_f64(1.5).unwrap_or_else(|| T::one()) * self.velocity * self.velocity;
        
        (pressure_term - kinetic_term) / self.radius
    }
    
    /// Time step the bubble dynamics
    pub fn step(&mut self, dt: T, ambient_pressure: T) {
        let accel = self.acceleration(ambient_pressure);
        
        // Update velocity and radius using explicit Euler
        self.velocity = self.velocity + accel * dt;
        self.radius = self.radius + self.velocity * dt;
        
        // Ensure minimum radius for stability
        let min_radius = T::from_f64(constants::MIN_BUBBLE_RADIUS).unwrap_or_else(|| T::one());
        if self.radius < min_radius {
            self.radius = min_radius;
            self.velocity = T::zero();
        }
    }
    
    /// Calculate collapse time (Rayleigh collapse time)
    pub fn rayleigh_collapse_time(&self, pressure_difference: T) -> T {
        // t_c = 0.915 * R_0 * sqrt(ρ / Δp)
        T::from_f64(0.915).unwrap_or_else(|| T::one()) * self.initial_radius * 
            (self.liquid_density / pressure_difference).sqrt()
    }
}

/// Cavitation model types
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CavitationModel<T: RealField> {
    /// Kunz model (2000)
    Kunz {
        /// Vaporization coefficient
        vaporization_coeff: T,
        /// Condensation coefficient
        condensation_coeff: T,
    },
    /// Schnerr-Sauer model (2001)
    SchnerrSauer {
        /// Bubble number density (#/m³)
        bubble_density: T,
        /// Initial bubble radius (m)
        initial_radius: T,
    },
    /// Zwart-Gerber-Belamri model (2004)
    ZGB {
        /// Nucleation site volume fraction
        nucleation_fraction: T,
        /// Bubble radius (m)
        bubble_radius: T,
        /// Vaporization coefficient
        f_vap: T,
        /// Condensation coefficient
        f_cond: T,
    },
}

impl<T: RealField + FromPrimitive + Copy> CavitationModel<T> {
    /// Calculate mass transfer rate (kg/m³/s)
    pub fn mass_transfer_rate(
        &self,
        pressure: T,
        vapor_pressure: T,
        void_fraction: T,
        density_liquid: T,
        density_vapor: T,
    ) -> T {
        match self {
            CavitationModel::Kunz { vaporization_coeff, condensation_coeff } => {
                let pressure_diff = pressure - vapor_pressure;
                
                if pressure_diff < T::zero() {
                    // Vaporization
                    let rate = *vaporization_coeff * density_vapor * 
                        (T::one() - void_fraction) * pressure_diff.abs() / 
                        (T::from_f64(0.5).unwrap_or_else(|| T::one()) * density_liquid);
                    rate
                } else {
                    // Condensation
                    let rate = *condensation_coeff * density_vapor * 
                        void_fraction * pressure_diff / 
                        (T::from_f64(0.5).unwrap_or_else(|| T::one()) * density_liquid);
                    -rate
                }
            },
            
            CavitationModel::SchnerrSauer { bubble_density, initial_radius: _ } => {
                // Calculate bubble radius from void fraction
                let n_b = bubble_density;
                let alpha = void_fraction;
                
                // R_B = [(3α)/(4πn_B(1-α))]^(1/3)
                let denominator = T::from_f64(4.0 * std::f64::consts::PI).unwrap_or_else(|| T::one()) * 
                    n_b * (T::one() - alpha);
                
                if denominator > T::from_f64(1e-10).unwrap_or_else(|| T::one()) {
                    let radius_cubed = T::from_f64(3.0).unwrap_or_else(|| T::one()) * alpha / denominator;
                    let radius = radius_cubed.powf(T::from_f64(1.0/3.0).unwrap_or_else(|| T::one()));
                    
                    // Mass transfer rate
                    let pressure_diff = pressure - vapor_pressure;
                    let sign = if pressure_diff < T::zero() { T::one() } else { -T::one() };
                    
                    let rate = sign * T::from_f64(3.0).unwrap_or_else(|| T::one()) * alpha * (T::one() - alpha) * 
                        density_vapor * (T::from_f64(2.0/3.0).unwrap_or_else(|| T::one()) * 
                        pressure_diff.abs() / density_liquid).sqrt() / radius;
                    
                    rate
                } else {
                    T::zero()
                }
            },
            
            CavitationModel::ZGB { nucleation_fraction, bubble_radius, f_vap, f_cond } => {
                let pressure_diff = pressure - vapor_pressure;
                let r_b = *bubble_radius;
                
                if pressure_diff < T::zero() {
                    // Vaporization
                    let rate = *f_vap * T::from_f64(3.0).unwrap_or_else(|| T::one()) * 
                        *nucleation_fraction * (T::one() - void_fraction) * 
                        density_vapor * (T::from_f64(2.0/3.0).unwrap_or_else(|| T::one()) * 
                        pressure_diff.abs() / density_liquid).sqrt() / r_b;
                    rate
                } else {
                    // Condensation
                    let rate = *f_cond * T::from_f64(3.0).unwrap_or_else(|| T::one()) * 
                        void_fraction * density_vapor * 
                        (T::from_f64(2.0/3.0).unwrap_or_else(|| T::one()) * pressure_diff / 
                        density_liquid).sqrt() / r_b;
                    -rate
                }
            }
        }
    }
}

/// Venturi cavitation parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiCavitation<T: RealField> {
    /// Inlet diameter (m)
    pub inlet_diameter: T,
    /// Throat diameter (m)
    pub throat_diameter: T,
    /// Outlet diameter (m)
    pub outlet_diameter: T,
    /// Convergent angle (radians)
    pub convergent_angle: T,
    /// Divergent angle (radians)
    pub divergent_angle: T,
    /// Inlet pressure (Pa)
    pub inlet_pressure: T,
    /// Outlet pressure (Pa)
    pub outlet_pressure: T,
    /// Fluid properties
    pub fluid_density: T,
    /// Fluid dynamic viscosity (Pa·s)
    pub fluid_viscosity: T,
    /// Vapor pressure of fluid (Pa)
    pub vapor_pressure: T,
    /// Surface tension coefficient (N/m)
    pub surface_tension: T,
}

impl<T: RealField + FromPrimitive + Copy> VenturiCavitation<T> {
    /// Calculate throat velocity using continuity
    pub fn throat_velocity(&self, inlet_velocity: T) -> T {
        let area_ratio = (self.inlet_diameter / self.throat_diameter).powi(2);
        inlet_velocity * area_ratio
    }
    
    /// Calculate pressure at throat using Bernoulli
    pub fn throat_pressure(&self, inlet_velocity: T) -> T {
        let v_throat = self.throat_velocity(inlet_velocity);
        
        // p_throat = p_inlet - 0.5 * ρ * (v_throat² - v_inlet²)
        let dynamic_diff = T::from_f64(0.5).unwrap_or_else(|| T::one()) * self.fluid_density * 
            (v_throat * v_throat - inlet_velocity * inlet_velocity);
        
        self.inlet_pressure - dynamic_diff
    }
    
    /// Calculate cavitation number at throat
    pub fn throat_cavitation_number(&self, inlet_velocity: T) -> T {
        let p_throat = self.throat_pressure(inlet_velocity);
        let v_throat = self.throat_velocity(inlet_velocity);
        
        let cav = CavitationNumber {
            reference_pressure: p_throat,
            vapor_pressure: self.vapor_pressure,
            density: self.fluid_density,
            velocity: v_throat,
        };
        
        cav.calculate()
    }
    
    /// Calculate pressure recovery coefficient
    pub fn pressure_recovery_coefficient(&self) -> T {
        // C_pr = 1 - (A_throat/A_outlet)²
        let area_ratio = (self.throat_diameter / self.outlet_diameter).powi(2);
        T::one() - area_ratio
    }
    
    /// Estimate cavitation inception point
    pub fn inception_velocity(&self) -> T {
        // Find velocity where σ = σ_inception
        let sigma_inception = T::from_f64(constants::CAVITATION_INCEPTION_THRESHOLD).unwrap_or_else(|| T::one());
        
        // σ = (p_inlet - p_v - 0.5*ρ*(v_throat² - v_inlet²)) / (0.5*ρ*v_throat²)
        // Solving for v_inlet when σ = σ_inception
        let area_ratio = (self.inlet_diameter / self.throat_diameter).powi(2);
        let ar2 = area_ratio * area_ratio;
        
        let numerator = T::from_f64(2.0).unwrap_or_else(|| T::one()) * (self.inlet_pressure - self.vapor_pressure);
        let denominator = self.fluid_density * (ar2 * (T::one() + sigma_inception) - T::one());
        
        if denominator > T::zero() {
            (numerator / denominator).sqrt()
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one()) // No cavitation possible
        }
    }
}

/// Cavitation damage model (erosion potential)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CavitationDamage<T: RealField> {
    /// Material properties
    pub material_hardness: T,  // Pa
    /// Material resilience against cavitation damage (J/m³)
    pub material_resilience: T,
    /// Cavitation intensity
    pub collapse_pressure: T,   // Pa
    /// Frequency of bubble collapses (Hz)
    pub collapse_frequency: T,
    /// Area affected by cavitation damage (m²)
    pub affected_area: T,
}

impl<T: RealField + FromPrimitive + Copy> CavitationDamage<T> {
    /// Calculate erosion rate (Hammitt model)
    pub fn erosion_rate(&self) -> T {
        // E = K * (p_collapse - p_threshold)^n * f * A
        let threshold = self.material_hardness * T::from_f64(0.1).unwrap_or_else(|| T::one());
        
        if self.collapse_pressure > threshold {
            let intensity = (self.collapse_pressure - threshold).powf(T::from_f64(2.0).unwrap_or_else(|| T::one()));
            let rate = T::from_f64(1e-12).unwrap_or_else(|| T::one()) * intensity * 
                self.collapse_frequency * self.affected_area / 
                self.material_resilience;
            rate
        } else {
            T::zero()
        }
    }
    
    /// Calculate cavitation intensity parameter
    pub fn intensity_parameter(&self) -> T {
        self.collapse_pressure * self.collapse_frequency.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_cavitation_number() {
        let cav = CavitationNumber::<f64> {
            reference_pressure: 101325.0,
            vapor_pressure: 2339.0,
            density: 998.0,
            velocity: 10.0,
        };
        
        let sigma = cav.calculate();
        assert_relative_eq!(sigma, 1.98, epsilon = 0.01);
        assert!(!cav.is_cavitating(1.0));
    }
    
    #[test]
    fn test_rayleigh_plesset() {
        let mut bubble = RayleighPlesset::<f64>::new(
            1e-6,  // 1 micron initial radius
            998.0, // water density
            0.0728, // surface tension
            0.001, // viscosity
            2339.0, // vapor pressure
        );
        
        let ambient = 101325.0;
        let accel = bubble.acceleration(ambient);
        assert!(accel < 0.0); // Should collapse under high pressure
        
        bubble.step(1e-9, ambient);
        assert!(bubble.radius < 1e-6); // Should shrink
    }
    
    #[test]
    fn test_venturi_cavitation() {
        let venturi = VenturiCavitation::<f64> {
            inlet_diameter: 0.1,
            throat_diameter: 0.05,
            outlet_diameter: 0.08,
            convergent_angle: 0.35, // 20 degrees
            divergent_angle: 0.14,  // 8 degrees
            inlet_pressure: 200000.0,
            outlet_pressure: 101325.0,
            fluid_density: 998.0,
            fluid_viscosity: 0.001,
            vapor_pressure: 2339.0,
            surface_tension: 0.0728,
        };
        
        let v_inlet = 5.0;
        let v_throat = venturi.throat_velocity(v_inlet);
        assert_relative_eq!(v_throat, 20.0, epsilon = 0.01);
        
        let p_throat = venturi.throat_pressure(v_inlet);
        assert!(p_throat < venturi.inlet_pressure);
        
        let sigma = venturi.throat_cavitation_number(v_inlet);
        assert!(sigma > 0.0);
    }
}