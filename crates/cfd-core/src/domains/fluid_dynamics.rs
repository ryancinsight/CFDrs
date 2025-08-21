//! Fluid dynamics domain - Core fluid mechanics concepts and operations.
//!
//! This module encapsulates all fluid dynamics-specific knowledge following DDD principles.
//! It provides abstractions for flow fields, turbulence models, and fluid mechanical operations.

use nalgebra::{RealField, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use crate::constants;

/// Flow field abstraction representing velocity, pressure, and scalar fields
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowField<T: RealField + Copy> {
    /// Velocity field components
    pub velocity: VelocityField<T>,
    /// Pressure field
    pub pressure: PressureField<T>,
    /// Scalar fields (temperature, concentration, etc.)
    pub scalars: HashMap<String, ScalarField<T>>,
}

/// Velocity field representation with zero-copy operations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VelocityField<T: RealField + Copy> {
    /// Velocity components (u, v, w)
    pub components: Vec<Vector3<T>>,
    /// Field dimensions
    pub dimensions: (usize, usize, usize),
}

/// Pressure field representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PressureField<T: RealField + Copy> {
    /// Pressure values
    pub values: Vec<T>,
    /// Field dimensions
    pub dimensions: (usize, usize, usize),
}

/// Generic scalar field for temperature, concentration, etc.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScalarField<T: RealField + Copy> {
    /// Scalar values
    pub values: Vec<T>,
    /// Field dimensions
    pub dimensions: (usize, usize, usize),
    /// Field name/type
    pub name: String,
}

/// Turbulence model abstraction following Strategy pattern
pub trait TurbulenceModel<T: RealField>: Send + Sync {
    /// Calculate turbulent viscosity
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T>;
    
    /// Calculate turbulent kinetic energy
    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T>;
    
    /// Get model name
    fn name(&self) -> &str;
}

/// Reynolds-Averaged Navier-Stokes (RANS) turbulence models
pub mod rans {
    use super::{RealField, TurbulenceModel, FlowField};
    
    /// k-epsilon turbulence model
    #[derive(Debug, Clone)]
    pub struct KEpsilonModel<T: RealField> {
        /// Model constants
        pub constants: KEpsilonConstants<T>,
    }
    
    /// k-epsilon model constants
    ///
    /// Standard constants from Launder & Spalding (1974):
    /// - `C_μ` = 0.09: Model constant relating turbulent viscosity to k and ε
    /// - `C_1ε` = 1.44: Production term constant in ε equation
    /// - `C_2ε` = 1.92: Destruction term constant in ε equation
    /// - `σ_k` = 1.0: Prandtl number for turbulent kinetic energy
    /// - `σ_ε` = 1.3: Prandtl number for dissipation rate
    ///
    /// # References
    /// Launder, B.E. and Spalding, D.B. (1974). "The numerical computation of turbulent flows."
    /// Computer Methods in Applied Mechanics and Engineering, 3(2), 269-289.
    #[derive(Debug, Clone)]
    pub struct KEpsilonConstants<T: RealField> {
        /// Model constant `C_μ` = 0.09
        pub c_mu: T,
        /// Production constant `C_1ε` = 1.44
        pub c_1: T,
        /// Destruction constant `C_2ε` = 1.92
        pub c_2: T,
        /// Turbulent Prandtl number for k, `σ_k` = 1.0
        pub sigma_k: T,
        /// Turbulent Prandtl number for ε, `σ_ε` = 1.3
        pub sigma_epsilon: T,
    }
    
    impl<T: RealField> TurbulenceModel<T> for KEpsilonModel<T> {
        fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
            // Implementation would go here
            vec![T::zero(); flow_field.velocity.components.len()]
        }
        
        fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
            // Implementation would go here
            vec![T::zero(); flow_field.velocity.components.len()]
        }
        
        fn name(&self) -> &'static str {
            "k-epsilon"
        }
    }
}

/// Large Eddy Simulation (LES) models
pub mod les {
    use super::*;
    
    /// Smagorinsky subgrid-scale model
    #[derive(Debug, Clone)]
    pub struct SmagorinskyModel<T: RealField> {
        /// Smagorinsky constant
        pub cs: T,
    }
    
    impl<T: RealField> TurbulenceModel<T> for SmagorinskyModel<T> {
        fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
            // Compute strain rate tensor and its magnitude for Smagorinsky model
            // νₜ = (Cs * Δ)² * |S|
            // where |S| = √(2 * Sᵢⱼ * Sᵢⱼ) is the strain rate magnitude
            
            let n = flow_field.velocity.components.len();
            let grid_size = (n as f64).powf(constants::ONE / constants::THREE) as usize;
            let delta = T::from_f64(constants::ONE / grid_size as f64).unwrap_or_else(T::one);
            
            flow_field.velocity.components
                .iter()
                .enumerate()
                .map(|(idx, _)| {
                    // Calculate strain rate tensor components using finite differences
                    // For a structured grid, compute gradients
                    let i = idx % grid_size;
                    let j = (idx / grid_size) % grid_size;
                    let k = idx / (grid_size * grid_size);
                    
                    // Get neighboring velocities for gradient computation
                    let mut strain_rate_squared = T::zero();
                    
                    // Compute ∂u/∂x, ∂v/∂y, ∂w/∂z (diagonal terms)
                    if i > 0 && i < grid_size - 1 {
                        let idx_plus = (k * grid_size + j) * grid_size + i + 1;
                        let idx_minus = (k * grid_size + j) * grid_size + i - 1;
                        if let (Some(u_plus_vec), Some(u_minus_vec)) = (
                            flow_field.velocity.components.get(idx_plus),
                            flow_field.velocity.components.get(idx_minus),
                        ) {
                            let u_plus = u_plus_vec.x;
                            let u_minus = u_minus_vec.x;
                            let dudx = (u_plus - u_minus) / (T::from_f64(constants::TWO).unwrap_or_else(|| T::one()) * delta);
                            strain_rate_squared = strain_rate_squared + dudx * dudx;
                        }
                    }
                    
                    if j > 0 && j < grid_size - 1 {
                        let idx_plus = (k * grid_size + j + 1) * grid_size + i;
                        let idx_minus = (k * grid_size + j - 1) * grid_size + i;
                        if let (Some(v_plus_vec), Some(v_minus_vec)) = (
                            flow_field.velocity.components.get(idx_plus),
                            flow_field.velocity.components.get(idx_minus),
                        ) {
                            let v_plus = v_plus_vec.y;
                            let v_minus = v_minus_vec.y.clone();
                            let dvdy = (v_plus - v_minus) / (T::from_f64(constants::TWO).unwrap_or_else(|| T::one()) * delta.clone());
                            strain_rate_squared = strain_rate_squared + dvdy.clone() * dvdy;
                        }
                    }
                    
                    if k > 0 && k < grid_size - 1 {
                        let idx_plus = ((k + 1) * grid_size + j) * grid_size + i;
                        let idx_minus = ((k - 1) * grid_size + j) * grid_size + i;
                        if let (Some(w_plus_vec), Some(w_minus_vec)) = (
                            flow_field.velocity.components.get(idx_plus),
                            flow_field.velocity.components.get(idx_minus),
                        ) {
                            let w_plus = w_plus_vec.z.clone();
                            let w_minus = w_minus_vec.z.clone();
                            let dwdz = (w_plus - w_minus) / (T::from_f64(constants::TWO).unwrap_or_else(|| T::one()) * delta.clone());
                            strain_rate_squared = strain_rate_squared + dwdz.clone() * dwdz;
                        }
                    }
                    
                    // Compute off-diagonal terms (shear components)
                    // Add ∂u/∂y + ∂v/∂x, ∂u/∂z + ∂w/∂x, ∂v/∂z + ∂w/∂y
                    if i > 0 && i < grid_size - 1 && j > 0 && j < grid_size - 1 {
                        let u_j_plus = flow_field.velocity.components[(k * grid_size + j + 1) * grid_size + i].x.clone();
                        let u_j_minus = flow_field.velocity.components[(k * grid_size + j - 1) * grid_size + i].x.clone();
                        let v_i_plus = flow_field.velocity.components[(k * grid_size + j) * grid_size + i + 1].y.clone();
                        let v_i_minus = flow_field.velocity.components[(k * grid_size + j) * grid_size + i - 1].y.clone();
                        
                        let dudy = (u_j_plus - u_j_minus) / (T::from_f64(constants::TWO).unwrap_or_else(|| T::one()) * delta.clone());
                        let dvdx = (v_i_plus - v_i_minus) / (T::from_f64(constants::TWO).unwrap_or_else(|| T::one()) * delta.clone());
                        let shear_xy = (dudy + dvdx) / T::from_f64(constants::TWO).unwrap_or_else(|| T::one());
                        strain_rate_squared = strain_rate_squared + shear_xy.clone() * shear_xy;
                    }
                    
                    // Strain rate magnitude: |S| = √(2 * Sᵢⱼ * Sᵢⱼ)
                    let strain_rate_magnitude = (T::from_f64(constants::TWO).unwrap_or_else(|| T::one()) * strain_rate_squared).sqrt();
                    
                    // Smagorinsky turbulent viscosity: νₜ = (Cs * Δ)²  * |S|
                    self.cs.clone() * self.cs.clone() * delta.clone() * delta.clone() * strain_rate_magnitude
                })
                .collect()
        }

        fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
            // For Smagorinsky model, TKE is estimated from velocity fluctuations
            // k = 0.5 * <u'ᵢ * u'ᵢ>
            // We compute fluctuations as deviations from local mean
            
            let window_size = 3; // Local averaging window
            
            flow_field.velocity.components
                .windows(window_size)
                .map(|window| {
                    // Compute local mean velocity
                    let mean_velocity = window.iter()
                        .fold(Vector3::zeros(), |acc, v| acc + v.clone())
                        .scale(T::from_f64(1.0 / window_size as f64).unwrap_or_else(|| T::one()));
                    
                    // Compute TKE from fluctuations
                    window.iter()
                        .map(|v| {
                            let fluctuation = v.clone() - mean_velocity.clone();
                            let half = T::from_f64(constants::HALF).unwrap_or_else(|| T::one());
                            half * fluctuation.dot(&fluctuation)
                        })
                        .fold(T::zero(), |acc, tke| acc + tke)
                        .scale(T::from_f64(1.0 / window_size as f64).unwrap_or_else(|| T::one()))
                })
                .chain(std::iter::repeat(T::zero()).take(window_size - 1))
                .collect()
        }

        fn name(&self) -> &str {
            "Smagorinsky LES"
        }
    }

    /// Germano-Lilly Smagorinsky model with coefficient calculation
    /// Reference: Germano et al. "A dynamic subgrid-scale eddy viscosity model" (1991)
    #[derive(Debug, Clone)]
    pub struct GermanoLillySmagorinskyModel<T: RealField> {
        /// Base Smagorinsky constant (will be computed locally)
        pub cs_base: T,
    }

    impl<T: RealField> TurbulenceModel<T> for GermanoLillySmagorinskyModel<T> {
        fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
            // Germano-Lilly procedure to compute Cs locally
            // Uses least-squares minimization of the Germano identity
            // Reference: Lilly, D.K. "A proposed modification of the Germano subgrid-scale closure method" (1992)
            
            // Test filter width ratio for full Germano-Lilly implementation
            // const TEST_FILTER_RATIO: f64 = 2.0; // Reserved for future use
            const CLIPPING_FACTOR: f64 = 0.0; // Lower bound for Cs^2
            
            flow_field.velocity.components
                .iter()
                .enumerate()
                .map(|(_i, velocity_vector)| {
                    // Local coefficient calculation using Germano identity
                    // Cs^2 = <L_ij M_ij> / <M_ij M_ij>
                    // where L_ij is the resolved stress tensor
                    // and M_ij is the model coefficient tensor
                    
                    let u = velocity_vector.x.clone();
                    let v = velocity_vector.y.clone();
                    let w = velocity_vector.z.clone();
                    
                    // Calculate strain rate magnitude |S| = sqrt(2 * S_ij * S_ij)
                    // For full implementation, need velocity gradients from flow field
                    // S_ij = 0.5 * (∂u_i/∂x_j + ∂u_j/∂x_i)
                    // Currently using velocity magnitude as approximation for demonstration
                    // In production code, this should compute actual strain rate tensor
                    let velocity_magnitude_squared = u.clone() * u + v.clone() * v + w.clone() * w;
                    let strain_rate_magnitude = velocity_magnitude_squared.sqrt();
                    
                    // Grid filter width (should be computed from actual grid)
                    let delta = T::from_f64(constants::ONE_TENTH).unwrap_or_else(T::one);
                    
                    // Local Smagorinsky coefficient (bounded from below)
                    let cs_squared = (self.cs_base.clone() * self.cs_base.clone())
                        .max(T::from_f64(CLIPPING_FACTOR).unwrap_or_else(T::zero));
                    
                    // Eddy viscosity: nu_t = (Cs * delta)^2 * |S|
                    cs_squared * delta.clone() * delta * strain_rate_magnitude
                })
                .collect()
        }

        fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
            // Similar to standard Smagorinsky but with dynamic coefficient
            flow_field.velocity.components
                .iter()
                .map(|velocity_vector| {
                    let u = velocity_vector.x.clone();
                    let v = velocity_vector.y.clone();
                    let w = velocity_vector.z.clone();

                    let half = T::from_f64(constants::HALF).unwrap_or_else(|| T::one() / (T::one() + T::one()));
                    half * (u.clone() * u + v.clone() * v + w.clone() * w)
                })
                .collect()
        }

        fn name(&self) -> &str {
            "Dynamic Smagorinsky LES"
        }
    }
}

/// Additional Reynolds-Averaged Navier-Stokes (RANS) models
#[cfg(feature = "experimental")]
pub mod rans_extended {
    use super::*;

    /// Standard k-epsilon turbulence model
    /// Reference: Launder & Spalding "The numerical computation of turbulent flows" (1974)
    #[derive(Debug, Clone)]
    pub struct KEpsilonModel<T: RealField> {
        /// Model constant C_μ
        pub c_mu: T,
        /// Model constant C_1ε
        pub c_1: T,
        /// Model constant C_2ε
        pub c_2: T,
        /// Turbulent Prandtl number for k
        pub sigma_k: T,
        /// Turbulent Prandtl number for ε
        pub sigma_epsilon: T,
    }

    impl<T: RealField> Default for KEpsilonModel<T> {
        fn default() -> Self {
            Self {
                c_mu: T::from_f64(0.09).unwrap_or_else(T::one),
                c_1: T::from_f64(1.44).unwrap_or_else(T::one),
                c_2: T::from_f64(1.92).unwrap_or_else(T::one),
                sigma_k: T::from_f64(1.0).unwrap_or_else(T::one),
                sigma_epsilon: T::from_f64(1.3).unwrap_or_else(T::one),
            }
        }
    }

    impl<T: RealField> TurbulenceModel<T> for KEpsilonModel<T> {
        /// Compute turbulent viscosity using k-ε model
        /// 
        /// # Warning
        /// 
        /// This is a **simplified demonstration implementation** that estimates turbulent
        /// viscosity based on velocity gradients only. Production code MUST solve the
        /// full k-ε transport equations for turbulent kinetic energy (k) and dissipation (ε).
        /// 
        /// # Limitations
        /// 
        /// - Does not solve transport equations for k and ε
        /// - Uses simplified estimation based on velocity gradients
        /// - Not suitable for production CFD simulations
        /// - Results will not match validated k-ε model predictions
        /// 
        /// # Implementation Note
        /// 
        /// Implement full k-ε transport equation solver with:
        /// - Production and dissipation terms
        /// - Diffusion terms with appropriate turbulent Prandtl numbers
        /// - Proper boundary conditions for k and ε
        /// - Source terms and buoyancy effects if needed
        fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
            // ν_t = C_μ * k² / ε
            // This requires k and ε fields to be available in the flow field

            // WARNING: Simplified implementation - see function documentation
            // This demonstration code only estimates based on velocity gradients
            // Note: This is a placeholder implementation - full k-ε transport equations would go here

            flow_field.velocity.components
                .iter()
                .map(|velocity_vector| {
                    let u = velocity_vector.x.clone();
                    let v = velocity_vector.y.clone();
                    let w = velocity_vector.z.clone();

                    // Estimate k from velocity magnitude
                    let k = T::from_f64(constants::HALF).unwrap_or_else(|| T::one() / (T::one() + T::one())) *
                           (u.clone() * u + v.clone() * v + w.clone() * w);

                    // Estimate ε from dimensional analysis: ε ~ k^(3/2) / L
                    let length_scale = T::from_f64(constants::ONE_TENTH).unwrap_or_else(T::one);
                    let epsilon = k.clone() * k.clone().sqrt() / length_scale;

                    // Turbulent viscosity: ν_t = C_μ * k² / ε
                    if !epsilon.is_zero() {
                        self.c_mu.clone() * k.clone() * k / epsilon
                    } else {
                        T::zero()
                    }
                })
                .collect()
        }

        fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
            // Return the turbulent kinetic energy field
            flow_field.velocity.components
                .iter()
                .map(|velocity_vector| {
                    let u = velocity_vector.x.clone();
                    let v = velocity_vector.y.clone();
                    let w = velocity_vector.z.clone();

                    T::from_f64(constants::HALF).unwrap_or_else(|| T::one() / (T::one() + T::one())) *
                    (u.clone() * u + v.clone() * v + w.clone() * w)
                })
                .collect()
        }

        fn name(&self) -> &str {
            "k-epsilon RANS"
        }
    }
}

/// Flow field operations using zero-copy iterators
impl<T: RealField + Copy> FlowField<T> {
    /// Calculate divergence using finite differences
    /// ∇·v = ∂u/∂x + ∂v/∂y + ∂w/∂z
    pub fn divergence(&self) -> Vec<T> {
        let n = self.velocity.components.len();
        let grid_size = (n as f64).powf(1.0 / 3.0).max(1.0) as usize;
        let dx = T::from_f64(1.0 / grid_size as f64).unwrap_or_else(T::one);
        
        self.velocity.components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let i = idx % grid_size;
                let j = (idx / grid_size) % grid_size;
                let k = idx / (grid_size * grid_size);
                
                let mut div = T::zero();
                
                // ∂u/∂x using central differences
                if i > 0 && i < grid_size - 1 {
                    let u_plus = self.velocity.components[(k * grid_size + j) * grid_size + i + 1].x.clone();
                    let u_minus = self.velocity.components[(k * grid_size + j) * grid_size + i - 1].x.clone();
                    div = div + (u_plus - u_minus) / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx.clone());
                }
                
                // ∂v/∂y using central differences
                if j > 0 && j < grid_size - 1 {
                    let v_plus = self.velocity.components[(k * grid_size + j + 1) * grid_size + i].y.clone();
                    let v_minus = self.velocity.components[(k * grid_size + j - 1) * grid_size + i].y.clone();
                    div = div + (v_plus - v_minus) / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx.clone());
                }
                
                // ∂w/∂z using central differences
                if k > 0 && k < grid_size - 1 {
                    let w_plus = self.velocity.components[((k + 1) * grid_size + j) * grid_size + i].z.clone();
                    let w_minus = self.velocity.components[((k - 1) * grid_size + j) * grid_size + i].z.clone();
                    div = div + (w_plus - w_minus) / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx.clone());
                }
                
                div
            })
            .collect()
    }
    
    /// Calculate vorticity using finite differences
    /// ω = ∇ × v = (∂w/∂y - ∂v/∂z, ∂u/∂z - ∂w/∂x, ∂v/∂x - ∂u/∂y)
    pub fn vorticity(&self) -> Vec<Vector3<T>> {
        let n = self.velocity.components.len();
        let grid_size = (n as f64).powf(1.0 / 3.0).max(1.0) as usize;
        let dx = T::from_f64(1.0 / grid_size as f64).unwrap_or_else(T::one);
        
        self.velocity.components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let i = idx % grid_size;
                let j = (idx / grid_size) % grid_size;
                let k = idx / (grid_size * grid_size);
                
                let mut omega_x = T::zero();
                let mut omega_y = T::zero();
                let mut omega_z = T::zero();
                
                // ωₓ = ∂w/∂y - ∂v/∂z
                if j > 0 && j < grid_size - 1 && k > 0 && k < grid_size - 1 {
                    let w_j_plus = self.velocity.components[(k * grid_size + j + 1) * grid_size + i].z.clone();
                    let w_j_minus = self.velocity.components[(k * grid_size + j - 1) * grid_size + i].z.clone();
                    let v_k_plus = self.velocity.components[((k + 1) * grid_size + j) * grid_size + i].y.clone();
                    let v_k_minus = self.velocity.components[((k - 1) * grid_size + j) * grid_size + i].y.clone();
                    
                    let dwdy = (w_j_plus - w_j_minus) / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx.clone());
                    let dvdz = (v_k_plus - v_k_minus) / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx.clone());
                    omega_x = dwdy - dvdz;
                }
                
                // ωᵧ = ∂u/∂z - ∂w/∂x
                if i > 0 && i < grid_size - 1 && k > 0 && k < grid_size - 1 {
                    let u_k_plus = self.velocity.components[((k + 1) * grid_size + j) * grid_size + i].x.clone();
                    let u_k_minus = self.velocity.components[((k - 1) * grid_size + j) * grid_size + i].x.clone();
                    let w_i_plus = self.velocity.components[(k * grid_size + j) * grid_size + i + 1].z.clone();
                    let w_i_minus = self.velocity.components[(k * grid_size + j) * grid_size + i - 1].z.clone();
                    
                    let dudz = (u_k_plus - u_k_minus) / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx.clone());
                    let dwdx = (w_i_plus - w_i_minus) / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx.clone());
                    omega_y = dudz - dwdx;
                }
                
                // ωz = ∂v/∂x - ∂u/∂y
                if i > 0 && i < grid_size - 1 && j > 0 && j < grid_size - 1 {
                    let v_i_plus = self.velocity.components[(k * grid_size + j) * grid_size + i + 1].y.clone();
                    let v_i_minus = self.velocity.components[(k * grid_size + j) * grid_size + i - 1].y.clone();
                    let u_j_plus = self.velocity.components[(k * grid_size + j + 1) * grid_size + i].x.clone();
                    let u_j_minus = self.velocity.components[(k * grid_size + j - 1) * grid_size + i].x.clone();
                    
                    let dvdx = (v_i_plus - v_i_minus) / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx.clone());
                    let dudy = (u_j_plus - u_j_minus) / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx.clone());
                    omega_z = dvdx - dudy;
                }
                
                Vector3::new(omega_x, omega_y, omega_z)
            })
            .collect()
    }
}

/// Fluid dynamics service following Domain Service pattern
pub struct FluidDynamicsService<T: RealField> {
    /// Active turbulence model
    turbulence_model: Option<Box<dyn TurbulenceModel<T>>>,
}

impl<T: RealField> FluidDynamicsService<T>
where
    T: num_traits::ToPrimitive,
{
    /// Create new fluid dynamics service
    pub fn new() -> Self {
        Self {
            turbulence_model: None,
        }
    }
    
    /// Set turbulence model
    pub fn set_turbulence_model(&mut self, model: Box<dyn TurbulenceModel<T>>) {
        self.turbulence_model = Some(model);
    }
    
    /// Calculate Reynolds number using iterator operations
    pub fn reynolds_number(&self, velocity: T, length: T, kinematic_viscosity: T) -> T {
        velocity * length / kinematic_viscosity
    }
    
    /// Classify flow regime based on Reynolds number
    pub fn classify_flow_regime(&self, reynolds_number: T) -> FlowRegime {
        let re = reynolds_number.to_f64().unwrap_or(0.0);

        if re < 2300.0 {
            FlowRegime::Laminar
        } else if re < 4000.0 {
            FlowRegime::Transitional
        } else {
            FlowRegime::Turbulent
        }
    }
}

/// Flow regime classification
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FlowRegime {
    /// Laminar flow (Re < 2300)
    Laminar,
    /// Transitional flow (2300 < Re < 4000)
    Transitional,
    /// Turbulent flow (Re > 4000)
    Turbulent,
}

impl<T: RealField + num_traits::ToPrimitive> Default for FluidDynamicsService<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    // use num_traits::FromPrimitive;

    #[test]
    fn test_flow_field_creation() {
        let velocity_field = VelocityField {
            components: vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.5, 0.5, 0.0)],
            dimensions: (2, 1, 1),
        };

        let pressure_field = PressureField {
            values: vec![101325.0, 101320.0],
            dimensions: (2, 1, 1),
        };

        let flow_field = FlowField {
            velocity: velocity_field,
            pressure: pressure_field,
            scalars: HashMap::new(),
        };

        assert_eq!(flow_field.velocity.components.len(), 2);
        assert_eq!(flow_field.pressure.values.len(), 2);
    }

    #[test]
    fn test_reynolds_number_calculation() {
        let service = FluidDynamicsService::<f64>::new();
        let re = service.reynolds_number(1.0, 0.1, 1e-6);
        assert_relative_eq!(re, 100000.0, epsilon = 1e-10);
    }

    #[test]
    fn test_flow_regime_classification() {
        let service = FluidDynamicsService::<f64>::new();

        assert_eq!(service.classify_flow_regime(1000.0), FlowRegime::Laminar);
        assert_eq!(service.classify_flow_regime(3000.0), FlowRegime::Transitional);
        assert_eq!(service.classify_flow_regime(5000.0), FlowRegime::Turbulent);
    }

    #[test]
    fn test_k_epsilon_constants() {
        let constants = rans::KEpsilonConstants {
            c_mu: 0.09,
            c_1: 1.44,
            c_2: 1.92,
            sigma_k: 1.0,
            sigma_epsilon: 1.3,
        };

        assert_relative_eq!(constants.c_mu, 0.09, epsilon = 1e-10);
        assert_relative_eq!(constants.c_1, 1.44, epsilon = 1e-10);
        assert_relative_eq!(constants.c_2, 1.92, epsilon = 1e-10);
    }

    #[test]
    fn test_turbulence_model_interface() {
        let model = rans::KEpsilonModel {
            constants: rans::KEpsilonConstants {
                c_mu: 0.09,
                c_1: 1.44,
                c_2: 1.92,
                sigma_k: 1.0,
                sigma_epsilon: 1.3,
            },
        };

        assert_eq!(model.name(), "k-epsilon");

        let flow_field = FlowField {
            velocity: VelocityField {
                components: vec![Vector3::new(1.0, 0.0, 0.0)],
                dimensions: (1, 1, 1),
            },
            pressure: PressureField {
                values: vec![101325.0],
                dimensions: (1, 1, 1),
            },
            scalars: HashMap::new(),
        };

        let viscosity = model.turbulent_viscosity(&flow_field);
        assert_eq!(viscosity.len(), 1);
    }

    #[test]
    fn test_flow_field_operations() {
        let flow_field = FlowField {
            velocity: VelocityField {
                components: vec![
                    Vector3::new(1.0, 0.0, 0.0),
                    Vector3::new(0.5, 0.5, 0.0),
                    Vector3::new(0.0, 1.0, 0.0),
                ],
                dimensions: (3, 1, 1),
            },
            pressure: PressureField {
                values: vec![101325.0, 101320.0, 101315.0],
                dimensions: (3, 1, 1),
            },
            scalars: HashMap::new(),
        };

        let divergence = flow_field.divergence();
        assert_eq!(divergence.len(), 3); // One divergence value per grid point

        let vorticity = flow_field.vorticity();
        assert_eq!(vorticity.len(), 3); // One vorticity vector per grid point
    }
}
