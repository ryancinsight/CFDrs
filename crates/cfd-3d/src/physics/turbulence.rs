//! Turbulence modeling abstractions and implementations
//!
//! Provides trait-based turbulence model implementations including
//! LES (Smagorinsky) and RANS (k-epsilon, Mixing Length) models.
//!
//! Moved from cfd-core to cfd-3d as these models are specific to 3D flow fields.

use cfd_core::physics::fluid_dynamics::fields::{FlowField, VelocityField};
use cfd_core::physics::fluid_dynamics::{rans::RANSModel, turbulence::TurbulenceModel};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Smagorinsky subgrid-scale model for LES
#[derive(Debug, Clone)]
pub struct SmagorinskyModel<T: RealField + Copy> {
    /// Smagorinsky constant (typically 0.1-0.2)
    pub cs: T,
    /// Base Smagorinsky constant for dynamic model
    pub cs_base: T,
}

impl<T: RealField + Copy + FromPrimitive> SmagorinskyModel<T> {
    /// Create a new Smagorinsky model with standard constant
    pub fn new(cs: T) -> Self {
        Self { cs, cs_base: cs }
    }

    /// Calculate strain rate magnitude at a grid point
    fn calculate_strain_rate_at_point(
        &self,
        flow_field: &FlowField<T>,
        i: usize,
        j: usize,
        k: usize,
        delta: T,
    ) -> T {
        let (nx, ny, nz) = flow_field.velocity.dimensions;

        // Calculate strain rate tensor components
        let mut s11 = T::zero();
        let mut s22 = T::zero();
        let mut s33 = T::zero();

        // Central differences for velocity gradients
        if i > 0 && i < nx - 1 {
            if let (Some(u_plus), Some(u_minus)) = (
                flow_field.velocity.get(i + 1, j, k),
                flow_field.velocity.get(i - 1, j, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                s11 = (u_plus.x - u_minus.x) / (two * delta);
            }
        }

        if j > 0 && j < ny - 1 {
            if let (Some(v_plus), Some(v_minus)) = (
                flow_field.velocity.get(i, j + 1, k),
                flow_field.velocity.get(i, j - 1, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                s22 = (v_plus.y - v_minus.y) / (two * delta);
            }
        }

        if k > 0 && k < nz - 1 {
            if let (Some(w_plus), Some(w_minus)) = (
                flow_field.velocity.get(i, j, k + 1),
                flow_field.velocity.get(i, j, k - 1),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                s33 = (w_plus.z - w_minus.z) / (two * delta);
            }
        }

        // Off-diagonal components: Sij = 0.5 * (∂ui/∂xj + ∂uj/∂xi)
        let mut s12 = T::zero();
        let mut s13 = T::zero();
        let mut s23 = T::zero();

        // S12 = 0.5 * (∂u/∂y + ∂v/∂x)
        if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
            if let (Some(u_jp), Some(u_jm), Some(v_ip), Some(v_im)) = (
                flow_field.velocity.get(i, j + 1, k),
                flow_field.velocity.get(i, j - 1, k),
                flow_field.velocity.get(i + 1, j, k),
                flow_field.velocity.get(i - 1, j, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                let du_dy = (u_jp.x - u_jm.x) / (two * delta);
                let dv_dx = (v_ip.y - v_im.y) / (two * delta);
                s12 = (du_dy + dv_dx) / two;
            }
        }

        // S13 = 0.5 * (∂u/∂z + ∂w/∂x)
        if i > 0 && i < nx - 1 && k > 0 && k < nz - 1 {
            if let (Some(u_kp), Some(u_km), Some(w_ip), Some(w_im)) = (
                flow_field.velocity.get(i, j, k + 1),
                flow_field.velocity.get(i, j, k - 1),
                flow_field.velocity.get(i + 1, j, k),
                flow_field.velocity.get(i - 1, j, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                let du_dz = (u_kp.x - u_km.x) / (two * delta);
                let dw_dx = (w_ip.z - w_im.z) / (two * delta);
                s13 = (du_dz + dw_dx) / two;
            }
        }

        // S23 = 0.5 * (∂v/∂z + ∂w/∂y)
        if j > 0 && j < ny - 1 && k > 0 && k < nz - 1 {
            if let (Some(v_kp), Some(v_km), Some(w_jp), Some(w_jm)) = (
                flow_field.velocity.get(i, j, k + 1),
                flow_field.velocity.get(i, j, k - 1),
                flow_field.velocity.get(i, j + 1, k),
                flow_field.velocity.get(i, j - 1, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                let dv_dz = (v_kp.y - v_km.y) / (two * delta);
                let dw_dy = (w_jp.z - w_jm.z) / (two * delta);
                s23 = (dv_dz + dw_dy) / two;
            }
        }

        // Strain rate magnitude: |S| = sqrt(2 * Sij * Sij)
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        let s_mag_sq =
            two * (s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23));
        s_mag_sq.sqrt()
    }
}

impl<T: RealField + Copy + FromPrimitive> TurbulenceModel<T> for SmagorinskyModel<T> {
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let delta = T::from_f64(1.0 / nx as f64).unwrap_or_else(T::one);
        let prefactor = self.cs * self.cs * delta * delta;

        let mut viscosity = Vec::with_capacity(nx * ny * nz);

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let strain_rate =
                        self.calculate_strain_rate_at_point(flow_field, i, j, k, delta);
                    viscosity.push(prefactor * strain_rate);
                }
            }
        }
        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // For Smagorinsky model, estimate TKE from strain rate
        // k ≈ (Cs * Δ * |S|)²
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let delta = T::from_f64(1.0 / nx as f64).unwrap_or_else(T::one);
        let cs_delta = self.cs * delta;

        let mut tke = Vec::with_capacity(nx * ny * nz);

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let strain_rate =
                        self.calculate_strain_rate_at_point(flow_field, i, j, k, delta);
                    let cs_delta_s = cs_delta * strain_rate;
                    tke.push(cs_delta_s * cs_delta_s);
                }
            }
        }
        tke
    }

    fn name(&self) -> &'static str {
        "Smagorinsky"
    }
}

/// Mixing length turbulence model
#[derive(Debug, Clone)]
pub struct MixingLengthModel<T: RealField + Copy> {
    /// Mixing length scale
    pub length_scale: T,
    /// von Kármán constant
    pub kappa: T,
}

impl<T: RealField + Copy + FromPrimitive> MixingLengthModel<T> {
    /// Create a new mixing length model
    pub fn new(length_scale: T) -> Self {
        let kappa = T::from_f64(cfd_core::physics::constants::physics::fluid::VON_KARMAN)
            .unwrap_or_else(T::one);
        Self {
            length_scale,
            kappa,
        }
    }

    /// Calculate velocity gradient magnitude
    #[allow(clippy::similar_names)] // CFD derivatives use standard notation
    fn calculate_velocity_gradient(&self, velocity: &VelocityField<T>, idx: usize) -> T {
        let (nx, ny, nz) = velocity.dimensions;
        let i = idx % nx;
        let j = (idx / nx) % ny;
        let k = idx / (nx * ny);

        let delta = T::from_f64(1.0 / nx as f64).unwrap_or_else(T::one);
        let two = T::from_f64(2.0).unwrap_or_else(T::one);

        let mut grad_u_sq = T::zero();

        // Calculate ∂u/∂y for wall-bounded flows (primary gradient)
        if j > 0 && j < ny - 1 {
            if let (Some(u_jp), Some(u_jm)) = (velocity.get(i, j + 1, k), velocity.get(i, j - 1, k))
            {
                let dudy = (u_jp.x - u_jm.x) / (two * delta);
                grad_u_sq = dudy * dudy;
            }
        }

        // Add other gradient components for 3D flows
        if i > 0 && i < nx - 1 {
            if let (Some(u_ip), Some(u_im)) = (velocity.get(i + 1, j, k), velocity.get(i - 1, j, k))
            {
                #[allow(clippy::similar_names)] // CFD derivatives use standard notation
                {
                    let dudx = (u_ip.x - u_im.x) / (two * delta);
                    let dvdx = (u_ip.y - u_im.y) / (two * delta);
                    let dwdx = (u_ip.z - u_im.z) / (two * delta);
                    grad_u_sq = grad_u_sq + dudx * dudx + dvdx * dvdx + dwdx * dwdx;
                }
            }
        }

        if k > 0 && k < nz - 1 {
            if let (Some(u_kp), Some(u_km)) = (velocity.get(i, j, k + 1), velocity.get(i, j, k - 1))
            {
                #[allow(clippy::similar_names)] // CFD derivatives use standard notation
                {
                    let dudz = (u_kp.x - u_km.x) / (two * delta);
                    let dvdz = (u_kp.y - u_km.y) / (two * delta);
                    let dwdz = (u_kp.z - u_km.z) / (two * delta);
                    grad_u_sq = grad_u_sq + dudz * dudz + dvdz * dvdz + dwdz * dwdz;
                }
            }
        }

        grad_u_sq.sqrt()
    }
}

impl<T: RealField + Copy + FromPrimitive> TurbulenceModel<T> for MixingLengthModel<T> {
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // νₜ = l² * |∂u/∂y| (Prandtl's mixing length hypothesis)
        flow_field
            .velocity
            .components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let grad_u = self.calculate_velocity_gradient(&flow_field.velocity, idx);
                self.length_scale * self.length_scale * grad_u
            })
            .collect()
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // Estimate TKE from mixing length and velocity gradient
        // k ≈ (l * |∂u/∂y|)²
        flow_field
            .velocity
            .components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let grad_u = self.calculate_velocity_gradient(&flow_field.velocity, idx);
                let l_grad_u = self.length_scale * grad_u;
                l_grad_u * l_grad_u
            })
            .collect()
    }

    fn name(&self) -> &'static str {
        "MixingLength"
    }
}

/// k-epsilon turbulence model state
#[derive(Debug, Clone)]
pub struct KEpsilonState<T: RealField + Copy> {
    /// Turbulent kinetic energy
    pub k: Vec<T>,
    /// Dissipation rate
    pub epsilon: Vec<T>,
}

/// k-epsilon turbulence model
#[derive(Debug, Clone)]
pub struct KEpsilonModel<T: RealField + Copy> {
    /// Model constants
    pub constants: KEpsilonConstants<T>,
    /// Current state (k and epsilon fields)
    pub state: Option<KEpsilonState<T>>,
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
pub struct KEpsilonConstants<T: RealField + Copy> {
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

impl<T: RealField + Copy + FromPrimitive> KEpsilonConstants<T> {
    /// Create standard k-epsilon constants
    pub fn standard() -> Self {
        Self {
            c_mu: T::from_f64(cfd_core::physics::constants::physics::turbulence::K_EPSILON_C_MU)
                .unwrap_or_else(T::one),
            c_1: T::from_f64(cfd_core::physics::constants::physics::turbulence::K_EPSILON_C1)
                .unwrap_or_else(T::one),
            c_2: T::from_f64(cfd_core::physics::constants::physics::turbulence::K_EPSILON_C2)
                .unwrap_or_else(T::one),
            sigma_k: T::from_f64(cfd_core::physics::constants::physics::turbulence::K_EPSILON_SIGMA_K)
                .unwrap_or_else(T::one),
            sigma_epsilon: T::from_f64(
                cfd_core::physics::constants::physics::turbulence::K_EPSILON_SIGMA_EPSILON,
            )
            .unwrap_or_else(T::one),
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for KEpsilonModel<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive> KEpsilonModel<T> {
    /// Create a new k-epsilon model with standard constants
    #[must_use]
    pub fn new() -> Self {
        Self {
            constants: KEpsilonConstants::standard(),
            state: None,
        }
    }

    /// Create with custom constants
    pub fn with_constants(constants: KEpsilonConstants<T>) -> Self {
        Self {
            constants,
            state: None,
        }
    }

    /// Initialize state with high Reynolds number approximation
    pub fn initialize_state(&mut self, flow_field: &FlowField<T>) {
        let n = flow_field.velocity.components.len();

        // Initialize k based on turbulence intensity (typically 1-5% of mean flow)
        let turbulence_intensity = T::from_f64(0.05).unwrap_or_else(T::one);
        let mut k_field = Vec::with_capacity(n);

        for vel in &flow_field.velocity.components {
            let u_mag = vel.norm();
            // k = 3/2 * (U * I)^2 where I is turbulence intensity
            let three_half = T::from_f64(1.5).unwrap_or_else(T::one);
            let k = three_half * (u_mag * turbulence_intensity).powi(2);
            k_field.push(k);
        }

        // Initialize epsilon based on mixing length scale
        // ε = C_μ^(3/4) * k^(3/2) / l
        let mixing_length = T::from_f64(0.1).unwrap_or_else(T::one); // 10% of domain size
        let mut epsilon_field = Vec::with_capacity(n);

        let c_mu_34 = self
            .constants
            .c_mu
            .powf(T::from_f64(0.75).unwrap_or_else(T::one));
        for &k in &k_field {
            let epsilon = c_mu_34 * k.powf(T::from_f64(1.5).unwrap_or_else(T::one)) / mixing_length;
            epsilon_field.push(epsilon);
        }

        self.state = Some(KEpsilonState {
            k: k_field,
            epsilon: epsilon_field,
        });
    }
}

impl<T: RealField + Copy + FromPrimitive> TurbulenceModel<T> for KEpsilonModel<T> {
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // νₜ = C_μ * k² / ε
        match &self.state {
            Some(state) => state
                .k
                .iter()
                .zip(state.epsilon.iter())
                .map(|(&k, &eps)| {
                    if eps > T::from_f64(1e-10).unwrap_or_else(T::zero) {
                        self.constants.c_mu * k * k / eps
                    } else {
                        T::zero()
                    }
                })
                .collect(),
            None => {
                // If not initialized, return zero viscosity
                vec![T::zero(); flow_field.velocity.components.len()]
            }
        }
    }

    fn turbulent_kinetic_energy(&self, _flow_field: &FlowField<T>) -> Vec<T> {
        match &self.state {
            Some(state) => state.k.clone(),
            None => Vec::new(),
        }
    }

    fn name(&self) -> &'static str {
        "k-epsilon"
    }
}

impl<T: RealField + Copy + FromPrimitive> RANSModel<T> for KEpsilonModel<T> {
    fn dissipation_rate(&self, _flow_field: &FlowField<T>) -> Vec<T> {
        match &self.state {
            Some(state) => state.epsilon.clone(),
            None => Vec::new(),
        }
    }

    fn constants(&self) -> &dyn std::any::Any {
        &self.constants
    }
}
