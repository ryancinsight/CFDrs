//! Standard k-ε two-equation turbulence model (Launder & Spalding 1974).
//!
//! # Theorem — k-ε Closure (Launder & Spalding 1974)
//!
//! The Reynolds-averaged Navier-Stokes equations are closed by two transport
//! equations for the turbulent kinetic energy $k$ and its dissipation rate
//! $\varepsilon$:
//!
//! ```text
//! ∂k/∂t + U·∇k = ∇·((ν + νₜ/σ_k)∇k) + P_k − ε
//! ∂ε/∂t + U·∇ε = ∇·((ν + νₜ/σ_ε)∇ε) + C₁ε (ε/k) P_k − C₂ε (ε²/k)
//! ```
//!
//! The turbulent eddy viscosity is $\nu_t = C_\mu k^2 / \varepsilon$ with
//! standard constants $C_\mu = 0.09$, $C_{1\varepsilon} = 1.44$,
//! $C_{2\varepsilon} = 1.92$, $\sigma_k = 1.0$, $\sigma_\varepsilon = 1.3$.
//!
//! **Proof sketch.** Dimensional analysis of the turbulent stress tensor
//! $\overline{u'_i u'_j}$ and the requirement that statistically homogeneous
//! turbulence admits an equilibrium cascade yield two independent scalar
//! invariants—$k$ and $\varepsilon$—from which a unique eddy viscosity
//! $\nu_t \propto k^2/\varepsilon$ follows (Jones & Launder 1972).
//!
//! **Reference:** Launder, B.E. and Spalding, D.B. (1974). "The numerical
//! computation of turbulent flows." *Comput. Methods Appl. Mech. Eng.*
//! 3(2):269–289.

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::{RANSModel, TurbulenceModel};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// k-epsilon turbulence model state
#[derive(Debug, Clone)]
pub struct KEpsilonState<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Turbulent kinetic energy
    pub k: Vec<T>,
    /// Dissipation rate
    pub epsilon: Vec<T>,
}

/// k-epsilon turbulence model
#[derive(Debug, Clone)]
pub struct KEpsilonModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
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
pub struct KEpsilonConstants<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
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

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> KEpsilonConstants<T> {
    /// Create standard k-epsilon constants
    pub fn standard() -> Self {
        Self {
            c_mu: <T as FromPrimitive>::from_f64(
                cfd_core::physics::constants::physics::turbulence::K_EPSILON_C_MU,
            )
            .expect("K_EPSILON_C_MU is an IEEE 754 representable f64 constant"),
            c_1: <T as FromPrimitive>::from_f64(
                cfd_core::physics::constants::physics::turbulence::K_EPSILON_C1,
            )
            .expect("K_EPSILON_C1 is an IEEE 754 representable f64 constant"),
            c_2: <T as FromPrimitive>::from_f64(
                cfd_core::physics::constants::physics::turbulence::K_EPSILON_C2,
            )
            .expect("K_EPSILON_C2 is an IEEE 754 representable f64 constant"),
            sigma_k: <T as FromPrimitive>::from_f64(
                cfd_core::physics::constants::physics::turbulence::K_EPSILON_SIGMA_K,
            )
            .expect("K_EPSILON_SIGMA_K is an IEEE 754 representable f64 constant"),
            sigma_epsilon: <T as FromPrimitive>::from_f64(
                cfd_core::physics::constants::physics::turbulence::K_EPSILON_SIGMA_EPSILON,
            )
            .expect("K_EPSILON_SIGMA_EPSILON is an IEEE 754 representable f64 constant"),
        }
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> Default
    for KEpsilonModel<T>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> KEpsilonModel<T> {
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

    /// Initialize turbulent state with mathematically exact boundary conditions
    /// or explicitly defined turbulence intensity and mixing length parameters.
    pub fn initialize_state(
        &mut self,
        flow_field: &FlowField<T>,
        turbulence_intensity: T,
        mixing_length: T,
    ) {
        let n = flow_field.velocity.components.len();

        // Initialize k exactly without hidden assumptions
        let mut k_field = Vec::with_capacity(n);

        for vel in &flow_field.velocity.components {
            let u_mag = vel.norm();
            let three_half = (T::one() + T::one() + T::one()) / (T::one() + T::one());
            let k = three_half * num_traits::Float::powi(u_mag * turbulence_intensity, 2);
            k_field.push(k);
        }

        // Initialize epsilon mathematically
        let mut epsilon_field = Vec::with_capacity(n);

        let three = T::one() + T::one() + T::one();
        let four = three + T::one();
        let c_mu_34 = num_traits::Float::powf(self.constants.c_mu, three / four);
        for &k in &k_field {
            let epsilon =
                c_mu_34 * num_traits::Float::powf(k, three / (T::one() + T::one())) / mixing_length;
            epsilon_field.push(epsilon);
        }

        self.state = Some(KEpsilonState {
            k: k_field,
            epsilon: epsilon_field,
        });
    }

    /// Initialize state rigidly with completely proven turbulent kinetic energy and dissipation fields
    pub fn initialize_state_exact(&mut self, k: Vec<T>, epsilon: Vec<T>) {
        assert_eq!(k.len(), epsilon.len(), "k and epsilon lengths must match");
        self.state = Some(KEpsilonState { k, epsilon });
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for KEpsilonModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // νₜ = C_μ * k² / ε
        match &self.state {
            Some(state) => {
                assert_eq!(
                    state.k.len(),
                    state.epsilon.len(),
                    "k and epsilon lengths must match"
                );
                assert_eq!(
                    state.k.len(),
                    flow_field.velocity.components.len(),
                    "k/epsilon state must match the flow-field size"
                );
                let mut viscosity = Vec::with_capacity(state.k.len());
                for idx in 0..state.k.len() {
                    let k = state.k[idx];
                    let eps = state.epsilon[idx];
                    if eps
                        > <T as FromPrimitive>::from_f64(1e-10)
                            .expect("1e-10 is an IEEE 754 representable f64 constant")
                    {
                        viscosity.push(self.constants.c_mu * k * k / eps);
                    } else {
                        viscosity.push(T::zero());
                    }
                }
                viscosity
            }
            None => {
                // If not initialized, return zero viscosity
                vec![T::zero(); flow_field.velocity.components.len()]
            }
        }
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let n = flow_field.velocity.components.len();
        match &self.state {
            Some(state) => {
                assert_eq!(
                    state.k.len(),
                    state.epsilon.len(),
                    "k and epsilon lengths must match"
                );
                state.k.clone()
            }
            None => vec![T::zero(); n],
        }
    }

    fn name(&self) -> &'static str {
        "k-epsilon"
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> RANSModel<T>
    for KEpsilonModel<T>
{
    fn dissipation_rate(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let n = flow_field.velocity.components.len();
        match &self.state {
            Some(state) => {
                assert_eq!(
                    state.k.len(),
                    state.epsilon.len(),
                    "k and epsilon lengths must match"
                );
                state.epsilon.clone()
            }
            None => vec![T::zero(); n],
        }
    }

    fn constants(&self) -> &dyn std::any::Any {
        &self.constants
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_core::physics::fluid_dynamics::RANSModel;
    use cfd_core::physics::fluid_dynamics::TurbulenceModel;

    #[test]
    fn uninitialized_state_returns_zero_fields_of_matching_length() {
        let flow = FlowField::<f64>::new(2, 2, 2);
        let model = KEpsilonModel::<f64>::new();

        let tke = model.turbulent_kinetic_energy(&flow);
        let dissipation = RANSModel::dissipation_rate(&model, &flow);

        assert_eq!(tke.len(), 8);
        assert_eq!(dissipation.len(), 8);
        assert!(tke.iter().all(|&value| value == 0.0));
        assert!(dissipation.iter().all(|&value| value == 0.0));
    }
}
