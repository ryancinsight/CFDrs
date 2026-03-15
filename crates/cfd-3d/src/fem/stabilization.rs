//! SUPG/PSPG stabilization for FEM
//!
//! # Theorem — SUPG Stability (Brooks & Hughes 1982)
//!
//! The Streamline Upwind Petrov–Galerkin (SUPG) method modifies the test
//! functions by adding a streamline perturbation:
//!
//! ```text
//! w_h = N_i + τ_SUPG (u · ∇N_i)
//! ```
//!
//! This restores coercivity of the bilinear form for advection-dominated flows
//! ($Pe_h > 1$), ensuring nodal stability without introducing excessive
//! numerical diffusion.
//!
//! # Theorem — PSPG Circumvention of LBB (Tezduyar 1991)
//!
//! The Pressure-Stabilising Petrov–Galerkin (PSPG) method adds a pressure
//! stabilisation term:
//!
//! ```text
//! B_PSPG = Σ_e ∫_{Ω_e} τ_PSPG ∇q · R_m dΩ
//! ```
//!
//! where $R_m$ is the momentum residual. This permits equal-order
//! velocity-pressure interpolation by circumventing the Babuška–Brezzi (LBB)
//! inf-sup condition.
//!
//! # Theorem — Stabilisation Parameter (Tezduyar & Osawa 2000)
//!
//! The element-level stabilisation parameter is:
//!
//! ```text
//! τ = [(2/Δt)² + (2|u|/h)² + (4ν/h²)²]^{-1/2}
//! ```
//!
//! This reduces to the correct limits: $\tau \to h/(2|u|)$ for $Pe_h \gg 1$
//! and $\tau \to h^2/(4\nu)$ for $Pe_h \ll 1$.
//!
//! References:
//! - Brooks, A.N. and Hughes, T.J.R. (1982). "Streamline upwind/Petrov-Galerkin formulations
//!   for convection dominated flows with particular emphasis on the incompressible Navier-Stokes equations"
//! - Tezduyar, T.E. (1991). "Stabilized finite element formulations for incompressible flow computations"

use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

// Named constants for stabilization
const TWO: f64 = 2.0;
const FOUR: f64 = 4.0;

/// SUPG/PSPG stabilization parameters
pub struct StabilizationParameters<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Element size
    h: T,
    /// Kinematic viscosity
    nu: T,
    /// Velocity magnitude
    u_mag: T,
    /// Time step (for transient problems)
    dt: Option<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>
    StabilizationParameters<T>
{
    /// Create new stabilization parameters
    pub fn new(h: T, nu: T, velocity: Vector3<T>, dt: Option<T>) -> Self {
        let u_mag = velocity.norm();
        Self { h, nu, u_mag, dt }
    }

    /// Calculate SUPG stabilization parameter tau
    ///
    /// Based on Tezduyar (1991) formulation:
    /// τ = [(2/Δt)² + (2U/h)² + (4ν/h²)²]^(-1/2)
    ///
    /// For steady-state: τ = [(2U/h)² + (4ν/h²)²]^(-1/2)
    pub fn tau_supg(&self) -> T {
        let two = <T as FromPrimitive>::from_f64(TWO)
            .expect("TWO (2.0) is representable in all IEEE 754 types");
        let four = <T as FromPrimitive>::from_f64(FOUR)
            .expect("FOUR (4.0) is representable in all IEEE 754 types");

        // Advection term: (2U/h)²
        let advection_term = if self.u_mag > T::zero() {
            let term = (two * self.u_mag) / self.h;
            term * term
        } else {
            T::zero()
        };

        // Diffusion term: (4ν/h²)²
        let diffusion_term = {
            let term = (four * self.nu) / (self.h * self.h);
            term * term
        };

        // Time term: (2/Δt)² (only for transient)
        let time_term = if let Some(dt) = self.dt {
            let term = two / dt;
            term * term
        } else {
            T::zero()
        };

        // Combined tau
        let sum = time_term + advection_term + diffusion_term;
        if sum > T::zero() {
            T::one() / num_traits::Float::sqrt(sum)
        } else {
            T::zero()
        }
    }

    /// Calculate PSPG stabilization parameter tau
    ///
    /// For pressure stabilization, uses same formulation as SUPG
    pub fn tau_pspg(&self) -> T {
        self.tau_supg()
    }

    /// Calculate element Peclet number
    /// Pe = U*h/(2ν)
    pub fn peclet_number(&self) -> T {
        if self.nu > T::zero() {
            (self.u_mag * self.h)
                / (<T as FromPrimitive>::from_f64(TWO)
                    .expect("TWO (2.0) is representable in all IEEE 754 types")
                    * self.nu)
        } else {
            T::zero()
        }
    }

    /// Calculate element Reynolds number
    /// Re = U*h/ν
    pub fn element_reynolds(&self) -> T {
        if self.nu > T::zero() {
            (self.u_mag * self.h) / self.nu
        } else {
            T::zero()
        }
    }

    /// Get optimal stabilization based on flow regime
    pub fn optimal_tau(&self) -> T {
        let pe = self.peclet_number();
        let tau_supg = self.tau_supg();

        // Adjust based on Peclet number
        if pe < T::one() {
            // Diffusion-dominated: reduce stabilization
            tau_supg * pe
        } else if pe
            > <T as FromPrimitive>::from_f64(100.0)
                .expect("100.0 is representable in all IEEE 754 types")
        {
            // Highly advection-dominated: use full stabilization
            tau_supg
        } else {
            // Intermediate regime: smooth transition
            let factor = num_traits::Float::max(T::one() - T::one() / pe, T::zero());
            tau_supg * factor
        }
    }
}

/// Default grad-div stabilization constant γ (Olshanskii & Reusken 2004).
///
/// Typical range: γ ∈ [0.1, 10]. The default value of 1.0 provides a good
/// balance between mass conservation improvement and conditioning.
pub const GRAD_DIV_GAMMA_DEFAULT: f64 = 1.0;

/// Olshanskii & Reusken (2004) grad-div stabilization parameter.
///
/// ## Theorem — Grad-Div Stabilization for Incompressible FEM
///
/// For inf-sup stable velocity-pressure finite element pairs (e.g.,
/// Taylor-Hood P2/P1), the grad-div stabilization term:
///
/// ```text
/// τ_div · (∇·v_h, ∇·w_h)
/// ```
///
/// improves pointwise mass conservation from O(h^k) to O(h^{k+1}) without
/// degrading the velocity error estimate. The stabilization parameter is:
///
/// ```text
/// τ_div = γ · h²
/// ```
///
/// where γ ∈ [0.1, 10] is a user-tunable constant (default γ = 1.0)
/// and h is the local element diameter.
///
/// **Proof sketch**: The bilinear form B(v_h, w_h) = τ_div·(∇·v_h, ∇·w_h)
/// adds a penalty on divergence-free violation. By the Lax-Milgram theorem,
/// the augmented problem remains well-posed if τ_div > 0, and the optimal
/// error estimate follows from Brezzi's stability framework (1974).
///
/// **Reference**: Olshanskii, M.A. & Reusken, A. (2004). "Grad-div
/// stabilization for the Stokes equations", *Math. Comp.* 73(248):1699-1718.
pub fn grad_div_parameter(h_element: f64, gamma: f64) -> f64 {
    gamma * h_element * h_element
}

/// Calculate element size for different element types
pub fn calculate_element_size<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    vertices: &[Vector3<T>],
    velocity_direction: &Vector3<T>,
) -> T {
    // For tetrahedral elements (4 vertices)
    if vertices.len() == 4 {
        calculate_tetrahedral_size(vertices, velocity_direction)
    } else if vertices.len() == 8 {
        // For hexahedral elements
        calculate_hexahedral_size(vertices, velocity_direction)
    } else {
        // Default: use minimum edge length
        calculate_min_edge_length(vertices)
    }
}

/// Calculate size for tetrahedral element
fn calculate_tetrahedral_size<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    vertices: &[Vector3<T>],
    velocity_direction: &Vector3<T>,
) -> T {
    if velocity_direction.norm() > T::zero() {
        // Directional element size in flow direction
        let dir = velocity_direction.normalize();

        // Project element onto flow direction
        let mut min_proj = <T as RealField>::max_value().unwrap_or_else(T::one);
        let mut max_proj = <T as RealField>::min_value().unwrap_or_else(T::zero);

        for vertex in vertices {
            let proj = vertex.dot(&dir);
            min_proj = num_traits::Float::min(min_proj, proj);
            max_proj = num_traits::Float::max(max_proj, proj);
        }

        num_traits::Float::abs(max_proj - min_proj)
    } else {
        // No flow: use characteristic length
        calculate_min_edge_length(vertices)
    }
}

/// Calculate size for hexahedral element
fn calculate_hexahedral_size<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    vertices: &[Vector3<T>],
    velocity_direction: &Vector3<T>,
) -> T {
    // Similar to tetrahedral but for 8 vertices
    calculate_tetrahedral_size(vertices, velocity_direction)
}

/// Calculate minimum edge length of element
fn calculate_min_edge_length<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    vertices: &[Vector3<T>],
) -> T {
    let mut min_length = <T as RealField>::max_value().unwrap_or_else(T::one);

    for i in 0..vertices.len() {
        for j in i + 1..vertices.len() {
            let edge_length = (vertices[i] - vertices[j]).norm();
            min_length = num_traits::Float::min(min_length, edge_length);
        }
    }

    min_length
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_grad_div_parameter_scales_with_h_squared() {
        let gamma = GRAD_DIV_GAMMA_DEFAULT; // 1.0
        let h = 0.1;
        let tau = grad_div_parameter(h, gamma);
        // τ = γ · h² = 1.0 · 0.01 = 0.01
        assert_relative_eq!(tau, 0.01, epsilon = 1e-15);
    }

    #[test]
    fn test_grad_div_parameter_zero_gamma() {
        // γ = 0 disables stabilization entirely
        let tau = grad_div_parameter(0.5, 0.0);
        assert_eq!(tau, 0.0);
    }

    #[test]
    fn test_grad_div_parameter_positive() {
        // For any positive γ and positive h, τ must be non-negative
        for &gamma in &[0.0, 0.1, 1.0, 5.0, 10.0] {
            for &h in &[0.0, 0.001, 0.01, 0.1, 1.0] {
                let tau = grad_div_parameter(h, gamma);
                assert!(
                    tau >= 0.0,
                    "τ must be non-negative for γ={gamma}, h={h}, got {tau}"
                );
            }
        }
    }

    #[test]
    fn test_grad_div_parameter_different_gamma_values() {
        let h = 0.2;
        // τ = γ · 0.04
        assert_relative_eq!(grad_div_parameter(h, 0.1), 0.004, epsilon = 1e-15);
        assert_relative_eq!(grad_div_parameter(h, 1.0), 0.04, epsilon = 1e-15);
        assert_relative_eq!(grad_div_parameter(h, 10.0), 0.4, epsilon = 1e-15);
    }
}
