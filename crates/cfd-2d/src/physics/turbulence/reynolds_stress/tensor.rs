//! Reynolds stress tensor storage type.
//!
//! Defines [`ReynoldsStressTensor`], the 2D second-moment tensor that carries
//! all six independent Reynolds stress components together with the scalar
//! turbulent kinetic energy `k` and dissipation rate `ε`.
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \overline{u_i^\prime u_j^\prime}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\overline{u_i^\prime u_i^\prime} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

use nalgebra::{DMatrix, RealField};

/// Reynolds stress tensor storage (6 independent components for 2D: xx, xy, yy)
///
/// ## Realizability Invariants (Lumley, 1978)
/// - ⟨u'u'⟩ ≥ 0, ⟨v'v'⟩ ≥ 0  (positive semi-definite)
/// - |⟨u'v'⟩| ≤ √(⟨u'u'⟩ ⟨v'v'⟩)  (Cauchy-Schwarz)
/// - −2/3 ≤ b_ij ≤ 2/3  (Lumley triangle)
/// - k = (1/2)(⟨u'u'⟩ + ⟨v'v'⟩)  (definition of TKE)
#[derive(Debug, Clone)]
pub struct ReynoldsStressTensor<T: RealField + Copy> {
    /// ⟨u'u'⟩ — streamwise normal stress
    pub xx: DMatrix<T>,
    /// ⟨u'v'⟩ — shear stress
    pub xy: DMatrix<T>,
    /// ⟨v'v'⟩ — wall-normal normal stress
    pub yy: DMatrix<T>,
    /// Turbulent kinetic energy k = (1/2)(⟨u'u'⟩ + ⟨v'v'⟩)
    pub k: DMatrix<T>,
    /// Dissipation rate ε
    pub epsilon: DMatrix<T>,
    /// Optional anisotropic dissipation tensor ε_xx component
    pub epsilon_xx: Option<DMatrix<T>>,
    /// Optional anisotropic dissipation tensor ε_xy component
    pub epsilon_xy: Option<DMatrix<T>>,
    /// Optional anisotropic dissipation tensor ε_yy component
    pub epsilon_yy: Option<DMatrix<T>>,
}
