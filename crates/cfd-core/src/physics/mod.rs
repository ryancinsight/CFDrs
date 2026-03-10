//! Physics module consolidating fluid properties, boundary conditions, and constants.
//!
//! This module provides a unified interface to the physical models and properties
//! used in CFD simulations, following a deep vertical module structure.
//!
//! ## Theorem — Navier-Stokes Existence & Uniqueness (Leray 1934)
//!
//! **Theorem**: For incompressible Newtonian flow with ν > 0 and initial condition
//! u₀ ∈ L²(Ω), there exists a global weak solution u ∈ L²(0,T; H¹(Ω)) satisfying
//! the Navier-Stokes equations in the distributional sense. Uniqueness holds in 2D;
//! in 3D, uniqueness is an open Millennium Prize Problem for large data.
//!
//! **Operational Implication**: All solvers in this suite operate in the **proven**
//! existence regime (small Reynolds number or 2D flow). 3D turbulent solvers use
//! RANS/LES closures that circumvent the 3D uniqueness gap.
//!
//! ## Theorem — Dirichlet / Neumann Well-Posedness
//!
//! **Theorem (Lax 1954 / Lions-Magenes 1972)**: The elliptic boundary value problem:
//!
//! ```text
//! -∇²u = f  on Ω,   u = g_D  on Γ_D,   ∂u/∂n = g_N  on Γ_N
//! ```
//!
//! is well-posed (unique solution in H¹(Ω)) if:
//! - Γ_D has positive measure: |Γ_D| > 0 (prevents pressure indeterminacy).
//! - For pure Neumann (Γ_D = ∅): compatibility condition ∫_Ω f dx + ∫_∂Ω g_N ds = 0.
//!
//! **Implementation Invariant**: The boundary condition manager enforces that
//! every pressure field has at least one Dirichlet node to prevent singular systems.
//!
//! ## Theorem — Frame Invariance (Noll 1958 / Coleman-Noll)
//!
//! **Theorem**: For a constitutively admissible fluid, the Cauchy stress tensor σ
//! must satisfy the principle of material frame indifference (objectivity):
//!
//! ```text
//! σ* = Q·σ·Qᵀ  under change of observer Q ∈ SO(3)
//! ```
//!
//! **Implications**: Implemented constitutive models (Newtonian, Power-Law, Casson,
//! Herschel-Bulkley) are frame-invariant. The viscous stress depends only on the
//! symmetric rate-of-deformation tensor D = (∇u + ∇uᵀ)/2, not on frame-dependent
//! (antisymmetric) vorticity components.
//!
//! ## Theorem — Buckingham Π Theorem (Dimensional Analysis)
//!
//! **Theorem (Buckingham 1914)**: A physically consistent equation involving n
//! dimensional variables and k independent fundamental dimensions can be reduced to
//! an equation in n-k dimensionless Π groups.
//!
//! **Application**: Reynolds number Re = ρuL/μ, Weber number We = ρu²L/σ,
//! Womersley number Wo = L√(ω/ν), and Capillary number Ca = μu/σ are the
//! canonical Π groups governing microfluidic flow regimes implemented here.

pub mod boundary;
pub mod cavitation;
pub mod cell_interaction;
pub mod constants;
pub mod fluid;
pub mod fluid_dynamics;
pub mod hemolysis;
pub mod material;

pub mod values;

#[allow(missing_docs)]
pub mod api {
    pub use super::boundary::{
        BoundaryCondition, BoundaryConditionManager, BoundaryConditionSet, WallType,
    };
    pub use super::constants::{mathematical, physics};
    pub use super::fluid::{ConstantPropertyFluid, FluidProperties};
    pub use super::fluid_dynamics::{
        FlowClassifier, FlowField, FlowOperations, FlowRegime, FluidDynamicsService, PressureField,
        RANSModel, RhieChowInterpolation, ScalarField, TurbulenceModel, VelocityField,
    };
    pub use super::hemolysis::{
        BloodTrauma, BloodTraumaSeverity, HemolysisCalculator, HemolysisModel, PlateletActivation,
    };
    pub use super::material::{MaterialDatabase, SolidProperties};
    pub use super::values::{Pressure, ReynoldsNumber, Temperature, Velocity};
}
