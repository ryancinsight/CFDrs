//! Volume of Fluid (VOF) method for interface tracking in 3D multiphase flows.
//!
//! The VOF method tracks interfaces by advecting volume fractions,
//! providing excellent mass conservation properties.
//!
//! ## Theorem — Discrete Volume Conservation (PLIC/Geometric Advection)
//!
//! **Theorem**: Under face-flux balancing with exact geometric clipping,
//! the total phase volume is invariant to machine precision across all timesteps:
//!
//! ```text
//! Σᵢ αᵢⁿ⁺¹ Vᵢ = Σᵢ αᵢⁿ Vᵢ      (∀n ≥ 0)
//! ```
//!
//! where αᵢ ∈ [0,1] is the volume fraction in cell i and Vᵢ its volume.
//!
//! **Proof sketch**: The PLIC reconstruction computes a locally planar interface
//! whose orientation satisfies the Youngs or Swartz normal estimate. The advection
//! flux F_f through face f is computed geometrically as the volume of the polygon
//! swept during Δt. Face fluxes are balanced: what leaves one cell enters the
//! adjacent cell. Global conservation follows by telescoping.
//!
//! ## Theorem — PLIC Reconstruction Accuracy (Rider & Kothe 1998)
//!
//! **Theorem**: For a smooth interface with curvature κ, the PLIC plane reconstructed
//! from target volume fraction α achieves an interface position error of O(h²):
//!
//! ```text
//! |d_reconstructed - d_exact| = O(h²)
//! ```
//!
//! where h is the cell size. This gives a 2nd-order accurate interface representation.
//!
//! ## Theorem — Young-Laplace Pressure Jump (Surface Tension)
//!
//! **Theorem**: At a sharp interface between two fluids with surface tension σ,
//! the pressure jump satisfies the Young-Laplace equation:
//!
//! ```text
//! [p] = p₁ - p₂ = σ (κ₁ + κ₂) = σ · κ_total
//! ```
//!
//! where κ₁, κ₂ are the principal curvatures. In the Continuum Surface Force (CSF)
//! model, this is implemented as a body force f_s = σ κ δ_s n̂ where δ_s smears
//! the interface over O(h) cells.
//!
//! ## Theorem — Maximum Principle for Volume Fractions
//!
//! **Theorem**: For conservative VOF advection with bounded face fluxes, the
//! volume fraction satisfies the maximum principle:
//!
//! ```text
//! α ∈ [0, 1]  is preserved for all time steps
//! ```
//!
//! in the absence of source terms. Values outside [0,1] indicate numerical errors
//! and are clamped with warning.
//!
//! ## Invariants (Runtime)
//! - `VofSolver::step()` must produce αᵢ ∈ [0,1] for all cells i after each step.
//! - Total volume conservation error must be < 1e-12 per step (machine epsilon).
//! - Interface normals must satisfy |n̂| = 1 (unit normal, enforced after reconstruction).

mod advection;
mod bubble_dynamics;
mod cavitation_solver;
mod cavitation_types;
mod config;
mod initialization;
mod plic_geometry;
mod reconstruction;
mod solver;

pub use bubble_dynamics::BubbleDynamicsConfig;
pub use cavitation_solver::CavitationVofSolver;
pub use cavitation_types::{CavitationStatistics, CavitationVofConfig};
pub use config::{constants, VofConfig};
pub use solver::VofSolver;

// Re-export key types for convenience
pub use advection::AdvectionMethod;
pub use plic_geometry::volume_under_plane_3d;
pub use reconstruction::InterfaceReconstruction;
pub use reconstruction::{
    height_function_normal_2d, height_function_normal_3d, youngs_normal_2d,
};
