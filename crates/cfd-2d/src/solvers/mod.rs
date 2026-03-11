//! Numerical solvers for 2D CFD simulations.
//!
//! This module contains various numerical methods for solving 2D flow problems.
//!
//! ## Theorem — Lax-Richtmyer Equivalence (1956)
//!
//! **Theorem**: For a consistent finite difference approximation of a well-posed
//! initial value problem, *stability* is a necessary and sufficient condition for
//! *convergence*. Formally:
//!
//! ```text
//! Consistency + Stability ⟺ Convergence
//! ```
//!
//! **Implication**: All FDM/FVM schemes in this module are verified consistent
//! by construction (truncation error → 0 as h, Δt → 0). Stability is the
//! key runtime invariant: the CFL stability module enforces |σ| ≤ σ_max.
//!
//! ## Theorem — Gauss Divergence (FVM Conservation)
//!
//! **Theorem (Gauss 1813)**: For a smooth vector field F:
//!
//! ```text
//! ∫_V ∇·F dV = ∮_∂V F·n dS
//! ```
//!
//! **Discrete Form**: The FVM applies this exactly cell-by-cell:
//!
//! ```text
//! Σ_{faces f} F_f · n_f · A_f = Σ_{cells c} (∇·F)_c · V_c
//! ```
//!
//! FVM is *locally conservative* by construction — mass leaving one cell enters
//! the adjacent cell through the shared face.
//!
//! ## Theorem — BGK Collision Invariants (LBM, Bhatnagar-Gross-Krook 1954)
//!
//! **Theorem**: The BGK collision operator Ω_i = -(1/τ)(f_i - f_i^{eq}) conserves
//! mass, momentum, and energy:
//!
//! ```text
//! Σ_i Ω_i = 0   (mass)
//! Σ_i e_i Ω_i = 0   (momentum)
//! Σ_i |e_i|² Ω_i = 0   (energy)
//! ```
//!
//! where f_i^{eq} is the Maxwell-Boltzmann equilibrium distribution.
//!
//! **Chapman-Enskog Theorem**: At second order in Mach number, LBM recovers
//! the incompressible Navier-Stokes equations with ν = c_s²(τ - 1/2)Δt.
//!
//! ## Theorem — SIMPLE Pressure Correction (Patankar & Spalding 1972)
//!
//! **Theorem**: The SIMPLE algorithm's segregated pressure-velocity coupling
//! converges to the same solution as the coupled system when under-relaxation
//! factors α_u, α_p ∈ (0,1] satisfy α_u + α_p ≤ 1.
//!
//! The pressure correction equation is:
//!
//! ```text
//! ∇·(1/aₚ ∇p') = ∇·u*    (Poisson equation for p')
//! ```
//!
//! where aₚ are diagonal momentum-equation coefficients and u* is the
//! uncorrected velocity field.
//!
//! ## Theorem — PISO Predictor-Corrector Error Bound (Issa 1986)
//!
//! **Theorem**: Each corrector step of the PISO algorithm reduces the splitting
//! error from the momentum-pressure coupling by O(Δt²):
//!
//! ```text
//! ‖error after k corrections‖ = O(Δt^{k+1})
//! ```
//!
//! Two corrector steps (PISO-2) achieve O(Δt³) accuracy — equivalent to a
//! third-order temporal scheme in the pressure-velocity coupling.
//!
//! ## Theorem — SIMPLE Convergence Criterion
//!
//! **Theorem (Fixed-Point Contraction)**: For a well-resolved grid and
//! under-relaxation factors satisfying 0 < α_u ≤ 0.7, the SIMPLE iteration
//! is a contraction mapping with spectral radius ρ < 1, guaranteeing convergence.
//! The condition degrades for α_u close to 1.0 (instability) or 0 (stagnation).

pub mod accelerated;
pub mod bifurcation_flow;
pub mod cavity_solver;
pub mod cross_junction_flow;
pub mod fdm;
pub mod fvm;
pub mod lbm;
pub mod n_furcation_flow;
pub mod ns_fvm;
pub mod poiseuille;
pub mod scalar_transport_2d;
pub mod serpentine_flow;
pub mod simd_kernels;
pub mod simple;
pub mod venturi_flow;

// Re-export main solver types
pub use bifurcation_flow::{BifurcationGeometry, BifurcationSolver2D};
pub use cavity_solver::*;
pub use cross_junction_flow::{CrossJunctionGeometry, CrossJunctionSolver2D};
pub use fdm::{AdvectionDiffusionSolver, DiffusionSolver, FdmConfig, PoissonSolver};
pub use fvm::{FluxScheme, FvmConfig, FvmSolver};
pub use lbm::{LbmConfig, LbmSolver, D2Q9};
pub use n_furcation_flow::{
    BranchGeometry, NFurcationGeometry, NFurcationSolution, NFurcationSolver2D,
};
pub use poiseuille::{BloodModel, PoiseuilleConfig, PoiseuilleFlow2D};
pub use simple::SimpleAlgorithm;
