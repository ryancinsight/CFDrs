//! 3D bifurcation and trifurcation solvers with FEM validation
//!
//! This module implements full 3D CFD simulations for bifurcating networks with:
//! - Realistic conical/cylindrical branching geometries
//! - Finite Element Method (FEM) for Navier-Stokes equations
//! - Blood flow with non-Newtonian rheology
//! - Complete validation against analytical and literature solutions
//!
//! # Physics Background
//!
//! ## 3D Navier-Stokes Equations
//!
//! In a bifurcation, incompressible flow satisfies:
//!
//! **Momentum equation:**
//! ```text
//! ρ(∂u/∂t + (u·∇)u) = -∇p + ∇·τ + f
//! ```
//!
//! **Continuity equation:**
//! ```text
//! ∇·u = 0
//! ```
//!
//! where:
//! - ρ = fluid density
//! - u = velocity field
//! - p = pressure
//! - τ = stress tensor (depends on fluid model)
//! - f = body forces
//!
//! ## 3D Bifurcation Physics
//!
//! In 3D bifurcations, additional phenomena occur:
//!
//! 1. **Asymmetry**: Flow may not split equally even in symmetric geometry
//! 2. **Secondary flows**: Vorticity generation in curved sections
//! 3. **Wall shear stress**: Varies spatially due to 3D geometry
//! 4. **Blood rheology**: Non-Newtonian effects more pronounced in small vessels
//!
//! ## Branching Laws
//!
//! **Murray's Law (3D):**
//! ```text
//! D₀³ = D₁³ + D₂³
//! ```
//!
//! **Finet's Law (velocity ratio):**
//! ```text
//! u₁/u₂ = (D₂/D₁)²
//! ```
//!
//! # Geometry Modeling
//!
//! The bifurcation consists of:
//! - **Parent branch**: Cylinder with D₀ diameter
//! - **Transition zone**: Conical taper region (most complex)
//! - **Daughter branches**: Two cylinders with D₁, D₂ diameters
//!
//! The transition zone is critical because:
//! - Creates flow acceleration/deceleration zones
//! - Generates secondary flows and vorticity
//! - Causes pressure losses beyond simple Poiseuille
//! - Exhibits wall shear stress concentrations
//!
//! # Validation Strategy
//!
//! 1. **Mesh convergence**: h-refinement study (mesh independence)
//! 2. **Richardson extrapolation**: Observe convergence order
//! 3. **Conservation**: Verify mass conservation < 1e-10
//! 4. **Comparison with 1D**: Map 3D solution to centerline for 1D comparison
//! 5. **Literature**: Compare with published bifurcation data
//!
//! # References
//!
//! - Kassab, G.S. (2006). "Scaling laws of viscous flow in branching vessels"
//! - Alastruey, J., et al. (2012). "Arterial biomechanics: From physiology
//!   to mathematical modeling"
//! - Fung, Y.C. (1997). "Biomechanics: Circulation"

pub mod geometry;
pub mod solver;
pub mod validation;

pub use geometry::{BifurcationGeometry3D, ConicalTransition, BifurcationMesh};
pub use solver::{BifurcationSolver3D, BifurcationConfig3D};
pub use validation::{BifurcationValidator3D, BifurcationValidationResult3D, MeshRefinementConfig};
