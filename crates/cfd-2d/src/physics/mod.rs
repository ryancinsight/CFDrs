//! Physics models for 2D CFD simulations
//!
//! This module contains physical models including energy, momentum, turbulence, and vorticity.
//!
//! ## Mathematical Foundation
//!
//! ### Theorem: Navier-Stokes Equations (Incompressible Flow)
//!
//! **Statement**: The motion of incompressible Newtonian fluids is governed by the Navier-Stokes
//! equations, which express conservation of momentum and mass.
//!
//! **Continuity Equation (Mass Conservation)**:
//! ```text
//! ∇·u = 0
//! ```
//!
//! **Momentum Equation (Conservation of Momentum)**:
//! ```text
//! ∂u/∂t + (u·∇)u = -∇p/ρ + ν ∇²u + f
//! ```
//!
//! where:
//! - u is the velocity vector field [m/s]
//! - p is the pressure field [Pa]
//! - ρ is the fluid density [kg/m³]
//! - ν = μ/ρ is the kinematic viscosity [m²/s]
//! - f represents body forces [m/s²]
//!
//! **Derivation**: The Navier-Stokes equations are derived from:
//! 1. **Conservation of Momentum**: Newton's second law applied to fluid elements
//! 2. **Stress Tensor**: Cauchy stress tensor for Newtonian fluids
//! 3. **Continuity**: Conservation of mass for incompressible flow
//!
//! Starting from the general conservation of momentum:
//!
//! d(ρu)/dt = ∇·σ + ρf
//!
//! where σ is the Cauchy stress tensor. For Newtonian fluids:
//!
//! σ = -pI + μ[(∇u) + (∇u)ᵀ] - (2μ/3)(∇·u)I
//!
//! For incompressible flow (∇·u = 0), this simplifies to:
//!
//! ∂u/∂t + (u·∇)u + ∇p/ρ = ν ∇²u + f
//!
//! **Assumptions**:
//! 1. **Incompressible Flow**: ρ = constant, ∇·u = 0
//! 2. **Newtonian Fluid**: Linear stress-strain relationship
//! 3. **Continuum Hypothesis**: Fluid treated as continuous medium
//! 4. **No Chemical Reactions**: Constant composition
//! 5. **No External Forces**: f = 0 (or gravity, etc.)
//!
//! **Validity Conditions**: Applies to incompressible Newtonian fluids in continuum regime
//! (Knudsen number Kn << 1).
//!
//! **Literature**: Navier, C.L.M.H. (1823). "Mémoire sur les lois du mouvement des fluides".
//! Mémoires de l'Académie Royale des Sciences, 6, 389-440.
//!
//! Stokes, G.G. (1845). "On the theories of the internal friction of fluids in motion".
//! Transactions of the Cambridge Philosophical Society, 8, 287-319.
//!
//! ### Theorem: Continuity Equation (Mass Conservation Principle)
//!
//! **Statement**: For incompressible fluids, the density remains constant following fluid
//! particles, leading to the divergence-free velocity field condition.
//!
//! **Mathematical Formulation**:
//! ```text
//! Dρ/Dt = 0  ⇒  ∇·u = 0
//! ```
//!
//! where D/Dt is the material derivative.
//!
//! **Proof**: The continuity equation for compressible flow is:
//!
//! ∂ρ/∂t + ∇·(ρu) = 0
//!
//! For incompressible flow, ρ = constant, so:
//!
//! ∇·u = 0
//!
//! **Physical Interpretation**: The continuity equation ensures that fluid volume is conserved.
//! Incompressibility means that fluid elements cannot be compressed or expanded.
//!
//! **Boundary Conditions**: The no-slip condition u·n = u_wall on solid boundaries, combined
//! with incompressibility, ensures proper mass conservation.
//!
//! **Assumptions**:
//! 1. Incompressible fluid (Mach number Ma << 1)
//! 2. No mass sources or sinks
//! 3. Constant fluid properties
//!
//! **Validity Conditions**: Valid for low-speed flows where density variations are negligible.
//!
//! **Literature**: Landau, L.D., & Lifshitz, E.M. (1987). Fluid Mechanics (2nd ed.).
//! Pergamon Press. §1.1-1.2.
//!
//! ### Theorem: Momentum Equation (Newton's Second Law for Fluids)
//!
//! **Statement**: The momentum equation expresses Newton's second law applied to fluid
//! elements, including convective acceleration, pressure forces, viscous forces, and body forces.
//!
//! **Mathematical Formulation**:
//! ```text
//! ρ Du/Dt = -∇p + μ ∇²u + ρf
//! ```
//!
//! where Du/Dt = ∂u/∂t + (u·∇)u is the material derivative.
//!
//! **Physical Interpretation**:
//! - **∂u/∂t**: Local acceleration
//! - **(u·∇)u**: Convective acceleration (nonlinear term)
//! - **-∇p**: Pressure gradient force
//! - **μ ∇²u**: Viscous diffusion (momentum diffusion)
//! - **ρf**: Body forces (gravity, electromagnetic, etc.)
//!
//! **Derivation**: Apply Newton's second law to a fluid element:
//!
//! Force = mass × acceleration
//!
//! The forces acting on a fluid element are:
//! 1. **Pressure forces**: Surface forces from pressure differences
//! 2. **Viscous forces**: Surface forces from viscous stresses
//! 3. **Body forces**: Volume forces (gravity, etc.)
//!
//! **Assumptions**:
//! 1. Newtonian fluid with constant viscosity
//! 2. Incompressible flow
//! 3. Continuum hypothesis valid
//! 4. No chemical reactions or phase changes
//!
//! **Validity Conditions**: Applies to laminar and turbulent flows of Newtonian fluids
//! in continuum regime.
//!
//! **Literature**: Batchelor, G.K. (1967). An Introduction to Fluid Dynamics.
//! Cambridge University Press. Chapter 2.
//!
//! ### Theorem: Incompressibility Constraint (Pressure Poisson Equation)
//!
//! **Statement**: For incompressible flows, the pressure field is determined by solving
//! a Poisson equation derived from the continuity constraint.
//!
//! **Mathematical Formulation**: Taking the divergence of the momentum equation:
//!
//! ∇²p = -ρ ∇·[(u·∇)u] + μ ∇²(∇·u) + ∇·(ρf)
//!
//! For incompressible flow (∇·u = 0), this simplifies to:
//!
//! ```text
//! ∇²p = -ρ ∇·[(u·∇)u] + ∇·(ρf)
//! ```
//!
//! **Physical Interpretation**: The pressure field acts as a Lagrange multiplier to enforce
//! the incompressibility constraint ∇·u = 0.
//!
//! **Numerical Solution**: In fractional-step methods (SIMPLE, PISO), the pressure equation
//! is solved to correct the velocity field and enforce mass conservation.
//!
//! **Assumptions**:
//! 1. Incompressible flow (∇·u = 0 exactly)
//! 2. Constant density and viscosity
//! 3. Appropriate boundary conditions for pressure
//!
//! **Validity Conditions**: Valid for all incompressible flow solvers using pressure-based
//! methods.
//!
//! **Literature**: Chorin, A.J. (1968). "Numerical solution of the Navier-Stokes equations".
//! Mathematics of Computation, 22(104), 745-762.
//!
//! Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow. Hemisphere Publishing.

pub mod energy;
pub mod immersed_boundary;
pub mod momentum;
pub mod turbulence;
pub mod vorticity_stream;

// Re-export main physics types
pub use energy::EnergyEquationSolver;
pub use immersed_boundary::{BoundaryPoint, ImmersedBoundaryConfig, ImmersedBoundaryMethod};
pub use momentum::{MomentumCoefficients, MomentumComponent, MomentumSolver};
pub use turbulence::{KEpsilonModel, WallFunction};
pub use vorticity_stream::VorticityStreamSolver;
