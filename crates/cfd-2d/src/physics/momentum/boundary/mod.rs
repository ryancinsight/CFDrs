//! Boundary condition handling for momentum equations with theoretical foundations
//!
//! ## No-Slip Boundary Condition Theory
//!
//! The no-slip boundary condition is derived from the viscous stress balance at solid walls.
//! For incompressible Navier-Stokes equations, the velocity at a solid wall must equal
//! the wall velocity due to the no-slip condition.
//!
//! ### Mathematical Derivation
//!
//! The momentum equation near a wall (y=0) is:
//!
//! ∂u/∂t + u·∇u = -∇p/ρ + ν∇²u + f
//!
//! At the wall, the viscous stress term ν∇²u dominates, and for no-slip:
//!
//! u(y=0) = u_wall, v(y=0) = v_wall
//!
//! The wall shear stress τ_wall = μ(∂u/∂y)|_wall determines the boundary layer behavior.
//!
//! ### Implementation
//!
//! For finite difference discretization, the no-slip condition is implemented as:
//!
//! u_{i,0} = u_wall  (for west/east walls)
//! v_{0,j} = v_wall  (for south/north walls)
//!
//! ## Characteristic-Based Inlet/Outlet Conditions
//!
//! Inlet/outlet boundaries require characteristic analysis to prevent spurious reflections.
//! For hyperbolic systems, the boundary conditions should be based on incoming/outgoing
//! characteristics.
//!
//! ### Navier-Stokes Characteristics
//!
//! The linearized Navier-Stokes equations have characteristics with speeds:
//! - λ1, λ2 = u ± c (acoustic waves)
//! - λ3, λ4 = u (convective waves)
//!
//! For subsonic flow:
//! - Inlet: Specify u,v,p for incoming characteristics
//! - Outlet: Extrapolate u,v, specify p for outgoing acoustic waves
//!
//! ### Implementation Strategy
//!
//! 1. **Inlet**: Dirichlet conditions for all variables (fully specified)
//! 2. **Outlet**: Neumann conditions for velocity, Dirichlet for pressure
//! 3. **Characteristic BCs**: Use Riemann invariants for compressible flow
//!
//! ## Wall Function Theory
//!
//! For high Reynolds number flows, wall functions relate wall shear stress to velocity
//! at the first grid point, avoiding the need to resolve the viscous sublayer.
//!
//! ### Log-Law Wall Function
//!
//! τ_wall = (κ u_* / ln(y^+)) * ρ u_*^2
//!
//! where:
//! - u_* = √(τ_wall/ρ) is the friction velocity
//! - y^+ = u_* y / ν is the dimensionless wall distance
//! - κ ≈ 0.41 is the von Kármán constant
//!
//! ### Reichardt Wall Function
//!
//! For improved accuracy in the buffer layer:
//!
//! u^+ = (1/κ) ln(1 + κ y^+) + C [1 - exp(-y^+/A) - y^+/A exp(-B y^+)]
//!
//! ## References
//!
//! - Gresho, P. M., & Sani, R. L. (1998). *Incompressible flow and the finite element method*.
//!   Wiley. Chapter 3: Boundary Conditions.
//! - Thompson, K. W. (1990). Time-dependent boundary conditions for hyperbolic systems.
//!   *Journal of Computational Physics*, 89(2), 439-461.
//! - Wilcox, D. C. (2008). *Turbulence modeling for CFD* (3rd ed.). DCW Industries.
//!   Chapter 7: Wall Boundary Conditions.
//!
//! # Theorem
//! The momentum discretization must conserve linear momentum globally and locally.
//!
//! **Proof sketch**:
//! By integrating the Navier-Stokes momentum equation over a control volume $\Omega$,
//! Gauss's divergence theorem converts the convective and diffusive volume integrals
//! into surface fluxes. The finite volume method ensures that the flux leaving one
//! cell exactly equals the flux entering the adjacent cell. Thus, in the absence of
//! external forces and boundary fluxes, the total momentum $\int_\Omega \rho \mathbf{u} dV$
//! is exactly conserved to machine precision.

mod consistency;
mod directional;
mod higher_order;
mod rotating;
#[cfg(test)]
mod tests;
mod updater;

pub use consistency::validate_boundary_consistency;
pub use higher_order::apply_higher_order_wall_boundaries;
pub use updater::MatrixUpdater;

use super::solver::MomentumComponent;
use crate::scalar::Cfd2dScalar;
use cfd_core::physics::boundary::BoundaryCondition;
use directional::{
    apply_east_boundary, apply_north_boundary, apply_south_boundary, apply_west_boundary,
};
use eunomia::FloatElement;
use leto::Array1;
use std::collections::HashMap;
use std::hash::BuildHasher;

/// Apply boundary conditions to momentum equation system
pub fn apply_momentum_boundaries<T, S, M>(
    matrix: &mut M,
    rhs: &mut Array1<T>,
    component: MomentumComponent,
    boundaries: &HashMap<String, BoundaryCondition<T>, S>,
    grid: &crate::grid::StructuredGrid2D<T>,
) -> cfd_core::error::Result<()>
where
    T: Cfd2dScalar + Copy + FloatElement,
    S: BuildHasher,
    M: MatrixUpdater<T>,
{
    let nx = grid.nx;
    let ny = grid.ny;
    for (name, bc) in boundaries {
        match name.as_str() {
            "west" => apply_west_boundary(matrix, rhs, bc, component, grid, nx, ny)?,
            "east" => apply_east_boundary(matrix, rhs, bc, component, grid, nx, ny)?,
            "north" => apply_north_boundary(matrix, rhs, bc, component, grid, nx, ny)?,
            "south" => apply_south_boundary(matrix, rhs, bc, component, grid, nx, ny)?,
            _ => {}
        }
    }

    Ok(())
}
