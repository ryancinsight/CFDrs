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

mod directional;

use super::solver::MomentumComponent;
use crate::grid::traits::Grid2D;
use cfd_core::physics::boundary::{BoundaryCondition, BoundaryError};
use cfd_math::sparse::SparseMatrixBuilder;
use directional::{
    apply_east_boundary, apply_north_boundary, apply_south_boundary, apply_west_boundary,
};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;
use std::hash::BuildHasher;

/// Apply rotating wall boundary condition: u_wall = ω × r
/// where r is the position vector from center of rotation
fn apply_rotating_wall_bc<T: RealField + Copy + FromPrimitive>(
    component: MomentumComponent,
    omega: &nalgebra::Vector3<T>,
    center: &nalgebra::Vector3<T>,
    grid: &crate::grid::StructuredGrid2D<T>,
    idx: usize,
) -> T {
    let nx = grid.nx;
    let i = idx % nx;
    let j = idx / nx;

    let cell_center = grid.cell_center(i, j).unwrap();
    let x = cell_center.x;
    let y = cell_center.y;

    let r_x = x - center.x;
    let r_y = y - center.y;

    let omega_z = omega.z;

    match component {
        MomentumComponent::U => -omega_z * r_y,
        MomentumComponent::V => omega_z * r_x,
    }
}

/// Apply boundary conditions to momentum equation system
pub fn apply_momentum_boundaries<T, S>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    component: MomentumComponent,
    boundaries: &HashMap<String, BoundaryCondition<T>, S>,
    grid: &crate::grid::StructuredGrid2D<T>,
) -> cfd_core::error::Result<()>
where
    T: RealField + Copy + FromPrimitive,
    S: BuildHasher,
{
    let nx = grid.nx;
    let ny = grid.ny;
    for (name, bc) in boundaries {
        let name: &String = name;
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

/// Get the Dirichlet value for a boundary condition and component
fn get_dirichlet_value<T: RealField + Copy>(
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
) -> Option<T> {
    match bc {
        BoundaryCondition::Dirichlet {
            value,
            component_values,
        } => {
            if let Some(comps) = component_values {
                let idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                if idx < comps.len() {
                    comps[idx].or(Some(*value))
                } else {
                    Some(*value)
                }
            } else {
                Some(*value)
            }
        }
        BoundaryCondition::VelocityInlet { velocity } => {
            let idx = match component {
                MomentumComponent::U => 0,
                MomentumComponent::V => 1,
            };
            Some(velocity[idx])
        }
        BoundaryCondition::Wall { wall_type } => match wall_type {
            cfd_core::physics::boundary::WallType::NoSlip => Some(T::zero()),
            cfd_core::physics::boundary::WallType::Moving { velocity } => {
                let idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                Some(velocity[idx])
            }
            _ => None,
        },
        _ => None,
    }
}

/// Validate boundary condition consistency
pub fn validate_boundary_consistency<T, S>(
    boundaries: &HashMap<String, BoundaryCondition<T>, S>,
    _grid: &crate::grid::StructuredGrid2D<T>,
) -> Result<(), BoundaryError>
where
    T: RealField + Copy + FromPrimitive,
    S: BuildHasher,
{
    for direction in &["north", "south", "east", "west"] {
        if !boundaries.contains_key(*direction) {
            return Err(BoundaryError::InvalidRegion(format!(
                "Missing required boundary: {direction}"
            )));
        }
    }

    for (name, bc) in boundaries {
        if let BoundaryCondition::Periodic { partner } = bc {
            if !boundaries.contains_key(partner) {
                return Err(BoundaryError::InvalidRegion(format!(
                    "Periodic partner '{partner}' not found for boundary '{name}'"
                )));
            }
            if let Some(partner_bc) = boundaries.get(partner) {
                if let BoundaryCondition::Periodic { partner: p2 } = partner_bc {
                    if p2 != name {
                        return Err(BoundaryError::InvalidRegion(format!(
                            "Periodic mismatch: '{name}' points to '{partner}', but '{partner}' points to '{p2}'"
                        )));
                    }
                } else {
                    return Err(BoundaryError::InvalidRegion(format!(
                        "Periodic partner '{partner}' is not periodic"
                    )));
                }
            }
        }
    }

    let corners = [
        ("west", "south"),
        ("east", "south"),
        ("west", "north"),
        ("east", "north"),
    ];

    for (b1_name, b2_name) in corners {
        let b1 = boundaries.get(b1_name).unwrap();
        let b2 = boundaries.get(b2_name).unwrap();

        for component in [MomentumComponent::U, MomentumComponent::V] {
            if let (Some(v1), Some(v2)) = (
                get_dirichlet_value(b1, component),
                get_dirichlet_value(b2, component),
            ) {
                let diff = (v1 - v2).abs();
                let epsilon = T::default_epsilon()
                    * T::from_f64(100.0).unwrap_or_else(num_traits::Zero::zero);
                if diff > epsilon {
                    return Err(BoundaryError::InvalidRegion(format!(
                        "Corner conflict at {b1_name}-{b2_name}: Component {component:?} values mismatch ({v1:?} vs {v2:?})"
                    )));
                }
            }
        }
    }

    Ok(())
}

/// Apply higher-order boundary conditions for improved near-wall accuracy.
/// Implements quadratic extrapolation for better velocity gradients near walls.
pub fn apply_higher_order_wall_boundaries<T, S>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    component: MomentumComponent,
    boundaries: &HashMap<String, BoundaryCondition<T>, S>,
    grid: &crate::grid::StructuredGrid2D<T>,
) -> cfd_core::error::Result<()>
where
    T: RealField + Copy + FromPrimitive,
    S: BuildHasher,
{
    let nx = grid.nx;
    let ny = grid.ny;

    for (name, bc) in boundaries {
        let name: &String = name;
        if let BoundaryCondition::Wall {
            wall_type: cfd_core::physics::boundary::WallType::NoSlip,
        } = bc
        {
            match name.as_str() {
                "west" => apply_higher_order_west_wall(matrix, rhs, component, grid, nx, ny)?,
                "east" => apply_higher_order_east_wall(matrix, rhs, component, grid, nx, ny)?,
                "north" => apply_higher_order_north_wall(matrix, rhs, component, grid, nx, ny)?,
                "south" => apply_higher_order_south_wall(matrix, rhs, component, grid, nx, ny)?,
                _ => {}
            }
        }
    }

    Ok(())
}

/// Quadratic extrapolation: u_0 = (4*u_1 - u_2)/3 for west wall
fn apply_higher_order_west_wall<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    _component: MomentumComponent,
    _grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    let four = T::from_f64(4.0).unwrap_or_else(|| T::one() + T::one() + T::one() + T::one());
    let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());

    for j in 0..ny {
        let idx_0 = j * nx;
        let idx_1 = j * nx + 1;
        let idx_2 = j * nx + 2;

        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, T::one())?;
        rhs[idx_0] = T::zero();
    }

    Ok(())
}

/// Quadratic extrapolation for east wall
fn apply_higher_order_east_wall<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    _component: MomentumComponent,
    _grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    let four = T::from_f64(4.0).unwrap_or_else(|| T::one() + T::one() + T::one() + T::one());
    let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());

    for j in 0..ny {
        let idx_0 = j * nx + nx - 1;
        let idx_1 = j * nx + nx - 2;
        let idx_2 = j * nx + nx - 3;

        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, T::one())?;
        rhs[idx_0] = T::zero();
    }

    Ok(())
}

/// Quadratic extrapolation for north wall
fn apply_higher_order_north_wall<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    _component: MomentumComponent,
    _grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    let four = T::from_f64(4.0).unwrap_or_else(|| T::one() + T::one() + T::one() + T::one());
    let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());

    for i in 0..nx {
        let idx_0 = (ny - 1) * nx + i;
        let idx_1 = (ny - 2) * nx + i;
        let idx_2 = (ny - 3) * nx + i;

        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, T::one())?;
        rhs[idx_0] = T::zero();
    }

    Ok(())
}

/// Quadratic extrapolation for south wall
fn apply_higher_order_south_wall<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    _component: MomentumComponent,
    _grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    _ny: usize,
) -> cfd_core::error::Result<()> {
    let four = T::from_f64(4.0).unwrap_or_else(|| T::one() + T::one() + T::one() + T::one());
    let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());

    for i in 0..nx {
        let idx_0 = i;
        let idx_1 = nx + i;
        let idx_2 = 2 * nx + i;

        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, T::one())?;
        rhs[idx_0] = T::zero();
    }

    Ok(())
}
