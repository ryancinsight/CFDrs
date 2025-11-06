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

use super::solver::MomentumComponent;
use crate::grid::traits::Grid2D;
use cfd_core::boundary::BoundaryCondition;
use cfd_math::sparse::SparseMatrixBuilder;
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
    // Convert linear index to 2D coordinates (row-major order: idx = j * nx + i)
    let nx = grid.nx;
    let i = idx % nx;
    let j = idx / nx;

    // Get cell center position using grid coordinate methods
    let cell_center = grid.cell_center(i, j).unwrap();
    let x = cell_center.x;
    let y = cell_center.y;

    // Position vector from center of rotation (2D: only x,y components)
    let r_x = x - center.x;
    let r_y = y - center.y;

    // Angular velocity vector (2D: only z-component)
    let omega_z = omega.z;

    // Velocity = ω × r = ω_z * (-r_y, r_x) for 2D rotation
    match component {
        MomentumComponent::U => -omega_z * r_y,  // u = -ω_z * r_y
        MomentumComponent::V => omega_z * r_x,   // v = ω_z * r_x
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
    // Apply boundary conditions based on location
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

/// Apply higher-order boundary conditions for improved near-wall accuracy
/// This function implements quadratic extrapolation for better velocity gradients near walls
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

    // Apply higher-order boundary conditions for walls to improve near-wall gradients
    for (name, bc) in boundaries {
        if let BoundaryCondition::Wall { wall_type } = bc {
            match wall_type {
                cfd_core::boundary::WallType::NoSlip => {
                    match name.as_str() {
                        "west" => apply_higher_order_west_wall(matrix, rhs, component, grid, nx, ny)?,
                        "east" => apply_higher_order_east_wall(matrix, rhs, component, grid, nx, ny)?,
                        "north" => apply_higher_order_north_wall(matrix, rhs, component, grid, nx, ny)?,
                        "south" => apply_higher_order_south_wall(matrix, rhs, component, grid, nx, ny)?,
                        _ => {}
                    }
                }
                cfd_core::boundary::WallType::Moving { .. } => {
                    // For moving walls, use standard Dirichlet BC for now
                    // Could be extended to higher-order moving wall BCs
                }
                _ => {}
            }
        }
    }

    Ok(())
}

/// Apply higher-order no-slip boundary condition on west wall
/// Uses quadratic extrapolation: u_0 = (4*u_1 - u_2)/3 for better near-wall gradients
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
        let idx_0 = j * nx;           // Boundary cell (i=0)
        let idx_1 = j * nx + 1;       // First interior cell (i=1)
        let idx_2 = j * nx + 2;       // Second interior cell (i=2)

        // Quadratic extrapolation: u_0 = (4*u_1 - u_2)/3
        // Rearranged: 3*u_0 - 4*u_1 + u_2 = 0
        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, T::one())?;
        rhs[idx_0] = T::zero();
    }

    Ok(())
}

/// Apply higher-order no-slip boundary condition on east wall
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
        let idx_0 = j * nx + nx - 1;  // Boundary cell (i=nx-1)
        let idx_1 = j * nx + nx - 2;  // First interior cell (i=nx-2)
        let idx_2 = j * nx + nx - 3;  // Second interior cell (i=nx-3)

        // Quadratic extrapolation: u_0 = (4*u_1 - u_2)/3
        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, T::one())?;
        rhs[idx_0] = T::zero();
    }

    Ok(())
}

/// Apply higher-order no-slip boundary condition on north wall
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
        let idx_0 = (ny - 1) * nx + i;  // Boundary cell (j=ny-1)
        let idx_1 = (ny - 2) * nx + i;  // First interior cell (j=ny-2)
        let idx_2 = (ny - 3) * nx + i;  // Second interior cell (j=ny-3)

        // Quadratic extrapolation: u_0 = (4*u_1 - u_2)/3
        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, T::one())?;
        rhs[idx_0] = T::zero();
    }

    Ok(())
}

/// Apply higher-order no-slip boundary condition on south wall
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
        let idx_0 = i;                 // Boundary cell (j=0)
        let idx_1 = nx + i;            // First interior cell (j=1)
        let idx_2 = 2 * nx + i;        // Second interior cell (j=2)

        // Quadratic extrapolation: u_0 = (4*u_1 - u_2)/3
        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, T::one())?;
        rhs[idx_0] = T::zero();
    }

    Ok(())
}

fn apply_west_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
    grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    for j in 0..ny {
        let idx = j * nx;

        match bc {
            BoundaryCondition::Dirichlet { value } => {
                // For Dirichlet BC: diagonal is already 1 (set in assemble_system)
                // Just set RHS to boundary value
                rhs[idx] = *value;
            }
            BoundaryCondition::VelocityInlet { velocity } => {
                // For velocity inlet: set velocity to specified value
                let component_idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                rhs[idx] = velocity[component_idx];
            }
            BoundaryCondition::Wall { wall_type } => {
                match wall_type {
                    cfd_core::boundary::WallType::NoSlip => {
                        // For no-slip wall: set velocity to zero
                        rhs[idx] = T::zero();
                    }
                    cfd_core::boundary::WallType::Slip => {
                        // For slip wall: zero normal gradient (handled in assemble_system)
                        // RHS is already set correctly
                    }
                    cfd_core::boundary::WallType::Moving { velocity } => {
                        // For moving wall: set velocity to wall velocity
                        let component_idx = match component {
                            MomentumComponent::U => 0,
                            MomentumComponent::V => 1,
                        };
                        rhs[idx] = velocity[component_idx];
                    }
                    cfd_core::boundary::WallType::Rotating { omega, center } => {
                        rhs[idx] = apply_rotating_wall_bc(component, omega, center, grid, idx);
                    }
                }
            }
            BoundaryCondition::Neumann { gradient } => {
                // Zero gradient condition: u_i - u_{i+1} = 0
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx + 1, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::Periodic { partner: _ } => {
                // Periodic BC: u(0,j) = u(nx-1,j)
                // Implemented as: u[idx] - u[idx + nx - 1] = 0
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + nx - 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Symmetry => {
                // Symmetry BC: zero normal gradient
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::PressureInlet { .. } | BoundaryCondition::PressureOutlet { .. } => {
                // Pressure BC: zero gradient velocity (extrapolate from interior)
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Outflow => {
                // Outflow BC: zero gradient
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::CharacteristicInlet { velocity, .. } => {
                // Characteristic-based inlet: specify incoming waves
                // For incompressible flow: set velocity from specified inlet conditions
                if let Some(vel) = velocity {
                    let component_idx = match component {
                        MomentumComponent::U => 0,
                        MomentumComponent::V => 1,
                    };
                    rhs[idx] = vel[component_idx];
                } else {
                    // Fallback to zero velocity if not specified
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::CharacteristicOutlet { extrapolate_velocity, .. } => {
                // Characteristic-based outlet: extrapolate velocity if requested
                if *extrapolate_velocity {
                    // Extrapolate from interior: u_boundary = u_interior
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx + 1, -T::one())?;
                    rhs[idx] = T::zero();
                } else {
                    // Zero gradient if not extrapolating
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx + 1, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            _ => {}
        }
    }

    Ok(())
}

fn apply_east_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
    grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    for j in 0..ny {
        let idx = j * nx + nx - 1;

        match bc {
            BoundaryCondition::Dirichlet { value } => {
                // For Dirichlet BC: diagonal is already 1 (set in assemble_system)
                // Just set RHS to boundary value
                rhs[idx] = *value;
            }
            BoundaryCondition::VelocityInlet { velocity } => {
                // For velocity inlet: set velocity to specified value
                let component_idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                rhs[idx] = velocity[component_idx];
            }
            BoundaryCondition::Wall { wall_type } => {
                match wall_type {
                    cfd_core::boundary::WallType::NoSlip => {
                        // For no-slip wall: set velocity to zero
                        rhs[idx] = T::zero();
                    }
                    cfd_core::boundary::WallType::Slip => {
                        // For slip wall: zero normal gradient (handled in assemble_system)
                        // RHS is already set correctly
                    }
                    cfd_core::boundary::WallType::Moving { velocity } => {
                        // For moving wall: set velocity to wall velocity
                        let component_idx = match component {
                            MomentumComponent::U => 0,
                            MomentumComponent::V => 1,
                        };
                        rhs[idx] = velocity[component_idx];
                    }
                    cfd_core::boundary::WallType::Rotating { omega, center } => {
                        rhs[idx] = apply_rotating_wall_bc(component, omega, center, grid, idx);
                    }
                }
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx - 1, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::Periodic { partner: _ } => {
                // Periodic BC: u(nx-1,j) = u(0,j)
                // Implemented as: u[idx] - u[j * nx] = 0
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, j * nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Symmetry => {
                // Symmetry BC: zero normal gradient
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::PressureInlet { .. } | BoundaryCondition::PressureOutlet { .. } => {
                // Pressure BC: zero gradient velocity (extrapolate from interior)
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Outflow => {
                // Outflow BC: zero gradient
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            _ => {}
        }
    }

    Ok(())
}

fn apply_north_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
    grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    for i in 0..nx {
        let idx = (ny - 1) * nx + i;

        match bc {
            BoundaryCondition::Dirichlet { value } => {
                // For Dirichlet BC: diagonal is already 1 (set in assemble_system)
                // Just set RHS to boundary value
                rhs[idx] = *value;
            }
            BoundaryCondition::VelocityInlet { velocity } => {
                // For velocity inlet: set velocity to specified value
                let component_idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                rhs[idx] = velocity[component_idx];
            }
            BoundaryCondition::Wall { wall_type } => {
                match wall_type {
                    cfd_core::boundary::WallType::NoSlip => {
                        // For no-slip wall: set velocity to zero
                        rhs[idx] = T::zero();
                    }
                    cfd_core::boundary::WallType::Slip => {
                        // For slip wall: zero normal gradient (handled in assemble_system)
                        // RHS is already set correctly
                    }
                    cfd_core::boundary::WallType::Moving { velocity } => {
                        // For moving wall: set velocity to wall velocity
                        let component_idx = match component {
                            MomentumComponent::U => 0,
                            MomentumComponent::V => 1,
                        };
                        rhs[idx] = velocity[component_idx];
                    }
                    cfd_core::boundary::WallType::Rotating { omega, center } => {
                        rhs[idx] = apply_rotating_wall_bc(component, omega, center, grid, idx);
                    }
                }
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx - nx, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::Periodic { partner: _ } => {
                // Periodic BC: u(i,ny-1) = u(i,0)
                // Implemented as: u[idx] - u[i] = 0
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, i, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Symmetry => {
                // Symmetry BC: zero normal gradient
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::PressureInlet { .. } | BoundaryCondition::PressureOutlet { .. } => {
                // Pressure BC: zero gradient velocity (extrapolate from interior)
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Outflow => {
                // Outflow BC: zero gradient
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            _ => {}
        }
    }

    Ok(())
}

fn apply_south_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
    grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    for i in 0..nx {
        let idx = i;

        match bc {
            BoundaryCondition::Dirichlet { value } => {
                // For Dirichlet BC: diagonal is already 1 (set in assemble_system)
                // Just set RHS to boundary value
                rhs[idx] = *value;
            }
            BoundaryCondition::VelocityInlet { velocity } => {
                // For velocity inlet: set velocity to specified value
                let component_idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                rhs[idx] = velocity[component_idx];
            }
            BoundaryCondition::Wall { wall_type } => {
                match wall_type {
                    cfd_core::boundary::WallType::NoSlip => {
                        // For no-slip wall: set velocity to zero
                        rhs[idx] = T::zero();
                    }
                    cfd_core::boundary::WallType::Slip => {
                        // For slip wall: zero normal gradient (handled in assemble_system)
                        // RHS is already set correctly
                    }
                    cfd_core::boundary::WallType::Moving { velocity } => {
                        // For moving wall: set velocity to wall velocity
                        let component_idx = match component {
                            MomentumComponent::U => 0,
                            MomentumComponent::V => 1,
                        };
                        rhs[idx] = velocity[component_idx];
                    }
                    cfd_core::boundary::WallType::Rotating { omega, center } => {
                        rhs[idx] = apply_rotating_wall_bc(component, omega, center, grid, idx);
                    }
                }
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx + nx, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::Periodic { partner: _ } => {
                // Periodic BC: u(i,0) = u(i,ny-1)
                // Implemented as: u[idx] - u[(ny-1)*nx + i] = 0
                let partner_idx = (ny - 1) * nx + i; // Calculate partner index
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, partner_idx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Symmetry => {
                // Symmetry BC: zero normal gradient
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::PressureInlet { .. } | BoundaryCondition::PressureOutlet { .. } => {
                // Pressure BC: zero gradient velocity (extrapolate from interior)
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Outflow => {
                // Outflow BC: zero gradient
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            _ => {}
        }
    }

    Ok(())
}
