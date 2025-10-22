//! Boundary condition handling for momentum equations

use super::solver::MomentumComponent;
use cfd_core::boundary::BoundaryCondition;
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;
use std::hash::BuildHasher;

/// Apply boundary conditions to momentum equation system
pub fn apply_momentum_boundaries<T, S>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    component: MomentumComponent,
    boundaries: &HashMap<String, BoundaryCondition<T>, S>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()>
where
    T: RealField + Copy + FromPrimitive,
    S: BuildHasher,
{
    // Apply boundary conditions based on location
    for (name, bc) in boundaries {
        match name.as_str() {
            "west" => apply_west_boundary(matrix, rhs, bc, component, nx, ny)?,
            "east" => apply_east_boundary(matrix, rhs, bc, component, nx, ny)?,
            "north" => apply_north_boundary(matrix, rhs, bc, component, nx, ny)?,
            "south" => apply_south_boundary(matrix, rhs, bc, component, nx, ny)?,
            _ => {}
        }
    }

    Ok(())
}

fn apply_west_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    _component: MomentumComponent,
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
            _ => {}
        }
    }

    Ok(())
}

fn apply_east_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    _component: MomentumComponent,
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
    _component: MomentumComponent,
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
    _component: MomentumComponent,
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
