//! Boundary condition handling for momentum equations

use super::solver::MomentumComponent;
use cfd_core::boundary::BoundaryCondition;
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// Apply boundary conditions to momentum equation system
pub fn apply_momentum_boundaries<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    component: MomentumComponent,
    boundaries: &HashMap<String, BoundaryCondition<T>>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
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
                // For Dirichlet BC: use penalty method with moderate multiplier
                // Typical interior diagonal ~1e-3 to 1e-1, so 1e6 is sufficient
                let penalty = T::from_f64(1e6).unwrap_or(T::from_f64(1000.0).unwrap_or(T::one()));
                matrix.add_entry(idx, idx, penalty)?;
                rhs[idx] = *value * penalty;
            }
            BoundaryCondition::Neumann { gradient } => {
                // Zero gradient condition: u_i - u_{i+1} = 0
                if *gradient == T::zero() {
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
    _component: MomentumComponent,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    for j in 0..ny {
        let idx = j * nx + nx - 1;

        match bc {
            BoundaryCondition::Dirichlet { value } => {
                // For Dirichlet BC: use penalty method
                let penalty = T::from_f64(1e6).unwrap_or(T::from_f64(1000.0).unwrap_or(T::one()));
                matrix.add_entry(idx, idx, penalty)?;
                rhs[idx] = *value * penalty;
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx - 1, -T::one())?;
                    rhs[idx] = T::zero();
                }
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
                // For Dirichlet BC: use penalty method
                let penalty = T::from_f64(1e6).unwrap_or(T::from_f64(1000.0).unwrap_or(T::one()));
                matrix.add_entry(idx, idx, penalty)?;
                rhs[idx] = *value * penalty;
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx - nx, -T::one())?;
                    rhs[idx] = T::zero();
                }
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
    _ny: usize,
) -> cfd_core::error::Result<()> {
    for i in 0..nx {
        let idx = i;

        match bc {
            BoundaryCondition::Dirichlet { value } => {
                // For Dirichlet BC: use penalty method
                let penalty = T::from_f64(1e6).unwrap_or(T::from_f64(1000.0).unwrap_or(T::one()));
                matrix.add_entry(idx, idx, penalty)?;
                rhs[idx] = *value * penalty;
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx + nx, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            _ => {}
        }
    }

    Ok(())
}
