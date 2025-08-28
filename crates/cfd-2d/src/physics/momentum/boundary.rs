//! Boundary condition handling for momentum equations

use super::solver::MomentumComponent;
use cfd_core::boundary::BoundaryCondition;
use cfd_math::SparseMatrixBuilder;
use nalgebra::{DVector, RealField};
use std::collections::HashMap;

/// Apply boundary conditions to momentum equation system
pub fn apply_momentum_boundaries<T: RealField + Copy>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut DVector<T>,
    component: MomentumComponent,
    boundaries: &HashMap<String, BoundaryCondition<T>>,
    nx: usize,
    ny: usize,
) -> cfd_core::Result<()> {
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

fn apply_west_boundary<T: RealField + Copy>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut DVector<T>,
    bc: &BoundaryCondition<T>,
    _component: MomentumComponent,
    nx: usize,
    ny: usize,
) -> cfd_core::Result<()> {
    for j in 0..ny {
        let idx = j * nx;

        match bc {
            BoundaryCondition::Dirichlet(value) => {
                // Set velocity to specified value
                matrix.set_row(idx, idx, T::one());
                rhs[idx] = *value;
            }
            BoundaryCondition::Neumann(gradient) => {
                // Zero gradient condition
                if *gradient == T::zero() {
                    matrix.add(idx, idx + 1, T::one());
                    matrix.add(idx, idx, -T::one());
                }
            }
            _ => {}
        }
    }

    Ok(())
}

fn apply_east_boundary<T: RealField + Copy>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut DVector<T>,
    bc: &BoundaryCondition<T>,
    _component: MomentumComponent,
    nx: usize,
    ny: usize,
) -> cfd_core::Result<()> {
    for j in 0..ny {
        let idx = j * nx + nx - 1;

        match bc {
            BoundaryCondition::Dirichlet(value) => {
                matrix.set_row(idx, idx, T::one());
                rhs[idx] = *value;
            }
            BoundaryCondition::Neumann(gradient) => {
                if *gradient == T::zero() {
                    matrix.add(idx, idx - 1, T::one());
                    matrix.add(idx, idx, -T::one());
                }
            }
            _ => {}
        }
    }

    Ok(())
}

fn apply_north_boundary<T: RealField + Copy>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut DVector<T>,
    bc: &BoundaryCondition<T>,
    _component: MomentumComponent,
    nx: usize,
    ny: usize,
) -> cfd_core::Result<()> {
    for i in 0..nx {
        let idx = (ny - 1) * nx + i;

        match bc {
            BoundaryCondition::Dirichlet(value) => {
                matrix.set_row(idx, idx, T::one());
                rhs[idx] = *value;
            }
            BoundaryCondition::Neumann(gradient) => {
                if *gradient == T::zero() {
                    matrix.add(idx, idx - nx, T::one());
                    matrix.add(idx, idx, -T::one());
                }
            }
            _ => {}
        }
    }

    Ok(())
}

fn apply_south_boundary<T: RealField + Copy>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut DVector<T>,
    bc: &BoundaryCondition<T>,
    _component: MomentumComponent,
    nx: usize,
    _ny: usize,
) -> cfd_core::Result<()> {
    for i in 0..nx {
        let idx = i;

        match bc {
            BoundaryCondition::Dirichlet(value) => {
                matrix.set_row(idx, idx, T::one());
                rhs[idx] = *value;
            }
            BoundaryCondition::Neumann(gradient) => {
                if *gradient == T::zero() {
                    matrix.add(idx, idx + nx, T::one());
                    matrix.add(idx, idx, -T::one());
                }
            }
            _ => {}
        }
    }

    Ok(())
}
