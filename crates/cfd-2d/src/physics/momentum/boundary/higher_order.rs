use super::super::solver::MomentumComponent;
use super::MatrixUpdater;
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use cfd_core::physics::boundary::BoundaryCondition;
use eunomia::FloatElement;
use leto::Array1;
use std::collections::HashMap;
use std::hash::BuildHasher;

/// Apply higher-order boundary conditions for improved near-wall accuracy.
/// Implements quadratic extrapolation for better velocity gradients near walls.
pub fn apply_higher_order_wall_boundaries<T, S, M>(
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

    if matches!(component, MomentumComponent::V) {
        if let Some(BoundaryCondition::Wall {
            wall_type: cfd_core::physics::boundary::WallType::NoSlip,
        }) = boundaries.get("west")
        {
            apply_higher_order_west_wall(matrix, rhs, nx, ny)?;
        }
    }
    if matches!(component, MomentumComponent::V) {
        if let Some(BoundaryCondition::Wall {
            wall_type: cfd_core::physics::boundary::WallType::NoSlip,
        }) = boundaries.get("east")
        {
            apply_higher_order_east_wall(matrix, rhs, nx, ny)?;
        }
    }
    if matches!(component, MomentumComponent::U) {
        if let Some(BoundaryCondition::Wall {
            wall_type: cfd_core::physics::boundary::WallType::NoSlip,
        }) = boundaries.get("north")
        {
            apply_higher_order_north_wall(matrix, rhs, nx, ny)?;
        }
    }
    if matches!(component, MomentumComponent::U) {
        if let Some(BoundaryCondition::Wall {
            wall_type: cfd_core::physics::boundary::WallType::NoSlip,
        }) = boundaries.get("south")
        {
            apply_higher_order_south_wall(matrix, rhs, nx)?;
        }
    }

    Ok(())
}

/// Quadratic extrapolation: u_0 = (4*u_1 - u_2)/3 for west wall
fn apply_higher_order_west_wall<T: Cfd2dScalar + Copy + FloatElement, M: MatrixUpdater<T>>(
    matrix: &mut M,
    rhs: &mut Array1<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    let four: T = scalar::from_f64(4.0);
    let three: T = scalar::from_f64(3.0);

    for j in 1..ny.saturating_sub(1) {
        let idx_0 = j * nx;
        let idx_1 = j * nx + 1;
        let idx_2 = j * nx + 2;

        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, scalar::one())?;
        rhs[idx_0] = scalar::zero();
    }

    Ok(())
}

fn apply_higher_order_east_wall<T: Cfd2dScalar + Copy + FloatElement, M: MatrixUpdater<T>>(
    matrix: &mut M,
    rhs: &mut Array1<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    let four: T = scalar::from_f64(4.0);
    let three: T = scalar::from_f64(3.0);

    for j in 1..ny.saturating_sub(1) {
        let idx_0 = j * nx + nx - 1;
        let idx_1 = j * nx + nx - 2;
        let idx_2 = j * nx + nx - 3;

        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, scalar::one())?;
        rhs[idx_0] = scalar::zero();
    }

    Ok(())
}

fn apply_higher_order_north_wall<T: Cfd2dScalar + Copy + FloatElement, M: MatrixUpdater<T>>(
    matrix: &mut M,
    rhs: &mut Array1<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    let four: T = scalar::from_f64(4.0);
    let three: T = scalar::from_f64(3.0);

    for i in 0..nx {
        let idx_0 = (ny - 1) * nx + i;
        let idx_1 = (ny - 2) * nx + i;
        let idx_2 = (ny - 3) * nx + i;

        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, scalar::one())?;
        rhs[idx_0] = scalar::zero();
    }

    Ok(())
}

fn apply_higher_order_south_wall<T: Cfd2dScalar + Copy + FloatElement, M: MatrixUpdater<T>>(
    matrix: &mut M,
    rhs: &mut Array1<T>,
    nx: usize,
) -> cfd_core::error::Result<()> {
    let four: T = scalar::from_f64(4.0);
    let three: T = scalar::from_f64(3.0);

    for i in 0..nx {
        let idx_0 = i;
        let idx_1 = nx + i;
        let idx_2 = 2 * nx + i;

        matrix.add_entry(idx_0, idx_0, three)?;
        matrix.add_entry(idx_0, idx_1, -four)?;
        matrix.add_entry(idx_0, idx_2, scalar::one())?;
        rhs[idx_0] = scalar::zero();
    }

    Ok(())
}
