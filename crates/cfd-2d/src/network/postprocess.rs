use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::solvers::ns_fvm::NavierStokesSolver2D;

/// Extract maximum and mean wall shear stress from the solved 2D velocity field.
pub(crate) fn extract_field_wall_shear<T>(solver: &NavierStokesSolver2D<T>) -> (T, T)
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    let nx = solver.grid.nx;
    let ny = solver.grid.ny;
    let dx = solver.grid.dx;
    let dy = solver.grid.dy;
    let mut max_tau = T::zero();
    let mut sum_tau = T::zero();
    let mut count = 0u64;

    for i in 1..nx.saturating_sub(1) {
        for j in 1..ny.saturating_sub(1) {
            if !solver.field.mask[(i, j)] {
                continue;
            }
            let next_to_wall = !solver.field.mask[(i - 1, j)]
                || !solver.field.mask[(i + 1, j)]
                || !solver.field.mask[(i, j - 1)]
                || !solver.field.mask[(i, j + 1)]
                || i == 1
                || i == nx - 2
                || j == 1
                || j == ny - 2;
            if !next_to_wall {
                continue;
            }

            let gamma = solver.field.compute_shear_rate(i, j, dx, dy);
            let tau = solver.field.mu[(i, j)] * gamma;
            if tau > max_tau {
                max_tau = tau;
            }
            sum_tau += tau;
            count += 1;
        }
    }

    let mean_tau = if count > 0 {
        sum_tau / T::from_u64(count).expect("analytical constant conversion")
    } else {
        T::zero()
    };

    (max_tau, mean_tau)
}

/// Reconstruct outlet flow from the solved outlet velocity profile and the true channel area.
pub(crate) fn extract_field_outlet_flow_rate<T>(solver: &NavierStokesSolver2D<T>, area: T) -> T
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    if solver.grid.nx == 0 || solver.grid.ny == 0 || area <= T::zero() {
        return T::zero();
    }

    let Some(outlet_cell) = (0..solver.grid.nx)
        .rev()
        .find(|&i| (0..solver.grid.ny).any(|j| solver.field.mask[(i, j)]))
    else {
        return T::zero();
    };
    let outlet_face = (outlet_cell + 1).min(solver.grid.nx);
    let mut integrated_velocity = T::zero();
    let mut open_width = T::zero();

    for j in 0..solver.grid.ny {
        if !solver.field.mask[(outlet_cell, j)] {
            continue;
        }
        let dy = solver.grid.dy_at(j);
        integrated_velocity += solver.field.u[(outlet_face, j)] * dy;
        open_width += dy;
    }

    if open_width <= T::zero() {
        return T::zero();
    }

    let mean_outlet_velocity = integrated_velocity / open_width;
    mean_outlet_velocity * area
}

/// Extract mean inlet and outlet pressures from the solved 2D field.
///
/// The pressures are averaged across the first and last fluid columns to reduce
/// sensitivity to a single cell at a branched or venturi inlet/outlet.
pub(crate) fn extract_field_inlet_outlet_pressure<T>(solver: &NavierStokesSolver2D<T>) -> (T, T)
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    if solver.grid.nx == 0 || solver.grid.ny == 0 {
        return (T::zero(), T::zero());
    }

    let Some(inlet_cell) =
        (0..solver.grid.nx).find(|&i| (0..solver.grid.ny).any(|j| solver.field.mask[(i, j)]))
    else {
        return (T::zero(), T::zero());
    };
    let Some(outlet_cell) = (0..solver.grid.nx)
        .rev()
        .find(|&i| (0..solver.grid.ny).any(|j| solver.field.mask[(i, j)]))
    else {
        return (T::zero(), T::zero());
    };

    let mut inlet_pressure = T::zero();
    let mut inlet_count = 0u64;
    let mut outlet_pressure = T::zero();
    let mut outlet_count = 0u64;

    for j in 0..solver.grid.ny {
        if solver.field.mask[(inlet_cell, j)] {
            inlet_pressure += solver.field.p[(inlet_cell, j)];
            inlet_count += 1;
        }
        if solver.field.mask[(outlet_cell, j)] {
            outlet_pressure += solver.field.p[(outlet_cell, j)];
            outlet_count += 1;
        }
    }

    let inlet_pressure = if inlet_count > 0 {
        inlet_pressure / T::from_u64(inlet_count).expect("analytical constant conversion")
    } else {
        T::zero()
    };
    let outlet_pressure = if outlet_count > 0 {
        outlet_pressure / T::from_u64(outlet_count).expect("analytical constant conversion")
    } else {
        T::zero()
    };

    (inlet_pressure, outlet_pressure)
}
