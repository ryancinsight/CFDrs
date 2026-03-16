use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::solvers::ns_fvm::NavierStokesSolver2D;
use crate::solvers::venturi_flow::VenturiGeometry;

/// Populate the solver's fluid mask from a venturi geometry.
pub(crate) fn populate_venturi_mask<T>(
    solver: &mut NavierStokesSolver2D<T>,
    geom: &VenturiGeometry<T>,
    nx: usize,
    ny: usize,
) where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    let half_h = geom.w_inlet / T::from_f64(2.0).unwrap_or_else(T::one);
    for i in 0..nx {
        for j in 0..ny {
            let x = solver.grid.x_center(i);
            let y = solver.grid.y_center(j) - half_h;
            solver.field.mask[(i, j)] = geom.contains(x, y);
        }
    }
}

/// Populate the solver's fluid mask for a circular hydraulic aperture.
pub(crate) fn populate_circular_mask<T>(
    solver: &mut NavierStokesSolver2D<T>,
    diameter_m: f64,
    nx: usize,
    ny: usize,
) where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    let radius = T::from_f64(diameter_m / 2.0).unwrap_or_else(T::one);
    let center_y = radius;
    for i in 0..nx {
        for j in 0..ny {
            let y = solver.grid.y_center(j);
            solver.field.mask[(i, j)] = Float::abs(y - center_y) <= radius;
        }
    }
}

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
        sum_tau / T::from_u64(count).unwrap_or_else(T::one)
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

    let outlet_face = solver.grid.nx;
    let outlet_cell = solver.grid.nx - 1;
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
