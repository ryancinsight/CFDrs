//! Scalar-distribution boundary handling for LBM passive-scalar transport.
//!
//! The passive scalar uses a conservative boundary policy aligned with the flow
//! solver's domain classification:
//! - wall boundaries use exact D2Q9 bounce-back for no-flux conditions
//! - periodic boundaries wrap to the opposite edge
//! - all other boundaries use a zero-gradient copy from the nearest interior cell

use crate::solvers::lbm::lattice::D2Q9;
use crate::solvers::lbm::streaming::f_idx;
use cfd_core::physics::boundary::BoundaryCondition;
use nalgebra::RealField;
use std::collections::HashMap;

#[inline]
pub(crate) fn apply_scalar_boundaries<T: RealField + Copy>(
    g: &mut [T],
    boundaries: &HashMap<(usize, usize), BoundaryCondition<T>>,
    nx: usize,
    ny: usize,
) {
    for ((i, j), bc) in boundaries {
        match bc {
            BoundaryCondition::Wall { .. } => {
                apply_scalar_bounce_back(g, *i, *j, nx);
            }
            BoundaryCondition::Periodic { .. } => {
                let (src_i, src_j) = periodic_scalar_source(*i, *j, nx, ny);
                copy_scalar_populations(g, nx, src_i, src_j, *i, *j);
            }
            _ => {
                let (src_i, src_j) = zero_gradient_scalar_source(*i, *j, nx, ny);
                copy_scalar_populations(g, nx, src_i, src_j, *i, *j);
            }
        }
    }
}

#[inline]
fn apply_scalar_bounce_back<T: RealField + Copy>(g: &mut [T], i: usize, j: usize, nx: usize) {
    let mut values = [T::zero(); 9];
    for q in 0..9 {
        values[q] = g[f_idx(j, i, q, nx)];
    }
    for q in 0..9 {
        let q_opp = D2Q9::OPPOSITE[q];
        g[f_idx(j, i, q, nx)] = values[q_opp];
    }
}

#[inline]
fn copy_scalar_populations<T: RealField + Copy>(
    g: &mut [T],
    nx: usize,
    src_i: usize,
    src_j: usize,
    dst_i: usize,
    dst_j: usize,
) {
    let mut values = [T::zero(); 9];
    for q in 0..9 {
        values[q] = g[f_idx(src_j, src_i, q, nx)];
    }
    for q in 0..9 {
        g[f_idx(dst_j, dst_i, q, nx)] = values[q];
    }
}

#[inline]
fn zero_gradient_scalar_source(i: usize, j: usize, nx: usize, ny: usize) -> (usize, usize) {
    let src_i = if i == 0 {
        usize::from(nx > 1)
    } else if i + 1 == nx {
        nx.saturating_sub(2)
    } else {
        i
    };

    let src_j = if j == 0 {
        usize::from(ny > 1)
    } else if j + 1 == ny {
        ny.saturating_sub(2)
    } else {
        j
    };

    (src_i, src_j)
}

#[inline]
fn periodic_scalar_source(i: usize, j: usize, nx: usize, ny: usize) -> (usize, usize) {
    let src_i = if i == 0 {
        nx.saturating_sub(1)
    } else if i + 1 == nx {
        0
    } else {
        i
    };

    let src_j = if j == 0 {
        ny.saturating_sub(1)
    } else if j + 1 == ny {
        0
    } else {
        j
    };

    (src_i, src_j)
}

#[cfg(test)]
mod tests {
    use super::apply_scalar_boundaries;
    use crate::solvers::lbm::lattice::D2Q9;
    use crate::solvers::lbm::streaming::f_idx;
    use approx::assert_relative_eq;
    use cfd_core::physics::boundary::BoundaryCondition;
    use std::collections::HashMap;

    fn write_cell(g: &mut [f64], i: usize, j: usize, nx: usize, values: [f64; 9]) {
        for q in 0..9 {
            g[f_idx(j, i, q, nx)] = values[q];
        }
    }

    #[test]
    fn wall_boundary_bounce_back_swaps_antipodes() {
        let nx = 3_usize;
        let ny = 3_usize;
        let mut g = vec![0.0_f64; nx * ny * 9];
        let mut boundaries = HashMap::new();
        boundaries.insert((0, 1), BoundaryCondition::wall_no_slip());

        let values = [0.31, 0.12, 0.23, 0.47, 0.19, 0.15, 0.17, 0.21, 0.29];
        write_cell(&mut g, 0, 1, nx, values);

        apply_scalar_boundaries(&mut g, &boundaries, nx, ny);

        for q in 0..9 {
            let q_opp = D2Q9::OPPOSITE[q];
            assert_relative_eq!(g[f_idx(1, 0, q, nx)], values[q_opp], epsilon = 1e-15);
        }
    }

    #[test]
    fn non_wall_boundary_copies_adjacent_interior_cell() {
        let nx = 3_usize;
        let ny = 3_usize;
        let mut g = vec![0.0_f64; nx * ny * 9];
        let mut boundaries = HashMap::new();
        boundaries.insert((2, 1), BoundaryCondition::Outflow);

        let source = [0.41, 0.32, 0.27, 0.18, 0.55, 0.66, 0.73, 0.84, 0.95];
        let target = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09];
        write_cell(&mut g, 1, 1, nx, source);
        write_cell(&mut g, 2, 1, nx, target);

        apply_scalar_boundaries(&mut g, &boundaries, nx, ny);

        for q in 0..9 {
            assert_relative_eq!(g[f_idx(1, 2, q, nx)], source[q], epsilon = 1e-15);
        }
    }
}
