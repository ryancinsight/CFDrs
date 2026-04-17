//! Shared continuity-residual utilities for pressure-velocity solvers.
//!
//! # Theorem - Discrete Continuity Residual SSOT
//!
//! For incompressible flow on a Cartesian grid, the solver residual is the
//! maximum magnitude of a discrete divergence operator. Different coupling
//! schemes only change the sampling lattice:
//! - cell-centred fields use a central-difference divergence
//! - Rhie-Chow face fields use a face-flux divergence
//! - predictor/corrector assembly uses a forward-difference divergence
//!
//! **Proof sketch**: each routine evaluates the same conservation law with a
//! different discrete stencil. The maximum absolute divergence over all active
//! cells is therefore the common convergence measure, independent of the
//! solver family that produced the velocity field.

use nalgebra::RealField;
use num_traits::FromPrimitive;

#[inline]
pub(crate) fn pointwise_forward_continuity_residual<T, U, V>(
    i: usize,
    j: usize,
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    u_at: &mut U,
    v_at: &mut V,
) -> T
where
    T: RealField + Copy + FromPrimitive,
    U: FnMut(usize, usize) -> T,
    V: FnMut(usize, usize) -> T,
{
    if nx == 0 || ny == 0 {
        return T::zero();
    }

    let du_dx = if i + 1 < nx {
        (u_at(i + 1, j) - u_at(i, j)) / dx
    } else {
        T::zero()
    };
    let dv_dy = if j + 1 < ny {
        (v_at(i, j + 1) - v_at(i, j)) / dy
    } else {
        T::zero()
    };

    du_dx + dv_dy
}

#[inline]
pub(crate) fn max_forward_continuity_residual<T, U, V>(
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    mut u_at: U,
    mut v_at: V,
) -> T
where
    T: RealField + Copy + FromPrimitive,
    U: FnMut(usize, usize) -> T,
    V: FnMut(usize, usize) -> T,
{
    let mut max_divergence = T::zero();

    for i in 0..nx {
        for j in 0..ny {
            let residual = pointwise_forward_continuity_residual(
                i,
                j,
                nx,
                ny,
                dx,
                dy,
                &mut u_at,
                &mut v_at,
            );
            let abs_residual = residual.abs();
            if abs_residual > max_divergence {
                max_divergence = abs_residual;
            }
        }
    }

    max_divergence
}

#[inline]
pub(crate) fn max_central_continuity_residual<T, U, V>(
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    mut u_at: U,
    mut v_at: V,
) -> T
where
    T: RealField + Copy + FromPrimitive,
    U: FnMut(usize, usize) -> T,
    V: FnMut(usize, usize) -> T,
{
    if nx < 3 || ny < 3 {
        return T::zero();
    }

    let mut max_divergence = T::zero();
    let two = T::from_f64(2.0).expect("Exact mathematically representable f64");

    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let du_dx = (u_at(i + 1, j) - u_at(i - 1, j)) / (two * dx);
            let dv_dy = (v_at(i, j + 1) - v_at(i, j - 1)) / (two * dy);
            let abs_divergence = (du_dx + dv_dy).abs();
            if abs_divergence > max_divergence {
                max_divergence = abs_divergence;
            }
        }
    }

    max_divergence
}

#[inline]
pub(crate) fn max_face_continuity_residual<T, U, V>(
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    mut u_face_at: U,
    mut v_face_at: V,
) -> T
where
    T: RealField + Copy + FromPrimitive,
    U: FnMut(usize, usize) -> T,
    V: FnMut(usize, usize) -> T,
{
    if nx < 3 || ny < 3 {
        return T::zero();
    }

    let mut max_divergence = T::zero();

    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let du_dx = (u_face_at(i, j) - u_face_at(i - 1, j)) / dx;
            let dv_dy = (v_face_at(i, j) - v_face_at(i, j - 1)) / dy;
            let abs_divergence = (du_dx + dv_dy).abs();
            if abs_divergence > max_divergence {
                max_divergence = abs_divergence;
            }
        }
    }

    max_divergence
}

#[cfg(test)]
mod tests {
    use super::{
        max_central_continuity_residual, max_face_continuity_residual,
        max_forward_continuity_residual, pointwise_forward_continuity_residual,
    };
    use approx::assert_relative_eq;

    #[test]
    fn continuity_helpers_reproduce_affine_divergence() {
        let nx = 5_usize;
        let ny = 5_usize;

        let pointwise = pointwise_forward_continuity_residual(
            1,
            1,
            nx,
            ny,
            1.0_f64,
            1.0_f64,
            &mut |i, j| 2.0 * i as f64 + 3.0 * j as f64,
            &mut |i, j| -1.0 * i as f64 + 4.0 * j as f64,
        );
        let forward = max_forward_continuity_residual(
            nx,
            ny,
            1.0_f64,
            1.0_f64,
            |i, j| 2.0 * i as f64 + 3.0 * j as f64,
            |i, j| -1.0 * i as f64 + 4.0 * j as f64,
        );
        let central = max_central_continuity_residual(
            nx,
            ny,
            1.0_f64,
            1.0_f64,
            |i, j| 2.0 * i as f64 + 3.0 * j as f64,
            |i, j| -1.0 * i as f64 + 4.0 * j as f64,
        );
        let face = max_face_continuity_residual(
            nx,
            ny,
            1.0_f64,
            1.0_f64,
            |i, j| 2.0 * i as f64 + 3.0 * j as f64,
            |i, j| -1.0 * i as f64 + 4.0 * j as f64,
        );

        assert_relative_eq!(pointwise, 6.0, epsilon = 1e-12);
        assert_relative_eq!(forward, 6.0, epsilon = 1e-12);
        assert_relative_eq!(central, 6.0, epsilon = 1e-12);
        assert_relative_eq!(face, 6.0, epsilon = 1e-12);
    }
}