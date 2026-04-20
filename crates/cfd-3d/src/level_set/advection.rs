//! Level-set advection helpers.
//!
//! The production path uses WENO5-Z spatial reconstruction with an SSPRK3
//! (TVD-RK3) time integrator when the grid admits the required 7-point stencil.
//! Smaller domains fall back to first-order upwind so the solver remains
//! well-defined on the boundary-test meshes used throughout the crate.

use super::{config::LevelSetConfig, weno::weno5_derivative};
use cfd_core::error::Result;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

#[inline]
fn linear_index(nx: usize, ny: usize, i: usize, j: usize, k: usize) -> usize {
    k * ny * nx + j * nx + i
}

/// Advance the level set by one time step using the configured advection scheme.
#[allow(clippy::too_many_arguments)]
pub(super) fn advance<T>(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: T,
    dy: T,
    dz: T,
    config: &LevelSetConfig,
    dt: T,
    velocity: &[Vector3<T>],
    phi_previous: &[T],
    phi: &mut [T],
    phi_reinit: &mut [T],
) -> Result<()>
where
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
{
    let supports_weno = nx >= 7 && ny >= 7 && nz >= 7;
    if config.use_weno && supports_weno {
        advance_weno5_rk3(
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            dt,
            velocity,
            phi_previous,
            phi,
            phi_reinit,
        );
    } else {
        // First-order upwind is the mathematically admissible fallback when the
        // 7-point WENO stencil cannot be formed on at least one axis.
        advance_first_order_euler(nx, ny, nz, dx, dy, dz, dt, velocity, phi_previous, phi);
    }

    Ok(())
}

/// Copy the halo cells from one buffer into another.
#[allow(clippy::too_many_arguments)]
pub(super) fn copy_halo_between_buffers<T>(
    nx: usize,
    ny: usize,
    nz: usize,
    source: &[T],
    dest: &mut [T],
    halo: usize,
) where
    T: Copy,
{
    if halo == 0 || nx == 0 || ny == 0 || nz == 0 {
        return;
    }

    debug_assert_eq!(source.len(), dest.len());

    let halo_x = halo.min(nx);
    let halo_y = halo.min(ny);
    let halo_z = halo.min(nz);
    if halo_x == nx || halo_y == ny || halo_z == nz {
        dest.copy_from_slice(source);
        return;
    }

    let x_inner_start = halo_x;
    let x_inner_end = nx.saturating_sub(halo_x);
    let y_inner_start = halo_y;
    let y_inner_end = ny.saturating_sub(halo_y);
    let z_inner_start = halo_z;
    let z_inner_end = nz.saturating_sub(halo_z);
    let row_len = nx;
    let slab_len = nx * ny;

    for k in 0..z_inner_start {
        let base = k * slab_len;
        dest[base..base + slab_len].copy_from_slice(&source[base..base + slab_len]);
    }

    for k in z_inner_end..nz {
        let base = k * slab_len;
        dest[base..base + slab_len].copy_from_slice(&source[base..base + slab_len]);
    }

    for k in z_inner_start..z_inner_end {
        let slab_base = k * slab_len;

        for j in 0..y_inner_start {
            let row_base = slab_base + j * row_len;
            dest[row_base..row_base + row_len]
                .copy_from_slice(&source[row_base..row_base + row_len]);
        }

        for j in y_inner_end..ny {
            let row_base = slab_base + j * row_len;
            dest[row_base..row_base + row_len]
                .copy_from_slice(&source[row_base..row_base + row_len]);
        }

        for j in y_inner_start..y_inner_end {
            let row_base = slab_base + j * row_len;
            if halo_x > 0 {
                dest[row_base..row_base + x_inner_start]
                    .copy_from_slice(&source[row_base..row_base + x_inner_start]);
            }
            if x_inner_end < row_len {
                dest[row_base + x_inner_end..row_base + row_len]
                    .copy_from_slice(&source[row_base + x_inner_end..row_base + row_len]);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn advance_weno5_rk3<T>(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: T,
    dy: T,
    dz: T,
    dt: T,
    velocity: &[Vector3<T>],
    phi_previous: &[T],
    phi: &mut [T],
    phi_reinit: &mut [T],
) where
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
{
    let zero = T::zero();
    let one = T::one();
    let quarter = <T as FromPrimitive>::from_f64(0.25).expect("0.25 is representable in IEEE 754");
    let three_quarters =
        <T as FromPrimitive>::from_f64(0.75).expect("0.75 is representable in IEEE 754");
    let one_third =
        <T as FromPrimitive>::from_f64(1.0 / 3.0).expect("1/3 is representable in IEEE 754");
    let two_thirds =
        <T as FromPrimitive>::from_f64(2.0 / 3.0).expect("2/3 is representable in IEEE 754");

    rk3_stage_weno(
        nx,
        ny,
        nz,
        dx,
        dy,
        dz,
        dt,
        velocity,
        phi_previous,
        phi_previous,
        phi,
        zero,
        one,
    );
    copy_halo_between_buffers(nx, ny, nz, phi_previous, phi, 3);

    rk3_stage_weno(
        nx,
        ny,
        nz,
        dx,
        dy,
        dz,
        dt,
        velocity,
        phi_previous,
        phi,
        phi_reinit,
        three_quarters,
        quarter,
    );
    copy_halo_between_buffers(nx, ny, nz, phi_previous, phi_reinit, 3);

    rk3_stage_weno(
        nx,
        ny,
        nz,
        dx,
        dy,
        dz,
        dt,
        velocity,
        phi_previous,
        phi_reinit,
        phi,
        one_third,
        two_thirds,
    );
    copy_halo_between_buffers(nx, ny, nz, phi_previous, phi, 3);
}

#[allow(clippy::too_many_arguments)]
fn rk3_stage_weno<T>(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: T,
    dy: T,
    dz: T,
    dt: T,
    velocity: &[Vector3<T>],
    source: &[T],
    stage_input: &[T],
    dest: &mut [T],
    source_weight: T,
    stage_weight: T,
) where
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
{
    for k in 3..nz.saturating_sub(3) {
        for j in 3..ny.saturating_sub(3) {
            for i in 3..nx.saturating_sub(3) {
                let idx = linear_index(nx, ny, i, j, k);
                let vel = velocity[idx];
                let transport = weno_transport_term(nx, ny, dx, dy, dz, stage_input, vel, i, j, k);
                dest[idx] = source_weight * source[idx]
                    + stage_weight * (stage_input[idx] - dt * transport);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn advance_first_order_euler<T>(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: T,
    dy: T,
    dz: T,
    dt: T,
    velocity: &[Vector3<T>],
    phi_previous: &[T],
    phi: &mut [T],
) where
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
{
    for k in 1..nz.saturating_sub(1) {
        for j in 1..ny.saturating_sub(1) {
            for i in 1..nx.saturating_sub(1) {
                let idx = linear_index(nx, ny, i, j, k);
                let vel = velocity[idx];
                let transport =
                    first_order_transport_term(nx, ny, dx, dy, dz, phi_previous, vel, i, j, k);
                phi[idx] = phi_previous[idx] - dt * transport;
            }
        }
    }

    copy_halo_between_buffers(nx, ny, nz, phi_previous, phi, 1);
}

#[allow(clippy::too_many_arguments)]
fn weno_transport_term<T>(
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    dz: T,
    state: &[T],
    vel: Vector3<T>,
    i: usize,
    j: usize,
    k: usize,
) -> T
where
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
{
    let idx = linear_index(nx, ny, i, j, k);

    let vx = [
        state[linear_index(nx, ny, i - 3, j, k)],
        state[linear_index(nx, ny, i - 2, j, k)],
        state[linear_index(nx, ny, i - 1, j, k)],
        state[idx],
        state[linear_index(nx, ny, i + 1, j, k)],
        state[linear_index(nx, ny, i + 2, j, k)],
        state[linear_index(nx, ny, i + 3, j, k)],
    ];
    let vy = [
        state[linear_index(nx, ny, i, j - 3, k)],
        state[linear_index(nx, ny, i, j - 2, k)],
        state[linear_index(nx, ny, i, j - 1, k)],
        state[idx],
        state[linear_index(nx, ny, i, j + 1, k)],
        state[linear_index(nx, ny, i, j + 2, k)],
        state[linear_index(nx, ny, i, j + 3, k)],
    ];
    let vz = [
        state[linear_index(nx, ny, i, j, k - 3)],
        state[linear_index(nx, ny, i, j, k - 2)],
        state[linear_index(nx, ny, i, j, k - 1)],
        state[idx],
        state[linear_index(nx, ny, i, j, k + 1)],
        state[linear_index(nx, ny, i, j, k + 2)],
        state[linear_index(nx, ny, i, j, k + 3)],
    ];

    let dphi_dx = weno5_derivative(vx, dx, vel.x);
    let dphi_dy = weno5_derivative(vy, dy, vel.y);
    let dphi_dz = weno5_derivative(vz, dz, vel.z);
    vel.x * dphi_dx + vel.y * dphi_dy + vel.z * dphi_dz
}

#[allow(clippy::too_many_arguments)]
fn first_order_transport_term<T>(
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    dz: T,
    state: &[T],
    vel: Vector3<T>,
    i: usize,
    j: usize,
    k: usize,
) -> T
where
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
{
    let idx = linear_index(nx, ny, i, j, k);

    let dphi_dx = if vel.x > T::zero() {
        (state[idx] - state[linear_index(nx, ny, i - 1, j, k)]) / dx
    } else if vel.x < T::zero() {
        (state[linear_index(nx, ny, i + 1, j, k)] - state[idx]) / dx
    } else {
        T::zero()
    };

    let dphi_dy = if vel.y > T::zero() {
        (state[idx] - state[linear_index(nx, ny, i, j - 1, k)]) / dy
    } else if vel.y < T::zero() {
        (state[linear_index(nx, ny, i, j + 1, k)] - state[idx]) / dy
    } else {
        T::zero()
    };

    let dphi_dz = if vel.z > T::zero() {
        (state[idx] - state[linear_index(nx, ny, i, j, k - 1)]) / dz
    } else if vel.z < T::zero() {
        (state[linear_index(nx, ny, i, j, k + 1)] - state[idx]) / dz
    } else {
        T::zero()
    };

    vel.x * dphi_dx + vel.y * dphi_dy + vel.z * dphi_dz
}
