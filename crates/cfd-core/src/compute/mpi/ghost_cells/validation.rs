//! Ghost cell data validation and checksum verification.

use super::sync_exchange::GhostCellManager;
use crate::compute::mpi::decomposition::NeighborDirection;
use crate::compute::mpi::error::{MpiError, MpiResult};
use crate::compute::mpi::LocalSubdomain;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> GhostCellManager<T> {
    /// Validate ghost cell data consistency
    pub fn validate_ghost_cells(
        &self,
        velocity_u: &[Vec<Vector2<T>>],
        velocity_v: &[Vec<Vector2<T>>],
        pressure: &[Vec<T>],
        subdomain: &LocalSubdomain,
    ) -> MpiResult<bool> {
        self.ensure_dimensions(velocity_u, velocity_v, pressure, subdomain)?;
        let mut local_valid = true;

        for (neighbor_rank, neighbor_info) in &self.neighbors {
            match neighbor_info.direction {
                NeighborDirection::Left => {
                    let owned = self.checksum_column(
                        velocity_u,
                        velocity_v,
                        pressure,
                        subdomain.ghost_layers,
                        subdomain.total_ny(),
                    )?;
                    let ghost = self.checksum_column(
                        velocity_u,
                        velocity_v,
                        pressure,
                        subdomain.ghost_layers.saturating_sub(1),
                        subdomain.total_ny(),
                    )?;
                    let neighbor_checksum =
                        self.exchange_checksum(*neighbor_rank, true, owned, 900)?;
                    if !self.checksum_matches(ghost, neighbor_checksum) {
                        local_valid = false;
                    }
                }
                NeighborDirection::Right => {
                    let owned = self.checksum_column(
                        velocity_u,
                        velocity_v,
                        pressure,
                        subdomain.ghost_layers + subdomain.nx_local - 1,
                        subdomain.total_ny(),
                    )?;
                    let ghost = self.checksum_column(
                        velocity_u,
                        velocity_v,
                        pressure,
                        subdomain.ghost_layers + subdomain.nx_local,
                        subdomain.total_ny(),
                    )?;
                    let neighbor_checksum =
                        self.exchange_checksum(*neighbor_rank, false, owned, 901)?;
                    if !self.checksum_matches(ghost, neighbor_checksum) {
                        local_valid = false;
                    }
                }
                NeighborDirection::Bottom => {
                    let owned = self.checksum_row(
                        velocity_u,
                        velocity_v,
                        pressure,
                        subdomain.ghost_layers,
                        subdomain.total_nx(),
                    )?;
                    let ghost = self.checksum_row(
                        velocity_u,
                        velocity_v,
                        pressure,
                        subdomain.ghost_layers.saturating_sub(1),
                        subdomain.total_nx(),
                    )?;
                    let neighbor_checksum =
                        self.exchange_checksum(*neighbor_rank, true, owned, 902)?;
                    if !self.checksum_matches(ghost, neighbor_checksum) {
                        local_valid = false;
                    }
                }
                NeighborDirection::Top => {
                    let owned = self.checksum_row(
                        velocity_u,
                        velocity_v,
                        pressure,
                        subdomain.ghost_layers + subdomain.ny_local - 1,
                        subdomain.total_nx(),
                    )?;
                    let ghost = self.checksum_row(
                        velocity_u,
                        velocity_v,
                        pressure,
                        subdomain.ghost_layers + subdomain.ny_local,
                        subdomain.total_nx(),
                    )?;
                    let neighbor_checksum =
                        self.exchange_checksum(*neighbor_rank, false, owned, 903)?;
                    if !self.checksum_matches(ghost, neighbor_checksum) {
                        local_valid = false;
                    }
                }
                NeighborDirection::Front | NeighborDirection::Back => {
                    return Err(MpiError::NotAvailable(
                        "3D ghost-cell validation requires 3D field storage".to_string(),
                    ));
                }
            }
        }

        let mut reduced = if local_valid { 1i32 } else { 0i32 };
        self.communicator.all_reduce_min(&mut reduced);
        Ok(reduced == 1)
    }

    fn ensure_dimensions(
        &self,
        velocity_u: &[Vec<Vector2<T>>],
        velocity_v: &[Vec<Vector2<T>>],
        pressure: &[Vec<T>],
        subdomain: &LocalSubdomain,
    ) -> MpiResult<()> {
        let nx = subdomain.total_nx();
        let ny = subdomain.total_ny();
        if velocity_u.len() != nx || velocity_v.len() != nx || pressure.len() != nx {
            return Err(MpiError::GhostCellError(
                "Field dimensions do not match subdomain extents".to_string(),
            ));
        }
        for i in 0..nx {
            if velocity_u[i].len() != ny
                || velocity_v[i].len() != ny
                || pressure[i].len() != ny
            {
                return Err(MpiError::GhostCellError(
                    "Field dimensions do not match subdomain extents".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn checksum_column(
        &self,
        velocity_u: &[Vec<Vector2<T>>],
        velocity_v: &[Vec<Vector2<T>>],
        pressure: &[Vec<T>],
        i: usize,
        ny: usize,
    ) -> MpiResult<T> {
        if i >= velocity_u.len() || i >= velocity_v.len() || i >= pressure.len() {
            return Err(MpiError::GhostCellError(
                "Column index out of bounds".to_string(),
            ));
        }
        let mut sum = T::zero();
        for j in 0..ny {
            let u = velocity_u[i][j];
            let v = velocity_v[i][j];
            sum += u.x.abs() + u.y.abs() + v.x.abs() + v.y.abs() + pressure[i][j].abs();
        }
        Ok(sum)
    }

    fn checksum_row(
        &self,
        velocity_u: &[Vec<Vector2<T>>],
        velocity_v: &[Vec<Vector2<T>>],
        pressure: &[Vec<T>],
        j: usize,
        nx: usize,
    ) -> MpiResult<T> {
        let mut sum = T::zero();
        for i in 0..nx {
            let u = velocity_u[i][j];
            let v = velocity_v[i][j];
            sum += u.x.abs() + u.y.abs() + v.x.abs() + v.y.abs() + pressure[i][j].abs();
        }
        Ok(sum)
    }

    fn exchange_checksum(
        &self,
        neighbor_rank: i32,
        send_first: bool,
        checksum: T,
        tag: i32,
    ) -> MpiResult<T> {
        let send_buf = [checksum];
        let recv = if send_first {
            self.communicator.send(&send_buf, neighbor_rank, tag);
            self.communicator.receive::<T>(neighbor_rank, tag)
        } else {
            let received = self.communicator.receive::<T>(neighbor_rank, tag);
            self.communicator.send(&send_buf, neighbor_rank, tag);
            received
        };
        recv.into_iter().next().ok_or_else(|| {
            MpiError::GhostCellError("Checksum exchange failed".to_string())
        })
    }

    fn checksum_matches(&self, local: T, remote: T) -> bool {
        let abs_tol = T::from_f64_or_zero(1e-10);
        let rel_tol = T::from_f64_or_zero(1e-8);
        let diff = (local - remote).abs();
        let scale = local.abs().max(remote.abs());
        diff <= abs_tol + rel_tol * scale
    }
}
