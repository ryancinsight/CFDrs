//! Ghost cell management for MPI communication

use super::communicator::MpiCommunicator;
use super::decomposition::{NeighborDirection, NeighborInfo};
use super::error::{MpiError, MpiResult};
use super::LocalSubdomain;
use nalgebra::{RealField, Vector2};
use std::collections::HashMap;

/// Manages ghost cell communication between MPI processes
#[derive(Debug)]
pub struct GhostCellManager<T: RealField + Copy> {
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Neighbor information
    neighbors: HashMap<i32, NeighborInfo>,
    /// Number of ghost cell layers
    ghost_layers: usize,
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> GhostCellManager<T> {
    /// Create a new ghost cell manager
    pub fn new(
        communicator: MpiCommunicator,
        neighbors: HashMap<i32, NeighborInfo>,
        ghost_layers: usize,
    ) -> Self {
        Self {
            communicator,
            neighbors,
            ghost_layers,
        }
    }

    /// Update ghost cells for all fields
    pub fn update_ghost_cells(
        &self,
        velocity_u: &mut [Vec<Vector2<T>>],
        velocity_v: &mut [Vec<Vector2<T>>],
        pressure: &mut [Vec<T>],
        subdomain: &LocalSubdomain,
    ) -> MpiResult<()> {
        // For each neighbor, exchange ghost cells
        for (neighbor_rank, neighbor_info) in &self.neighbors {
            self.exchange_ghost_cells(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                *neighbor_rank,
                neighbor_info,
            )?;
        }

        Ok(())
    }

    /// Exchange ghost cells with a specific neighbor
    fn exchange_ghost_cells(
        &self,
        velocity_u: &[Vec<Vector2<T>>],
        velocity_v: &[Vec<Vector2<T>>],
        pressure: &[Vec<T>],
        subdomain: &LocalSubdomain,
        neighbor_rank: i32,
        neighbor_info: &NeighborInfo,
    ) -> MpiResult<()> {
        match neighbor_info.direction {
            NeighborDirection::Left => self.exchange_left_right(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                true,
            ),
            NeighborDirection::Right => self.exchange_left_right(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                false,
            ),
            NeighborDirection::Bottom => self.exchange_bottom_top(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                true,
            ),
            NeighborDirection::Top => self.exchange_bottom_top(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                false,
            ),
            _ => Ok(()), // TODO: Implement 3D ghost-cell exchange (front/back).
        }
    }

    /// Exchange ghost cells in left/right direction
    fn exchange_left_right(
        &self,
        velocity_u: &mut [Vec<Vector2<T>>],
        velocity_v: &mut [Vec<Vector2<T>>],
        pressure: &mut [Vec<T>],
        subdomain: &LocalSubdomain,
        neighbor_rank: i32,
        is_left: bool,
    ) -> MpiResult<()> {
        let ny = subdomain.total_ny();

        // Prepare send buffer (ghost layer from owned region)
        let send_i = if is_left {
            subdomain.ghost_layers
        } else {
            subdomain.ghost_layers + subdomain.nx_local - 1
        };

        // Prepare receive buffer position (ghost layer)
        let recv_i = if is_left {
            subdomain.ghost_layers - 1
        } else {
            subdomain.ghost_layers + subdomain.nx_local
        };

        // For Vector2<T>, we need to send individual components
        // Prepare send buffers for u velocity components
        let mut send_u_x = Vec::with_capacity(ny);
        let mut send_u_y = Vec::with_capacity(ny);
        let mut send_v_x = Vec::with_capacity(ny);
        let mut send_v_y = Vec::with_capacity(ny);
        let mut send_p = Vec::with_capacity(ny);

        for j in 0..ny {
            send_u_x.push(velocity_u[send_i][j].x);
            send_u_y.push(velocity_u[send_i][j].y);
            send_v_x.push(velocity_v[send_i][j].x);
            send_v_y.push(velocity_v[send_i][j].y);
            send_p.push(pressure[send_i][j]);
        }

        // Prepare receive buffers
        let mut recv_u_x = vec![T::zero(); ny];
        let mut recv_u_y = vec![T::zero(); ny];
        let mut recv_v_x = vec![T::zero(); ny];
        let mut recv_v_y = vec![T::zero(); ny];
        let mut recv_p = vec![T::zero(); ny];

        // Use unique tags for each field to avoid message conflicts
        let tag_base = if is_left { 100 } else { 200 };

        // Exchange data with proper ordering to avoid deadlocks
        if is_left {
            // Send to left, receive from left (standard MPI ordering)
            self.communicator.send(&send_u_x, neighbor_rank, tag_base);
            self.communicator
                .send(&send_u_y, neighbor_rank, tag_base + 1);
            self.communicator
                .send(&send_v_x, neighbor_rank, tag_base + 2);
            self.communicator
                .send(&send_v_y, neighbor_rank, tag_base + 3);
            self.communicator.send(&send_p, neighbor_rank, tag_base + 4);

            recv_u_x = self.communicator.receive(neighbor_rank, tag_base);
            recv_u_y = self.communicator.receive(neighbor_rank, tag_base + 1);
            recv_v_x = self.communicator.receive(neighbor_rank, tag_base + 2);
            recv_v_y = self.communicator.receive(neighbor_rank, tag_base + 3);
            recv_p = self.communicator.receive(neighbor_rank, tag_base + 4);
        } else {
            // Receive from right, send to right (reverse ordering to match left neighbor)
            recv_u_x = self.communicator.receive(neighbor_rank, tag_base);
            recv_u_y = self.communicator.receive(neighbor_rank, tag_base + 1);
            recv_v_x = self.communicator.receive(neighbor_rank, tag_base + 2);
            recv_v_y = self.communicator.receive(neighbor_rank, tag_base + 3);
            recv_p = self.communicator.receive(neighbor_rank, tag_base + 4);

            self.communicator.send(&send_u_x, neighbor_rank, tag_base);
            self.communicator
                .send(&send_u_y, neighbor_rank, tag_base + 1);
            self.communicator
                .send(&send_v_x, neighbor_rank, tag_base + 2);
            self.communicator
                .send(&send_v_y, neighbor_rank, tag_base + 3);
            self.communicator.send(&send_p, neighbor_rank, tag_base + 4);
        }

        // Copy received data to ghost cells
        for j in 0..ny {
            velocity_u[recv_i][j] = Vector2::new(recv_u_x[j], recv_u_y[j]);
            velocity_v[recv_i][j] = Vector2::new(recv_v_x[j], recv_v_y[j]);
            pressure[recv_i][j] = recv_p[j];
        }

        Ok(())
    }

    /// Exchange ghost cells in bottom/top direction
    fn exchange_bottom_top(
        &self,
        velocity_u: &mut [Vec<Vector2<T>>],
        velocity_v: &mut [Vec<Vector2<T>>],
        pressure: &mut [Vec<T>],
        subdomain: &LocalSubdomain,
        neighbor_rank: i32,
        is_bottom: bool,
    ) -> MpiResult<()> {
        let nx = subdomain.total_nx();

        // Prepare send buffer (ghost layer from owned region)
        let send_j = if is_bottom {
            subdomain.ghost_layers
        } else {
            subdomain.ghost_layers + subdomain.ny_local - 1
        };

        // Prepare receive buffer position (ghost layer)
        let recv_j = if is_bottom {
            subdomain.ghost_layers - 1
        } else {
            subdomain.ghost_layers + subdomain.ny_local
        };

        // Prepare send buffers for all i positions at this j layer
        let mut send_u_x = Vec::with_capacity(nx);
        let mut send_u_y = Vec::with_capacity(nx);
        let mut send_v_x = Vec::with_capacity(nx);
        let mut send_v_y = Vec::with_capacity(nx);
        let mut send_p = Vec::with_capacity(nx);

        for i in 0..nx {
            send_u_x.push(velocity_u[i][send_j].x);
            send_u_y.push(velocity_u[i][send_j].y);
            send_v_x.push(velocity_v[i][send_j].x);
            send_v_y.push(velocity_v[i][send_j].y);
            send_p.push(pressure[i][send_j]);
        }

        // Prepare receive buffers
        let mut recv_u_x = vec![T::zero(); nx];
        let mut recv_u_y = vec![T::zero(); nx];
        let mut recv_v_x = vec![T::zero(); nx];
        let mut recv_v_y = vec![T::zero(); nx];
        let mut recv_p = vec![T::zero(); nx];

        // Use unique tags for bottom/top communication (different from left/right)
        let tag_base = if is_bottom { 300 } else { 400 };

        // Exchange data with proper ordering to avoid deadlocks
        if is_bottom {
            // Send to bottom, receive from bottom
            self.communicator.send(&send_u_x, neighbor_rank, tag_base);
            self.communicator
                .send(&send_u_y, neighbor_rank, tag_base + 1);
            self.communicator
                .send(&send_v_x, neighbor_rank, tag_base + 2);
            self.communicator
                .send(&send_v_y, neighbor_rank, tag_base + 3);
            self.communicator.send(&send_p, neighbor_rank, tag_base + 4);

            recv_u_x = self.communicator.receive(neighbor_rank, tag_base);
            recv_u_y = self.communicator.receive(neighbor_rank, tag_base + 1);
            recv_v_x = self.communicator.receive(neighbor_rank, tag_base + 2);
            recv_v_y = self.communicator.receive(neighbor_rank, tag_base + 3);
            recv_p = self.communicator.receive(neighbor_rank, tag_base + 4);
        } else {
            // Receive from top, send to top
            recv_u_x = self.communicator.receive(neighbor_rank, tag_base);
            recv_u_y = self.communicator.receive(neighbor_rank, tag_base + 1);
            recv_v_x = self.communicator.receive(neighbor_rank, tag_base + 2);
            recv_v_y = self.communicator.receive(neighbor_rank, tag_base + 3);
            recv_p = self.communicator.receive(neighbor_rank, tag_base + 4);

            self.communicator.send(&send_u_x, neighbor_rank, tag_base);
            self.communicator
                .send(&send_u_y, neighbor_rank, tag_base + 1);
            self.communicator
                .send(&send_v_x, neighbor_rank, tag_base + 2);
            self.communicator
                .send(&send_v_y, neighbor_rank, tag_base + 3);
            self.communicator.send(&send_p, neighbor_rank, tag_base + 4);
        }

        // Copy received data to ghost cells
        for i in 0..nx {
            velocity_u[i][recv_j] = Vector2::new(recv_u_x[i], recv_u_y[i]);
            velocity_v[i][recv_j] = Vector2::new(recv_v_x[i], recv_v_y[i]);
            pressure[i][recv_j] = recv_p[i];
        }

        Ok(())
    }

    /// Validate ghost cell data consistency
    pub fn validate_ghost_cells(
        &self,
        _velocity_u: &[Vec<Vector2<T>>],
        _velocity_v: &[Vec<Vector2<T>>],
        _pressure: &[Vec<T>],
        _subdomain: &LocalSubdomain,
    ) -> MpiResult<bool> {
        // TODO: Validate ghost cell consistency across process boundaries (checksums + neighbor compare).
        Ok(true)
    }
}

/// Configuration for ghost cell updates
#[derive(Debug, Clone)]
pub struct GhostCellUpdate {
    /// Maximum iterations for convergence
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: f64,
    /// Communication pattern (blocking vs non-blocking)
    pub blocking: bool,
}

/// Asynchronous ghost cell exchange for communication overlap
impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> GhostCellManager<T> {
    /// Update ghost cells asynchronously for communication overlap with computation
    pub fn update_ghost_cells_async(
        &self,
        velocity_u: &mut [Vec<Vector2<T>>],
        velocity_v: &mut [Vec<Vector2<T>>],
        pressure: &mut [Vec<T>],
        subdomain: &LocalSubdomain,
    ) -> MpiResult<Vec<mpi::request::Request<'_, Vec<T>>>> {
        let mut requests = Vec::new();

        // Start asynchronous exchanges for each neighbor
        for (neighbor_rank, neighbor_info) in &self.neighbors {
            let reqs = self.exchange_ghost_cells_async(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                *neighbor_rank,
                neighbor_info,
            )?;
            requests.extend(reqs);
        }

        Ok(requests)
    }

    /// Asynchronous exchange with a specific neighbor
    fn exchange_ghost_cells_async(
        &self,
        velocity_u: &mut [Vec<Vector2<T>>],
        velocity_v: &mut [Vec<Vector2<T>>],
        pressure: &mut [Vec<T>],
        subdomain: &LocalSubdomain,
        neighbor_rank: i32,
        neighbor_info: &NeighborInfo,
    ) -> MpiResult<Vec<mpi::request::Request<'_, Vec<T>>>> {
        match neighbor_info.direction {
            NeighborDirection::Left => self.exchange_left_right_async(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                true,
            ),
            NeighborDirection::Right => self.exchange_left_right_async(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                false,
            ),
            NeighborDirection::Bottom => self.exchange_bottom_top_async(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                true,
            ),
            NeighborDirection::Top => self.exchange_bottom_top_async(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                false,
            ),
            _ => Ok(Vec::new()), // TODO: Implement 3D async ghost-cell exchange (front/back).
        }
    }

    /// Asynchronous left/right exchange
    fn exchange_left_right_async(
        &self,
        velocity_u: &mut [Vec<Vector2<T>>],
        velocity_v: &mut [Vec<Vector2<T>>],
        pressure: &mut [Vec<T>],
        subdomain: &LocalSubdomain,
        neighbor_rank: i32,
        is_left: bool,
    ) -> MpiResult<Vec<mpi::request::Request<'_, Vec<T>>>> {
        let ny = subdomain.total_ny();

        // Prepare send buffer (ghost layer from owned region)
        let send_i = if is_left {
            subdomain.ghost_layers
        } else {
            subdomain.ghost_layers + subdomain.nx_local - 1
        };

        // Prepare receive buffer position (ghost layer)
        let recv_i = if is_left {
            subdomain.ghost_layers - 1
        } else {
            subdomain.ghost_layers + subdomain.nx_local
        };

        // Prepare send buffers for u velocity components
        let mut send_u_x = Vec::with_capacity(ny);
        let mut send_u_y = Vec::with_capacity(ny);
        let mut send_v_x = Vec::with_capacity(ny);
        let mut send_v_y = Vec::with_capacity(ny);
        let mut send_p = Vec::with_capacity(ny);

        for j in 0..ny {
            send_u_x.push(velocity_u[send_i][j].x);
            send_u_y.push(velocity_u[send_i][j].y);
            send_v_x.push(velocity_v[send_i][j].x);
            send_v_y.push(velocity_v[send_i][j].y);
            send_p.push(pressure[send_i][j]);
        }

        // Use unique tags for async communication
        let tag_base = if is_left { 500 } else { 600 };
        let mut requests = Vec::new();

        // Start asynchronous sends and receives
        if is_left {
            // Send to left, receive from left
            requests.push(
                self.communicator
                    .send_async(&send_u_x, neighbor_rank, tag_base),
            );
            requests.push(
                self.communicator
                    .send_async(&send_u_y, neighbor_rank, tag_base + 1),
            );
            requests.push(
                self.communicator
                    .send_async(&send_v_x, neighbor_rank, tag_base + 2),
            );
            requests.push(
                self.communicator
                    .send_async(&send_v_y, neighbor_rank, tag_base + 3),
            );
            requests.push(
                self.communicator
                    .send_async(&send_p, neighbor_rank, tag_base + 4),
            );

            // Store receive requests for later completion
            let recv_u_x_req = self.communicator.receive_async(neighbor_rank, tag_base);
            let recv_u_y_req = self.communicator.receive_async(neighbor_rank, tag_base + 1);
            let recv_v_x_req = self.communicator.receive_async(neighbor_rank, tag_base + 2);
            let recv_v_y_req = self.communicator.receive_async(neighbor_rank, tag_base + 3);
            let recv_p_req = self.communicator.receive_async(neighbor_rank, tag_base + 4);

            requests.push(recv_u_x_req);
            requests.push(recv_u_y_req);
            requests.push(recv_v_x_req);
            requests.push(recv_v_y_req);
            requests.push(recv_p_req);
        } else {
            // Receive from right, send to right
            let recv_u_x_req = self.communicator.receive_async(neighbor_rank, tag_base);
            let recv_u_y_req = self.communicator.receive_async(neighbor_rank, tag_base + 1);
            let recv_v_x_req = self.communicator.receive_async(neighbor_rank, tag_base + 2);
            let recv_v_y_req = self.communicator.receive_async(neighbor_rank, tag_base + 3);
            let recv_p_req = self.communicator.receive_async(neighbor_rank, tag_base + 4);

            requests.push(recv_u_x_req);
            requests.push(recv_u_y_req);
            requests.push(recv_v_x_req);
            requests.push(recv_v_y_req);
            requests.push(recv_p_req);

            requests.push(
                self.communicator
                    .send_async(&send_u_x, neighbor_rank, tag_base),
            );
            requests.push(
                self.communicator
                    .send_async(&send_u_y, neighbor_rank, tag_base + 1),
            );
            requests.push(
                self.communicator
                    .send_async(&send_v_x, neighbor_rank, tag_base + 2),
            );
            requests.push(
                self.communicator
                    .send_async(&send_v_y, neighbor_rank, tag_base + 3),
            );
            requests.push(
                self.communicator
                    .send_async(&send_p, neighbor_rank, tag_base + 4),
            );
        }

        Ok(requests)
    }

    /// Asynchronous bottom/top exchange
    fn exchange_bottom_top_async(
        &self,
        velocity_u: &mut [Vec<Vector2<T>>],
        velocity_v: &mut [Vec<Vector2<T>>],
        pressure: &mut [Vec<T>],
        subdomain: &LocalSubdomain,
        neighbor_rank: i32,
        is_bottom: bool,
    ) -> MpiResult<Vec<mpi::request::Request<'_, Vec<T>>>> {
        let nx = subdomain.total_nx();

        // Prepare send buffer (ghost layer from owned region)
        let send_j = if is_bottom {
            subdomain.ghost_layers
        } else {
            subdomain.ghost_layers + subdomain.ny_local - 1
        };

        // Prepare send buffers for all i positions at this j layer
        let mut send_u_x = Vec::with_capacity(nx);
        let mut send_u_y = Vec::with_capacity(nx);
        let mut send_v_x = Vec::with_capacity(nx);
        let mut send_v_y = Vec::with_capacity(nx);
        let mut send_p = Vec::with_capacity(nx);

        for i in 0..nx {
            send_u_x.push(velocity_u[i][send_j].x);
            send_u_y.push(velocity_u[i][send_j].y);
            send_v_x.push(velocity_v[i][send_j].x);
            send_v_y.push(velocity_v[i][send_j].y);
            send_p.push(pressure[i][send_j]);
        }

        // Use unique tags for async bottom/top communication
        let tag_base = if is_bottom { 700 } else { 800 };
        let mut requests = Vec::new();

        // Start asynchronous sends and receives
        if is_bottom {
            // Send to bottom, receive from bottom
            requests.push(
                self.communicator
                    .send_async(&send_u_x, neighbor_rank, tag_base),
            );
            requests.push(
                self.communicator
                    .send_async(&send_u_y, neighbor_rank, tag_base + 1),
            );
            requests.push(
                self.communicator
                    .send_async(&send_v_x, neighbor_rank, tag_base + 2),
            );
            requests.push(
                self.communicator
                    .send_async(&send_v_y, neighbor_rank, tag_base + 3),
            );
            requests.push(
                self.communicator
                    .send_async(&send_p, neighbor_rank, tag_base + 4),
            );

            let recv_u_x_req = self.communicator.receive_async(neighbor_rank, tag_base);
            let recv_u_y_req = self.communicator.receive_async(neighbor_rank, tag_base + 1);
            let recv_v_x_req = self.communicator.receive_async(neighbor_rank, tag_base + 2);
            let recv_v_y_req = self.communicator.receive_async(neighbor_rank, tag_base + 3);
            let recv_p_req = self.communicator.receive_async(neighbor_rank, tag_base + 4);

            requests.push(recv_u_x_req);
            requests.push(recv_u_y_req);
            requests.push(recv_v_x_req);
            requests.push(recv_v_y_req);
            requests.push(recv_p_req);
        } else {
            // Receive from top, send to top
            let recv_u_x_req = self.communicator.receive_async(neighbor_rank, tag_base);
            let recv_u_y_req = self.communicator.receive_async(neighbor_rank, tag_base + 1);
            let recv_v_x_req = self.communicator.receive_async(neighbor_rank, tag_base + 2);
            let recv_v_y_req = self.communicator.receive_async(neighbor_rank, tag_base + 3);
            let recv_p_req = self.communicator.receive_async(neighbor_rank, tag_base + 4);

            requests.push(recv_u_x_req);
            requests.push(recv_u_y_req);
            requests.push(recv_v_x_req);
            requests.push(recv_v_y_req);
            requests.push(recv_p_req);

            requests.push(
                self.communicator
                    .send_async(&send_u_x, neighbor_rank, tag_base),
            );
            requests.push(
                self.communicator
                    .send_async(&send_u_y, neighbor_rank, tag_base + 1),
            );
            requests.push(
                self.communicator
                    .send_async(&send_v_x, neighbor_rank, tag_base + 2),
            );
            requests.push(
                self.communicator
                    .send_async(&send_v_y, neighbor_rank, tag_base + 3),
            );
            requests.push(
                self.communicator
                    .send_async(&send_p, neighbor_rank, tag_base + 4),
            );
        }

        Ok(requests)
    }

    /// Complete asynchronous ghost cell updates
    pub fn complete_async_updates(
        &self,
        velocity_u: &mut [Vec<Vector2<T>>],
        velocity_v: &mut [Vec<Vector2<T>>],
        pressure: &mut [Vec<T>],
        subdomain: &LocalSubdomain,
        requests: Vec<mpi::request::Request<'_, Vec<T>>>,
    ) -> MpiResult<()> {
        // TODO: Track request-to-field mapping and unpack received buffers into ghost layers.
        for request in requests {
            let _data = MpiCommunicator::wait(request);
        }

        Ok(())
    }
}

/// Statistics about ghost cell communication
#[derive(Debug, Clone)]
pub struct GhostCellStats {
    /// Total bytes sent
    pub bytes_sent: usize,
    /// Total bytes received
    pub bytes_received: usize,
    /// Number of messages sent
    pub messages_sent: usize,
    /// Number of messages received
    pub messages_received: usize,
    /// Communication time
    pub comm_time: std::time::Duration,
}

impl Default for GhostCellStats {
    fn default() -> Self {
        Self {
            bytes_sent: 0,
            bytes_received: 0,
            messages_sent: 0,
            messages_received: 0,
            comm_time: std::time::Duration::default(),
        }
    }
}
