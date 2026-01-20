//! Ghost cell management for MPI communication

use super::communicator::MpiCommunicator;
use super::decomposition::{NeighborDirection, NeighborInfo};
use super::error::{MpiError, MpiResult};
use super::LocalSubdomain;
use nalgebra::{RealField, Vector2};
use std::collections::HashMap;
use num_traits::FromPrimitive;

#[derive(Debug, Clone, Copy)]
enum FieldType {
    Ux,
    Uy,
    Vx,
    Vy,
    P,
}

#[derive(Debug, Clone)]
struct RecvOp {
    buffer_index: usize,
    neighbor_rank: i32,
    direction: NeighborDirection,
    field_type: FieldType,
    tag: i32,
}

#[derive(Debug, Clone)]
struct SendOp {
    buffer_index: usize,
    neighbor_rank: i32,
    tag: i32,
}

/// Context for asynchronous ghost cell exchange to ensure buffer lifetime
#[derive(Debug)]
pub struct GhostExchangeContext<T: RealField + Copy> {
    send_buffers: Vec<Vec<T>>,
    recv_buffers: Vec<Vec<T>>,
    recv_ops: Vec<RecvOp>,
    send_ops: Vec<SendOp>,
}

impl<T: RealField + Copy> GhostExchangeContext<T> {
    pub fn new() -> Self {
        Self {
            send_buffers: Vec::new(),
            recv_buffers: Vec::new(),
            recv_ops: Vec::new(),
            send_ops: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.send_buffers.clear();
        self.recv_buffers.clear();
        self.recv_ops.clear();
        self.send_ops.clear();
    }
}

impl<T: RealField + Copy> Default for GhostExchangeContext<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Manages ghost cell communication between MPI processes
#[derive(Debug)]
pub struct GhostCellManager<T: RealField + Copy> {
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Neighbor information
    neighbors: HashMap<i32, NeighborInfo>,
    /// Number of ghost cell layers
    ghost_layers: usize,
    _marker: std::marker::PhantomData<T>,
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
            _marker: std::marker::PhantomData,
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
    pub fn start_async_update<'a>(
        &self,
        velocity_u: &[Vec<Vector2<T>>],
        velocity_v: &[Vec<Vector2<T>>],
        pressure: &[Vec<T>],
        subdomain: &LocalSubdomain,
        context: &'a mut GhostExchangeContext<T>,
    ) -> MpiResult<Vec<mpi::request::Request<'a, Vec<T>>>> {
        context.clear();

        // Phase 1: Prepare all buffers and operations
        for (neighbor_rank, neighbor_info) in &self.neighbors {
            self.prepare_ghost_cells(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                *neighbor_rank,
                neighbor_info,
                context,
            )?;
        }

        let mut requests = Vec::new();
        let comm = &self.communicator;

        // Phase 2: Execute operations

        // Send operations - strict 1:1 mapping between buffers and ops
        for (buffer, op) in context.send_buffers.iter().zip(context.send_ops.iter()) {
            requests.push(comm.send_async(buffer, op.neighbor_rank, op.tag));
        }

        // Receive operations - strict 1:1 mapping between buffers and ops
        for (buffer, op) in context.recv_buffers.iter_mut().zip(context.recv_ops.iter()) {
            requests.push(comm.receive_async(op.neighbor_rank, op.tag, buffer));
        }

        Ok(requests)
    }

    /// Prepare ghost cell exchange with a specific neighbor
    fn prepare_ghost_cells(
        &self,
        velocity_u: &[Vec<Vector2<T>>],
        velocity_v: &[Vec<Vector2<T>>],
        pressure: &[Vec<T>],
        subdomain: &LocalSubdomain,
        neighbor_rank: i32,
        neighbor_info: &NeighborInfo,
        context: &mut GhostExchangeContext<T>,
    ) -> MpiResult<()> {
        match neighbor_info.direction {
            NeighborDirection::Left => self.prepare_left_right(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                true,
                context,
            ),
            NeighborDirection::Right => self.prepare_left_right(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                false,
                context,
            ),
            NeighborDirection::Bottom => self.prepare_bottom_top(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                true,
                context,
            ),
            NeighborDirection::Top => self.prepare_bottom_top(
                velocity_u,
                velocity_v,
                pressure,
                subdomain,
                neighbor_rank,
                false,
                context,
            ),
            _ => Ok(()), // TODO: Implement 3D async ghost-cell exchange (front/back).
        }
    }

    /// Prepare left/right exchange
    fn prepare_left_right(
        &self,
        velocity_u: &[Vec<Vector2<T>>],
        velocity_v: &[Vec<Vector2<T>>],
        pressure: &[Vec<T>],
        subdomain: &LocalSubdomain,
        neighbor_rank: i32,
        is_left: bool,
        context: &mut GhostExchangeContext<T>,
    ) -> MpiResult<()> {
        let ny = subdomain.total_ny();

        // Prepare send buffer (ghost layer from owned region)
        let send_i = if is_left {
            subdomain.ghost_layers
        } else {
            subdomain.ghost_layers + subdomain.nx_local - 1
        };

        // For unpacking
        let direction = if is_left {
            NeighborDirection::Left
        } else {
            NeighborDirection::Right
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

        // Push send buffers to context
        let send_idx = context.send_buffers.len();
        context.send_buffers.push(send_u_x);
        context.send_buffers.push(send_u_y);
        context.send_buffers.push(send_v_x);
        context.send_buffers.push(send_v_y);
        context.send_buffers.push(send_p);

        // Create receive buffers
        let recv_idx = context.recv_buffers.len();
        for _ in 0..5 {
            context.recv_buffers.push(vec![T::zero(); ny]);
        }

        let tag_base = if is_left { 500 } else { 600 };

        // Register ops for unpacking and request creation
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx,
            neighbor_rank,
            direction,
            field_type: FieldType::Ux,
            tag: tag_base,
        });
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx + 1,
            neighbor_rank,
            direction,
            field_type: FieldType::Uy,
            tag: tag_base + 1,
        });
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx + 2,
            neighbor_rank,
            direction,
            field_type: FieldType::Vx,
            tag: tag_base + 2,
        });
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx + 3,
            neighbor_rank,
            direction,
            field_type: FieldType::Vy,
            tag: tag_base + 3,
        });
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx + 4,
            neighbor_rank,
            direction,
            field_type: FieldType::P,
            tag: tag_base + 4,
        });

        for i in 0..5 {
            context.send_ops.push(SendOp {
                buffer_index: send_idx + i,
                neighbor_rank,
                tag: tag_base + i as i32,
            });
        }

        Ok(())
    }

    /// Prepare bottom/top exchange
    fn prepare_bottom_top(
        &self,
        velocity_u: &[Vec<Vector2<T>>],
        velocity_v: &[Vec<Vector2<T>>],
        pressure: &[Vec<T>],
        subdomain: &LocalSubdomain,
        neighbor_rank: i32,
        is_bottom: bool,
        context: &mut GhostExchangeContext<T>,
    ) -> MpiResult<()> {
        let nx = subdomain.total_nx();

        // Prepare send buffer (ghost layer from owned region)
        let send_j = if is_bottom {
            subdomain.ghost_layers
        } else {
            subdomain.ghost_layers + subdomain.ny_local - 1
        };

        let direction = if is_bottom {
            NeighborDirection::Bottom
        } else {
            NeighborDirection::Top
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

        // Push send buffers
        let send_idx = context.send_buffers.len();
        context.send_buffers.push(send_u_x);
        context.send_buffers.push(send_u_y);
        context.send_buffers.push(send_v_x);
        context.send_buffers.push(send_v_y);
        context.send_buffers.push(send_p);

        // Create receive buffers
        let recv_idx = context.recv_buffers.len();
        for _ in 0..5 {
            context.recv_buffers.push(vec![T::zero(); nx]);
        }

        let tag_base = if is_bottom { 700 } else { 800 };

        // Register ops
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx,
            neighbor_rank,
            direction,
            field_type: FieldType::Ux,
            tag: tag_base,
        });
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx + 1,
            neighbor_rank,
            direction,
            field_type: FieldType::Uy,
            tag: tag_base + 1,
        });
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx + 2,
            neighbor_rank,
            direction,
            field_type: FieldType::Vx,
            tag: tag_base + 2,
        });
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx + 3,
            neighbor_rank,
            direction,
            field_type: FieldType::Vy,
            tag: tag_base + 3,
        });
        context.recv_ops.push(RecvOp {
            buffer_index: recv_idx + 4,
            neighbor_rank,
            direction,
            field_type: FieldType::P,
            tag: tag_base + 4,
        });

        for i in 0..5 {
            context.send_ops.push(SendOp {
                buffer_index: send_idx + i,
                neighbor_rank,
                tag: tag_base + i as i32,
            });
        }

        Ok(())
    }

    /// Complete asynchronous ghost cell updates
    pub fn complete_async_updates(
        &self,
        velocity_u: &mut [Vec<Vector2<T>>],
        velocity_v: &mut [Vec<Vector2<T>>],
        pressure: &mut [Vec<T>],
        subdomain: &LocalSubdomain,
        requests: Vec<mpi::request::Request<'_, Vec<T>>>,
        context: &GhostExchangeContext<T>,
    ) -> MpiResult<()> {
        // Wait for all requests to complete
        for request in requests {
            let _ = MpiCommunicator::wait(request);
        }

        // Unpack received buffers into ghost layers
        for op in &context.recv_ops {
            let buffer = &context.recv_buffers[op.buffer_index];

            match op.direction {
                NeighborDirection::Left | NeighborDirection::Right => {
                    let is_left = op.direction == NeighborDirection::Left;
                    let recv_i = if is_left {
                        subdomain.ghost_layers - 1
                    } else {
                        subdomain.ghost_layers + subdomain.nx_local
                    };

                    let ny = subdomain.total_ny();
                    if buffer.len() != ny {
                        continue;
                    }

                    for j in 0..ny {
                        match op.field_type {
                            FieldType::Ux => velocity_u[recv_i][j].x = buffer[j],
                            FieldType::Uy => velocity_u[recv_i][j].y = buffer[j],
                            FieldType::Vx => velocity_v[recv_i][j].x = buffer[j],
                            FieldType::Vy => velocity_v[recv_i][j].y = buffer[j],
                            FieldType::P => pressure[recv_i][j] = buffer[j],
                        }
                    }
                }
                NeighborDirection::Bottom | NeighborDirection::Top => {
                    let is_bottom = op.direction == NeighborDirection::Bottom;
                    let recv_j = if is_bottom {
                        subdomain.ghost_layers - 1
                    } else {
                        subdomain.ghost_layers + subdomain.ny_local
                    };

                    let nx = subdomain.total_nx();
                    if buffer.len() != nx {
                        continue;
                    }

                    for i in 0..nx {
                        match op.field_type {
                            FieldType::Ux => velocity_u[i][recv_j].x = buffer[i],
                            FieldType::Uy => velocity_u[i][recv_j].y = buffer[i],
                            FieldType::Vx => velocity_v[i][recv_j].x = buffer[i],
                            FieldType::Vy => velocity_v[i][recv_j].y = buffer[i],
                            FieldType::P => pressure[i][recv_j] = buffer[i],
                        }
                    }
                }
                _ => {}
            }
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
