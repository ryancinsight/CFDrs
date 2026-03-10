//! Asynchronous ghost cell exchange for communication overlap with computation.

use super::sync_exchange::GhostCellManager;
use super::types::{FieldType, GhostExchangeContext, RecvOp, SendOp};
use crate::compute::mpi::communicator::MpiCommunicator;
use crate::compute::mpi::decomposition::{NeighborDirection, NeighborInfo};
use crate::compute::mpi::error::{MpiError, MpiResult};
use crate::compute::mpi::LocalSubdomain;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

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
            NeighborDirection::Front | NeighborDirection::Back => Err(MpiError::NotAvailable(
                "3D ghost-cell exchange requires 3D field storage".to_string(),
            )),
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
