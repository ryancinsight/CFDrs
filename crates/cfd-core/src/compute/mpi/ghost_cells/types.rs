//! Internal types for ghost cell exchange operations.

use super::super::decomposition::NeighborDirection;
use nalgebra::RealField;

#[derive(Debug, Clone, Copy)]
pub(super) enum FieldType {
    Ux,
    Uy,
    Vx,
    Vy,
    P,
}

#[derive(Debug, Clone)]
pub(super) struct RecvOp {
    pub(super) buffer_index: usize,
    pub(super) neighbor_rank: i32,
    pub(super) direction: NeighborDirection,
    pub(super) field_type: FieldType,
    pub(super) tag: i32,
}

#[derive(Debug, Clone)]
pub(super) struct SendOp {
    pub(super) buffer_index: usize,
    pub(super) neighbor_rank: i32,
    pub(super) tag: i32,
}

/// Context for asynchronous ghost cell exchange to ensure buffer lifetime
#[derive(Debug)]
pub struct GhostExchangeContext<T: RealField + Copy> {
    pub(super) send_buffers: Vec<Vec<T>>,
    pub(super) recv_buffers: Vec<Vec<T>>,
    pub(super) recv_ops: Vec<RecvOp>,
    pub(super) send_ops: Vec<SendOp>,
}

impl<T: RealField + Copy> GhostExchangeContext<T> {
    /// Create a new empty exchange context
    pub fn new() -> Self {
        Self {
            send_buffers: Vec::new(),
            recv_buffers: Vec::new(),
            recv_ops: Vec::new(),
            send_ops: Vec::new(),
        }
    }

    /// Clear all buffers and operations
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
