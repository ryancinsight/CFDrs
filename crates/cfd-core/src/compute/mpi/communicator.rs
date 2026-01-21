//! MPI communicator management and basic operations

use super::error::{MpiError, MpiResult};
use super::{Rank, Size};
use mpi::environment;
use mpi::traits::*;
use std::sync::Once;

/// Global MPI initialization flag
static MPI_INIT: Once = Once::new();
static mut MPI_INITIALIZED: bool = false;

/// MPI universe wrapper for safe initialization
#[derive(Debug)]
pub struct MpiUniverse {
    _guard: (), // Ensures this struct owns the MPI environment
}

impl MpiUniverse {
    /// Initialize MPI environment
    ///
    /// This must be called before any other MPI operations.
    /// Returns an error if MPI is already initialized or initialization fails.
    pub fn new() -> MpiResult<Self> {
        let mut initialized = false;

        MPI_INIT.call_once(|| match environment::initialize() {
            Ok(_) => {
                unsafe { MPI_INITIALIZED = true };
                initialized = true;
            }
            Err(e) => {
                initialized = false;
            }
        });

        if !initialized {
            return Err(MpiError::InitializationError(
                "Failed to initialize MPI environment".to_string(),
            ));
        }

        Ok(Self { _guard: () })
    }

    /// Get the world communicator
    pub fn world(&self) -> MpiCommunicator {
        MpiCommunicator {
            comm: environment::Universe::world(),
        }
    }
}

impl Drop for MpiUniverse {
    fn drop(&mut self) {
        // MPI will be finalized when the program exits
        // We don't explicitly finalize here as it might interfere with other parts of the program
    }
}

/// Wrapper around MPI communicator
#[derive(Debug, Clone)]
pub struct MpiCommunicator {
    comm: mpi::topology::SimpleCommunicator,
}

impl MpiCommunicator {
    /// Get the rank of this process
    pub fn rank(&self) -> Rank {
        self.comm.rank()
    }

    /// Get the total number of processes
    pub fn size(&self) -> Size {
        self.comm.size()
    }

    /// Check if this is the root process (rank 0)
    pub fn is_root(&self) -> bool {
        self.rank() == 0
    }

    /// Barrier synchronization
    pub fn barrier(&self) {
        self.comm.barrier();
    }

    /// Broadcast a value from root to all processes
    pub fn broadcast<T: Equivalence>(&self, value: &mut T) {
        self.comm.process_at_rank(0).broadcast_into(value);
    }

    /// All-reduce operation (sum all values across processes)
    pub fn all_reduce_sum<T: Equivalence + std::ops::AddAssign>(&self, value: &mut T) {
        let mut result = *value;
        self.comm
            .all_reduce_into(&result, value, mpi::collective::SystemOperation::sum());
    }

    /// All-reduce operation (find maximum across processes)
    pub fn all_reduce_max<T: Equivalence + PartialOrd + Copy>(&self, value: &mut T) {
        let mut result = *value;
        self.comm
            .all_reduce_into(&result, value, mpi::collective::SystemOperation::max());
    }

    /// All-reduce operation (find minimum across processes)
    pub fn all_reduce_min<T: Equivalence + PartialOrd + Copy>(&self, value: &mut T) {
        let mut result = *value;
        self.comm
            .all_reduce_into(&result, value, mpi::collective::SystemOperation::min());
    }

    /// Send data to another rank
    pub fn send<T: Equivalence>(&self, data: &[T], destination: Rank, tag: i32) {
        self.comm.process_at_rank(destination).send(data);
    }

    /// Receive data from another rank
    pub fn receive<T: Equivalence>(&self, source: Rank, tag: i32) -> Vec<T> {
        let (data, _status) = self.comm.process_at_rank(source).receive_vec();
        data
    }

    /// Non-blocking send
    pub fn send_async<T: Equivalence>(
        &self,
        data: &[T],
        destination: Rank,
        tag: i32,
    ) -> mpi::request::Request<'_, Vec<T>> {
        self.comm.process_at_rank(destination).immediate_send(data)
    }

    /// Non-blocking receive
    pub fn receive_async<T: Equivalence>(
        &self,
        source: Rank,
        tag: i32,
        buffer: &mut [T],
    ) -> mpi::request::Request<'_, Vec<T>> {
        self.comm
            .process_at_rank(source)
            .immediate_receive_into(buffer)
    }

    /// All-gather operation (gather values from all processes to all processes)
    pub fn all_gather<T: Equivalence>(&self, send_buf: &T, recv_buf: &mut [T]) {
        self.comm.all_gather_into(send_buf, recv_buf);
    }

    /// Wait for async operation to complete
    pub fn wait<T>(request: mpi::request::Request<'_, T>) -> T {
        request.wait().0
    }

    /// Check if MPI is available
    pub fn is_available() -> bool {
        unsafe { MPI_INITIALIZED }
    }
}

/// Utility functions for MPI operations
pub mod utils {
    use super::*;

    /// Get MPI rank and size, panicking if MPI not initialized
    pub fn rank_size() -> (Rank, Size) {
        let world = environment::Universe::world();
        (world.rank(), world.size())
    }

    /// Check if current process is root (rank 0)
    pub fn is_root() -> bool {
        rank_size().0 == 0
    }

    /// Synchronize all processes
    pub fn barrier() {
        environment::Universe::world().barrier();
    }

    /// Print from root process only
    pub fn print_root(message: &str) {
        if is_root() {
            println!("{}", message);
        }
    }

    /// Print with rank prefix
    pub fn print_rank(message: &str) {
        let (rank, size) = rank_size();
        println!("[{}/{}] {}", rank, size, message);
    }
}
