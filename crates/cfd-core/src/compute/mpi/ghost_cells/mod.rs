//! Ghost cell management for MPI communication.
//!
//! Provides synchronous and asynchronous ghost cell exchange between MPI
//! processes, including validation via checksums and communication statistics.

/// Asynchronous ghost cell exchange for overlapped communication.
mod async_exchange;
/// Ghost cell manager and synchronous exchange operations.
mod sync_exchange;
/// Internal types for ghost cell exchange operations.
mod types;
/// Ghost cell data validation and checksum verification.
mod validation;

pub use sync_exchange::GhostCellManager;
pub use types::{GhostCellStats, GhostCellUpdate, GhostExchangeContext};
