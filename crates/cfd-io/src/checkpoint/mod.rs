//! Checkpoint and restart functionality for CFD simulations.
//!
//! Provides robust checkpointing with validation and versioning.

mod compression;
mod data;
mod manager;
mod metadata;
mod validator;

pub use compression::CompressionStrategy;
pub use data::Checkpoint;
pub use manager::CheckpointManager;
pub use metadata::{CheckpointMetadata, CHECKPOINT_VERSION};
pub use validator::CheckpointValidator;

// Convenience functions removed - users must explicitly create CheckpointManager
// to ensure proper configuration for compression handling
