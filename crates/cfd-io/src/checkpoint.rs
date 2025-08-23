//! Checkpoint and restart functionality.

use cfd_core::error::Result;
use nalgebra::RealField;

/// Checkpoint data
pub struct Checkpoint<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

/// Checkpoint manager
pub struct CheckpointManager<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> CheckpointManager<T> {
    /// Create a new checkpoint manager
    #[must_use] pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Save checkpoint
    pub fn save(&self, _checkpoint: &Checkpoint<T>) -> Result<()> {
        Ok(())
    }

    /// Load checkpoint
    pub fn load(&self, _path: &std::path::Path) -> Result<Checkpoint<T>> {
        Ok(Checkpoint {
            _phantom: std::marker::PhantomData,
        })
    }
}

impl<T: RealField + Copy> Default for CheckpointManager<T> {
    fn default() -> Self {
        Self::new()
    }
}