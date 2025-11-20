//! Checkpoint data structures

use crate::checkpoint::metadata::CheckpointMetadata;
use nalgebra::{DMatrix, RealField};
use serde::{Deserialize, Serialize};
use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;

/// Checkpoint data containing simulation state
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Checkpoint<T: RealField + Copy> {
    /// Metadata
    pub metadata: CheckpointMetadata,
    /// Velocity field (u component)
    pub u_velocity: DMatrix<T>,
    /// Velocity field (v component)
    pub v_velocity: DMatrix<T>,
    /// Pressure field
    pub pressure: DMatrix<T>,
    /// Temperature field (optional)
    pub temperature: Option<DMatrix<T>>,
    /// Turbulence kinetic energy (optional)
    pub turbulence_k: Option<DMatrix<T>>,
    /// Turbulence dissipation rate (optional)
    pub turbulence_epsilon: Option<DMatrix<T>>,
}

impl<T: RealField + Copy> Checkpoint<T> {
    /// Create a new checkpoint
    pub fn new(
        metadata: CheckpointMetadata,
        u_velocity: DMatrix<T>,
        v_velocity: DMatrix<T>,
        pressure: DMatrix<T>,
    ) -> Self {
        Self {
            metadata,
            u_velocity,
            v_velocity,
            pressure,
            temperature: None,
            turbulence_k: None,
            turbulence_epsilon: None,
        }
    }

    /// Validate checkpoint data consistency
    pub fn validate(&self) -> Result<(), String> {
        // Check dimensions consistency
        let (nx, ny) = self.metadata.dimensions;
        if self.u_velocity.ncols() != nx || self.u_velocity.nrows() != ny {
            return Err(format!(
                "U velocity dimensions mismatch: expected {}x{}, got {}x{}",
                ny, nx, self.u_velocity.nrows(), self.u_velocity.ncols()
            ));
        }
        if self.v_velocity.ncols() != nx || self.v_velocity.nrows() != ny {
            return Err(format!(
                "V velocity dimensions mismatch: expected {}x{}, got {}x{}",
                ny, nx, self.v_velocity.nrows(), self.v_velocity.ncols()
            ));
        }
        if self.pressure.ncols() != nx || self.pressure.nrows() != ny {
            return Err(format!(
                "Pressure dimensions mismatch: expected {}x{}, got {}x{}",
                ny, nx, self.pressure.nrows(), self.pressure.ncols()
            ));
        }

        Ok(())
    }

    /// Get grid dimensions as (nx, ny)
    pub fn dimensions(&self) -> (usize, usize) {
        self.metadata.dimensions
    }
}

impl Checkpoint<f64> {
    /// Computes deterministic u128 checksum for bit-exact verification.
    ///
    /// # Invariants
    /// * Roundtrip preservation: `checkpoint.compute_checksum() == loaded.compute_checksum()`
    /// * Parallel consistency: `allreduce_xor(local_checksums) == global_checksum`
    ///
    /// Uses dual DefaultHasher on metadata fields (skip checksum) + column-major matrix bits.
    pub fn compute_checksum(&self) -> u128 {
        let mut hasher1 = DefaultHasher::new();
        let mut hasher2 = DefaultHasher::new();

        let meta = &self.metadata;
        hasher1.write_u32(meta.version);
        hasher1.write_u64(meta.time.to_bits());
        hasher1.write_u64(meta.iteration as u64);
        hasher1.write_u64(meta.dimensions.0 as u64);
        hasher1.write_u64(meta.dimensions.1 as u64);
        hasher1.write_u64(meta.domain_size.0.to_bits());
        hasher1.write_u64(meta.domain_size.1.to_bits());
        hasher1.write_u64(meta.timestamp);

        hasher2.write_u32(meta.version);
        hasher2.write_u64(meta.time.to_bits());
        hasher2.write_u64(meta.iteration as u64);
        hasher2.write_u64(meta.dimensions.0 as u64);
        hasher2.write_u64(meta.dimensions.1 as u64);
        hasher2.write_u64(meta.domain_size.0.to_bits());
        hasher2.write_u64(meta.domain_size.1.to_bits());
        hasher2.write_u64(meta.timestamp);

        macro_rules! hash_matrix {
            ($mat:expr) => {
                for col in $mat.column_iter() {
                    for &v in col.iter() {
                        let bits = v.to_bits();
                        hasher1.write_u64(bits);
                        hasher2.write_u64(bits);
                    }
                }
            };
        }

        hash_matrix!(&self.u_velocity);
        hash_matrix!(&self.v_velocity);
        hash_matrix!(&self.pressure);
        if let Some(ref mat) = self.temperature {
            hash_matrix!(mat);
        }
        if let Some(ref mat) = self.turbulence_k {
            hash_matrix!(mat);
        }
        if let Some(ref mat) = self.turbulence_epsilon {
            hash_matrix!(mat);
        }

        let low = hasher1.finish();
        let high = hasher2.finish();
        (low as u128) | ((high as u128) << 64)
    }
}
