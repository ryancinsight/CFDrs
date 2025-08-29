//! Checkpoint data structures

use crate::checkpoint::metadata::CheckpointMetadata;
use nalgebra::{DMatrix, RealField};
use serde::{Deserialize, Serialize};

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

    /// Add temperature field
    pub fn with_temperature(mut self, temperature: DMatrix<T>) -> Self {
        self.temperature = Some(temperature);
        self
    }

    /// Add turbulence fields
    pub fn with_turbulence(mut self, k: DMatrix<T>, epsilon: DMatrix<T>) -> Self {
        self.turbulence_k = Some(k);
        self.turbulence_epsilon = Some(epsilon);
        self
    }

    /// Validate checkpoint data consistency
    pub fn validate(&self) -> Result<(), String> {
        self.metadata.validate()?;

        let (ny, nx) = self.metadata.dimensions;

        if self.u_velocity.nrows() != ny || self.u_velocity.ncols() != nx {
            return Err("U velocity field dimension mismatch".to_string());
        }

        if self.v_velocity.nrows() != ny || self.v_velocity.ncols() != nx {
            return Err("V velocity field dimension mismatch".to_string());
        }

        if self.pressure.nrows() != ny || self.pressure.ncols() != nx {
            return Err("Pressure field dimension mismatch".to_string());
        }

        if let Some(ref temp) = self.temperature {
            if temp.nrows() != ny || temp.ncols() != nx {
                return Err("Temperature field dimension mismatch".to_string());
            }
        }

        if let Some(ref k) = self.turbulence_k {
            if k.nrows() != ny || k.ncols() != nx {
                return Err("Turbulence k field dimension mismatch".to_string());
            }
        }

        if let Some(ref eps) = self.turbulence_epsilon {
            if eps.nrows() != ny || eps.ncols() != nx {
                return Err("Turbulence epsilon field dimension mismatch".to_string());
            }
        }

        Ok(())
    }

    /// Get grid dimensions
    pub fn dimensions(&self) -> (usize, usize) {
        self.metadata.dimensions
    }

    /// Get simulation time
    pub fn time(&self) -> f64 {
        self.metadata.time
    }

    /// Get iteration number
    pub fn iteration(&self) -> usize {
        self.metadata.iteration
    }
}
