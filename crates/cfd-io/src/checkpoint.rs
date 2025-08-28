//! Checkpoint and restart functionality for CFD simulations.
//!
//! Provides robust checkpointing with validation and versioning.

use cfd_core::error::{Error, Result};
use nalgebra::{DMatrix, DVector, RealField};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

/// Version for checkpoint format compatibility
const CHECKPOINT_VERSION: u32 = 1;

/// Checkpoint metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CheckpointMetadata {
    /// Version of checkpoint format
    pub version: u32,
    /// Simulation time
    pub time: f64,
    /// Iteration number
    pub iteration: usize,
    /// Grid dimensions
    pub dimensions: (usize, usize),
    /// Physical domain size
    pub domain_size: (f64, f64),
    /// Timestamp of checkpoint creation
    pub timestamp: u64,
}

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

/// Checkpoint manager for saving and loading simulation states
pub struct CheckpointManager {
    /// Base directory for checkpoints
    checkpoint_dir: PathBuf,
    /// Maximum number of checkpoints to keep
    max_checkpoints: usize,
    /// Enable compression
    use_compression: bool,
}

impl CheckpointManager {
    /// Create a new checkpoint manager
    pub fn new(checkpoint_dir: impl AsRef<Path>) -> Result<Self> {
        let checkpoint_dir = checkpoint_dir.as_ref().to_path_buf();

        // Create directory if it doesn't exist
        std::fs::create_dir_all(&checkpoint_dir).map_err(|e| {
            Error::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to create checkpoint directory: {}", e),
            ))
        })?;

        Ok(Self {
            checkpoint_dir,
            max_checkpoints: 10,
            use_compression: true,
        })
    }

    /// Set maximum number of checkpoints to keep
    pub fn set_max_checkpoints(&mut self, max: usize) {
        self.max_checkpoints = max;
    }

    /// Set compression usage
    pub fn set_compression(&mut self, use_compression: bool) {
        self.use_compression = use_compression;
    }

    /// Save checkpoint to file
    pub fn save<T>(&self, checkpoint: &Checkpoint<T>, name: &str) -> Result<PathBuf>
    where
        T: RealField + Copy + Serialize,
    {
        // Validate checkpoint data
        self.validate_checkpoint(checkpoint)?;

        // Generate filename with timestamp
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_err(|e| {
                Error::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("System time error: {}", e),
                ))
            })?
            .as_secs();

        let filename = format!(
            "{}_iter{:06}_{}.json",
            name, checkpoint.metadata.iteration, timestamp
        );
        let filepath = self.checkpoint_dir.join(&filename);

        // Write checkpoint
        let file = File::create(&filepath)?;
        let writer = BufWriter::new(file);

        if self.use_compression {
            // Use compressed JSON with zstd
            let encoder = zstd::stream::write::Encoder::new(writer, 3).map_err(|e| {
                Error::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Compression error: {}", e),
                ))
            })?;
            serde_json::to_writer(encoder, checkpoint).map_err(|e| {
                Error::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Serialization error: {}", e),
                ))
            })?;
        } else {
            // Regular JSON
            serde_json::to_writer_pretty(writer, checkpoint).map_err(|e| {
                Error::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Serialization error: {}", e),
                ))
            })?;
        }

        // Clean up old checkpoints
        self.cleanup_old_checkpoints(name)?;

        Ok(filepath)
    }

    /// Load checkpoint from file
    pub fn load<T>(&self, filepath: &Path) -> Result<Checkpoint<T>>
    where
        T: RealField + Copy + for<'de> Deserialize<'de>,
    {
        if !filepath.exists() {
            return Err(Error::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("Checkpoint file not found: {:?}", filepath),
            )));
        }

        let file = File::open(filepath)?;
        let reader = BufReader::new(file);

        let checkpoint = if filepath.extension().and_then(|s| s.to_str()) == Some("zst") {
            // Compressed checkpoint
            let decoder = zstd::stream::read::Decoder::new(reader).map_err(|e| {
                Error::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Decompression error: {}", e),
                ))
            })?;
            serde_json::from_reader(decoder).map_err(|e| {
                Error::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Deserialization error: {}", e),
                ))
            })?
        } else {
            // Regular JSON
            serde_json::from_reader(reader).map_err(|e| {
                Error::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Deserialization error: {}", e),
                ))
            })?
        };

        // Validate loaded checkpoint
        self.validate_checkpoint(&checkpoint)?;

        // Check version compatibility
        if checkpoint.metadata.version > CHECKPOINT_VERSION {
            return Err(Error::InvalidConfiguration(format!(
                "Checkpoint version {} is newer than supported version {}",
                checkpoint.metadata.version, CHECKPOINT_VERSION
            )));
        }

        Ok(checkpoint)
    }

    /// Find latest checkpoint for a given simulation name
    pub fn find_latest(&self, name: &str) -> Result<Option<PathBuf>> {
        let entries = std::fs::read_dir(&self.checkpoint_dir)?;

        let mut checkpoints: Vec<_> = entries
            .filter_map(|entry| entry.ok())
            .filter(|entry| {
                entry
                    .file_name()
                    .to_str()
                    .map(|s| s.starts_with(name))
                    .unwrap_or(false)
            })
            .collect();

        if checkpoints.is_empty() {
            return Ok(None);
        }

        // Sort by modification time
        checkpoints.sort_by_key(|entry| {
            entry
                .metadata()
                .and_then(|m| m.modified())
                .unwrap_or(std::time::SystemTime::UNIX_EPOCH)
        });

        Ok(checkpoints.last().map(|e| e.path()))
    }

    /// Validate checkpoint data
    fn validate_checkpoint<T>(&self, checkpoint: &Checkpoint<T>) -> Result<()>
    where
        T: RealField + Copy,
    {
        let (nx, ny) = checkpoint.metadata.dimensions;

        // Check dimensions consistency
        if checkpoint.u_velocity.nrows() != nx || checkpoint.u_velocity.ncols() != ny {
            return Err(Error::InvalidConfiguration(
                "U velocity field dimensions mismatch".to_string(),
            ));
        }

        if checkpoint.v_velocity.nrows() != nx || checkpoint.v_velocity.ncols() != ny {
            return Err(Error::InvalidConfiguration(
                "V velocity field dimensions mismatch".to_string(),
            ));
        }

        if checkpoint.pressure.nrows() != nx || checkpoint.pressure.ncols() != ny {
            return Err(Error::InvalidConfiguration(
                "Pressure field dimensions mismatch".to_string(),
            ));
        }

        // Check for NaN or infinite values
        let has_invalid = |matrix: &DMatrix<T>| matrix.iter().any(|&x| !x.is_finite());

        if has_invalid(&checkpoint.u_velocity) {
            return Err(Error::Numerical(
                cfd_core::error::NumericalErrorKind::InvalidValue {
                    value: "NaN or Inf in velocity field".to_string(),
                },
            ));
        }

        if has_invalid(&checkpoint.v_velocity) {
            return Err(Error::Numerical(
                cfd_core::error::NumericalErrorKind::InvalidValue {
                    value: "NaN or Inf in velocity field".to_string(),
                },
            ));
        }

        if has_invalid(&checkpoint.pressure) {
            return Err(Error::Numerical(
                cfd_core::error::NumericalErrorKind::InvalidValue {
                    value: "NaN or Inf in velocity field".to_string(),
                },
            ));
        }

        Ok(())
    }

    /// Clean up old checkpoints keeping only the most recent
    fn cleanup_old_checkpoints(&self, name: &str) -> Result<()> {
        let entries = std::fs::read_dir(&self.checkpoint_dir)?;

        let mut checkpoints: Vec<_> = entries
            .filter_map(|entry| entry.ok())
            .filter(|entry| {
                entry
                    .file_name()
                    .to_str()
                    .map(|s| s.starts_with(name))
                    .unwrap_or(false)
            })
            .collect();

        if checkpoints.len() <= self.max_checkpoints {
            return Ok(());
        }

        // Sort by modification time
        checkpoints.sort_by_key(|entry| {
            entry
                .metadata()
                .and_then(|m| m.modified())
                .unwrap_or(std::time::SystemTime::UNIX_EPOCH)
        });

        // Remove oldest checkpoints
        let to_remove = checkpoints.len() - self.max_checkpoints;
        for entry in checkpoints.iter().take(to_remove) {
            std::fs::remove_file(entry.path())?;
        }

        Ok(())
    }
}

impl Default for CheckpointManager {
    fn default() -> Self {
        Self::new("checkpoints").expect("Failed to create default checkpoint manager")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_checkpoint_save_load() {
        let temp_dir = tempfile::tempdir().unwrap();
        let manager = CheckpointManager::new(temp_dir.path()).unwrap();

        // Create test checkpoint
        let checkpoint = Checkpoint::<f64> {
            metadata: CheckpointMetadata {
                version: CHECKPOINT_VERSION,
                time: 1.0,
                iteration: 100,
                dimensions: (10, 10),
                domain_size: (1.0, 1.0),
                timestamp: 0,
            },
            u_velocity: DMatrix::zeros(10, 10),
            v_velocity: DMatrix::zeros(10, 10),
            pressure: DMatrix::from_element(10, 10, 1.0),
            temperature: None,
            turbulence_k: None,
            turbulence_epsilon: None,
        };

        // Save checkpoint
        let filepath = manager.save(&checkpoint, "test").unwrap();
        assert!(filepath.exists());

        // Load checkpoint
        let loaded = manager.load::<f64>(&filepath).unwrap();
        assert_eq!(loaded.metadata.iteration, 100);
        assert_eq!(loaded.pressure[(5, 5)], 1.0);
    }

    #[test]
    fn test_find_latest_checkpoint() {
        let temp_dir = tempfile::tempdir().unwrap();
        let manager = CheckpointManager::new(temp_dir.path()).unwrap();

        // No checkpoints initially
        assert!(manager.find_latest("test").unwrap().is_none());

        // Create checkpoint
        let checkpoint = Checkpoint::<f64> {
            metadata: CheckpointMetadata {
                version: CHECKPOINT_VERSION,
                time: 1.0,
                iteration: 100,
                dimensions: (10, 10),
                domain_size: (1.0, 1.0),
                timestamp: 0,
            },
            u_velocity: DMatrix::zeros(10, 10),
            v_velocity: DMatrix::zeros(10, 10),
            pressure: DMatrix::zeros(10, 10),
            temperature: None,
            turbulence_k: None,
            turbulence_epsilon: None,
        };

        manager.save(&checkpoint, "test").unwrap();

        // Should find the checkpoint
        assert!(manager.find_latest("test").unwrap().is_some());
    }
}
