//! Checkpoint management and I/O operations

use crate::checkpoint::{Checkpoint, CompressionStrategy};
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

/// Checkpoint manager for handling save/load operations
pub struct CheckpointManager {
    /// Base directory for checkpoints
    base_dir: PathBuf,
    /// Compression strategy
    compression: CompressionStrategy,
    /// Maximum number of checkpoints to keep
    max_checkpoints: Option<usize>,
}

impl CheckpointManager {
    /// Create a new checkpoint manager
    pub fn new(base_dir: impl AsRef<Path>) -> Result<Self> {
        let base_dir = base_dir.as_ref().to_path_buf();

        // Create directory if it doesn't exist
        fs::create_dir_all(&base_dir)?;

        Ok(Self {
            base_dir,
            compression: CompressionStrategy::None,
            max_checkpoints: None,
        })
    }

    /// Set compression strategy
    pub fn set_compression(&mut self, compression: CompressionStrategy) {
        self.compression = compression;
    }

    /// Set maximum number of checkpoints to keep
    pub fn set_max_checkpoints(&mut self, max: usize) {
        self.max_checkpoints = Some(max);
    }

    /// Save a checkpoint
    pub fn save<T>(&self, checkpoint: &Checkpoint<T>) -> Result<PathBuf>
    where
        T: RealField + Copy + Serialize,
        for<'de> T: Deserialize<'de>,
    {
        let filename = format!(
            "checkpoint_iter_{:08}_t_{:.6}.bin",
            checkpoint.metadata.iteration, checkpoint.metadata.time
        );
        let path = self.base_dir.join(filename);

        let file = File::create(&path)?;
        let writer = BufWriter::new(file);

        match self.compression {
            CompressionStrategy::None => {
                bincode::serialize_into(writer, checkpoint).map_err(|e| {
                    Error::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Serialization error: {}", e),
                    ))
                })?;
            }
            CompressionStrategy::Zstd(level) => {
                let encoder = zstd::stream::Encoder::new(writer, level as i32).map_err(|e| {
                    Error::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Compression error: {}", e),
                    ))
                })?;

                bincode::serialize_into(encoder.auto_finish(), checkpoint).map_err(|e| {
                    Error::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Serialization error: {}", e),
                    ))
                })?;
            }
        }

        // Clean up old checkpoints if needed
        if let Some(max) = self.max_checkpoints {
            self.cleanup_old_checkpoints(max)?;
        }

        Ok(path)
    }

    /// Load a checkpoint
    pub fn load<T>(&self, path: impl AsRef<Path>) -> Result<Checkpoint<T>>
    where
        T: RealField + Copy + for<'de> Deserialize<'de>,
    {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let checkpoint = match self.compression {
            CompressionStrategy::None => bincode::deserialize_from(reader).map_err(|e| {
                Error::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Deserialization error: {}", e),
                ))
            })?,
            CompressionStrategy::Zstd(_) => {
                let decoder = zstd::stream::Decoder::new(reader).map_err(|e| {
                    Error::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Decompression error: {}", e),
                    ))
                })?;

                bincode::deserialize_from(decoder).map_err(|e| {
                    Error::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Deserialization error: {}", e),
                    ))
                })?
            }
        };

        Ok(checkpoint)
    }

    /// Find the latest checkpoint
    pub fn find_latest(&self) -> Result<Option<PathBuf>> {
        let mut checkpoints = Vec::new();

        for entry in fs::read_dir(&self.base_dir)? {
            let entry = entry?;
            let path = entry.path();

            if path.extension().and_then(|s| s.to_str()) == Some("bin") {
                if let Some(name) = path.file_stem().and_then(|s| s.to_str()) {
                    if name.starts_with("checkpoint_") {
                        checkpoints.push(path);
                    }
                }
            }
        }

        checkpoints.sort();
        Ok(checkpoints.into_iter().last())
    }

    /// Clean up old checkpoints
    fn cleanup_old_checkpoints(&self, max_to_keep: usize) -> Result<()> {
        let mut checkpoints = Vec::new();

        for entry in fs::read_dir(&self.base_dir)? {
            let entry = entry?;
            let path = entry.path();

            if path.extension().and_then(|s| s.to_str()) == Some("bin") {
                if let Some(name) = path.file_stem().and_then(|s| s.to_str()) {
                    if name.starts_with("checkpoint_") {
                        checkpoints.push(path);
                    }
                }
            }
        }

        checkpoints.sort();

        // Remove oldest checkpoints if we exceed the limit
        while checkpoints.len() > max_to_keep {
            let path = checkpoints.remove(0);
            fs::remove_file(path)?;
        }

        Ok(())
    }
}

/// Convenience function to save a checkpoint
pub fn save_checkpoint<T>(path: impl AsRef<Path>, checkpoint: &Checkpoint<T>) -> Result<()>
where
    T: RealField + Copy + Serialize,
    for<'de> T: Deserialize<'de>,
{
    let manager = CheckpointManager::new(path.as_ref().parent().unwrap_or(Path::new(".")))?;
    manager.save(checkpoint)?;
    Ok(())
}

/// Convenience function to load a checkpoint
pub fn load_checkpoint<T>(path: impl AsRef<Path>) -> Result<Checkpoint<T>>
where
    T: RealField + Copy + for<'de> Deserialize<'de>,
{
    let manager = CheckpointManager::new(path.as_ref().parent().unwrap_or(Path::new(".")))?;
    manager.load(path)
}
