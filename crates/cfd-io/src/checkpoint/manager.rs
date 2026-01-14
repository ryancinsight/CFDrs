//! Checkpoint management and I/O operations

use crate::checkpoint::{Checkpoint, CompressionStrategy};
use cfd_core::error::{Error, Result};
use std::fs::{self, File};

use std::io::{BufReader, BufWriter, ErrorKind, Read, Write};
use std::path::{Path, PathBuf};

/// Magic number for checkpoint files
const CHECKPOINT_MAGIC: &[u8; 8] = b"CFDCHKPT";

/// Compression type identifiers
const COMPRESSION_NONE: u8 = 0;
const COMPRESSION_ZSTD: u8 = 1;

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

    /// Save a checkpoint with self-describing file format
    ///
    /// # Invariants
    /// * Bit-exact preservation: `checksum(stored_data) == metadata.checksum`
    pub fn save(&self, checkpoint: &Checkpoint<f64>) -> Result<PathBuf> {
        let filename = format!(
            "checkpoint_iter_{:08}_tbits_{:016x}.bin",
            checkpoint.metadata.iteration,
            checkpoint.metadata.time.to_bits()
        );
        let path = self.base_dir.join(filename);

        let mut cp = checkpoint.clone();
        cp.metadata.checksum = cp.compute_checksum();
        cp.metadata.validate()?;

        let encoded = bincode::serde::encode_to_vec(&cp, bincode::config::standard())
            .map_err(|e| Error::Serialization(format!("bincode encode error: {e}")))?;

        let file = File::create(&path)?;
        let mut writer = BufWriter::new(file);

        // Write self-describing header
        writer.write_all(CHECKPOINT_MAGIC)?;
        let compression_type = match self.compression {
            CompressionStrategy::None => COMPRESSION_NONE,
            CompressionStrategy::Zstd(_) => COMPRESSION_ZSTD,
        };
        writer.write_all(&[compression_type])?;

        // Write checkpoint data with appropriate compression
        match self.compression {
            CompressionStrategy::None => {
                writer.write_all(&encoded)?;
            }
            CompressionStrategy::Zstd(level) => {
                let mut encoder =
                    zstd::stream::Encoder::new(writer, i32::from(level)).map_err(|e| {
                        Error::Io(std::io::Error::other(format!("Compression error: {e}")))
                    })?;
                encoder.write_all(&encoded)?;
                encoder.finish()?;
            }
        }

        // Clean up old checkpoints if needed
        if let Some(max) = self.max_checkpoints {
            self.cleanup_old_checkpoints(max)?;
        }

        Ok(path)
    }

    /// Load a checkpoint from self-describing file format
    ///
    /// # Invariants
    /// * Bit-exact verification: `compute_checksum(loaded) == metadata.checksum`
    pub fn load(&self, path: impl AsRef<Path>) -> Result<Checkpoint<f64>> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read and verify magic number
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != CHECKPOINT_MAGIC {
            return Err(Error::Io(std::io::Error::new(
                ErrorKind::InvalidData,
                "Not a valid checkpoint file: incorrect magic number",
            )));
        }

        // Read compression type
        let mut compression_type = [0u8; 1];
        reader.read_exact(&mut compression_type)?;

        let mut buffer = Vec::new();
        reader.read_to_end(&mut buffer)?;

        let data = match compression_type[0] {
            COMPRESSION_NONE => buffer,
            COMPRESSION_ZSTD => zstd::decode_all(&buffer[..]).map_err(|e| {
                Error::Io(std::io::Error::other(format!(
                    "Zstd decompression error: {e}"
                )))
            })?,
            _ => {
                return Err(Error::Io(std::io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Unsupported compression type: {}", compression_type[0]),
                )));
            }
        };

        let (checkpoint, _read_bytes): (Checkpoint<f64>, usize) =
            bincode::serde::decode_from_slice(&data, bincode::config::standard()).map_err(|e| {
                Error::Io(std::io::Error::other(format!(
                    "Bincode deserialization error: {e}"
                )))
            })?;

        checkpoint.metadata.validate()?;

        let recomputed = checkpoint.compute_checksum();
        if checkpoint.metadata.checksum != recomputed {
            return Err(Error::Io(std::io::Error::new(
                ErrorKind::InvalidData,
                "Checkpoint checksum mismatch: data corrupted or tampered",
            )));
        }

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

// REMOVED: These convenience functions were unsafe and misleading.
// They created a default CheckpointManager with no compression,
// which would fail to load compressed checkpoints.
// Users must explicitly create and configure a CheckpointManager.
