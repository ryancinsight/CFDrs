//! Binary file I/O operations with zero-copy abstractions.
//!
//! This module provides efficient binary serialization/deserialization
//! using iterator-based streaming and zero-copy operations.

use crate::error::{Error, Result};
use crate::leto_arrays::try_for_each_row_major;
use eunomia::RealField;
use leto::{Array1, Array2};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Binary writer with streaming capabilities and zero-copy operations
pub struct BinaryWriter<W: Write> {
    writer: BufWriter<W>,
}

impl<W: Write> BinaryWriter<W> {
    /// Create a new binary writer
    pub fn new(writer: W) -> Self {
        Self {
            writer: BufWriter::new(writer),
        }
    }

    /// Write serializable data using serde
    pub fn write<T: Serialize>(&mut self, data: &T) -> Result<()> {
        let encoded = serde_json::to_vec(data)
            .map_err(|e| Error::Io(std::io::Error::other(format!("Serialization error: {e}"))))?;
        self.writer
            .write_all(&encoded)
            .map_err(|e| Error::Io(std::io::Error::other(format!("Binary write error: {e}"))))?;
        Ok(())
    }

    /// Write vector data using iterator-based streaming
    pub fn write_vector<T: RealField + Serialize>(&mut self, vector: &Array1<T>) -> Result<()> {
        // Write length first
        self.write(&vector.size())?;

        // Stream vector data using iterator
        let [len] = vector.shape();
        for index in 0..len {
            let value = vector
                .get([index])
                .expect("invariant: generated Leto vector index is in bounds");
            self.write(value)?;
        }
        Ok(())
    }

    /// Write matrix data with zero-copy slicing
    pub fn write_matrix<T: RealField + Serialize>(&mut self, matrix: &Array2<T>) -> Result<()> {
        // Write dimensions
        self.write(&matrix.shape())?;

        // Stream matrix data in Leto's logical row-major order.
        try_for_each_row_major(matrix, |value| self.write(&value))
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().map_err(Error::Io)
    }
}

impl BinaryWriter<File> {
    /// Create a binary writer for a file
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path).map_err(Error::Io)?;
        Ok(Self::new(file))
    }
}

/// Binary reader with streaming capabilities and zero-copy operations
pub struct BinaryReader<R: Read> {
    reader: BufReader<R>,
}

impl<R: Read> BinaryReader<R> {
    /// Create a new binary reader
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
        }
    }

    /// Read deserializable data using serde
    pub fn read<T: for<'de> Deserialize<'de>>(&mut self) -> Result<T> {
        let mut buffer = Vec::new();
        self.reader
            .read_to_end(&mut buffer)
            .map_err(|e| Error::Io(std::io::Error::other(format!("Binary read error: {e}"))))?;
        serde_json::from_slice(&buffer)
            .map_err(|e| Error::Io(std::io::Error::other(format!("Deserialization error: {e}"))))
    }

    /// Read vector data using iterator-based construction
    pub fn read_vector<T: RealField + for<'de> Deserialize<'de>>(&mut self) -> Result<Array1<T>> {
        let len: usize = self.read()?;

        // Use iterator to collect vector elements efficiently
        let data: Result<Vec<T>> = (0..len).map(|_| self.read()).collect();

        Array1::from_shape_vec([len], data?).map_err(|error| {
            Error::InvalidInput(format!("Invalid Leto vector payload shape: {error}"))
        })
    }

    /// Read matrix data with efficient allocation
    pub fn read_matrix<T: RealField + for<'de> Deserialize<'de>>(&mut self) -> Result<Array2<T>> {
        let shape: [usize; 2] = self.read()?;

        // Use iterator to collect matrix elements efficiently
        let element_count = shape[0]
            .checked_mul(shape[1])
            .ok_or_else(|| Error::InvalidInput("Matrix shape overflows usize".to_string()))?;
        let data: Result<Vec<T>> = (0..element_count).map(|_| self.read()).collect();

        Array2::from_shape_vec(shape, data?).map_err(|error| {
            Error::InvalidInput(format!("Invalid Leto matrix payload shape: {error}"))
        })
    }
}

impl BinaryReader<File> {
    /// Create a binary reader for a file
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path).map_err(Error::Io)?;
        Ok(Self::new(file))
    }
}
