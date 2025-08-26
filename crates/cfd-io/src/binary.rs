//! Binary file I/O operations with zero-copy abstractions.
//!
//! This module provides efficient binary serialization/deserialization
//! using iterator-based streaming and zero-copy operations.

use cfd_core::error::{Error, Result};
use nalgebra::{RealField, DMatrix, DVector};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use bincode;
use serde::{Serialize, Deserialize};
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
    /// Write serializable data using bincode
    pub fn write<T: Serialize>(&mut self, data: &T) -> Result<()> {
        bincode::serialize_into(&mut self.writer, data)
            .map_err(|e| Error::SerializationError(format!("Binary write error: {e}")))?;
        Ok(())
    /// Write vector data using iterator-based streaming
    pub fn write_vector<T: RealField + Copy + Serialize>(&mut self, vector: &DVector<T>) -> Result<()> {
        // Write length first
        self.write(&vector.len())?;
        // Stream vector data using iterator
        vector.iter()
            .try_for_each(|value| self.write(value))
    /// Write matrix data with zero-copy slicing
    pub fn write_matrix<T: RealField + Copy + Serialize>(&mut self, matrix: &DMatrix<T>) -> Result<()> {
        // Write dimensions
        self.write(&(matrix.nrows(), matrix.ncols()))?;
        // Stream matrix data row by row using iterators
        matrix.row_iter()
            .try_for_each(|row| {
                row.iter().try_for_each(|value| self.write(value))
            })
    /// Flush the writer
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()
            .map_err(Error::Io)
impl BinaryWriter<File> {
    /// Create a binary writer for a file
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path)
            .map_err(Error::Io)?;
        Ok(Self::new(file))
/// Binary reader with streaming capabilities and zero-copy operations
pub struct BinaryReader<R: Read> {
    reader: BufReader<R>,
impl<R: Read> BinaryReader<R> {
    /// Create a new binary reader
    pub fn new(reader: R) -> Self {
            reader: BufReader::new(reader),
    /// Read deserializable data using bincode
    pub fn read<T: for<'de> Deserialize<'de>>(&mut self) -> Result<T> {
        bincode::deserialize_from(&mut self.reader)
            .map_err(|e| Error::SerializationError(format!("Binary read error: {e}")))
    /// Read vector data using iterator-based construction
    pub fn read_vector<T: RealField + Copy + for<'de> Deserialize<'de>>(&mut self) -> Result<DVector<T>> {
        let len: usize = self.read()?;
        // Use iterator to collect vector elements efficiently
        let data: Result<Vec<T>> = (0..len)
            .map(|_| self.read())
            .collect();
        Ok(DVector::from_vec(data?))
    /// Read matrix data with efficient allocation
    pub fn read_matrix<T: RealField + Copy + for<'de> Deserialize<'de>>(&mut self) -> Result<DMatrix<T>> {
        let (nrows, ncols): (usize, usize) = self.read()?;
        // Use iterator to collect matrix elements efficiently
        let data: Result<Vec<T>> = (0..nrows * ncols)
        Ok(DMatrix::from_vec(nrows, ncols, data?))
impl BinaryReader<File> {
    /// Create a binary reader for a file
        let file = File::open(path)
