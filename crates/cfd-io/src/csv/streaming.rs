//! Streaming CSV I/O for large datasets

use cfd_core::error::{Error, Result};
use csv::{Reader as CsvReaderImpl, Writer as CsvWriterImpl};
use nalgebra::RealField;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

/// Streaming CSV reader for memory-efficient processing
pub struct StreamingReader<T: RealField + Copy> {
    reader: Option<CsvReaderImpl<BufReader<File>>>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> StreamingReader<T> {
    /// Open a CSV file for streaming
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let reader = CsvReaderImpl::from_reader(BufReader::new(file));

        Ok(Self {
            reader: Some(reader),
            _phantom: std::marker::PhantomData,
        })
    }

    /// Read next row from the stream
    pub fn read_row(&mut self) -> Result<Option<Vec<T>>>
    where
        T: FromStr,
        <T as FromStr>::Err: std::fmt::Display,
    {
        let reader = self
            .reader
            .as_mut()
            .ok_or_else(|| Error::InvalidInput("Reader is closed".to_string()))?;

        match reader.records().next() {
            Some(result) => {
                let record = result.map_err(|e| {
                    Error::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("CSV error: {}", e),
                    ))
                })?;

                let row: Result<Vec<T>> = record
                    .iter()
                    .map(|field| {
                        field.parse::<T>().map_err(|e| {
                            Error::Io(std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                format!("Parse error: {}", e),
                            ))
                        })
                    })
                    .collect();

                Ok(Some(row?))
            }
            None => Ok(None),
        }
    }

    /// Process rows with a callback function
    pub fn process_rows<F>(&mut self, mut processor: F) -> Result<usize>
    where
        T: FromStr,
        <T as FromStr>::Err: std::fmt::Display,
        F: FnMut(Vec<T>) -> Result<()>,
    {
        let mut count = 0;

        while let Some(row) = self.read_row()? {
            processor(row)?;
            count += 1;
        }

        Ok(count)
    }

    /// Close the reader
    pub fn close(&mut self) {
        self.reader = None;
    }
}

/// Streaming CSV writer for memory-efficient output
pub struct StreamingWriter<T: RealField + Copy> {
    writer: Option<CsvWriterImpl<BufWriter<File>>>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> StreamingWriter<T> {
    /// Create a new streaming CSV writer
    pub fn create(path: &Path, headers: &[&str]) -> Result<Self> {
        let file = File::create(path)?;
        let mut writer = CsvWriterImpl::from_writer(BufWriter::new(file));

        // Write headers
        writer.write_record(headers).map_err(|e| {
            Error::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("CSV error: {}", e),
            ))
        })?;

        Ok(Self {
            writer: Some(writer),
            _phantom: std::marker::PhantomData,
        })
    }

    /// Write a single row
    pub fn write_row(&mut self, row: &[T]) -> Result<()>
    where
        T: std::fmt::Display,
    {
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| Error::InvalidInput("Writer is closed".to_string()))?;

        let string_row: Vec<String> = row.iter().map(|v| v.to_string()).collect();

        writer.write_record(&string_row).map_err(|e| {
            Error::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("CSV error: {}", e),
            ))
        })?;

        Ok(())
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<()> {
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| Error::InvalidInput("Writer is closed".to_string()))?;

        writer.flush().map_err(|e| {
            Error::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("CSV flush error: {}", e),
            ))
        })?;

        Ok(())
    }

    /// Close the writer
    pub fn close(&mut self) -> Result<()> {
        self.flush()?;
        self.writer = None;
        Ok(())
    }
}

impl<T: RealField + Copy> Drop for StreamingWriter<T> {
    fn drop(&mut self) {
        // Best effort flush on drop
        let _ = self.flush();
    }
}
