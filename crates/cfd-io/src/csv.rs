//! CSV file I/O operations.
//!
//! This module provides CSV file support for CFD simulation data
//! with streaming and iterator-based processing for efficiency.

use cfd_core::error::{Error, Result};
use csv::{Reader as CsvReaderImpl, Writer as CsvWriterImpl};
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

/// CSV writer for simulation data
pub struct CsvWriter<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> CsvWriter<T> {
    /// Create a new CSV writer
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Write time series data to CSV
    pub fn write_time_series<I>(&self, path: &Path, headers: &[&str], data: I) -> Result<()>
    where
        I: IntoIterator<Item = Vec<T>>,
        T: std::fmt::Display,
    {
        let file = File::create(path)?;
        let mut writer = CsvWriterImpl::from_writer(BufWriter::new(file));

        // Write headers
        writer
            .write_record(headers)
            .map_err(|e| Error::CsvError(e.to_string()))?;

        // Write data rows using iterator
        for row in data {
            let string_row: Vec<String> =
                row.iter().map(std::string::ToString::to_string).collect();
            writer
                .write_record(&string_row)
                .map_err(|e| Error::CsvError(e.to_string()))?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Write structured data using serde
    pub fn write_records<I, R>(&self, path: &Path, records: I) -> Result<()>
    where
        I: IntoIterator<Item = R>,
        R: Serialize,
    {
        let file = File::create(path)?;
        let mut writer = CsvWriterImpl::from_writer(BufWriter::new(file));

        for record in records {
            writer
                .serialize(record)
                .map_err(|e| Error::CsvError(e.to_string()))?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Create a streaming writer for large datasets
    pub fn create_stream_writer(&self, path: &Path) -> Result<StreamWriter> {
        let file = File::create(path)?;
        let writer = CsvWriterImpl::from_writer(BufWriter::new(file));
        Ok(StreamWriter { writer })
    }
}

impl<T: RealField + Copy> Default for CsvWriter<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Streaming CSV writer for large datasets
pub struct StreamWriter {
    writer: CsvWriterImpl<BufWriter<File>>,
}

impl StreamWriter {
    /// Write headers
    pub fn write_headers(&mut self, headers: &[&str]) -> Result<()> {
        self.writer
            .write_record(headers)
            .map_err(|e| Error::CsvError(e.to_string()))?;
        Ok(())
    }

    /// Write a single row
    pub fn write_row<T: std::fmt::Display>(&mut self, row: &[T]) -> Result<()> {
        let string_row: Vec<String> = row.iter().map(std::string::ToString::to_string).collect();
        self.writer
            .write_record(&string_row)
            .map_err(|e| Error::CsvError(e.to_string()))?;
        Ok(())
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

/// CSV reader for simulation data
pub struct CsvReader<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> CsvReader<T> {
    /// Create a new CSV reader
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Read time series data from CSV
    pub fn read_time_series(&self, path: &Path) -> Result<TimeSeriesData<T>>
    where
        T: std::str::FromStr,
        T::Err: std::fmt::Display,
    {
        let file = File::open(path)?;
        let mut reader = CsvReaderImpl::from_reader(BufReader::new(file));

        // Read headers
        let headers = reader
            .headers()
            .map_err(|e| Error::CsvError(e.to_string()))?;
        let header_vec: Vec<String> = headers
            .iter()
            .map(std::string::ToString::to_string)
            .collect();

        // Read data using iterator
        let mut data = Vec::new();
        for result in reader.records() {
            let record = result.map_err(|e| Error::CsvError(e.to_string()))?;
            let row: Result<Vec<T>> = record
                .iter()
                .map(|field| {
                    field
                        .parse::<T>()
                        .map_err(|e| Error::SerializationError(format!("Parse error: {e}")))
                })
                .collect();
            data.push(row?);
        }

        Ok(TimeSeriesData {
            headers: header_vec,
            data,
        })
    }

    /// Read structured records using serde
    pub fn read_records<R>(&self, path: &Path) -> Result<Vec<R>>
    where
        R: for<'de> Deserialize<'de>,
    {
        let file = File::open(path)?;
        let mut reader = CsvReaderImpl::from_reader(BufReader::new(file));

        let mut records = Vec::new();
        for result in reader.deserialize() {
            let record: R = result.map_err(|e| Error::CsvError(e.to_string()))?;
            records.push(record);
        }

        Ok(records)
    }

    /// Create a streaming reader for large datasets
    pub fn create_stream_reader(&self, path: &Path) -> Result<StreamReader> {
        let file = File::open(path)?;
        let reader = CsvReaderImpl::from_reader(BufReader::new(file));
        Ok(StreamReader { reader })
    }
}

impl<T: RealField + Copy> Default for CsvReader<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Streaming CSV reader for large datasets
pub struct StreamReader {
    reader: CsvReaderImpl<BufReader<File>>,
}

impl StreamReader {
    /// Get headers
    pub fn headers(&mut self) -> Result<Vec<String>> {
        let headers = self
            .reader
            .headers()
            .map_err(|e| Error::CsvError(e.to_string()))?;
        Ok(headers
            .iter()
            .map(std::string::ToString::to_string)
            .collect())
    }

    /// Iterator over records
    pub fn records(&mut self) -> impl Iterator<Item = Result<csv::StringRecord>> + '_ {
        self.reader
            .records()
            .map(|r| r.map_err(|e| Error::CsvError(e.to_string())))
    }

    /// Iterator over typed records
    pub fn typed_records<R>(&mut self) -> impl Iterator<Item = Result<R>> + '_
    where
        R: for<'de> Deserialize<'de> + 'static,
    {
        self.reader
            .deserialize()
            .map(|r| r.map_err(|e| Error::CsvError(e.to_string())))
    }
}

/// Time series data structure
#[derive(Debug, Clone)]
pub struct TimeSeriesData<T> {
    /// Column headers
    pub headers: Vec<String>,
    /// Data rows
    pub data: Vec<Vec<T>>,
}

impl<T: RealField + Copy> TimeSeriesData<T> {
    /// Get number of rows
    #[must_use]
    pub fn num_rows(&self) -> usize {
        self.data.len()
    }

    /// Get number of columns
    #[must_use]
    pub fn num_cols(&self) -> usize {
        self.headers.len()
    }

    /// Get column by name
    #[must_use]
    pub fn get_column(&self, name: &str) -> Option<Vec<T>>
    where
        T: Clone,
    {
        let index = self.headers.iter().position(|h| h == name)?;
        Some(self.data.iter().map(|row| row[index]).collect())
    }

    /// Get column by index
    #[must_use]
    pub fn get_column_by_index(&self, index: usize) -> Option<Vec<T>>
    where
        T: Clone,
    {
        if index >= self.num_cols() {
            return None;
        }
        Some(self.data.iter().map(|row| row[index]).collect())
    }

    /// Iterator over rows
    pub fn rows(&self) -> impl Iterator<Item = &Vec<T>> {
        self.data.iter()
    }

    /// Iterator over columns
    pub fn columns(&self) -> impl Iterator<Item = Vec<&T>> + '_ {
        (0..self.num_cols()).map(move |i| self.data.iter().map(|row| &row[i]).collect())
    }
}

/// CSV configuration for simulation output
#[derive(Debug, Clone)]
pub struct CsvConfig {
    /// Delimiter character
    pub delimiter: u8,
    /// Whether to quote fields
    pub quote_style: csv::QuoteStyle,
    /// Whether to write headers
    pub has_headers: bool,
}

impl Default for CsvConfig {
    fn default() -> Self {
        Self {
            delimiter: b',',
            quote_style: csv::QuoteStyle::Necessary,
            has_headers: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_csv_write_read() {
        let test_file = NamedTempFile::new().expect("Failed to create test file");
        let path = test_file.path();

        // Write data
        let writer = CsvWriter::<f64>::new();
        let headers = vec!["time", "value1", "value2"];
        let data = vec![
            vec![0.0, 1.0, 2.0],
            vec![0.1, 1.1, 2.1],
            vec![0.2, 1.2, 2.2],
        ];
        writer
            .write_time_series(path, &headers, data.clone())
            .expect("Failed to write time series data");

        // Read data
        let reader = CsvReader::<f64>::new();
        let result = reader
            .read_time_series(path)
            .expect("Failed to read time series data");

        assert_eq!(result.headers, headers);
        assert_eq!(result.data, data);
    }

    #[test]
    fn test_streaming() {
        let test_file = NamedTempFile::new().expect("Failed to create test file");
        let path = test_file.path();

        // Write using streaming
        let writer = CsvWriter::<f64>::new();
        let mut stream = writer
            .create_stream_writer(path)
            .expect("CRITICAL: Add proper error handling");
        stream
            .write_headers(&["x", "y"])
            .expect("CRITICAL: Add proper error handling");
        stream
            .write_row(&[1.0, 2.0])
            .expect("CRITICAL: Add proper error handling");
        stream
            .write_row(&[3.0, 4.0])
            .expect("CRITICAL: Add proper error handling");
        stream.flush().expect("CRITICAL: Add proper error handling");

        // Read using streaming
        let reader = CsvReader::<f64>::new();
        let mut stream = reader
            .create_stream_reader(path)
            .expect("CRITICAL: Add proper error handling");
        let headers = stream
            .headers()
            .expect("CRITICAL: Add proper error handling");
        assert_eq!(headers, vec!["x", "y"]);

        let records: Vec<_> = stream.records().collect();
        assert_eq!(records.len(), 2);
    }
}
