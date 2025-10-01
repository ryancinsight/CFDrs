//! CSV writing functionality

use cfd_core::error::{Error, Result};
use csv::Writer as CsvWriterImpl;
use nalgebra::RealField;
use std::fs::File;
use std::io::BufWriter;
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
        writer.write_record(headers).map_err(|e| {
            Error::Io(std::io::Error::other(
                format!("CSV error: {e}"),
            ))
        })?;

        // Write data rows using iterator
        for row in data {
            let string_row: Vec<String> =
                row.iter().map(std::string::ToString::to_string).collect();
            writer.write_record(&string_row).map_err(|e| {
                Error::Io(std::io::Error::other(
                    format!("CSV error: {e}"),
                ))
            })?;
        }

        writer.flush().map_err(|e| {
            Error::Io(std::io::Error::other(
                format!("CSV flush error: {e}"),
            ))
        })?;

        Ok(())
    }

    /// Write field data to CSV
    pub fn write_field_data(
        &self,
        path: &Path,
        x: &[T],
        y: &[T],
        field: &[T],
        field_name: &str,
    ) -> Result<()>
    where
        T: std::fmt::Display,
    {
        if x.len() != y.len() || x.len() != field.len() {
            return Err(Error::InvalidInput("Mismatched array lengths".to_string()));
        }

        let headers = vec!["x", "y", field_name];
        let data = x
            .iter()
            .zip(y.iter())
            .zip(field.iter())
            .map(|((x_val, y_val), f_val)| {
                vec![x_val.to_string(), y_val.to_string(), f_val.to_string()]
            });

        let file = File::create(path)?;
        let mut writer = CsvWriterImpl::from_writer(BufWriter::new(file));

        writer.write_record(&headers).map_err(|e| {
            Error::Io(std::io::Error::other(
                format!("CSV error: {e}"),
            ))
        })?;

        for row in data {
            writer.write_record(&row).map_err(|e| {
                Error::Io(std::io::Error::other(
                    format!("CSV error: {e}"),
                ))
            })?;
        }

        Ok(())
    }
}

impl<T: RealField + Copy> Default for CsvWriter<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Convenience function to write time series data
pub fn write_time_series<T, I>(path: &Path, headers: &[&str], data: I) -> Result<()>
where
    T: RealField + Copy + std::fmt::Display,
    I: IntoIterator<Item = Vec<T>>,
{
    let writer = CsvWriter::<T>::new();
    writer.write_time_series(path, headers, data)
}
