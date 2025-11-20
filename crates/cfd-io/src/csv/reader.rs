//! CSV reading functionality

use cfd_core::error::{Error, Result};
use csv::Reader as CsvReaderImpl;
use nalgebra::RealField;
use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::str::FromStr;

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
    pub fn read_time_series(&self, path: &Path) -> Result<(Vec<String>, Vec<Vec<T>>)>
    where
        T: FromStr,
        <T as FromStr>::Err: std::fmt::Display,
    {
        let file = File::open(path)?;
        let mut reader = CsvReaderImpl::from_reader(BufReader::new(file));

        // Read headers
        let headers = reader
            .headers()
            .map_err(|e| Error::Io(std::io::Error::other(format!("CSV error: {e}"))))?
            .iter()
            .map(String::from)
            .collect();

        // Read data rows
        let mut data = Vec::new();
        for result in reader.records() {
            let record =
                result.map_err(|e| Error::Io(std::io::Error::other(format!("CSV error: {e}"))))?;

            let row: Result<Vec<T>> = record
                .iter()
                .map(|field| {
                    field.parse::<T>().map_err(|e| {
                        Error::Io(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!("Parse error: {e}"),
                        ))
                    })
                })
                .collect();

            data.push(row?);
        }

        Ok((headers, data))
    }

    /// Read field data from CSV
    pub fn read_field_data(&self, path: &Path) -> Result<(Vec<T>, Vec<T>, Vec<T>)>
    where
        T: FromStr,
        <T as FromStr>::Err: std::fmt::Display,
    {
        let (_headers, data) = self.read_time_series(path)?;

        if data.is_empty() {
            return Err(Error::InvalidInput("Empty CSV file".to_string()));
        }

        if data[0].len() < 3 {
            return Err(Error::InvalidInput(
                "CSV must have at least 3 columns".to_string(),
            ));
        }

        let mut x = Vec::with_capacity(data.len());
        let mut y = Vec::with_capacity(data.len());
        let mut field = Vec::with_capacity(data.len());

        for row in data {
            x.push(row[0]);
            y.push(row[1]);
            field.push(row[2]);
        }

        Ok((x, y, field))
    }

    /// Create an iterator over CSV records
    pub fn iter_records<'a, R: Deserialize<'a>>(&self, path: &Path) -> Result<CsvIterator<'a, R>> {
        CsvIterator::new(path)
    }
}

impl<T: RealField + Copy> Default for CsvReader<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Iterator over CSV records
pub struct CsvIterator<'a, R: Deserialize<'a>> {
    reader: CsvReaderImpl<BufReader<File>>,
    _phantom: std::marker::PhantomData<&'a R>,
}

impl<'a, R: Deserialize<'a>> CsvIterator<'a, R> {
    /// Create a new CSV iterator
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let reader = CsvReaderImpl::from_reader(BufReader::new(file));
        Ok(Self {
            reader,
            _phantom: std::marker::PhantomData,
        })
    }
}

impl<R: for<'de> Deserialize<'de>> Iterator for CsvIterator<'_, R> {
    type Item = Result<R>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.deserialize().next().map(|result| {
            result.map_err(|e| {
                Error::Io(std::io::Error::other(format!(
                    "CSV deserialization error: {e}"
                )))
            })
        })
    }
}

/// Convenience function to read field data
pub fn read_field_data<T>(path: &Path) -> Result<(Vec<T>, Vec<T>, Vec<T>)>
where
    T: RealField + Copy + FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    let reader = CsvReader::<T>::new();
    reader.read_field_data(path)
}
