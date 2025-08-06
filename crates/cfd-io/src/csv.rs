//! CSV file format support.

use cfd_core::Result;
use nalgebra::RealField;

/// CSV file writer
pub struct CsvWriter<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> CsvWriter<T> {
    /// Create a new CSV writer
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField> Default for CsvWriter<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// CSV file reader
pub struct CsvReader<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> CsvReader<T> {
    /// Create a new CSV reader
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField> Default for CsvReader<T> {
    fn default() -> Self {
        Self::new()
    }
}