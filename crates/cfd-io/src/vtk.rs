//! VTK file format support.

use cfd_core::Result;
use nalgebra::RealField;

/// VTK file writer
pub struct VtkWriter<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> VtkWriter<T> {
    /// Create a new VTK writer
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Write to file
    pub fn write(&self, _path: &std::path::Path) -> Result<()> {
        Ok(())
    }
}

impl<T: RealField> Default for VtkWriter<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// VTK file reader
pub struct VtkReader<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> VtkReader<T> {
    /// Create a new VTK reader
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Read from file
    pub fn read(&self, _path: &std::path::Path) -> Result<()> {
        Ok(())
    }
}

impl<T: RealField> Default for VtkReader<T> {
    fn default() -> Self {
        Self::new()
    }
}