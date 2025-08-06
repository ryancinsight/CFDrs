//! Binary file format support.

use cfd_core::Result;
use nalgebra::RealField;

/// Binary file writer
pub struct BinaryWriter<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> BinaryWriter<T> {
    /// Create a new binary writer
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField> Default for BinaryWriter<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Binary file reader
pub struct BinaryReader<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> BinaryReader<T> {
    /// Create a new binary reader
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField> Default for BinaryReader<T> {
    fn default() -> Self {
        Self::new()
    }
}