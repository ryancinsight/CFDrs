//! JSON file format support.

use cfd_core::Result;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// JSON configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsonConfig {
    /// Configuration name
    pub name: String,
}

/// JSON file writer
pub struct JsonWriter<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> JsonWriter<T> {
    /// Create a new JSON writer
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField> Default for JsonWriter<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// JSON file reader
pub struct JsonReader<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> JsonReader<T> {
    /// Create a new JSON reader
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField> Default for JsonReader<T> {
    fn default() -> Self {
        Self::new()
    }
}