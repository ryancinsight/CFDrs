//! CSV data types and structures

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Generic CSV record
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CsvRecord<T: RealField + Copy> {
    /// Time or iteration number
    pub time: T,
    /// Data values
    pub values: Vec<T>,
}

/// Field data structure for 2D/3D fields
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldData<T: RealField + Copy> {
    /// X coordinates
    pub x: Vec<T>,
    /// Y coordinates
    pub y: Vec<T>,
    /// Z coordinates (optional)
    pub z: Option<Vec<T>>,
    /// Field values
    pub field: Vec<T>,
    /// Field name
    pub name: String,
}

/// Time series data structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeSeriesData<T: RealField + Copy> {
    /// Time values
    pub time: Vec<T>,
    /// Data columns
    pub data: Vec<Vec<T>>,
    /// Column headers
    pub headers: Vec<String>,
}

impl<T: RealField + Copy> FieldData<T> {
    /// Create new 2D field data
    pub fn new_2d(x: Vec<T>, y: Vec<T>, field: Vec<T>, name: String) -> Self {
        Self {
            x,
            y,
            z: None,
            field,
            name,
        }
    }

    /// Create new 3D field data
    pub fn new_3d(x: Vec<T>, y: Vec<T>, z: Vec<T>, field: Vec<T>, name: String) -> Self {
        Self {
            x,
            y,
            z: Some(z),
            field,
            name,
        }
    }

    /// Check if data is 3D
    pub fn is_3d(&self) -> bool {
        self.z.is_some()
    }

    /// Get number of points
    pub fn num_points(&self) -> usize {
        self.field.len()
    }
}

impl<T: RealField + Copy> TimeSeriesData<T> {
    /// Create new time series data
    pub fn new(time: Vec<T>, data: Vec<Vec<T>>, headers: Vec<String>) -> Self {
        Self {
            time,
            data,
            headers,
        }
    }

    /// Get number of time steps
    pub fn num_timesteps(&self) -> usize {
        self.time.len()
    }

    /// Get number of variables
    pub fn num_variables(&self) -> usize {
        self.data.len()
    }

    /// Get data for a specific variable
    pub fn get_variable(&self, index: usize) -> Option<&Vec<T>> {
        self.data.get(index)
    }

    /// Get data for a specific time step
    pub fn get_timestep(&self, index: usize) -> Option<Vec<T>> {
        if index >= self.time.len() {
            return None;
        }

        Some(self.data.iter().map(|col| col[index]).collect())
    }
}
