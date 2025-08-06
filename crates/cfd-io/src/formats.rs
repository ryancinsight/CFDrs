//! File format definitions and utilities.

/// Supported file formats
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    /// VTK format
    Vtk,
    /// CSV format
    Csv,
    /// JSON format
    Json,
    /// Binary format
    Binary,
}