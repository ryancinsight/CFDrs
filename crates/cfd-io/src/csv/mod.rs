//! CSV file I/O operations.
//!
//! This module provides CSV file support for CFD simulation data
//! with streaming and iterator-based processing for efficiency.

mod reader;
mod streaming;
mod types;
mod writer;

pub use reader::{CsvIterator, CsvReader};
pub use streaming::{StreamingReader, StreamingWriter};
pub use types::{CsvRecord, FieldData, TimeSeriesData};
pub use writer::CsvWriter;

// Re-export for backward compatibility (will deprecate in future)
pub use reader::read_field_data;
pub use writer::write_time_series;
